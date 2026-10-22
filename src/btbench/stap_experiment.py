"""Frozen CPU study of angle-time STAP on explicit, non-acoustic IQ surrogates.

The independent calibration packet is available to both STAP and a demodulated
SVD control. No filter is passed blood/tissue truth or prescribed test velocity.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import itertools
import json
import platform
import time
import numpy as np
import scipy
from scipy import ndimage,stats
from threadpoolctl import threadpool_limits
from .core import Acquisition
from .filters import sliding_filter,temporal_highpass
from .stap import ViewFilterResult,local_stap_views,estimate_tissue_velocity,lcmv_velocity_spectrum,steering
from .measurements import vector_doppler,relative_error,wall_shear_scalar


@dataclass
class AngleScene:
    acquisition: Acquisition
    blood_views: np.ndarray
    tissue_views: np.ndarray
    noise_views: np.ndarray
    calibration_views: np.ndarray
    training_roi: np.ndarray
    geometry: dict
    velocity_true_zx_mps: np.ndarray


def angle_scene(seed=0,n_frames=96,shape=(8,25),prf_hz=1000.,f0_hz=5e6,
        umax_relative_mps=.02,clutter_db=30.,snr_db=20.,
        tissue_velocity_zx_mps=(.001,0.),tissue_acceleration_mps2=0.,
        tissue_velocity_spread_mps=0.,decorrelation=.995,angle_error_deg=0.,
        gain_phase_error_rad=0.,training_contamination=0.,training_radius=None,
        uniform_flow=False,phase_sign=1):
    """Common speckle envelope, exact bistatic Doppler phase at sequential views.

    This is a specified statistical-IQ forward model. It omits RF propagation,
    angle decorrelation, scatterer transport, attenuation, aberration and PSF
    mixing at the wall. Native-blood scattering is not claimed to be simulated.
    """
    if n_frames<16 or not 0<=decorrelation<=1 or not 0<=training_contamination<=1:
        raise ValueError("Invalid scene dimensions or decorrelation/training fraction")
    rng=np.random.default_rng(seed)
    z,x=np.linspace(0,.003,shape[0]),np.linspace(-.0015,.0015,shape[1])
    rr=np.broadcast_to(abs(x)[None,:],shape)
    flow=rr<.001-1e-12
    training=rr>.00115
    profile=np.ones(shape) if uniform_flow else np.maximum(1-(rr/.001)**2,0)
    base=np.asarray(tissue_velocity_zx_mps,float)
    velocity=np.zeros((*shape,2));velocity[...,0]=base[0]+umax_relative_mps*profile
    # Axial wall translation leaves the x-normal vessel geometry unchanged.
    tx=np.deg2rad(np.array([-15.,0.,15.]));rx=np.zeros(3)
    offsets=np.arange(3)/(3*prf_hz)
    geom=np.stack([np.cos(tx)+np.cos(rx),np.sin(tx)+np.sin(rx)],axis=1)
    gains=np.array([1.,.92*np.exp(.1j),1.08*np.exp(-.06j)])
    true_gains=gains*np.exp(1j*np.array([0,gain_phase_error_rad,-gain_phase_error_rad]))
    times=np.arange(n_frames)/prf_hz
    tj=times[:,None]+offsets[None,:]
    scale=phase_sign*2*np.pi*f0_hz/1540.
    def field():
        a=ndimage.gaussian_filter(rng.normal(size=shape)+1j*rng.normal(size=shape),.6)
        return a/np.sqrt(np.mean(abs(a)**2))
    state=field();envelope=[]
    for _ in times:
        state=decorrelation*state+np.sqrt(1-decorrelation**2)*field()
        envelope.append(state.copy())
    envelope=np.asarray(envelope)*flow
    doppler=np.einsum("vd,zxd->vzx",geom,velocity)
    blood=envelope[:,None]*np.exp(1j*scale*tj[:,:,None,None]*doppler[None])*true_gains[None,:,None,None]
    tissue_v=np.broadcast_to(base,(*shape,2)).copy()
    tissue_v[...,0]+=tissue_velocity_spread_mps*np.linspace(-.5,.5,shape[0])[:,None]
    tissue_fd=np.einsum("vd,zxd->vzx",geom,tissue_v)
    tissue_phase=scale*(tj[:,:,None,None]*tissue_fd[None]
        +.5*tissue_acceleration_mps2*tj[:,:,None,None]**2*geom[None,:,0,None,None])
    texture=field()*(.15+.85*(~flow))
    tissue=10**(clutter_db/20)*texture[None,None]*np.exp(1j*tissue_phase)*true_gains[None,:,None,None]
    noise=(rng.normal(size=blood.shape)+1j*rng.normal(size=blood.shape))/np.sqrt(2)*10**(-snr_db/20)
    # An independent calibration packet: no test mixture component is reused.
    cal_times=np.arange(-48,0)/prf_hz
    ct=cal_times[:,None]+offsets[None,:]
    cphase=scale*(ct[:,:,None,None]*tissue_fd[None]
        +.5*tissue_acceleration_mps2*ct[:,:,None,None]**2*geom[None,:,0,None,None])
    cal_noise=(rng.normal(size=(48,3,*shape))+1j*rng.normal(size=(48,3,*shape)))/np.sqrt(2)*10**(-snr_db/20)
    calibration=10**(clutter_db/20)*field()[None,None]*np.exp(1j*cphase)*true_gains[None,:,None,None]+cal_noise
    if training_contamination:
        # A deliberate annotation error, selected before filtering, to test
        # sensitivity to blood-bearing training pixels. Evaluation uses flow.
        lumen=np.flatnonzero(flow)
        n=int(np.ceil(training.sum()*training_contamination/max(1-training_contamination,1e-6)))
        chosen=rng.choice(lumen,min(n,len(lumen)),replace=False)
        training=training.copy();training.ravel()[chosen]=True
    wall_z=z[1:-1:2]
    wall={"positions_zx_m":[[float(zz),xx] for zz in wall_z for xx in (-.001,.001)],
        "inward_normals_zx":[[0,nn] for zz in wall_z for nn in (1,-1)],
        "distances_m":np.linspace(.0002,.0006,5).tolist(),"velocity_mps":float(base[0]),
        "viscosity_pa_s":.0035}
    actual=tx.copy();assumed=tx.copy();assumed[-1]+=np.deg2rad(angle_error_deg)
    geometry={"prf_hz":prf_hz,"f0_hz":f0_hz,"sound_speed_mps":1540.,"phase_sign":phase_sign,
        "tx_rad":assumed.tolist(),"rx_rad":rx.tolist(),"view_offsets_s":offsets.tolist(),
        "gains_real":gains.real.tolist(),"gains_imag":gains.imag.tolist(),
        "actual_tx_rad":actual.tolist()}
    a=Acquisition((blood+tissue+noise).mean(1),prf_hz,f0_hz,phase_sign=phase_sign,
        iq_views=blood+tissue+noise,z_m=z,x_m=x,tissue_roi=~flow,flow_mask=flow,
        near_wall_mask=flow&(rr>.65e-3),blood_true=blood.mean(1),tissue_true=tissue.mean(1),
        velocity_true_mps=velocity[...,0],wall=wall,
        metadata={"generator":"angle_time_statistical_IQ_surrogate","seed":seed,"contrast_agent":False,
            "not_acoustic_ground_truth":True,"absolute_time_origin_s":0.,"independent_calibration_frames":48,
            "umax_relative_mps":umax_relative_mps,"clutter_db":clutter_db,"snr_db":snr_db,
            "decorrelation":decorrelation,"tissue_acceleration_mps2":tissue_acceleration_mps2,
            "tissue_velocity_spread_mps":tissue_velocity_spread_mps,"angle_error_deg":angle_error_deg,
            "gain_phase_error_rad":gain_phase_error_rad,"training_contamination":training_contamination,
            "uniform_flow":uniform_flow})
    return AngleScene(a,blood,tissue,noise,calibration,training,geometry,velocity)


def view_svd(views,window=64,hop=32,selector="fixed",discard=1,upper=None,per_view=False,modulation=None):
    """All angular samples retained. Optional calibrated phase demodulation.

    Joint SVD folds views into the spatial dimension, not into slow time.
    The per-view arm fits an independent spatial/temporal SVD for each view.
    """
    t,a,z,x=views.shape
    mod=np.ones((t,a,1,1),complex) if modulation is None else modulation[:,:,None,None]
    aligned=views*mod.conj()
    if per_view:
        fits=[sliding_filter(aligned[:,j],window=window,hop=hop,selector=selector,discard=discard,upper=upper) for j in range(a)]
        out=np.stack([r.iq for r in fits],axis=1)*mod
        def apply(other):
            b=other*mod.conj()
            return np.stack([r.apply_frozen(b[:,j]) for j,r in enumerate(fits)],axis=1)*mod
        diag={"per_view":True,"fits":[r.diagnostics for r in fits]}
    else:
        fit=sliding_filter(aligned.reshape(t,a*z,x),window=window,hop=hop,selector=selector,discard=discard,upper=upper)
        out=fit.iq.reshape(views.shape)*mod
        def apply(other):
            b=(other*mod.conj()).reshape(t,a*z,x)
            return fit.apply_frozen(b).reshape(views.shape)*mod
        diag={"per_view":False,**fit.diagnostics}
    return ViewFilterResult(out,{**diag,"retained_input_views":a,"calibration_demodulated":modulation is not None},apply)


def _auc(power,flow):
    a=np.asarray(power)[flow];b=np.asarray(power)[~flow]
    if not len(a) or not len(b) or not np.isfinite(power).all():return None
    ranks=stats.rankdata(np.r_[a,b])
    return float((ranks[:len(a)].sum()-len(a)*(len(a)+1)/2)/(len(a)*len(b)))


def view_metrics(result,scene):
    a=scene.acquisition;g=scene.geometry;v=result.views
    velocity,vd=vector_doppler(v,g["tx_rad"],g["rx_rad"],a.prf_hz,a.f0_hz,
        a.sound_speed_mps,a.phase_sign)
    flow=a.flow_mask
    b=result.apply_frozen(scene.blood_views);t=result.apply_frozen(scene.tissue_views)
    n=result.apply_frozen(scene.noise_views)
    pb=max(float(np.sum(abs(b)**2)),1e-30)
    pt=max(float(np.sum(abs(scene.tissue_views)**2)),1e-30)
    original_b=max(float(np.sum(abs(scene.blood_views)**2)),1e-30)
    power=np.mean(abs(v)**2,axis=(0,1))
    axial=velocity[...,0]
    w=wall_shear_scalar(axial,a.z_m,a.x_m,a.wall)
    wr=wall_shear_scalar(a.velocity_true_mps,a.z_m,a.x_m,a.wall)
    m={"vector_velocity_nrmse":relative_error(velocity,scene.velocity_true_zx_mps,np.repeat(flow[...,None],2,axis=-1)),
        "axial_velocity_nrmse":relative_error(axial,a.velocity_true_mps,flow),
        "near_wall_axial_nrmse":relative_error(axial,a.velocity_true_mps,a.near_wall_mask),
        "velocity_invalid_fraction":float(np.mean(~np.isfinite(velocity[flow]))),
        "directional_wss_nrmse":relative_error(w,wr) if np.linalg.norm(wr)>1e-10 else None,
        "wss_invalid_fraction":float(np.mean(~np.isfinite(w))),
        "blood_iq_nrmse":float(np.linalg.norm(v-scene.blood_views)/np.sqrt(original_b)),
        "frozen_blood_power_retention":float(np.sum(abs(b)**2)/original_b),
        "frozen_tissue_power_retention":float(np.sum(abs(t)**2)/pt),
        "frozen_tissue_to_blood_power_ratio":float(np.sum(abs(t)**2)/pb),
        "frozen_noise_to_blood_power_ratio":float(np.sum(abs(n)**2)/pb),
        "frozen_mixture_linearity_error":float(np.linalg.norm(v-b-t-n)/max(np.linalg.norm(v),1e-30)),
        "power_doppler_auc":_auc(power,flow),"geometry_condition_number":vd["geometry_condition_number"]}
    return m


def default_stap_config():
    return {"version":1,"blas_threads":1,"seeds":[0,1,2],
        "generator":{"n_frames":96,"shape":[8,25],"prf_hz":1000.,"f0_hz":5e6,"decorrelation":.995},
        "factorial":{"umax_relative_mps":[.003,.02],"clutter_db":[20.,40.],"tissue_acceleration_mps2":[0.,.03]},
        "stress_cases":[
            {"label":"view_angle_error","angle_error_deg":3.},
            {"label":"view_phase_error","gain_phase_error_rad":.35},
            {"label":"training_contamination","training_contamination":.25},
            {"label":"velocity_spread","tissue_velocity_spread_mps":.003},
            {"label":"identical_velocity","uniform_flow":True,"umax_relative_mps":0.,"decorrelation":1.},
            {"label":"opposed_flow","umax_relative_mps":-.006}],
        "filter":{"taps":8,"hop":4,"guard":1,"training_radius":32,"loading_rel":.01},
        "baselines":{"long_window":64,"long_hop":32,"svd_discard":1,"butterworth_cutoff_hz":30.,
            "butterworth_order":3,"butterworth_phase_mode":"zero_phase","stap_long_taps":16,"stap_long_hop":8},
        "methods":["unfiltered","butterworth_30Hz","svd_joint_matched","svd_per_view_matched",
            "svd_joint_long","baranger_joint_long","calibrated_demodulated_svd","stap_orthogonal",
            "stap_covariance","stap_covariance_16taps","stap_training_starved_control"],
        "lcmv_candidates_zx_mps":[[float(z),0.] for z in np.linspace(-.005,.045,51)],
        "oracle_svd_grid":[{"discard":k,"upper":h} for k in range(9) for h in range(k,9)],
        "oracle_objective":"blood_iq_nrmse",
        "interpretation":"frozen synthetic study; calibration-assisted STAP and demodulated-SVD controls"}


def _clean(v):
    if isinstance(v,dict):return {k:_clean(w) for k,w in v.items()}
    if isinstance(v,(list,tuple)):return [_clean(w) for w in v]
    if isinstance(v,np.ndarray):return _clean(v.tolist())
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v


def run_stap_study(config,output):
    """Run every frozen case; preserve failures, controls and all candidate scores."""
    if config.get("oracle_objective","blood_iq_nrmse")!="blood_iq_nrmse":
        raise ValueError("STAP study oracle is explicitly defined for blood_iq_nrmse")
    bands=[{"discard":k,"upper":None} if isinstance(k,int) else k for k in config.get("oracle_svd_grid",[])]
    if any(int(b["discard"])!=b["discard"] or b["discard"]<0 or b["discard"]>config["filter"]["taps"]
           or (b["upper"] is not None and not b["discard"]<=b["upper"]<=config["filter"]["taps"]) for b in bands):
        raise ValueError("Oracle ranks must be within the matched snapshot length")
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError("Choose a new STAP output directory")
    out.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(_clean(config),sort_keys=True,indent=2,allow_nan=False)
    (out/"frozen_config.json").write_text(encoded)
    config_hash=hashlib.sha256(encoded.encode()).hexdigest()
    keys=list(config["factorial"])
    conditions=[dict(zip(keys,v)) for v in itertools.product(*(config["factorial"][k] for k in keys))]
    conditions=[{"label":"factorial",**c} for c in conditions]+config.get("stress_cases",[])
    records=[];begin=time.perf_counter()
    with threadpool_limits(limits=int(config.get("blas_threads",1))):
        for cond,seed in itertools.product(conditions,config["seeds"]):
            p={**config["generator"],**{k:v for k,v in cond.items() if k!="label"},"seed":seed}
            scene=angle_scene(**p)
            a=scene.acquisition;views=a.iq_views;g=scene.geometry;f=config["filter"];base=config["baselines"]
            # Calibration ROI is defined from an independently annotated wall
            # band; do not use the deliberately contaminated training mask here.
            cal_roi=np.broadcast_to(abs(a.x_m)[None,:]>.00115,a.iq.shape[1:])
            tv=estimate_tissue_velocity(scene.calibration_views,cal_roi,g["tx_rad"],g["rx_rad"],
                a.prf_hz,a.f0_hz,a.sound_speed_mps,a.phase_sign)
            gains=np.asarray(g["gains_real"])+1j*np.asarray(g["gains_imag"])
            skw={k:g[k] for k in ("prf_hz","f0_hz","tx_rad","rx_rad","sound_speed_mps","phase_sign","view_offsets_s")}
            skw.update(gains=gains,tissue_velocities_zx_mps=[tv.tolist()],training_roi=scene.training_roi,**f)
            mod=steering([tv],g["tx_rad"],g["rx_rad"],np.arange(len(views))/a.prf_hz,a.f0_hz,
                a.sound_speed_mps,a.phase_sign,g["view_offsets_s"]).reshape(len(views),views.shape[1])
            methods=dict([
                ("unfiltered",lambda:ViewFilterResult(views,{},lambda x:x)),
                ("butterworth_30Hz",lambda:(lambda r:ViewFilterResult(r.iq,r.diagnostics,r.apply_frozen))(
                    temporal_highpass(views,base["butterworth_cutoff_hz"],a.prf_hz,base["butterworth_order"],base["butterworth_phase_mode"]))),
                ("svd_joint_matched",lambda:view_svd(views,window=f["taps"],hop=f["hop"],discard=base["svd_discard"])),
                ("svd_per_view_matched",lambda:view_svd(views,window=f["taps"],hop=f["hop"],discard=base["svd_discard"],per_view=True)),
                ("svd_joint_long",lambda:view_svd(views,window=base["long_window"],hop=base["long_hop"],discard=base["svd_discard"])),
                ("baranger_joint_long",lambda:view_svd(views,window=base["long_window"],hop=base["long_hop"],selector="baranger")),
                ("calibrated_demodulated_svd",lambda:view_svd(views,window=f["taps"],hop=f["hop"],discard=base["svd_discard"],modulation=mod)),
                ("stap_orthogonal",lambda:local_stap_views(views,**skw)),
                ("stap_covariance",lambda:local_stap_views(views,adaptive=True,**skw)),
                ("stap_covariance_16taps",lambda:local_stap_views(views,adaptive=True,
                    **{**skw,"taps":base["stap_long_taps"],"hop":base["stap_long_hop"]})),
                ("stap_training_starved_control",lambda:local_stap_views(views,adaptive=True,
                    **{**skw,"training_radius":2,"guard":1}))])
            rows=[]
            for name in config["methods"]:
                fn=methods[name]
                tic=time.perf_counter();row={"method":name,"uses_test_truth_for_selection":False}
                try:
                    result=fn()
                    row.update(status="ok",metrics=view_metrics(result,scene),diagnostics=result.diagnostics)
                except (ValueError,MemoryError,np.linalg.LinAlgError) as exc:
                    row.update(status="failed",reason=str(exc))
                row["wall_time_s"]=time.perf_counter()-tic
                rows.append(row)
            # This positive control estimates the finite-ensemble readout floor.
            br=ViewFilterResult(scene.blood_views,{},lambda x:x)
            vref,_=vector_doppler(br.views,g["tx_rad"],g["rx_rad"],a.prf_hz,a.f0_hz,phase_sign=a.phase_sign)
            reference={"evaluation_only":True,"axial_velocity_nrmse":relative_error(vref[...,0],a.velocity_true_mps,a.flow_mask)}
            oracle_curve=[];best=None
            for band in bands:
                result=view_svd(views,window=f["taps"],hop=f["hop"],**band)
                loss=float(np.linalg.norm(result.views-scene.blood_views)/max(np.linalg.norm(scene.blood_views),1e-30))
                oracle_curve.append({**band,"loss":loss})
                if best is None or loss<best[0]:best=(loss,band,result)
            if best is not None:
                rows.append({"method":"oracle_svd_joint_matched","status":"ok","uses_test_truth_for_selection":True,
                    "objective":"blood_iq_nrmse","selected_band":best[1],"candidate_curve":oracle_curve,
                    "metrics":view_metrics(best[2],scene),"diagnostics":best[2].diagnostics})
            center=(len(a.z_m)//2,len(a.x_m)//2)
            det=lcmv_velocity_spectrum(views,config["lcmv_candidates_zx_mps"],[tv.tolist()],
                g["tx_rad"],g["rx_rad"],a.prf_hz,a.f0_hz,scene.training_roi,center,
                taps=f["taps"],hop=f["hop"],guard=f["guard"],loading_rel=f["loading_rel"],
                view_offsets_s=g["view_offsets_s"],gains=gains,phase_sign=a.phase_sign)
            target=scene.velocity_true_zx_mps[center]
            selected=det["selected_velocity_zx_mps"]
            det["evaluation_velocity_error_mps"]=None if selected is None else float(np.linalg.norm(np.asarray(selected)-target))
            det["grid_quantization_floor_mps"]=float(np.min(np.linalg.norm(np.asarray(config["lcmv_candidates_zx_mps"])-target,axis=1)))
            record={"case":cond,"seed":seed,"scene_parameters":p,"geometry":g,
                "views_sha256":hashlib.sha256(views.tobytes()).hexdigest(),
                "calibration_sha256":hashlib.sha256(scene.calibration_views.tobytes()).hexdigest(),
                "estimated_tissue_velocity_zx_mps":tv.tolist(),"independent_calibration":True,
                "calibration_control":reference,"results":rows,"lcmv_detector":det}
            records.append(record)
            print(json.dumps({"stap_cases_completed":len(records),"case":cond["label"],"seed":seed}),flush=True)
    code=hashlib.sha256()
    for file in sorted(Path(__file__).parent.glob("*.py")):
        code.update(file.name.encode());code.update(file.read_bytes())
    report={"config_sha256":config_hash,"code_sha256":code.hexdigest(),"records":records,
        "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
            "blas_threads":config.get("blas_threads",1),"elapsed_s":time.perf_counter()-begin},
        "scope":"Angle-time statistical IQ and fixed-wall 2-D directional shear; no acoustic/experimental validation",
        "fairness":"All filters receive all views; calibrated demodulated SVD receives the same independent tissue calibration as STAP"}
    (out/"results.json").write_text(json.dumps(_clean(report),indent=2,allow_nan=False))
    lines=["# First angle-time STAP study","",
        "Frozen statistical-IQ study. All methods retain the same views. Calibration is independently acquired.",
        "Oracle SVD uses same-case blood truth and is evaluation-only. LCMV is a separate single-target detector.","",
        "| Method | Completed | Median blood retention | Median tissue retention | Median axial NRMSE | Median WSS NRMSE |","|---|---:|---:|---:|---:|---:|"]
    ids=list(dict.fromkeys(r["method"] for c in records for r in c["results"]))
    for name in ids:
        rows=[r for c in records for r in c["results"] if r["method"]==name and r["status"]=="ok"]
        def median(key):
            vals=[r["metrics"].get(key) for r in rows]
            vals=[v for v in vals if v is not None]
            return "n/a" if not vals else f"{np.median(vals):.4g}"
        lines.append(f"| {name} | {len(rows)}/{len(records)} | {median('frozen_blood_power_retention')} | {median('frozen_tissue_power_retention')} | {median('axial_velocity_nrmse')} | {median('directional_wss_nrmse')} |")
    lines.extend(["","Medians pool deliberately different conditions, including failure controls. Use the condition records for interpretation.",
        "Velocity spread, acceleration, calibration error, equal blood/tissue velocity and insufficient training are explicit stress cases."])
    (out/"SUMMARY.md").write_text("\n".join(lines)+"\n")
    return report
