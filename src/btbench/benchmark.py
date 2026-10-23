"""Explicit method registry; filters receive IQ/settings, never reference fields."""
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

from .filters import svd_filter,aer_filter,local_filter,rpca_filter,temporal_highpass,sliding_filter
from .axes import motion_packet_filter,multilag_filter,ica_filter
from .measurements import evaluate,frozen_component_metrics,evaluate_windows,injection_response
from .svd_literature import randomized_svd_filter,huang_evd_filter,gaussian_local_evd,song_local_filter
from .core import FilterResult


def run_filter(iq,spec,tissue_roi=None,views=None,prf_hz=None,view_geometry=None):
    p=dict(spec.get("params",{}))
    method=spec["method"]
    if prf_hz is not None:
        if method=="highpass":
            if "prf_hz" in p and not np.isclose(p["prf_hz"],prf_hz,rtol=1e-9,atol=0):
                raise ValueError("Filter PRF differs from the acquisition revisit rate")
            p.setdefault("prf_hz",prf_hz)
        if method=="song_local" or p.get("selector") in ("song","temporal_frequency"):
            options=dict(p.get("selector_options",{}))
            if "prf_hz" in options and not np.isclose(options["prf_hz"],prf_hz,rtol=1e-9,atol=0):
                raise ValueError("Spectral selector PRF differs from the acquisition revisit rate")
            options.setdefault("prf_hz",prf_hz)
            p["selector_options"]=options
    if method=="stap":
        from .stap import local_stap
        if views is None:
            raise ValueError("STAP requires uncompounded IQ views")
        if not np.allclose(iq,np.mean(views,axis=1)):
            raise ValueError("STAP comparison needs iq to be the same coherent view mean")
        for key,value in (view_geometry or {}).items():
            if key in p and not np.allclose(p[key],value,rtol=1e-5,atol=1e-10):
                raise ValueError("STAP geometry differs from acquisition: "+key)
            p.setdefault(key,value)
        return local_stap(views,**p)
    if method=="view_svd":
        from .stap_experiment import view_svd
        from .stap import steering
        if views is None or not np.allclose(iq,np.mean(views,axis=1)):
            raise ValueError("View SVD requires views and their coherent-mean baseline")
        calibration=p.pop("calibration_velocity_zx_mps",None)
        modulation=None
        if calibration is not None:
            g=dict(view_geometry or {})
            required=("tx_rad","rx_rad","prf_hz","f0_hz")
            if any(k not in g for k in required):raise ValueError("Calibrated view SVD needs acquisition geometry")
            modulation=steering(calibration,g["tx_rad"],g["rx_rad"],
                g.get("start_time_s",0.)+np.arange(len(views))/g["prf_hz"],g["f0_hz"],
                g.get("sound_speed_mps",1540.),g.get("phase_sign",1),g.get("view_offsets_s")).reshape(len(views),views.shape[1])
        r=view_svd(views,modulation=modulation,**p)
        return FilterResult(r.views.mean(axis=1),r.diagnostics,
            apply_frozen_views=lambda x:r.apply_frozen(x).mean(axis=1))
    if method=="identity":
        return FilterResult(np.array(iq,copy=True),{"control":"unfiltered mixture"},lambda x:np.array(x,copy=True))
    registry={"svd":svd_filter,"aer":aer_filter,"local":local_filter,"rpca":rpca_filter,
              "axis1":motion_packet_filter,"axis2":multilag_filter,"ica":ica_filter,
              "highpass":temporal_highpass,"sliding":sliding_filter,
              "rsvd":randomized_svd_filter,"song_local":song_local_filter,
              "huang_evd":huang_evd_filter,"gaussian_evd":gaussian_local_evd}
    if method not in registry:
        raise ValueError("Unknown IQ method "+method)
    if method in ("axis1","axis2","ica"):
        p["tissue_roi"]=tissue_roi
    return registry[method](iq,**p)


def _run_acquisition(a,spec):
    geometry={"prf_hz":a.prf_hz,"f0_hz":a.f0_hz,"sound_speed_mps":a.sound_speed_mps,
              "phase_sign":a.phase_sign,"start_time_s":0. if a.timestamps_s is None else float(a.timestamps_s[0])}
    for name,key in (("tx_angles_rad","tx_rad"),("rx_angles_rad","rx_rad"),("view_offsets_s","view_offsets_s")):
        value=getattr(a,name)
        if value is not None:geometry[key]=value
    return run_filter(a.iq,spec,a.tissue_roi,a.iq_views,a.prf_hz,geometry)


def grid_specs(base,grid):
    keys=list(grid)
    for values in itertools.product(*(grid[k] for k in keys)):
        spec={**base,"params":{**base.get("params",{}),**dict(zip(keys,values))}}
        yield spec


class OracleInfeasible(ValueError):
    def __init__(self,diagnostics):
        super().__init__("No finite oracle candidate satisfies the declared preservation constraints")
        self.diagnostics=diagnostics


def pareto_front(rows,x="frozen_blood_iq_nrmse",y="frozen_tissue_power_retention"):
    """Indices of nondominated finite points (both losses minimized)."""
    points=[]
    for i,row in enumerate(rows):
        m=row.get("metrics",{})
        v=[m.get(x),m.get(y)]
        if all(a is not None and np.isfinite(a) for a in v):
            points.append((i,np.asarray(v)))
    return [i for i,p in points if not any(
        np.all(q<=p+1e-12) and np.any(q<p-1e-12) for j,q in points if i!=j)]


def oracle(acquisition,base,grid,objective="blood_iq_nrmse",geometry_factor=2.,constraints=None):
    """Same-case truth-assisted bound over the DECLARED candidate grid only.

    Never used as a deployable setting. Every candidate is logged. Separate
    objectives are needed for velocity, IQ, and WSS; no implicit equivalence.
    """
    frozen=objective.startswith("frozen_") or bool(constraints)
    if frozen and (acquisition.blood_true is None or acquisition.tissue_true is None):
        raise ValueError("Separation oracle requires independent blood and tissue components")
    if not frozen and acquisition.velocity_true_mps is None and objective!="blood_iq_nrmse":
        raise ValueError("Oracle requires independent velocity/WSS truth")
    if objective=="blood_iq_nrmse" and acquisition.blood_true is None:
        raise ValueError("IQ oracle requires independent blood-IQ truth")
    grid=dict(grid)
    if grid.get("discard")=="all":
        p=int(np.prod(acquisition.iq.shape[1:]))
        t=min(base.get("params",{}).get("window") or len(acquisition.iq),len(acquisition.iq))
        grid["discard"]=list(range(min(p,t)+1))
    constraints=constraints or {}
    for limits in constraints.values():
        if not isinstance(limits,dict) or not limits or set(limits)-{"min","max"} or not all(np.isfinite(v) for v in limits.values()):
            raise ValueError("Constraints must contain finite min/max metric bounds")
    best=None;curve=[]
    for spec in grid_specs(base,grid):
        try:
            result=_run_acquisition(acquisition,spec)
        except (ValueError,MemoryError,np.linalg.LinAlgError) as exc:
            curve.append({"params":spec["params"],"loss":None,"status":"invalid","reason":str(exc)})
            continue
        metrics=evaluate(result.iq,acquisition,geometry_factor)
        if acquisition.blood_true is not None and acquisition.tissue_true is not None:
            metrics.update(frozen_component_metrics(result,acquisition))
        loss=metrics.get(objective)
        if loss is None:
            raise ValueError("Oracle objective unavailable: "+objective)
        feasible=bool(np.isfinite(loss))
        for key,limits in constraints.items():
            v=metrics.get(key)
            feasible=feasible and v is not None and np.isfinite(v)
            if feasible:
                feasible=(v>=limits.get("min",-np.inf) and v<=limits.get("max",np.inf))
        curve.append({"params":spec["params"],"loss":float(loss),"feasible":bool(feasible),"metrics":metrics})
        if feasible and (best is None or loss<best[0]):
            best=(loss,result,spec,metrics)
    diag={"oracle":True,"selection_uses_test_truth":True,"objective":objective,
          "constraints":constraints,"selected":None if best is None else best[2],"curve":curve,
          "feasible_candidates":sum(c.get("feasible",False) for c in curve),
          "pareto_candidate_indices":pareto_front(curve)}
    if best is None:
        raise OracleInfeasible(diag)
    return best[1],diag


def standard_methods(packet=64,long_window=256,prf_hz=1000.):
    """Factorial controls; short window and long window report distinct latency."""
    methods=[]
    for centered in (False,True):
        tag="centered" if centered else "uncentered"
        for sel in ("fixed","energy","baranger"):
            methods.append({"id":f"svd_{sel}_{tag}","method":"svd","params":{
                "centered":centered,"window":packet,"selector":sel,"discard":1}})
        methods.append({"id":f"svd_long_{tag}","method":"sliding","params":{
            "method":"svd","centered":centered,"window":long_window,"hop":packet,"discard":1}})
        methods.append({"id":f"svd_local_{tag}","method":"local","params":{
            "method":"svd","centered":centered,"patch_shape":[6,6],"stride":[3,3],"window":packet,"discard":1}})
        methods.append({"id":f"song_local_reconstruction_{tag}","method":"song_local","params":{
            "centered":centered,"patch_shape":[6,6],"stride":[3,3],"window":packet,
            "selector_options":{"prf_hz":prf_hz,"slope_threshold":0.1,"tail_tolerance":0.1}}})
        methods.append({"id":f"svd_temporal_frequency_{tag}","method":"svd","params":{
            "centered":centered,"window":packet,"selector":"temporal_frequency",
            "selector_options":{"prf_hz":prf_hz,"cutoff_hz":.03*prf_hz}}})
        methods.append({"id":f"evd_fixed_{tag}","method":"svd","params":{
            "centered":centered,"window":packet,"discard":1,"backend":"gram"}})
        methods.append({"id":f"song_rsvd_{tag}","method":"rsvd","params":{
            "centered":centered,"window":packet,"discard":1,"seed":0}})
        methods.append({"id":f"song_rsvd_downsampled_{tag}","method":"rsvd","params":{
            "centered":centered,"window":packet,"discard":1,"seed":0,"spatial_fraction":.5}})
        for rule,beta in (("median",1.43),("ncc",1.43),("sigmoid",2.)):
            methods.append({"id":f"huang_evd_{rule}_{tag}","method":"huang_evd","params":{
                "centered":centered,"window":packet,"rule":rule,"beta":beta}})
        methods.append({"id":f"huang_gaussian_iq_extension_{tag}","method":"gaussian_evd","params":{
            "centered":centered,"window":packet,"grid_shape":[2,2],"rule":"ncc"}})
        for domain in ("temporal","spatial"):
            for alpha in (0.5,0.25):
                methods.append({"id":f"aer_{domain}_alpha{alpha}_{tag}","method":"aer",
                    "params":{"domain":domain,"centered":centered,"epsilon_rel":1e-3,"alpha":alpha,"window":long_window}})
        methods.append({"id":f"aer_noise_floor_{tag}","method":"aer","params":{
            "centered":centered,"domain":"temporal","epsilon_rule":"noise_floor","window":long_window}})
        methods.append({"id":f"axis1_{tag}","method":"axis1","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3]}})
        methods.append({"id":f"axis1_rigid_{tag}","method":"axis1","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],"rigid":True}})
        methods.append({"id":f"axis2_{tag}","method":"axis2","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],
            "rank":6,"lags":[1,2,4,8,16],"use_memory":False}})
        methods.append({"id":f"axis2_with_memory_{tag}","method":"axis2","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],
            "rank":6,"lags":[1,2,4,8,16],"use_memory":True}})
    methods += [
        {"id":"ica_circular_complex","method":"ica","params":{"packet":packet,"patch_shape":[6,6],"stride":[3,3],
            "rank":6,"centered":True,"ica_nonconvergence":"svd_fallback"}},
        {"id":"rpca_pcp","method":"rpca","params":{}},
    ]
    for fraction in (.01,.03,.08):
        methods.append({"id":f"butterworth_hp_{fraction:g}_prf","method":"highpass","params":{
            "cutoff_hz":fraction*prf_hz,"prf_hz":prf_hz,"order":3,"phase_mode":"zero_phase"}})
    methods.append({"id":"butterworth_hp_causal_0.03_prf","method":"highpass","params":{
        "cutoff_hz":.03*prf_hz,"prf_hz":prf_hz,"order":3,"phase_mode":"causal"}})
    return methods


def run_benchmark(acquisition,config,output):
    out=Path(output)
    if (out/"results.json").exists():
        raise FileExistsError("Output already contains results; choose a new directory")
    out.mkdir(parents=True,exist_ok=True)
    a=acquisition
    results=[]
    start_all=time.perf_counter()
    threads=int(config.get("blas_threads",1))
    with threadpool_limits(limits=threads):
        for spec in config["methods"]:
            t0=time.perf_counter()
            row={"id":spec["id"],"method":spec["method"],"oracle":False,"spec":spec}
            try:
                if "oracle_grid" in spec:
                    result,diag=oracle(a,spec,spec["oracle_grid"],spec.get("objective","blood_iq_nrmse"),
                        config.get("geometry_factor",2.),spec.get("constraints"))
                    row.update(diag)
                else:
                    result=_run_acquisition(a,spec)
                metrics=evaluate(result.iq,a,config.get("geometry_factor",2.))
                if getattr(result,"power_doppler",None) is not None:
                    pd=result.power_doppler
                    row["paper_power_doppler"]={"mean":float(pd.mean()),"max":float(pd.max()),
                        "interpretation":"separate incoherent power aggregation; no phase/velocity output",
                        "flow_mean":float(pd[a.flow_mask].mean()) if a.flow_mask is not None else None,
                        "tissue_mean":float(pd[a.tissue_roi].mean()) if a.tissue_roi is not None else None}
                metrics.update(frozen_component_metrics(result,a))
                row.update(status="ok",metrics=metrics,diagnostics=result.diagnostics)
                if spec["method"]=="rpca" and a.blood_true is not None and config.get("injection_amplitudes"):
                    fn=lambda x:run_filter(x,spec,a.tissue_roi,prf_hz=a.prf_hz)
                    row["counterfactual_injection"]=injection_response(fn,a.iq,a.blood_true,config["injection_amplitudes"])
                row["readout"]=evaluate_windows(result.iq,a,config.get("measurement_window",64),
                    config.get("measurement_hop",64),config.get("geometry_factor",2.))
                if config.get("save_iq",False):
                    # IDs are labels, never paths supplied by external datasets.
                    name=hashlib.sha256(spec["id"].encode()).hexdigest()[:16]
                    np.save(out/f"filtered_{name}.npy",result.iq)
            except (ValueError,MemoryError,np.linalg.LinAlgError) as exc:
                row.update(status="infeasible" if isinstance(exc,OracleInfeasible) else "skipped",reason=str(exc))
                if isinstance(exc,OracleInfeasible):row.update(exc.diagnostics)
            row["wall_time_s"]=time.perf_counter()-t0
            try:
                import resource,sys
                peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                row["process_high_water_rss_bytes"]=int(peak if sys.platform=="darwin" else peak*1024)
            except ImportError:
                row["process_high_water_rss_bytes"]=None
            results.append(row)
    code=hashlib.sha256()
    for module in sorted(Path(__file__).parent.glob("*.py")):
        code.update(module.name.encode());code.update(module.read_bytes())
    report={"config":config,"results":results,"acquisition":a.metadata,"code_sha256":code.hexdigest(),
            "iq_shape":list(a.iq.shape),"input_sha256":hashlib.sha256(a.iq.tobytes()).hexdigest(),
            "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
            "platform":platform.platform(),"blas_threads":threads,"elapsed_s":time.perf_counter()-start_all},
            "primary_metrics":["frozen_blood_iq_nrmse","frozen_tissue_power_retention",
                "frozen_blood_spectral_tv","frozen_noise_power_gain"],
            "data_only_pareto_ids":[results[i]["id"] for i in pareto_front([
                r if not r.get("oracle") else {} for r in results])],
            "interpretation":"Component fidelity requires independent source truth; downstream readouts are separate endpoints"}
    # Strict JSON: invalid numerical estimates become null with explicit fractions
    # and candidate selection records, not nonstandard Infinity tokens.
    def clean(v):
        if isinstance(v,dict): return {k:clean(w) for k,w in v.items()}
        if isinstance(v,(tuple,list)): return [clean(w) for w in v]
        if isinstance(v,(float,np.floating)) and not np.isfinite(v): return None
        if isinstance(v,np.generic): return v.item()
        return v
    (out/"results.json").write_text(json.dumps(clean(report),indent=2,allow_nan=False))
    lines=["# Benchmark run", "", "Component metrics use the mixture-fitted operator, with no fitted gain or phase correction.",
           "Experimental anatomy and processed velocity maps do not supply independent component truth.","",
           "| Method | Status | Blood IQ distortion | Tissue power retention | Blood power retention | Spectral TV | Noise gain | Time (s) |",
           "|---|---|---:|---:|---:|---:|---:|---:|"]
    for row in results:
        m=row.get("metrics",{})
        def fmt(key):
            v=m.get(key)
            return "—" if v is None else f"{v:.3g}"
        lines.append(f"| {row['id']} | {row['status']} | {fmt('frozen_blood_iq_nrmse')} | {fmt('frozen_tissue_power_retention')} | {fmt('frozen_blood_power_retention')} | {fmt('frozen_blood_spectral_tv')} | {fmt('frozen_noise_power_gain')} | {row['wall_time_s']:.3g} |")
    lines.extend(["","Data-only Pareto set (this case and declared settings only): "+", ".join(report["data_only_pareto_ids"]),
        "","## Downstream Doppler and directional WSS","",
        "| Method | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE |","|---|---:|---:|---:|"])
    for row in results:
        m=row.get("metrics",{})
        values=[m.get(k) for k in ("velocity_nrmse","near_wall_velocity_nrmse","wss_nrmse")]
        if any(v is not None for v in values):
            lines.append("| "+row["id"]+" | "+" | ".join("—" if v is None else f"{v:.3g}" for v in values)+" |")
    (out/"SUMMARY.md").write_text("\n".join(lines)+"\n")
    return report
