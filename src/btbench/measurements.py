"""Shared measurement chain, with explicit geometry and invalid-estimate handling."""
import numpy as np
from dataclasses import replace
from scipy.interpolate import RegularGridInterpolator


def lag_one(iq, min_coherence=0.0):
    cross = np.mean(iq[1:]*iq[:-1].conj(),axis=0)
    p1 = np.mean(np.abs(iq[1:])**2,axis=0)
    p0 = np.mean(np.abs(iq[:-1])**2,axis=0)
    coherence = np.abs(cross)/np.maximum(np.sqrt(p1*p0),1e-30)
    phase = np.angle(cross)
    phase[(p1*p0 <= 1e-50) | (coherence < min_coherence)] = np.nan
    return phase, coherence


def doppler_velocity(iq, prf_hz, f0_hz, sound_speed_mps=1540.0, phase_sign=1,
                     geometry_factor=2.0, min_coherence=0.0):
    """Scalar velocity along known direction d: geometry_factor=(n_tx+n_rx).d.

    Default 2 measures axial velocity for normal pulse-echo. No angle inferred
    from intensity, no de-aliasing, no post-hoc gain correction.
    """
    if phase_sign not in (-1,1) or abs(geometry_factor) < 1e-3:
        raise ValueError("Invalid phase sign or unobservable flow direction")
    phase, coherence = lag_one(iq,min_coherence)
    velocity = phase_sign*phase*sound_speed_mps*prf_hz/(2*np.pi*f0_hz*geometry_factor)
    return velocity, coherence


def vector_doppler(iq_views, tx_angles_rad, rx_angles_rad, prf_hz, f0_hz,
                   sound_speed_mps=1540., phase_sign=1):
    """iq_views: (time,view,z,x), angles measured from z; output (z,x,2)=(vz,vx).

    Each view must have its own correct revisit rate and constant time offset.
    This compact function requires a common uniform revisit PRF across views.
    No motion compensation of sequential transmit offsets is presumed.
    """
    tx, rx = np.asarray(tx_angles_rad),np.asarray(rx_angles_rad)
    if iq_views.ndim != 4 or len(tx) != iq_views.shape[1] or tx.shape != rx.shape:
        raise ValueError("View/angle shape mismatch")
    a = np.stack((np.cos(tx)+np.cos(rx),np.sin(tx)+np.sin(rx)),axis=1)
    if np.linalg.matrix_rank(a) < 2 or np.linalg.cond(a) > 100:
        raise ValueError("Vector velocity not observable with these angles")
    b = np.stack([doppler_velocity(iq_views[:,j],prf_hz,f0_hz,sound_speed_mps,
                                  phase_sign,geometry_factor=1)[0]
                  for j in range(len(tx))])
    v = np.einsum("dv,vzx->zxd",np.linalg.pinv(a),b)
    return v, {"geometry_condition_number":float(np.linalg.cond(a))}


def wall_shear_scalar(velocity, z_m, x_m, wall, fit="free", degree=2):
    """Known-direction scalar WSS at supplied walls, fluid-side local polynomial.

    wall: positions_zx_m (W,2), inward_normals_zx (W,2), distances_m (K,),
    viscosity_pa_s, optional velocity_mps (W,). Positive n enters fluid.
    This is a 2-D directional shear comparator, NOT full vector/tensor WSS.
    A moving/deforming boundary requires appropriately tracked geometry.
    """
    pos = np.asarray(wall["positions_zx_m"],float)
    normal = np.asarray(wall["inward_normals_zx"],float)
    d = np.asarray(wall["distances_m"],float)
    if pos.shape != normal.shape or pos.ndim != 2 or pos.shape[1] != 2:
        raise ValueError("Wall positions/normals must be (W,2), (z,x)")
    if np.any(d <= 0) or not np.allclose(np.linalg.norm(normal,axis=1),1):
        raise ValueError("Need strictly fluid-sided distances and unit normals")
    if wall["viscosity_pa_s"] <= 0 or degree < 1 or fit not in ("free","no_slip"):
        raise ValueError("Invalid WSS fit")
    interp = RegularGridInterpolator((z_m,x_m),velocity,bounds_error=False,fill_value=np.nan)
    vals = interp(pos[:,None,:]+normal[:,None,:]*d[None,:,None])
    wall_v = np.broadcast_to(np.asarray(wall.get("velocity_mps",0.)),(len(pos),))
    # Scale distances for numerical conditioning; restore SI gradient at end.
    scale = d.max()
    columns = range(1,degree+1) if fit == "no_slip" else range(degree+1)
    design = np.stack([(d/scale)**j for j in columns],axis=1)
    out = np.full(len(pos),np.nan)
    for i,row in enumerate(vals):
        ok = np.isfinite(row)
        if ok.sum() < design.shape[1] or np.linalg.matrix_rank(design[ok]) < design.shape[1]:
            continue
        target = row-wall_v[i] if fit == "no_slip" else row
        coef = np.linalg.lstsq(design[ok],target[ok],rcond=None)[0]
        gradient = coef[0 if fit == "no_slip" else 1]/scale
        out[i] = wall["viscosity_pa_s"]*gradient
    return out


def relative_error(a,b,mask=None):
    """Invalid predictions inside valid truth mask incur failure, not silent deletion."""
    mask = np.isfinite(b) if mask is None else np.asarray(mask,bool)&np.isfinite(b)
    if not mask.any():
        return None
    if not np.isfinite(np.asarray(a)[mask]).all():
        return float("inf")
    return float(np.linalg.norm((np.asarray(a)-b)[mask])/max(np.linalg.norm(b[mask]),1e-30))


def evaluate(iq, acquisition, geometry_factor=2., fit="free"):
    a = acquisition
    if a.velocity_true_mps is not None and a.velocity_true_mps.ndim==3:
        # Readout compares a window's Doppler estimate with its prescribed
        # temporal mean velocity; amplitude-dependent Doppler bias remains visible.
        a=replace(a,velocity_true_mps=np.mean(a.velocity_true_mps,axis=0))
    v,coh = doppler_velocity(iq,a.prf_hz,a.f0_hz,a.sound_speed_mps,a.phase_sign,geometry_factor)
    mask = a.flow_mask if a.flow_mask is not None else np.ones(v.shape,bool)
    out = {"coherence_flow_median":float(np.median(coh[mask])) if mask.any() else None,
           "velocity_invalid_fraction":float(np.mean(~np.isfinite(v[mask]))) if mask.any() else None,
           "power_output":float(np.mean(np.abs(iq)**2))}
    if a.velocity_true_mps is not None:
        out["velocity_nrmse"] = relative_error(v,a.velocity_true_mps,mask)
        if a.near_wall_mask is not None:
            out["near_wall_velocity_nrmse"] = relative_error(v,a.velocity_true_mps,a.near_wall_mask)
        out["velocity_mean_bias_mps"] = float(np.mean(v[mask]-a.velocity_true_mps[mask]))
        # Peak bias is based on each map's own maximum, with no chosen reference pixel.
        out["velocity_peak_bias_mps"] = float(np.max(v[mask])-np.max(a.velocity_true_mps[mask]))
    if a.blood_true is not None:
        out["blood_iq_nrmse"] = float(np.linalg.norm(iq-a.blood_true)/max(np.linalg.norm(a.blood_true),1e-30))
        out["output_to_blood_power_ratio"] = float(np.sum(abs(iq)**2)/max(np.sum(abs(a.blood_true)**2),1e-30))
        # This ratio is NOT pure blood retention: leakage contributes to it.
    if a.tissue_roi is not None and a.tissue_roi.any():
        out["tissue_roi_output_power_ratio"] = float(np.sum(abs(iq[:,a.tissue_roi])**2)/max(np.sum(abs(a.iq[:,a.tissue_roi])**2),1e-30))
    if a.wall is not None and a.z_m is not None and a.x_m is not None:
        wss = wall_shear_scalar(v,a.z_m,a.x_m,a.wall,fit)
        out["wss_invalid_fraction"] = float(np.mean(~np.isfinite(wss)))
        if a.velocity_true_mps is not None:
            ref = wall_shear_scalar(a.velocity_true_mps,a.z_m,a.x_m,a.wall,fit)
            out["wss_nrmse"] = relative_error(wss,ref)
        if "true_pa" in a.wall:
            out["wss_analytic_nrmse"] = relative_error(wss,np.asarray(a.wall["true_pa"]))
    return out


def evaluate_windows(iq,acquisition,window=64,hop=64,geometry_factor=2.):
    """One shared readout support for every filter; no cross-window lag pairs.

    Dynamic truth is averaged only within the matching measurement interval.
    Unknown moving wall geometry is not inferred from the blood field.
    """
    a=acquisition
    if window<3 or hop<1:
        raise ValueError("Invalid measurement window/hop")
    window=min(window,len(iq))
    starts=list(range(0,len(iq)-window+1,hop))
    if starts[-1]!=len(iq)-window:starts.append(len(iq)-window)
    records=[];velocities=[]
    for start in starts:
        sl=slice(start,start+window)
        truth=a.velocity_true_mps
        if truth is not None and truth.ndim==3:truth=truth[sl]
        sub=replace(a,iq=a.iq[sl],iq_views=None,timestamps_s=None,
            tx_angles_rad=None,rx_angles_rad=None,view_offsets_s=None,
            blood_views_true=None,tissue_views_true=None,
            blood_true=None if a.blood_true is None else a.blood_true[sl],
            tissue_true=None if a.tissue_true is None else a.tissue_true[sl],velocity_true_mps=truth)
        metrics=evaluate(iq[sl],sub,geometry_factor)
        records.append({"start":start,"stop":start+window,"metrics":metrics})
        velocities.append(doppler_velocity(iq[sl],a.prf_hz,a.f0_hz,a.sound_speed_mps,a.phase_sign,geometry_factor)[0])
    mask=a.flow_mask if a.flow_mask is not None else np.ones(iq.shape[1:],bool)
    stack=np.stack(velocities)
    temporal_std=float(np.mean(np.std(stack[:,mask],axis=0)))
    # No nanmean: failures must not improve an aggregate by dropping pixels.
    return {"window_frames":window,"hop_frames":hop,"records":records,
        "unattributed_velocity_temporal_std_mps":temporal_std,
        "interpretation":"temporal variability includes real pulsatility and artifacts"}


def frozen_component_metrics(result,acquisition):
    """Evaluate known components with the SAME operator fitted to their mixture.

    Independently refitting a filter to isolated blood/tissue would give a
    different data-adaptive operator and an invalid decomposition diagnostic.
    """
    a=acquisition
    unavailable={"frozen_component_metrics_available":False,
        "frozen_component_unavailable_reason":"requires a frozen operator and independent component truth"}
    if a.blood_true is None or a.tissue_true is None:
        return unavailable
    if result.apply_frozen_views is not None and a.blood_views_true is not None and a.tissue_views_true is not None:
        op=result.apply_frozen_views
        b_input,t_input,total=a.blood_views_true,a.tissue_views_true,a.iq_views
        domain="uncompounded_views_to_coherent_mean"
    elif result.apply_frozen is not None:
        op=result.apply_frozen
        b_input,t_input,total=a.blood_true,a.tissue_true,a.iq
        domain="compounded_iq"
    else:
        return unavailable
    b,t=op(b_input),op(t_input)
    n=op(total-b_input-t_input)
    power=lambda x:float(np.sum(abs(x)**2))
    ratio=lambda num,den:float(num/den) if den>0 else None
    out={"frozen_component_metrics_available":True,"frozen_component_input_domain":domain,
        "frozen_blood_power_retention":ratio(power(b),power(a.blood_true)),
        "frozen_tissue_power_retention":ratio(power(t),power(a.tissue_true)),
        "frozen_tissue_to_blood_power_ratio":ratio(power(t),power(b)),
        "frozen_noise_power_gain":ratio(power(n),power(a.iq-a.blood_true-a.tissue_true)),
        "frozen_noise_to_blood_power_ratio":ratio(power(n),power(b)),
        "frozen_component_closure_nrmse":relative_norm(b+t+n-result.iq,a.iq),
        "frozen_blood_iq_nrmse":relative_norm(b-a.blood_true,a.blood_true),
        "frozen_blood_amplitude_nrmse":relative_norm(abs(b)-abs(a.blood_true),abs(a.blood_true))}
    # A zero output has undefined leakage/output-blood and spectral shape;
    # explicit flags prevent a zero denominator from looking like success.
    out["frozen_blood_output_zero"]=power(b)<=1e-28*max(power(a.blood_true),1e-30)
    out.update(blood_fidelity(b,a.blood_true,a.flow_mask))
    if a.near_wall_mask is not None:
        out["frozen_near_wall_blood_iq_nrmse"]=relative_norm(
            (b-a.blood_true)[:,a.near_wall_mask],a.blood_true[:,a.near_wall_mask])
    return out


def relative_norm(error,reference):
    """Uncalibrated complex error. A zero reference has no relative metric."""
    den=np.linalg.norm(reference)
    return float(np.linalg.norm(error)/den) if den>0 else None


def blood_fidelity(filtered,reference,flow_mask=None):
    """Isolated blood diagnostics, with physical amplitude and phase retained.

    No fitted gain, phase, alignment, or frequency-dependent equalization.
    Normalized spectral shape is reported separately from absolute power.
    Weak support is the lowest truth-power quartile of nonzero flow pixels;
    it is evaluation-only and is never passed to a filter.
    """
    mask=np.ones(reference.shape[1:],bool) if flow_mask is None else np.asarray(flow_mask,bool)
    truth_power=np.mean(abs(reference)**2,axis=0)
    active=mask&(truth_power>0)
    if not active.any():
        return {}
    r,f=reference[:,active],filtered[:,active]
    out_power=np.mean(abs(f)**2,axis=0)
    out={"frozen_blood_spatial_power_nrmse":relative_norm(out_power-truth_power[active],truth_power[active]),
         "frozen_flow_blood_iq_nrmse":relative_norm(f-r,r)}
    p_ref=np.sum(abs(np.fft.fft(r,axis=0))**2,axis=1)
    p_out=np.sum(abs(np.fft.fft(f,axis=0))**2,axis=1)
    spectral_valid=bool(p_out.sum()>1e-28*max(p_ref.sum(),1e-30))
    out["frozen_blood_spectral_shape_valid"]=spectral_valid
    out["frozen_blood_spectral_tv"]=(float(.5*np.sum(abs(p_out/p_out.sum()-p_ref/p_ref.sum())))
        if spectral_valid else None)
    r_cov=np.mean(r[1:]*r[:-1].conj(),axis=0)
    f_cov=np.mean(f[1:]*f[:-1].conj(),axis=0)
    out["frozen_blood_lag1_covariance_nrmse"]=relative_norm(f_cov-r_cov,r_cov)
    local=np.sqrt(np.sum(abs(f-r)**2,axis=0)/np.sum(abs(r)**2,axis=0))
    out["frozen_blood_pixel_nrmse_p90"]=float(np.quantile(local,.9))
    retention=np.sum(abs(f)**2,axis=0)/np.sum(abs(r)**2,axis=0)
    out["frozen_blood_pixel_power_retention_p10"]=float(np.quantile(retention,.1))
    weak=truth_power[active]<=np.quantile(truth_power[active],.25)
    out["frozen_weak_blood_iq_nrmse"]=relative_norm((f-r)[:,weak],r[:,weak])
    return out


def injection_response(filter_fn,mixture,blood,amplitudes=(.1,1.)):
    """Refitted counterfactual blood injection, suitable for nonlinear filters.

    D_alpha=[F(T+N+alpha B)-F(T+N)]/alpha. This is a finite incremental
    response, not an additive decomposition of F(X) or a frozen linear map.
    The injected source is independent simulation truth; all fits see IQ only.
    """
    if not amplitudes or any(not np.isfinite(a) or a<=0 for a in amplitudes):
        raise ValueError("Injection amplitudes must be finite and positive")
    def fit(x):
        r=filter_fn(x)
        return (r.iq,r.diagnostics) if hasattr(r,"iq") else (r,{})
    background=mixture-blood
    baseline,background_diagnostics=fit(background)
    rows=[]
    for alpha in amplitudes:
        perturbed,diagnostics=fit(background+alpha*blood)
        incremental=(perturbed-baseline)/alpha
        rows.append({"amplitude":float(alpha),"incremental_blood_iq_nrmse":relative_norm(incremental-blood,blood),
            "incremental_power_ratio":float(np.sum(abs(incremental)**2)/max(np.sum(abs(blood)**2),1e-30)),
            "fit_diagnostics":diagnostics})
    return {"interpretation":"refitted counterfactual response; no additive source attribution",
            "background_output_power":float(np.mean(abs(baseline)**2)),
            "background_fit_diagnostics":background_diagnostics,"responses":rows}
