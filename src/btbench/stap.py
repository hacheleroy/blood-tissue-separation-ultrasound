"""Axis 3 prototype: physical space/angle/slow-time steering and local nulling.

Requires uncompounded complex views and supplied propagation geometry.
No receive-channel information is invented from compounded IQ.
"""
import numpy as np
from dataclasses import dataclass,field
from scipy import linalg
from .core import FilterResult


def steering(velocities_zx_mps, tx_rad, rx_rad, frame_times_s, f0_hz,
             sound_speed_mps=1540., phase_sign=1, view_offsets_s=None, gains=None):
    v = np.atleast_2d(velocities_zx_mps).astype(float)
    tx,rx = np.asarray(tx_rad,float),np.asarray(rx_rad,float)
    times0=np.asarray(frame_times_s,float)
    if (tx.ndim!=1 or not len(tx) or tx.shape != rx.shape or v.shape[1] != 2
            or phase_sign not in (-1,1) or times0.ndim!=1 or not len(times0)
            or not all(np.isfinite(a).all() for a in (v,tx,rx,times0))
            or min(f0_hz,sound_speed_mps)<=0):
        raise ValueError("Steering geometry mismatch")
    geom = np.stack([np.cos(tx)+np.cos(rx),np.sin(tx)+np.sin(rx)],axis=1)
    offsets = np.zeros(len(tx)) if view_offsets_s is None else np.asarray(view_offsets_s)
    gain = np.ones(len(tx)) if gains is None else np.asarray(gains)
    if (offsets.shape!=tx.shape or gain.shape!=tx.shape or not np.isfinite(offsets).all()
            or not np.isfinite(gain).all() or np.any(abs(gain)==0)):
        raise ValueError("View timing/gain shape or values invalid")
    times = times0[:,None]+offsets[None,:]
    phase = phase_sign*2*np.pi*f0_hz/sound_speed_mps*times[:,:,None]*(geom@v.T)[None,:,:]
    return (gain[None,:,None]*np.exp(1j*phase)).reshape(-1,len(v))


def null_projector(st, covariance=None, loading_rel=1e-3):
    """P ST=0. With covariance this is a loaded oblique clutter projector."""
    st = np.asarray(st,complex)
    if st.ndim==1:st=st[:,None]
    if st.ndim!=2 or not len(st) or not np.isfinite(st).all():
        raise ValueError("Tissue steering must be a finite (snapshot,component) matrix")
    if covariance is None:
        return np.eye(len(st))-st@np.linalg.pinv(st)
    if loading_rel <= 0:
        raise ValueError("Adaptive covariance requires positive diagonal loading")
    cov = (covariance+covariance.conj().T)/2
    scale = max(np.trace(cov).real/len(cov),1e-30)
    ri_st = linalg.solve(cov+loading_rel*scale*np.eye(len(cov)),st,assume_a="pos")
    return np.eye(len(st))-st@np.linalg.pinv(st.conj().T@ri_st)@ri_st.conj().T


def lcmv_weights(sb,st,covariance,loading_rel=1e-3):
    """Unit response to supplied blood steering, zero to supplied tissue steering."""
    st=np.asarray(st)
    if st.ndim==1:st=st[:,None]
    s = np.column_stack((np.ravel(sb),st))
    if loading_rel<=0:raise ValueError("LCMV requires positive diagonal loading")
    if np.linalg.matrix_rank(s)<s.shape[1]:
        raise ValueError("Blood and tissue steering constraints overlap; distortionless null is unidentifiable")
    c = (covariance+covariance.conj().T)/2
    scale = max(np.trace(c).real/len(c),1e-30)
    ris = linalg.solve(c+loading_rel*scale*np.eye(len(c)),s,assume_a="pos")
    f = np.zeros(s.shape[1]);f[0]=1
    gram=s.conj().T@ris
    gram=(gram+gram.conj().T)/2
    if np.linalg.cond(gram)>1e12:
        raise ValueError("Loaded LCMV constraints are numerically unidentifiable")
    return ris@linalg.solve(gram,f,assume_a="pos")


@dataclass
class ViewFilterResult:
    views: np.ndarray
    diagnostics: dict=field(default_factory=dict)
    apply_frozen: object=field(default=None,repr=False)


def local_stap_views(iq_views,prf_hz,f0_hz,tx_rad,rx_rad,tissue_velocities_zx_mps,
               sound_speed_mps=1540.,phase_sign=1,view_offsets_s=None,gains=None,
               taps=4,guard=1,training_radius=3,training_roi=None,adaptive=False,
               loading_rel=1e-2,hop=1,start_time_s=0.,min_training=None,
               max_operator_bytes=512*1024**2):
    """Short-snapshot overlap/add with guard-excluded spatial training.

    roi must be independently defined tissue-only support for adaptive mode.
    This is not a calibrated full-array STAP detector. Unknown array gains,
    near-field steering, and deformation require an acquisition-specific model.
    """
    iq_views=np.asarray(iq_views)
    if (iq_views.ndim != 4 or len(iq_views)<taps or taps<2 or not np.iscomplexobj(iq_views)
            or not np.isfinite(iq_views).all() or prf_hz<=0):
        raise ValueError("Need (time,view,z,x) data with >=taps frames")
    t,a,z,x = iq_views.shape
    if (len(tx_rad)!=a or not 1<=hop<=taps or guard<0 or training_radius<=guard
            or loading_rel<=0):
        raise ValueError("Invalid view count, hop, guard, radius or loading")
    if adaptive and (training_roi is None or np.shape(training_roi)!=(z,x)):
        raise ValueError("Adaptive STAP requires tissue-only training ROI")
    starts=list(range(0,t-taps+1,hop))
    if starts[-1]!=t-taps:starts.append(t-taps)
    d=taps*a
    operator_bytes=len(starts)*(z*x if adaptive else 1)*d*d*16
    if operator_bytes>max_operator_bytes:
        raise MemoryError("Frozen STAP operators exceed memory cap; use smaller packets/patches")
    y = np.zeros_like(iq_views,dtype=complex)
    mass = np.zeros(t)
    window = np.exp(-2*np.linspace(-1,1,taps)**2)
    insufficient=0;sample_counts=[];low_support=0;operators=[];ranks=[]
    zz,xx=np.indices((z,x))
    masks=[]
    if adaptive:
        for iz in range(z):
            for ix in range(x):
                dist=np.maximum(abs(zz-iz),abs(xx-ix))
                masks.append(((dist>guard)&(dist<=training_radius)&np.asarray(training_roi,bool)).ravel())
    for start in starts:
        # Angular phase evolves with absolute snapshot time. Resetting every
        # view's time to zero is incorrect for nonzero transverse velocity.
        st = steering(tissue_velocities_zx_mps,tx_rad,rx_rad,
            start_time_s+(np.arange(taps)+start)/prf_hz,f0_hz,sound_speed_mps,
            phase_sign,view_offsets_s,gains)
        ranks.append(int(np.linalg.matrix_rank(st)))
        required=max(4,ranks[-1]+2) if min_training is None else int(min_training)
        block = iq_views[start:start+taps].reshape(d,z*x)
        out = np.empty_like(block)
        base_p = null_projector(st)
        op=np.repeat(base_p[None],z*x,axis=0) if adaptive else base_p
        for iz in range(z):
            for ix in range(x):
                p = base_p
                if adaptive:
                    train=masks[iz*x+ix]
                    sample_counts.append(int(train.sum()))
                    if train.sum()<2*d:low_support+=1
                    if train.sum()>=required:
                        samples=block[:,train]
                        cov=samples@samples.conj().T/samples.shape[1]
                        p=null_projector(st,cov,loading_rel)
                    else:
                        insufficient+=1
                    op[iz*x+ix]=p
                out[:,iz*x+ix]=p@block[:,iz*x+ix]
        operators.append((start,op))
        y[start:start+taps] += out.reshape(taps,a,z,x)*window[:,None,None,None]
        mass[start:start+taps] += window
    def apply(other):
        if other.shape!=iq_views.shape:raise ValueError("Frozen STAP view shape mismatch")
        out=np.zeros_like(other,dtype=complex)
        for start,op in operators:
            block=other[start:start+taps].reshape(d,z*x)
            b=np.einsum("pij,jp->ip",op,block) if adaptive else op@block
            out[start:start+taps]+=b.reshape(taps,a,z,x)*window[:,None,None,None]
        return out/mass[:,None,None,None]
    return ViewFilterResult(y/mass[:,None,None,None],
        {"taps":taps,"adaptive":adaptive,"training_fallbacks":insufficient,
         "input_views":a,"snapshot_dimension":d,"hop":hop,"snapshots":len(starts),
         "loading_rel":loading_rel,"training_samples_min":min(sample_counts) if sample_counts else None,
         "training_samples_max":max(sample_counts) if sample_counts else None,
         "fits_below_twice_dimension":low_support,"tissue_rank_max":max(ranks),
         "operator_bytes":operator_bytes,"absolute_snapshot_times":True,
         "temporal_window":"Gaussian exp(-2 u^2), u in [-1,1]",
         "coherent_views_preserved":True,"model":"local far-field angle-time prototype"},apply)


def local_stap(iq_views,*args,**kwargs):
    """Legacy compounded-IQ API. View-level truth diagnostics use local_stap_views."""
    r=local_stap_views(iq_views,*args,**kwargs)
    return FilterResult(r.views.mean(axis=1),r.diagnostics,
        apply_frozen_views=lambda other:r.apply_frozen(other).mean(axis=1))


def estimate_tissue_velocity(calibration_views,roi,tx_rad,rx_rad,prf_hz,f0_hz,
                             sound_speed_mps=1540.,phase_sign=1):
    """Vector Doppler from an independent calibration packet and annotated ROI.

    No test blood/tissue component arrays or prescribed test velocity enter.
    """
    a=np.asarray(calibration_views)
    roi=np.asarray(roi,bool)
    if a.ndim!=4 or roi.shape!=a.shape[2:] or not roi.any():
        raise ValueError("Need calibration views and a nonempty matching tissue ROI")
    cross=np.mean(a[1:,:,roi]*a[:-1,:,roi].conj(),axis=(0,2))
    tx,rx=np.asarray(tx_rad),np.asarray(rx_rad)
    g=np.stack([np.cos(tx)+np.cos(rx),np.sin(tx)+np.sin(rx)],axis=1)
    if np.linalg.matrix_rank(g)<2 or np.linalg.cond(g)>100:
        raise ValueError("Calibration geometry cannot identify 2-D tissue velocity")
    freq=phase_sign*np.angle(cross)*sound_speed_mps*prf_hz/(2*np.pi*f0_hz)
    return np.linalg.lstsq(g,freq,rcond=None)[0]


def lcmv_velocity_spectrum(views,candidates_zx_mps,tissue_velocities_zx_mps,
        tx_rad,rx_rad,prf_hz,f0_hz,training_roi,target_zx,taps=8,hop=8,
        guard=1,loading_rel=1e-2,view_offsets_s=None,gains=None,
        sound_speed_mps=1540.,phase_sign=1,min_steering_sine=1e-3):
    """Single-target velocity-bank detection; this is not a blood-IQ reconstruction.

    Constraints and scores use only data, geometry and a fixed candidate grid.
    Normalize each output power by its loaded training-covariance response.
    Unresolvable blood/tissue constraints are rejected and counted explicitly.
    """
    t,a,z,x=views.shape
    iz,ix=target_zx
    if not (0<=iz<z and 0<=ix<x) or not 2<=taps<=t or not 1<=hop<=taps:
        raise ValueError("Invalid target or snapshot support")
    zz,xx=np.indices((z,x))
    train=np.asarray(training_roi,bool)&(np.maximum(abs(zz-iz),abs(xx-ix))>guard)
    if train.sum()<4:raise ValueError("Too few guard-excluded training pixels")
    candidates=np.asarray(candidates_zx_mps,float)
    scores=np.zeros(len(candidates));raw=np.zeros(len(candidates));counts=np.zeros(len(candidates),int)
    rejected=0;max_response_error=0.;max_noise_gain=0.
    for start in range(0,t-taps+1,hop):
        time=(np.arange(taps)+start)/prf_hz
        kw=dict(tx_rad=tx_rad,rx_rad=rx_rad,frame_times_s=time,f0_hz=f0_hz,
            sound_speed_mps=sound_speed_mps,phase_sign=phase_sign,view_offsets_s=view_offsets_s,gains=gains)
        st=steering(tissue_velocities_zx_mps,**kw)
        sb=steering(candidates,**kw)
        block=views[start:start+taps].reshape(taps*a,z,x)
        samples=block[:,train]
        c=samples@samples.conj().T/samples.shape[1]
        loaded=c+loading_rel*max(np.trace(c).real/len(c),1e-30)*np.eye(len(c))
        orth=null_projector(st)
        for j in range(len(candidates)):
            sine=np.linalg.norm(orth@sb[:,j])/np.linalg.norm(sb[:,j])
            if sine<min_steering_sine:
                rejected+=1;continue
            try:w=lcmv_weights(sb[:,j],st,c,loading_rel)
            except (ValueError,np.linalg.LinAlgError):
                rejected+=1;continue
            response_error=max(float(abs(w.conj()@sb[:,j]-1)),float(np.max(abs(w.conj()@st))))
            if response_error>1e-7:
                rejected+=1;continue
            power=float(abs(w.conj()@block[:,iz,ix])**2)
            raw[j]+=power
            scores[j]+=power/max(float((w.conj()@loaded@w).real),1e-30)
            counts[j]+=1
            max_response_error=max(max_response_error,response_error)
            max_noise_gain=max(max_noise_gain,float(np.vdot(w,w).real))
    score=np.divide(scores,counts,out=np.full(len(scores),-np.inf),where=counts>0)
    best=int(np.argmax(score)) if counts.any() else None
    return {"candidates_zx_mps":candidates.tolist(),
        "normalized_power":np.where(np.isfinite(score),score,np.nan).tolist(),
        "raw_power":np.divide(raw,counts,out=np.full(len(raw),np.nan),where=counts>0).tolist(),
        "selected_velocity_zx_mps":None if best is None else candidates[best].tolist(),
        "rejected_snapshot_constraints":rejected,"constraint_response_error_max":max_response_error,
        "white_noise_gain_max":max_noise_gain,"training_pixels":int(train.sum()),
        "output":"velocity-grid detection at one target; no WSS or reconstructed IQ",
        "grid_selection_uses_test_truth":False}
