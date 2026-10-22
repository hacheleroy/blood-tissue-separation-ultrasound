"""Axis 3 prototype: physical space/angle/slow-time steering and local nulling.

Requires uncompounded complex views and supplied propagation geometry.
No receive-channel information is invented from compounded IQ.
"""
import numpy as np
from scipy import linalg
from .core import FilterResult


def steering(velocities_zx_mps, tx_rad, rx_rad, frame_times_s, f0_hz,
             sound_speed_mps=1540., phase_sign=1, view_offsets_s=None, gains=None):
    v = np.atleast_2d(velocities_zx_mps)
    tx,rx = np.asarray(tx_rad),np.asarray(rx_rad)
    if tx.shape != rx.shape or v.shape[1] != 2 or phase_sign not in (-1,1):
        raise ValueError("Steering geometry mismatch")
    geom = np.stack([np.cos(tx)+np.cos(rx),np.sin(tx)+np.sin(rx)],axis=1)
    offsets = np.zeros(len(tx)) if view_offsets_s is None else np.asarray(view_offsets_s)
    gain = np.ones(len(tx)) if gains is None else np.asarray(gains)
    times = np.asarray(frame_times_s)[:,None]+offsets[None,:]
    phase = phase_sign*2*np.pi*f0_hz/sound_speed_mps*times[:,:,None]*(geom@v.T)[None,:,:]
    return (gain[None,:,None]*np.exp(1j*phase)).reshape(-1,len(v))


def null_projector(st, covariance=None, loading_rel=1e-3):
    """P ST=0. With covariance this is a loaded oblique clutter projector."""
    st = np.atleast_2d(st)
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
    s = np.column_stack((np.ravel(sb),st))
    if np.linalg.matrix_rank(s)<s.shape[1]:
        raise ValueError("Blood and tissue steering constraints overlap; distortionless null is unidentifiable")
    c = (covariance+covariance.conj().T)/2
    scale = max(np.trace(c).real/len(c),1e-30)
    ris = linalg.solve(c+loading_rel*scale*np.eye(len(c)),s,assume_a="pos")
    f = np.zeros(s.shape[1]);f[0]=1
    return ris@np.linalg.pinv(s.conj().T@ris)@f


def local_stap(iq_views,prf_hz,f0_hz,tx_rad,rx_rad,tissue_velocities_zx_mps,
               sound_speed_mps=1540.,phase_sign=1,view_offsets_s=None,
               taps=4,guard=1,training_radius=3,training_roi=None,adaptive=False):
    """Short-snapshot overlap/add with guard-excluded spatial training.

    roi must be independently defined tissue-only support for adaptive mode.
    This is not a calibrated full-array STAP detector. Unknown array gains,
    near-field steering, and deformation require an acquisition-specific model.
    """
    if iq_views.ndim != 4 or len(iq_views)<taps or taps<2:
        raise ValueError("Need (time,view,z,x) data with >=taps frames")
    t,a,z,x = iq_views.shape
    if adaptive and (training_roi is None or np.shape(training_roi)!=(z,x)):
        raise ValueError("Adaptive STAP requires tissue-only training ROI")
    st = steering(tissue_velocities_zx_mps,tx_rad,rx_rad,np.arange(taps)/prf_hz,
                  f0_hz,sound_speed_mps,phase_sign,view_offsets_s)
    y = np.zeros_like(iq_views,dtype=complex)
    mass = np.zeros(t)
    window = np.hanning(taps+2)[1:-1]
    insufficient=0
    for start in range(t-taps+1):
        block = iq_views[start:start+taps].reshape(taps*a,z,x)
        out = np.empty_like(block)
        base_p = null_projector(st)
        for iz in range(z):
            for ix in range(x):
                p = base_p
                if adaptive:
                    zz,xx=np.ogrid[:z,:x]
                    dist=np.maximum(abs(zz-iz),abs(xx-ix))
                    train=(dist>guard)&(dist<=training_radius)&training_roi
                    if train.sum()>=st.shape[1]+2:
                        samples=block[:,train]
                        cov=samples@samples.conj().T/samples.shape[1]
                        p=null_projector(st,cov)
                    else:
                        insufficient+=1
                out[:,iz,ix]=p@block[:,iz,ix]
        y[start:start+taps] += out.reshape(taps,a,z,x)*window[:,None,None,None]
        mass[start:start+taps] += window
    return FilterResult(np.mean(y/mass[:,None,None,None],axis=1),
        {"taps":taps,"adaptive":adaptive,"training_fallbacks":insufficient,
         "input_views":a,"model":"local far-field angle-time prototype"})
