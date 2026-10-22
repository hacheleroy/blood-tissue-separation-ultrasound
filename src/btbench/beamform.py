"""Bounded CPU reference plane-wave DAS, preserving complex slow-time phase.

Acquisition-specific time references must be supplied explicitly. This module
does not guess multiplexing, t0 conventions, or turn a B-mode recipe into IQ.
Compare its B-mode/point-target result with the provider's reconstruction before
quantitative use. Production RF reconstruction should use the provider/zea stack.
"""
import numpy as np
from scipy.signal import hilbert


def plane_wave_das(rf,element_zx_m,z_m,x_m,tx_angles_rad,fs_hz,f0_hz,
                   initial_times_s,tx_reference_times_s,tx_origins_zx_m=None,
                   sound_speed_mps=1540.,f_number=1.5):
    """rf: real (T,A,S,E); return unfiltered (T,A,Z,X) complex IQ.

    t_tx(r)=t_ref[a]+(r-origin[a]).n_tx[a]/c;
    sample index=(t_tx+||r-r_e||/c-initial_time[a])*fs.
    No compounding across unknown angle timestamps. No implicit TGC reversal.
    """
    rf=np.asarray(rf)
    if rf.ndim!=4 or np.iscomplexobj(rf):
        raise ValueError("RF must be real (time,transmit,sample,element)")
    t,a,s,e=rf.shape
    el=np.asarray(element_zx_m)
    theta=np.asarray(tx_angles_rad)
    ti=np.broadcast_to(initial_times_s,(a,))
    tr=np.broadcast_to(tx_reference_times_s,(a,))
    origin=np.zeros((a,2)) if tx_origins_zx_m is None else np.asarray(tx_origins_zx_m)
    if el.shape!=(e,2) or theta.shape!=(a,) or origin.shape!=(a,2):
        raise ValueError("DAS geometry mismatch")
    zz,xx=np.meshgrid(z_m,x_m,indexing="ij")
    pos=np.stack((zz.ravel(),xx.ravel()),axis=1)
    rx=np.linalg.norm(pos[:,None,:]-el[None,:,:],axis=2)/sound_speed_mps
    aperture=abs(pos[:,None,1]-el[None,:,1])<=np.maximum(pos[:,None,0],1e-12)/(2*f_number)
    # Smooth finite aperture, with a common amplitude normalization across frames.
    aperture=aperture.astype(float)
    mass=np.maximum(aperture.sum(1),1)
    out=np.empty((t,a,len(z_m),len(x_m)),complex)
    analytic=hilbert(rf,axis=2)
    for j,angle in enumerate(theta):
        direction=np.array([np.cos(angle),np.sin(angle)])
        tx=(pos-origin[j])@direction/sound_speed_mps+tr[j]
        idx=(tx[:,None]+rx-ti[j])*fs_hz
        i=np.floor(idx).astype(int)
        f=idx-i
        valid=(i>=0)&(i<s-1)
        i=np.clip(i,0,s-2)
        chan=np.arange(e)[None,:]
        for frame in range(t):
            vals=(1-f)*analytic[frame,j,i,chan]+f*analytic[frame,j,i+1,chan]
            iq=np.sum(vals*aperture*valid,axis=1)/mass
            # Fixed, spatial demodulation has no effect on slow-time phase slope.
            iq*=np.exp(-1j*2*np.pi*f0_hz*(2*pos[:,0]/sound_speed_mps))
            out[frame,j]=iq.reshape(len(z_m),len(x_m))
    return out
