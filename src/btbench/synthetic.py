"""Reproducible non-microbubble IQ surrogates, not acoustic/RBC simulations."""
import numpy as np
from scipy import ndimage
from .core import Acquisition


def poiseuille(seed=0,n_frames=256,shape=(12,16),prf_hz=1000.,f0_hz=5e6,
               umax_mps=0.025,radius_m=1e-3,clutter_db=30.,snr_db=20.,
               tissue_velocity_mps=0.001,tissue_acceleration_mps2=0.005,
               radial_amplitude_exponent=0.,decorrelation=0.98,motion_event=False):
    """Vertical cylindrical Poiseuille central slice, axial flow and x-normal walls.

    Dense correlated complex Gaussian speckle + local Doppler phase; tissue
    has a coherent texture and accelerating phase. RBC aggregation, cell-free
    layer, acoustic scattering, PSF transport, and microstreaming are NOT modeled.
    radial_amplitude_exponent explores amplitude gradients without bubbles.
    Spatial correlation is imposed on innovations; phase varies with x.
    """
    if not 0 <= decorrelation < 1:
        raise ValueError("decorrelation must be in [0,1)")
    shape=tuple(shape)
    rng=np.random.default_rng(seed)
    z=np.linspace(0,3e-3,shape[0])
    x=np.linspace(-1.5*radius_m,1.5*radius_m,shape[1])
    rr=np.broadcast_to(abs(x)[None,:],shape)
    flow=rr<radius_m
    velocity=umax_mps*np.maximum(1-(rr/radius_m)**2,0)
    dt=1/prf_hz
    phase_step=4*np.pi*f0_hz/1540*velocity*dt
    amp=np.maximum(1-(rr/radius_m)**2,0)**(radial_amplitude_exponent/2)*flow
    def random_field():
        a=ndimage.gaussian_filter(rng.normal(size=shape)+1j*rng.normal(size=shape),0.7)
        return a/np.sqrt(np.mean(abs(a)**2))
    b=np.empty((n_frames,)+shape,complex)
    state=random_field()
    for t in range(n_frames):
        state=decorrelation*state*np.exp(1j*phase_step)+np.sqrt(1-decorrelation**2)*random_field()
        b[t]=state*amp
    texture=random_field()*(0.15+0.85*(~flow))
    times=np.arange(n_frames)*dt
    disp=tissue_velocity_mps*times+0.5*tissue_acceleration_mps2*times**2
    phi=4*np.pi*f0_hz/1540*disp
    tissue=10**(clutter_db/20)*texture[None]*np.exp(1j*phi)[:,None,None]
    if motion_event:
        # A sudden texture displacement, explicitly a registration stress test.
        i=n_frames//2
        moved=ndimage.shift(texture,(1,1),order=1,mode="nearest")
        tissue[i:]=10**(clutter_db/20)*moved[None]*np.exp(1j*phi[i:,None,None])
    noise=(rng.normal(size=b.shape)+1j*rng.normal(size=b.shape))/np.sqrt(2)*10**(-snr_db/20)
    wall_z=z[1:-1:2]
    positions=np.array([(zz,xx) for zz in wall_z for xx in (-radius_m,radius_m)])
    normals=np.array([(0,nn) for zz in wall_z for nn in (1,-1)])
    wall={"positions_zx_m":positions.tolist(),"inward_normals_zx":normals.tolist(),
          "distances_m":np.linspace(radius_m*0.08,radius_m*0.45,5).tolist(),
          "viscosity_pa_s":0.0035,"velocity_mps":0.,
          "true_pa":np.full(len(positions),2*0.0035*umax_mps/radius_m).tolist()}
    return Acquisition(b+tissue+noise,prf_hz,f0_hz,z_m=z,x_m=x,
        tissue_roi=~flow,flow_mask=flow,near_wall_mask=flow&(rr>0.65*radius_m),
        blood_true=b,tissue_true=tissue,velocity_true_mps=velocity,wall=wall,
        metadata={"generator":"statistical_IQ_surrogate","seed":seed,
        "not_acoustic_ground_truth":True,"contrast_agent":False,
        "clutter_db":clutter_db,"snr_db":snr_db})
