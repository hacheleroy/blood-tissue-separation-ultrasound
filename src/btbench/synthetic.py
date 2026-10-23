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


def advected_scatterers(seed=0,n_frames=96,shape=(10,17),prf_hz=1000.,f0_hz=5e6,
        umax_mps=.02,radius_m=.001,clutter_db=30.,snr_db=20.,
        tissue_velocity_mps=.001,tissue_acceleration_mps2=0.,
        lateral_motion_m=0.,strain_motion_m=0.,motion_frequency_hz=4.,
        psf_sigma_zx_m=(.00008,.00016),n_blood=700,n_tissue=500,phase_sign=1,uniform_flow=False):
    """Linear Gaussian-PSF scatterer transport, with separate source truth.

    Blood reflectors advect through a parabolic axial velocity field; fixed
    signed reflectivities are rendered at their moving positions. Tissue may
    translate and deform. Finite PSF support mixes vessel and wall echoes.
    A padded finite reservoir avoids periodic wrapping and phase resets.

    This is a 2-D baseband forward model, not an RF/channel or RBC simulation.
    It omits 3-D transport, attenuation, multiple scattering, aggregation and
    channel aberration. Prescribed fluid velocity is distinct from the biased
    finite-speckle, PSF-mixed Doppler readout.
    """
    shape=tuple(shape)
    if (n_frames<3 or len(shape)!=2 or min(shape)<2 or not np.isfinite([
            prf_hz,f0_hz,radius_m,umax_mps,tissue_velocity_mps,tissue_acceleration_mps2,
            lateral_motion_m,strain_motion_m,motion_frequency_hz,clutter_db,snr_db,*psf_sigma_zx_m]).all() or
            min(prf_hz,f0_hz,radius_m)<=0 or phase_sign not in (-1,1) or
            min(n_blood,n_tissue)<1 or min(psf_sigma_zx_m)<=0):
        raise ValueError("Invalid scatterer scene parameters")
    rng=np.random.default_rng(seed);c=1540.
    z=np.linspace(0,.003,shape[0]);x=np.linspace(-1.6*radius_m,1.6*radius_m,shape[1])
    zz,xx=np.meshgrid(z,x,indexing="ij")
    grid=np.stack((zz.ravel(),xx.ravel()),axis=1)
    times=np.arange(n_frames)/prf_hz
    duration=times[-1]
    pad=abs(umax_mps)*duration+abs(tissue_velocity_mps)*duration+.5*abs(tissue_acceleration_mps2)*duration**2+6*max(psf_sigma_zx_m)+abs(strain_motion_m)
    blood_pos=np.column_stack((rng.uniform(-pad,z[-1]+pad,n_blood),rng.uniform(-radius_m,radius_m,n_blood)))
    tissue_x=rng.uniform(radius_m,1.6*radius_m,n_tissue)*rng.choice([-1,1],n_tissue)
    tissue_pos=np.column_stack((rng.uniform(-pad,z[-1]+pad,n_tissue),tissue_x))
    b_amp=rng.choice([-1.,1.],n_blood);t_amp=rng.choice([-1.,1.],n_tissue)
    speeds=np.full(n_blood,umax_mps) if uniform_flow else umax_mps*(1-(blood_pos[:,1]/radius_m)**2)
    wave_number=phase_sign*4*np.pi*f0_hz/c
    grid_phase=np.exp(-1j*wave_number*grid[:,0])
    sigma=np.asarray(psf_sigma_zx_m)
    def render(pos,amp):
        d=(grid[:,None,:]-pos[None,:,:])/sigma
        envelope=np.exp(-.5*np.sum(d*d,axis=2))
        return ((envelope@(amp*np.exp(1j*wave_number*pos[:,0])))*grid_phase).reshape(shape)
    b=np.empty((n_frames,)+shape,complex);tissue=np.empty_like(b)
    for i,time in enumerate(times):
        bp=blood_pos.copy();bp[:,0]+=speeds*time
        tp=tissue_pos.copy()
        periodic=np.sin(2*np.pi*motion_frequency_hz*time)
        tp[:,0]+=tissue_velocity_mps*time+.5*tissue_acceleration_mps2*time**2
        tp[:,0]+=strain_motion_m*(tissue_pos[:,0]-z[-1]/2)/z[-1]*periodic
        tp[:,1]+=lateral_motion_m*periodic
        b[i]=render(bp,b_amp);tissue[i]=render(tp,t_amp)
    flow=np.broadcast_to(abs(x)[None,:]<radius_m,shape)
    # Fixed forward-model source strengths: no per-frame normalization.
    b/=max(np.sqrt(np.mean(abs(b[0,flow])**2)),1e-30)
    tissue*=10**(clutter_db/20)/max(np.sqrt(np.mean(abs(tissue[0,~flow])**2)),1e-30)
    noise=(rng.normal(size=b.shape)+1j*rng.normal(size=b.shape))/np.sqrt(2)*10**(-snr_db/20)
    profile=np.ones_like(x) if uniform_flow else np.maximum(1-(x/radius_m)**2,0)
    fluid_velocity=np.broadcast_to(umax_mps*profile,shape).copy()
    return Acquisition(b+tissue+noise,prf_hz,f0_hz,phase_sign=phase_sign,z_m=z,x_m=x,
        flow_mask=flow,tissue_roi=~flow,near_wall_mask=flow&(abs(xx)>.65*radius_m),
        blood_true=b,tissue_true=tissue,velocity_true_mps=fluid_velocity,
        metadata={"generator":"linear_PSF_advected_scatterers","seed":seed,"contrast_agent":False,
            "not_acoustic_ground_truth":True,"psf_sigma_zx_m":list(psf_sigma_zx_m),
            "n_blood":n_blood,"n_tissue":n_tissue,"source_strengths":"fixed from first frame",
            "flow_model":"finite padded reservoir; parabolic axial transport",
            "lateral_motion_m":lateral_motion_m,"strain_motion_m":strain_motion_m,
            "clutter_db":clutter_db,"snr_db":snr_db,"uniform_flow":uniform_flow,
            "fluid_structure_coupling":False})
