"""Shared-view, bounded-history operators for the three research axes.

PublicInput contains observations and calibration only. Source arrays and known
warps belong to the evaluator. One small spatial tile is shared by every arm;
this is an IQ software experiment, not an acoustic or experimental validation.
"""
from dataclasses import dataclass, field
import hashlib
import json
import numpy as np
from scipy import linalg, ndimage, optimize
from .axes import joint_diagonalize, component_weights
from .stap import estimate_tissue_velocity, null_projector
from .svd_literature import baranger_band
from .registration import MotionGate, register_displacement


AXIS_NAMES = {
    'axis1': 'Axis 1: Motion-compensated virtual long ensemble',
    'axis2': 'Axis 2: Motion-demodulated multi-lag separation',
    'axis3': 'Axis 3: Local ultrasound STAP',
}


@dataclass(frozen=True)
class PublicInput:
    views: np.ndarray
    calibration_views: np.ndarray
    training_roi: np.ndarray
    geometry: dict

    def __post_init__(self):
        x, c = self.views, self.calibration_views
        if (x.ndim != 4 or c.ndim != 4 or x.shape[1:] != c.shape[1:]
                or min(len(x), len(c)) < 3 or not np.iscomplexobj(x)
                or not np.iscomplexobj(c) or not np.isfinite(x).all()
                or not np.isfinite(c).all()):
            raise ValueError('Need finite complex (time,view,z,x) observations and calibration')
        if self.training_roi.shape != x.shape[2:] or not self.training_roi.any():
            raise ValueError('Need a nonempty, shared annotation ROI')
        g = self.geometry
        if min(g['prf_hz'], g['f0_hz'], g['sound_speed_mps']) <= 0 or g['phase_sign'] not in (-1, 1):
            raise ValueError('Invalid acquisition constants')
        for key in ('tx_rad', 'rx_rad', 'view_offsets_s'):
            if np.shape(g[key]) != (x.shape[1],) or not np.isfinite(g[key]).all():
                raise ValueError('One geometry entry per view is required')


@dataclass
class SharedScene:
    public: PublicInput
    blood: np.ndarray
    tissue: np.ndarray
    noise: np.ndarray
    known_displacement: np.ndarray
    metadata: dict


def array_hash(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(str(a.shape).encode() + a.dtype.str.encode() + a.tobytes()).hexdigest()


def motion_basis(shape):
    z, x = np.meshgrid(np.linspace(-1, 1, shape[0]), np.linspace(-1, 1, shape[1]), indexing='ij')
    return np.stack([np.ones(shape), z, x, np.sin(np.pi*z)*np.sin(np.pi*x)])


def inverse_coordinates(displacement, iterations=40):
    """Invert reference->lab coordinates by fixed-point iteration; no wrapping."""
    grid = np.indices(displacement.shape[1:], dtype=float)
    q = grid.copy()
    for _ in range(iterations):
        q = grid - np.stack([ndimage.map_coordinates(d, q, order=1, mode='nearest') for d in displacement])
    return q


def sample_frame(frame, coordinates, zero_extension=False):
    out = np.empty_like(frame, dtype=complex)
    mode='grid-constant' if zero_extension else 'constant'
    for j in range(len(frame)):
        out[j] = (ndimage.map_coordinates(frame[j].real, coordinates, order=1, mode=mode, cval=0.)
                  + 1j*ndimage.map_coordinates(frame[j].imag, coordinates, order=1, mode=mode, cval=0.))
    return out


def warp_series(views, displacement, inverse=False):
    grid = np.indices(views.shape[2:], dtype=float)
    out = np.empty_like(views, dtype=complex)
    for t in range(len(views)):
        coordinates = inverse_coordinates(displacement[t]) if inverse else grid + displacement[t]
        out[t] = sample_frame(views[t], coordinates)
    return out


def shared_scene(seed=0, n_frames=128, calibration_frames=64, shape=(8, 16),
                 prf_hz=1000., f0_hz=5e6, flow_mps=.006, clutter_db=30.,
                 snr_db=20., decorrelation=.985, motion='none', amplitude_pixels=.8,
                 phase_sign=1, acceleration_mps2=0., view_phase_drift_rad=0.,
                 low_rank=False, identical_statistics=False):
    """Specified angular-IQ fields followed by an imposed spatial coordinate map.

    Calibration and evaluation use disjoint times and independent blood/noise
    draws, sharing a fixed tissue texture. Calibration has no imposed spatial
    deformation. Envelope transport and Doppler are specified separately; this
    is deliberately not an RF/particle-scattering model. Resampling and finite
    boundaries remain part of the scored output.
    """
    if min(shape) < 4 or n_frames < 16 or calibration_frames < 8 or not 0 <= decorrelation <= 1:
        raise ValueError('Invalid shared scene dimensions or correlation')
    if motion not in ('none', 'translation', 'smooth', 'unmodeled'):
        raise ValueError('Unknown imposed motion')
    rng = np.random.default_rng(seed)
    z, x = np.meshgrid(np.linspace(-1, 1, shape[0]), np.linspace(-.0015, .0015, shape[1]), indexing='ij')
    flow = abs(x) < .00085
    roi = abs(x) > .00105
    profile = np.maximum(1 - (x/.00085)**2, 0.) * flow
    tx = np.deg2rad([-18., 0., 18.]); rx = np.zeros(3)
    offsets = np.arange(3)/(3*prf_hz)
    geom = np.stack([np.cos(tx)+np.cos(rx), np.sin(tx)+np.sin(rx)], axis=1)
    velocity = np.array([.0015, .0005])
    gains = np.array([1., .92*np.exp(.1j), 1.08*np.exp(-.06j)])
    scale = phase_sign*2*np.pi*f0_hz/1540.
    def texture():
        a = ndimage.gaussian_filter(rng.normal(size=shape)+1j*rng.normal(size=shape), .65)
        return a/max(np.sqrt(np.mean(abs(a)**2)), 1e-30)
    tissue_texture = texture()*(.2+.8*(~flow))
    blood_texture = texture()*flow
    def fields(times, test=False):
        tj = times[:, None]+offsets[None]
        relative = np.zeros(shape) if identical_statistics else flow_mps*(flow if low_rank else profile)
        doppler = np.einsum('vd,zxd->vzx', geom,
                            np.broadcast_to(velocity, (*shape, 2)) + np.stack([relative, np.zeros(shape)], -1))
        state = texture()*flow; scalar = 1.+0j; env = []
        for _ in times:
            if identical_statistics:
                env.append(blood_texture)
            elif low_rank:
                scalar = decorrelation*scalar + np.sqrt(1-decorrelation**2)*(rng.normal()+1j*rng.normal())/np.sqrt(2)
                env.append(blood_texture*scalar)
            else:
                state = decorrelation*state + np.sqrt(1-decorrelation**2)*texture()*flow
                env.append(state.copy())
        phase_t = scale*(tj[:,:,None,None]*(geom@velocity)[None,:,None,None]
                         + (.5*acceleration_mps2 if test else 0.)*tj[:,:,None,None]**2*geom[None,:,0,None,None])
        phase_b = scale*tj[:,:,None,None]*doppler[None]
        gain = gains*np.exp(1j*np.array([0., view_phase_drift_rad, -view_phase_drift_rad])) if test else gains
        b = np.asarray(env)[:,None]*np.exp(1j*phase_b)*gain[None,:,None,None]
        t = 10**(clutter_db/20)*tissue_texture[None,None]*np.exp(1j*phase_t)*gain[None,:,None,None]
        return b, t
    times = np.arange(n_frames)/prf_hz
    b, t = fields(times, True)
    cb, ct = fields(np.arange(-calibration_frames, 0)/prf_hz)
    basis = motion_basis(shape); displacement = np.zeros((n_frames, 2, *shape))
    phase = 2*np.pi*np.arange(n_frames)/n_frames
    if motion != 'none':
        displacement[:,0] = amplitude_pixels*np.sin(phase)[:,None,None]
        displacement[:,1] = .7*amplitude_pixels*np.sin(2*phase)[:,None,None]
    if motion in ('smooth', 'unmodeled'):
        displacement[:,0] += .35*amplitude_pixels*np.sin(phase)[:,None,None]*basis[3]
        displacement[:,1] += .25*amplitude_pixels*np.sin(phase)[:,None,None]*basis[1]
    if motion == 'unmodeled':
        extra = np.cos(2*np.pi*z)*np.cos(2*np.pi*x/.0015)
        displacement[:,0] += .25*amplitude_pixels*np.sin(phase)[:,None,None]*extra
    b = warp_series(b, displacement, inverse=True)
    t = warp_series(t, displacement, inverse=True)
    def noise(a):
        return (rng.normal(size=a.shape)+1j*rng.normal(size=a.shape))/np.sqrt(2)*10**(-snr_db/20)
    n = noise(b); calibration = cb+ct+noise(cb)
    g = {'prf_hz':prf_hz, 'f0_hz':f0_hz, 'sound_speed_mps':1540., 'phase_sign':phase_sign,
         'tx_rad':tx.tolist(), 'rx_rad':rx.tolist(), 'view_offsets_s':offsets.tolist()}
    public = PublicInput(b+t+n, calibration, roi, g)
    return SharedScene(public, b, t, n, displacement,
        {'generator':'specified angular IQ with imposed image-coordinate motion; not acoustic truth',
         'seed':seed, 'motion':motion, 'shared_tissue_texture_in_disjoint_calibration':True,
         'calibration_time_stop_exclusive_s':0., 'evaluation_start_s':0.,
         'known_motion_is_evaluation_only':True, 'contrast_agent':False})


def calibrate(public, return_complex_reference=False):
    """One common, data-only phase/gain calibration for every arm."""
    g = public.geometry; cal = public.calibration_views; roi = public.training_roi
    v = estimate_tissue_velocity(cal, roi, g['tx_rad'], g['rx_rad'], g['prf_hz'],
                                 g['f0_hz'], g['sound_speed_mps'], g['phase_sign'])
    tx, rx = np.asarray(g['tx_rad']), np.asarray(g['rx_rad'])
    geom = np.stack([np.cos(tx)+np.cos(rx), np.sin(tx)+np.sin(rx)], 1)
    rate = g['phase_sign']*2*np.pi*g['f0_hz']/g['sound_speed_mps']*(geom@v)
    ct = np.arange(-len(cal), 0)[:,None]/g['prf_hz'] + np.asarray(g['view_offsets_s'])[None]
    aligned = cal*np.exp(-1j*ct*rate[None])[:,:,None,None]
    ref = aligned[:,0,roi]
    gains = np.array([np.vdot(ref, aligned[:,j,roi])/max(np.vdot(ref,ref).real,1e-30)
                      for j in range(cal.shape[1])])
    if np.any(abs(gains)<1e-6):
        raise ValueError('Unusable view calibration')
    times = np.arange(len(public.views))[:,None]/g['prf_hz'] + np.asarray(g['view_offsets_s'])[None]
    modulation = gains[None]*np.exp(1j*times*rate[None])
    reference = np.mean(abs(aligned/gains[None,:,None,None]), axis=(0,1))
    diagnostics={'estimated_velocity_zx_mps':v.tolist(),
        'relative_gains_real':gains.real.tolist(), 'relative_gains_imag':gains.imag.tolist(),
        'calibration_sha256':array_hash(cal), 'training_roi_sha256':array_hash(roi)}
    result=(modulation[:,:,None,None],reference,diagnostics)
    if return_complex_reference:return (*result,np.mean(aligned/gains[None,:,None,None],axis=0))
    return result


def estimate_displacement(views, reference, roi, mode, iterations=30, reference_views=None, options=None):
    """Forward registration with independent spatial validation and zero fallback.

    See registration.py for the complex-reference model and explicit ungated
    diagnostic arms. The envelope-only call remains available for simple inputs.
    """
    options=options or MotionGate(max_evaluations=iterations)
    return register_displacement(views,reference,roi,mode,options,reference_views)


def bounded_moments(x, packet, horizon, rho=1., lags=(0,)):
    """Exact exponentially weighted packet moments with expiration at horizon."""
    if (packet < 1 or horizon < packet or not lags or min(lags) < 0
            or horizon % packet or x.shape[1] % packet or not 0 < rho <= 1 or max(lags)>=packet):
        raise ValueError('History and record must contain whole packets; lags < packet')
    keep = horizon//packet; queue = []; state = {}; snapshots = []
    for stop in range(packet,x.shape[1]+1,packet):
        block=x[:,stop-packet:stop]; moments={}
        for lag in lags:
            a,b=block[:,lag:],block[:,:packet-lag]
            moments[lag]=(a@b.conj().T/a.shape[1],a.mean(1,keepdims=True),b.mean(1,keepdims=True))
            if lag not in state:state[lag]=[0.,np.zeros((len(x),len(x)),complex),np.zeros((len(x),1),complex),np.zeros((len(x),1),complex)]
            mass,c,mp,mm=state[lag];v,p,m=moments[lag]
            state[lag]=[rho*mass+1,rho*c+v,rho*mp+p,rho*mm+m]
        queue.append(moments)
        if len(queue)>keep:
            old=queue.pop(0)
            for lag in lags:
                v,p,m=old[lag];s=state[lag];s[0]-=rho**keep;s[1]-=rho**keep*v;s[2]-=rho**keep*p;s[3]-=rho**keep*m
        snapshots.append({lag:tuple(a.copy() if isinstance(a,np.ndarray) else a for a in values)
                          for lag,values in state.items()})
    return snapshots


@dataclass
class SharedResult:
    views: np.ndarray
    diagnostics: dict
    apply_frozen: object = field(repr=False)
    source_diagnostic: object = field(default=None,repr=False)


class SharedContext:
    """Cached data-only fits. Known-motion controls use an explicit separate path."""
    def __init__(self, public, packet=16, score_start=64, motion='none', known_motion=None, motion_options=None):
        n=public.views.shape[0]
        if packet < 1 or n%packet or score_start%packet or not 0<=score_start<n:
            raise ValueError('Whole packets and a nonempty scoring interval are required')
        if motion=='known':
            if known_motion is None or np.shape(known_motion)!=(n,2,*public.views.shape[2:]):
                raise ValueError('Known motion is an explicitly evaluation-only control')
        elif known_motion is not None:
            raise ValueError('Known motion cannot enter a data-only arm')
        self.public=public;self.packet=packet;self.score_start=score_start
        self.modulation,reference,cal,reference_views=calibrate(public,return_complex_reference=True)
        if motion=='known':
            self.displacement=np.asarray(known_motion);md={'registration':'known evaluation-only warp'}
        else:self.displacement,md=estimate_displacement(public.views/self.modulation,reference,public.training_roi,motion,
            reference_views=reference_views,options=MotionGate(**motion_options) if motion_options else None)
        self.motion=motion;self.calibration=cal;self.motion_diagnostics=md
        if 'audit_rows' in md:
            md['audit_sha256']=hashlib.sha256(json.dumps(md['audit_rows'],separators=(',',':'),allow_nan=False).encode()).hexdigest()
        self.motion_summary={k:v for k,v in md.items() if k!='audit_rows'}
        self.has_displacement=bool(np.any(self.displacement))
        grid=np.indices(public.views.shape[2:],dtype=float)
        self.alignment_coordinates=grid[None]+self.displacement
        self.restoration_coordinates=(np.stack([inverse_coordinates(d) for d in self.displacement])
                                      if self.has_displacement else self.alignment_coordinates)
        self.aligned=self.align(public.views);self.x=self.aligned.reshape(n,-1).T
        self.cache={};self.stops=list(range(score_start+packet,n+1,packet))

    def align(self, other):
        value=other/self.modulation
        return (np.stack([sample_frame(f,q,zero_extension=True) for f,q in zip(value,self.alignment_coordinates)])
                if self.has_displacement else value)

    def restore_block(self, aligned, start):
        stop=start+len(aligned)
        value=(np.stack([sample_frame(f,q,zero_extension=True) for f,q in zip(aligned,self.restoration_coordinates[start:stop])])
               if self.has_displacement else aligned)
        return value*self.modulation[start:stop]

    def history(self, stop, horizon, rho=1., centered=False):
        start=max(0,stop-horizon);x=self.x[:,start:stop]
        weights=rho**((stop-1-np.arange(start,stop))//self.packet)
        mean=(x@weights[:,None])/weights.sum() if centered else np.zeros((len(x),1),complex)
        return x-mean,weights,mean

    def decomposition(self, stop, horizon, rho=1., centered=False):
        key=('svd',stop,horizon,rho,centered)
        if key not in self.cache:
            x,w,mean=self.history(stop,horizon,rho,centered)
            u,s,_=linalg.svd(x*np.sqrt(w)[None],full_matrices=False,lapack_driver='gesdd')
            self.cache[key]=(u,s,mean)
        return self.cache[key]

    def moments(self, horizon, rho, lags):
        key=('moments',horizon,rho,tuple(lags))
        if key not in self.cache:self.cache[key]=bounded_moments(self.x,self.packet,horizon,rho,lags)
        return self.cache[key]


def _lag_basis(context, stop, spec):
    horizon=spec['horizon'];rho=spec.get('rho',.75);rank=spec.get('rank',4)
    centered=spec.get('centered',False);lags=tuple(spec.get('lags',[1,2,4,8]))
    ridge=spec.get('ridge_rel',1e-3);threshold=spec.get('epsilon_rel',1e-6)
    key=('lag_basis',stop,horizon,rho,rank,centered,lags,ridge,threshold)
    if key in context.cache:return context.cache[key]
    state=context.moments(horizon,rho,(0,)+lags)[stop//context.packet-1]
    covs={}
    for lag,(mass,c,mp,mm) in state.items():
        covs[lag]=c/mass-(mp@mm.conj().T/mass**2 if centered else 0.)
    ev,u=linalg.eigh(covs[0],subset_by_index=[max(0,len(covs[0])-rank),len(covs[0])-1])
    top=max(float(ev[-1]),1e-30);valid=ev>top*threshold;ev=ev[valid];u=u[:,valid]
    mean=state[0][2]/state[0][0] if centered else np.zeros((len(context.x),1),complex)
    if not len(ev):raise ValueError('Empty retained lag basis')
    whitening=u/np.sqrt(ev+ridge*top)[None]
    stack=[]
    for lag in lags:
        c=whitening.conj().T@covs[lag]@whitening
        stack.extend([(c+c.conj().T)/2,(c-c.conj().T)/(2j)])
    q,info=joint_diagonalize(stack,max_sweeps=40,tol=1e-6)
    demix=q.conj().T@whitening.conj().T;mix=np.linalg.pinv(demix)
    hist,_,_=context.history(stop,horizon,rho,centered)
    sources=demix@hist
    _,score=component_weights(sources,mix,lags,phase_width=spec.get('phase_width',.1))
    info.update({'rank':len(ev),'whitening_eigenvalue_condition':float(ev[-1]/ev[0]),
                 'minimum_pairwise_lag_signature_distance':score['minimum_pairwise_lag_signature_distance']})
    value=(demix,mix,mean,info,sources,lags)
    context.cache[key]=value
    return value


def fit_shared(context, spec):
    """Fit using public observations only; return all views on the common interval."""
    method=spec['method'];horizon=spec.get('horizon',64);packet=context.packet
    if horizon%packet or horizon<packet:raise ValueError('Horizon must be a positive packet multiple')
    if method=='identity' and context.motion!='none':
        raise ValueError('Identity uses no registration; use the explicit roundtrip control')
    shape=context.public.views.shape;plans=[];records=[];pieces=[];p=context.x.shape[0]
    for stop in context.stops:
        start=stop-packet;current=context.x[:,start:stop];centered=spec.get('centered',False)
        mean=np.zeros((p,1),complex);detail={};basis=None
        if method in ('identity','roundtrip'):
            operator=np.eye(p)
        elif method in ('svd','weighted_svd','baranger','lee_aer'):
            rho=spec.get('rho',.75) if method=='weighted_svd' else 1.
            if method=='lee_aer':centered=True
            u,s,mean=context.decomposition(stop,horizon,rho,centered)
            if method=='baranger':
                k,h,detail=baranger_band(u);operator=u[:,k:h]@u[:,k:h].conj().T
            elif method=='lee_aer':
                domain=spec.get('domain','temporal');den=p if domain=='temporal' else min(stop,horizon)-1
                ev=s*s/den
                eps=spec.get('epsilon_abs')
                if eps is None:eps=spec.get('epsilon_rel',1e-3)*max(float(ev[0]),1e-30)
                if eps<=0:raise ValueError('Shared AER requires positive regularization')
                base=eps**-.5;g=(ev+eps)**-.5
                scale=1/base if spec.get('unit_operator',False) else 1.
                operator=scale*(base*np.eye(p)+(u*(g-base)[None])@u.conj().T)
                detail={'domain':domain,'covariance_denominator':den,'epsilon_absolute':eps,
                        'unit_operator_extension':spec.get('unit_operator',False),'operator_gain_max':scale*base}
            else:
                discard=spec.get('discard',1)
                if discard=='all':operator=np.zeros((p,p),complex);k=len(s)
                else:k=min(int(discard),len(s));operator=np.eye(p)-u[:,:k]@u[:,:k].conj().T
                detail={'discard':k,'rho':rho}
        elif method=='axis1':
            rho=spec.get('rho',.75);state=context.moments(horizon,rho,(0,))[stop//packet-1][0]
            mass,c,mp,_=state;mean=mp/mass if centered else mean
            cov=c/mass-(mean@mean.conj().T if centered else 0.)
            discard=spec.get('discard',1)
            if discard=='all':operator=np.zeros((p,p),complex);k=min(stop,horizon)
            else:
                k=min(int(discard),min(stop,horizon))
                if k < 0:raise ValueError('Discard must be nonnegative')
                if k:
                    _,u=linalg.eigh(cov,subset_by_index=[p-k,p-1])
                    operator=np.eye(p)-u@u.conj().T
                else:operator=np.eye(p)
            detail={'discard':k,'rho':rho,'moment_mass':float(mass),'bounded_expiration':True}
        elif method=='axis2':
            demix,mix,mean,info,sources,lags=_lag_basis(context,stop,spec)
            weights,score=component_weights(sources,mix,lags,phase_width=spec.get('phase_width',.1),
                                            coherence_floor=spec.get('coherence_floor',.1))
            # Exact small-rank SVD: no full spatial SVD is needed to cap gain.
            qm,rm=linalg.qr(mix,mode='economic');qd,rd=linalg.qr(demix.conj().T,mode='economic')
            a,d,b=linalg.svd((rm*weights[None])@rd.conj().T,full_matrices=False)
            gain=float(d[0]);cap=spec.get('gain_cap',1.)
            operator=(qm@a*(np.minimum(d,cap) if cap is not None else d)[None])@(b@qd.conj().T)
            detail={**info,'blood_weights':weights.tolist(),'unconstrained_operator_gain_max':gain,
                    'operator_gain_max':min(gain,cap) if cap is not None else gain}
            basis=(demix,mix,weights)
        elif method=='stap':
            # In the common calibrated frame, the nominal bulk-tissue manifold
            # is constant in time and view. Drift/mismatch remains in the data.
            d=packet*shape[1];steer=np.ones((d,1),complex);orth=null_projector(steer)
            blocks=context.aligned[max(0,stop-horizon):stop].reshape(-1,packet,shape[1],*shape[2:])
            matrices=blocks.reshape(len(blocks),d,-1)
            mask=context.public.training_roi.ravel();zz,xx=np.indices(shape[2:]);ops=[];counts=[];fallback=0;gains=[]
            for iz in range(shape[2]):
                for ix in range(shape[3]):
                    eligible=mask&(np.maximum(abs(zz-iz),abs(xx-ix)).ravel()>spec.get('guard',1))
                    samples=matrices[:,:,eligible].transpose(1,0,2).reshape(d,-1)
                    counts.append(samples.shape[1]);op=orth
                    if spec.get('adaptive',False) and samples.shape[1]>=4:
                        op=null_projector(steer,samples@samples.conj().T/samples.shape[1],spec.get('loading_rel',.01))
                    elif spec.get('adaptive',False):fallback+=1
                    ops.append(op);gains.append(float(linalg.svdvals(op)[0]))
            operator=np.stack(ops)
            detail={'snapshot_dimension':d,'training_snapshots_min':min(counts),'training_snapshots_max':max(counts),
                    'fits_below_twice_dimension':sum(k<2*d for k in counts),'fallbacks':fallback,
                    'training_snapshots_not_assumed_independent':True,'operator_gain_max':max(gains),
                    'history_packets':len(blocks),'adaptive':spec.get('adaptive',False)}
        else:raise ValueError('Unknown shared method: '+method)
        if method=='stap':
            block=context.aligned[start:stop].reshape(packet*shape[1],-1)
            out=np.einsum('pij,jp->ip',operator,block).reshape(packet,*shape[1:])
        else:out=(operator@(current-mean)).T.reshape(packet,*shape[1:])
        pieces.append(context.restore_block(out,start))
        plans.append((start,stop,operator,mean,basis,centered))
        records.append({'start':start,'stop':stop,'history_start':max(0,stop-horizon),**detail})
    result=np.concatenate(pieces)
    # Identity does not undergo registration; roundtrip is its explicit loss control.
    def frozen(other):
        aligned=context.align(other);x=aligned.reshape(len(aligned),-1).T;out=[]
        for start,stop,operator,_,_,centered in plans:
            if method=='stap':
                block=aligned[start:stop].reshape(packet*shape[1],-1)
                value=np.einsum('pij,jp->ip',operator,block).reshape(packet,*shape[1:])
            else:
                rho=spec.get('rho',.75) if method in ('axis1','axis2','weighted_svd') else 1.
                left=max(0,stop-horizon);weights=rho**((stop-1-np.arange(left,stop))//packet)
                m=x[:,left:stop]@weights[:,None]/weights.sum() if centered else np.zeros((p,1),complex)
                value=(operator@(x[:,start:stop]-m)).T.reshape(packet,*shape[1:])
            out.append(context.restore_block(value,start))
        return np.concatenate(out)
    def source_diagnostic(blood):
        if method!='axis2':return None
        candidates=[];heuristics=[];success=[]
        for start,stop,_,mean,basis,_ in plans:
            demix,mix,w=basis;sources=demix@(context.x[:,start:stop]-mean)
            columns=[]
            for j in range(len(w)):
                contribution=(mix[:,j,None]*sources[j,None,:]).T.reshape(packet,*shape[1:])
                columns.append(context.restore_block(contribution,start).ravel())
            design=np.stack(columns,axis=1);target=blood[start:stop].ravel()
            fit=optimize.lsq_linear(np.vstack([design.real,design.imag]),np.r_[target.real,target.imag],
                                    bounds=(0.,1.),tol=1e-8,max_iter=100)
            candidates.append((design@fit.x).reshape(packet,*shape[1:]));heuristics.append((design@w).reshape(packet,*shape[1:]));success.append(bool(fit.success))
        truth=blood[context.score_start:];norm=max(np.linalg.norm(truth),1e-30)
        return {'evaluation_only':True,'selection_use':False,'solver_success':all(success),
                'best_diagonal_uncapped_nrmse':float(np.linalg.norm(np.concatenate(candidates)-truth)/norm),
                'heuristic_uncapped_nrmse':float(np.linalg.norm(np.concatenate(heuristics)-truth)/norm),
                'interpretation':'truth-fitted real weights in [0,1] in the fixed uncapped basis; diagnostic only, not an upper bound on every capped operator'}
    return SharedResult(result,{'method':method,'motion':context.motion,'horizon_frames':horizon,
        'packet_frames':packet,'latency_max_frames':packet,'spatial_tile':list(shape[2:]),'input_views':shape[1],
        'calibration':context.calibration,'motion_fit':context.motion_summary,'blocks':records,
        'resampling':'continuous complex bilinear zero extension; generator retains its historical boundary rule',
        'frozen_attribution':'conditional on this operator representation'},frozen,source_diagnostic)
