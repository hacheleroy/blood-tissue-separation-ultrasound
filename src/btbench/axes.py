"""Axes 1 and 2: rigid/phase compensation, confidence memory, complex SOBI.

Heuristic component scores and rigid registration are research prototypes.
Nonrigid motion and a trained/validated component classifier are future work.
"""
import numpy as np
from scipy import linalg, ndimage

from .core import FilterResult,casorati,frames,patches,chunks,prepare
from .filters import svd_filter


def shift_complex(a, shift):
    return ndimage.shift(a.real,shift,order=1,mode="constant",cval=0,prefilter=False)+1j*ndimage.shift(a.imag,shift,order=1,mode="constant",cval=0,prefilter=False)


def compensate(iq,tissue_roi=None,rigid=False,max_shift=3):
    """Bright-tissue phase demodulation; optional integer 2-D envelope alignment.

    Phase estimates use only tissue_roi (or brightest 25% of time-mean power).
    Output is in tissue coordinates; inverse restores lab-frame IQ phase.
    """
    if tissue_roi is None:
        power = np.mean(abs(iq)**2,axis=0)
        tissue_roi = power >= np.percentile(power,75)
    tissue_roi = np.asarray(tissue_roi,bool)
    if tissue_roi.shape != iq.shape[1:] or not tissue_roi.any():
        raise ValueError("Need nonempty tissue ROI with matching grid")
    y = np.array(iq,dtype=complex)
    shifts = np.zeros((len(iq),2),int)
    ref = abs(iq[0])*tissue_roi
    if rigid:
        for t in range(1,len(iq)):
            cur = abs(iq[t])*tissue_roi
            corr = np.fft.ifftn(np.fft.fftn(ref)*np.fft.fftn(cur).conj()).real
            loc = np.asarray(np.unravel_index(corr.argmax(),corr.shape))
            loc = np.where(loc>np.array(corr.shape)//2,loc-np.array(corr.shape),loc)
            if np.max(abs(loc)) <= max_shift:
                shifts[t] = loc
                y[t] = shift_complex(iq[t],loc)
    phase = np.zeros(len(iq))
    q = np.ones(len(iq))
    for t in range(1,len(iq)):
        a,b = y[t,tissue_roi],y[t-1,tissue_roi]
        cross = np.vdot(b,a)
        phase[t] = phase[t-1]+np.angle(cross)
        q[t] = np.clip(abs(cross)/max(np.linalg.norm(a)*np.linalg.norm(b),1e-30),0,1)
    y *= np.exp(-1j*phase)[:,None,None]
    return y,{"phase":phase,"shifts":shifts,"confidence":q,"rigid":rigid}


def restore(iq,motion):
    y = iq*np.exp(1j*motion["phase"])[:,None,None]
    valid = np.ones(iq.shape,bool)
    if motion["rigid"]:
        for t,s in enumerate(motion["shifts"]):
            y[t] = shift_complex(y[t],-s)
            # Round trip loses boundary samples; expose validity, never wrap.
            a = shift_complex(np.ones(iq.shape[1:],complex),s)
            valid[t] = shift_complex(a,-s).real > 0.999
    return y,valid


class LagMemory:
    """Exponentially weighted raw pair moments; center after merging.

    Separate confidence mass normalizes q-weighting. For tau>0 retain both
    endpoint means, preserving the complex non-Hermitian lag information.
    """
    def __init__(self,lags,rho=0.9):
        if not 0 <= rho < 1:
            raise ValueError("rho must lie in [0,1)")
        self.lags = tuple(lags)
        self.rho = rho
        self.state = {}

    def update(self,x,q=1.,centered=False):
        for lag in self.lags:
            if lag >= x.shape[1]:
                raise ValueError("Every lag must be smaller than every packet")
            later = x[:,lag:]
            earlier = x[:,:x.shape[1]-lag]
            c = later@earlier.conj().T/later.shape[1]
            mp,mm = later.mean(1,keepdims=True),earlier.mean(1,keepdims=True)
            if lag not in self.state:
                self.state[lag] = [0.,np.zeros_like(c),np.zeros_like(mp),np.zeros_like(mm)]
            mass,cs,ps,ms = self.state[lag]
            f = (1-self.rho)*q
            self.state[lag] = [self.rho*mass+f,self.rho*cs+f*c,
                               self.rho*ps+f*mp,self.rho*ms+f*mm]
        return self.get(centered)

    def get(self,centered=False):
        out = {}
        for lag,(mass,c,mp,mm) in self.state.items():
            if mass <= 0:
                raise ValueError("No reliable covariance state")
            out[lag] = c/mass-(mp@mm.conj().T/mass**2 if centered else 0)
        return out


def joint_diagonalize(matrices,max_sweeps=60,tol=1e-8):
    """Complex unitary Jacobi AJD of a stack of Hermitian matrices.

    Each non-Hermitian lag covariance is supplied as Hermitian real AND
    Hermitian imaginary parts; dropping the second part loses Doppler sign.
    Returned Q obeys Q^H M_j Q approximately diagonal.
    """
    mats = np.array(matrices,dtype=complex,copy=True)
    r = mats.shape[1]
    q = np.eye(r,dtype=complex)
    for sweep in range(max_sweeps):
        largest = 0.
        for p in range(r-1):
            for j in range(p+1,r):
                h = np.stack(((mats[:,p,p]-mats[:,j,j]).real,
                              2*mats[:,p,j].real,2*mats[:,p,j].imag))
                _,v = np.linalg.eigh(h@h.T)
                n = v[:,-1]
                if n[0] < 0:
                    n = -n
                c = np.sqrt(max((1+n[0])/2,0))
                s = (n[1]-1j*n[2])/max(2*c,1e-30)
                largest = max(largest,abs(s))
                g = np.eye(r,dtype=complex)
                g[p,p]=g[j,j]=c
                g[j,p]=s
                g[p,j]=-s.conjugate()
                mats = g.conj().T@mats@g
                q = q@g
        if largest < tol:
            break
    off = mats.copy()
    off[:,np.arange(r),np.arange(r)] = 0
    return q,{"sweeps":sweep+1,"ajd_relative_offdiag":float(np.linalg.norm(off)/max(np.linalg.norm(mats),1e-30))}


def component_weights(sources,mixing,lags,phase_width=0.10,coherence_floor=0.1):
    """Data-only heuristic: coherent near-DC broad modes are tissue; incoherent
    modes are noise. Slow blood sharing tissue statistics can still be lost.
    Scores are logged and are NOT a validated tissue/blood classifier.
    """
    power = np.mean(abs(sources)**2,axis=1)
    correlations = []
    for lag in lags:
        correlations.append(np.mean(sources[:,lag:]*sources[:,:-lag].conj(),axis=1)/np.maximum(power,1e-30))
    corr = np.stack(correlations,axis=1)
    r1 = corr[:,0]
    phase = abs(np.angle(r1))/lags[0]
    coh = np.clip(abs(r1),0,1)
    a = abs(mixing)**2
    extent = (a.sum(0)**2)/(np.maximum((a*a).sum(0),1e-30)*len(a))
    dc_tissue = np.exp(-(phase/phase_width)**2)*coh*np.sqrt(extent)
    flow_coherence = np.clip((coh-coherence_floor)/(1-coherence_floor),0,1)
    weights = np.clip((1-dc_tissue)*flow_coherence,0,1)
    return weights,{"blood_weights":weights.tolist(),"phase_rad_per_frame":phase.tolist(),
                    "lag1_coherence":coh.tolist(),"spatial_extent":extent.tolist()}


def _separate(x,covs,lags,rank=8,epsilon_rel=1e-6,ica=False,seed=0):
    c0 = (covs[0]+covs[0].conj().T)/2
    ev,u = linalg.eigh(c0)
    idx = np.flatnonzero(ev>max(ev[-1],1e-30)*epsilon_rel)[-rank:]
    if not len(idx):
        return np.zeros_like(x),{"rank":0}
    b = u[:,idx]/np.sqrt(ev[idx])[None,:]
    z = b.conj().T@x
    if ica:
        q,info = complex_fastica_rotation(z,seed=seed)
    else:
        stack = []
        for lag in lags:
            a = b.conj().T@covs[lag]@b
            stack += [(a+a.conj().T)/2,(a-a.conj().T)/(2j)]
        q,info = joint_diagonalize(stack)
    demix = q.conj().T@b.conj().T
    sources = demix@x
    mixing = np.linalg.pinv(demix)
    weights,score = component_weights(sources,mixing,lags)
    return mixing@(weights[:,None]*sources),{"rank":len(idx),**info,**score}


def complex_fastica_rotation(z,seed=0,max_iter=300,tol=1e-7):
    """Symmetric circular-complex kurtosis ICA in a whitened reduced space."""
    rng = np.random.default_rng(seed)
    w = rng.normal(size=(len(z),len(z)))+1j*rng.normal(size=(len(z),len(z)))
    def orth(a):
        d,u = linalg.eigh(a@a.conj().T)
        return (u/np.sqrt(np.maximum(d,1e-15)))@u.conj().T@a
    w = orth(w)
    converged = False
    for it in range(max_iter):
        y = w@z
        new = (y*abs(y)**2)@z.conj().T/z.shape[1]-2*np.mean(abs(y)**2,axis=1)[:,None]*w
        new = orth(new)
        change = np.max(abs(abs(np.diag(new@w.conj().T))-1))
        w = new
        if change < tol:
            converged=True
            break
    return w.conj().T,{"ica_iterations":it+1,"ica_converged":converged}


def _packet_pipeline(iq,kind,packet=64,centered=False,patch_shape=(8,8),stride=(4,4),
                     tissue_roi=None,rigid=False,rho=0.9,q_min=0.3,rank=8,
                     discard=1,lags=(1,2,4,8,16),use_memory=True,seed=0):
    aligned,motion = compensate(iq,tissue_roi,rigid)
    if packet <= max(lags) or len(iq) % packet and len(iq) % packet <= max(lags):
        # chunks merges only 1-2 frames; avoid silent unequal lag supports.
        if kind in ("sobi","ica"):
            raise ValueError("Packets, including last packet, must exceed the largest lag")
    y = np.zeros_like(aligned)
    den = np.zeros(iq.shape[1:])
    records = []
    for patch_id,(sl,w) in enumerate(patches(iq.shape[1:],patch_shape,stride)):
        allx = casorati(aligned[(slice(None),)+sl])
        out = np.empty_like(allx)
        memory = LagMemory((0,)+tuple(lags) if kind != "svd_memory" else (0,),rho)
        for ts in chunks(len(iq),packet):
            raw = allx[:,ts]
            confidence = float(np.min(motion["confidence"][ts]))
            reset = confidence < q_min
            if not use_memory or reset:
                memory = LagMemory(memory.lags,rho)
            if reset:
                # A low-confidence motion event uses a short local SVD fallback;
                # no unreliable registration enters covariance memory.
                fb = svd_filter(frames(raw,(raw.shape[1],)+w.shape),centered,discard=discard)
                out[:,ts] = casorati(fb.iq)
                records.append({"patch":patch_id,"start":ts.start,"q":confidence,"reset":True,"fallback":"local_svd"})
                continue
            covs = memory.update(raw,confidence,centered)
            x,mean = prepare(raw,centered)
            if centered:
                # Same pooled mean for covariance and transform (lag-zero mean).
                mass,_,mp,_ = memory.state[0]
                x = raw-mp/mass
            if kind == "svd_memory":
                ev,u = linalg.eigh((covs[0]+covs[0].conj().T)/2)
                if discard < 0 or discard > len(ev):
                    raise ValueError("Invalid rejection rank")
                tissue = u[:,-discard:] if discard else u[:,:0]
                out[:,ts] = x-tissue@(tissue.conj().T@x)
                info = {"discard":discard}
            else:
                out[:,ts],info = _separate(x,covs,lags,rank,ica=kind=="ica",seed=seed)
            records.append({"patch":patch_id,"start":ts.start,"q":confidence,"reset":False,**info})
        y[(slice(None),)+sl] += frames(out,(len(iq),)+w.shape)*w
        den[sl] += w
    y,valid = restore(y/den,motion)
    return FilterResult(y,{"centered":centered,"packets":records,
        "motion_phase_rad":motion["phase"].tolist(),"motion_shifts_pixels":motion["shifts"].tolist(),
        "valid_fraction":float(valid.mean()),"lab_frame_phase_restored":True,
        "registration":"integer rigid 2D" if rigid else "bulk IQ phase only",
        "component_classifier":"data-only heuristic" if kind != "svd_memory" else None})


def motion_packet_filter(iq,**kwargs):
    return _packet_pipeline(iq,"svd_memory",**kwargs)


def multilag_filter(iq,**kwargs):
    return _packet_pipeline(iq,"sobi",**kwargs)


def ica_filter(iq,**kwargs):
    return _packet_pipeline(iq,"ica",use_memory=False,**kwargs)
