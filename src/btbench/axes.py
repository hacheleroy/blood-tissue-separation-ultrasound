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


def apply_motion(iq,motion):
    """Replay fitted registration/demodulation, with no new motion estimation."""
    out=np.array(iq,dtype=complex)
    if len(out)!=len(motion["phase"]):
        raise ValueError("Frozen motion length mismatch")
    if motion["rigid"]:
        for t,s in enumerate(motion["shifts"]):
            out[t]=shift_complex(out[t],s)
    return out*np.exp(-1j*motion["phase"])[:,None,None]


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
    return q,{"sweeps":sweep+1,"ajd_converged":bool(largest<tol),
              "ajd_largest_rotation":float(largest),
              "ajd_relative_offdiag":float(np.linalg.norm(off)/max(np.linalg.norm(mats),1e-30))}


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
    separation=np.linalg.norm(corr[:,None,:]-corr[None,:,:],axis=2)/np.sqrt(len(lags))
    np.fill_diagonal(separation,np.inf)
    return weights,{"blood_weights":weights.tolist(),"phase_rad_per_frame":phase.tolist(),
                    "normalized_lag_signatures_real":corr.real.tolist(),
                    "normalized_lag_signatures_imag":corr.imag.tolist(),
                    "minimum_pairwise_lag_signature_distance":float(separation.min()) if len(corr)>1 else None,
                    "lag1_coherence":coh.tolist(),"spatial_extent":extent.tolist()}


def _separate(x,covs,lags,rank=8,epsilon_rel=1e-6,ica=False,seed=0,ica_options=None,
              whitening_regularization_rel=0.,operator_gain_cap=None):
    if (not np.isfinite([epsilon_rel,whitening_regularization_rel]).all() or
            not 0<=epsilon_rel<1 or whitening_regularization_rel<0 or
            (operator_gain_cap is not None and (not np.isfinite(operator_gain_cap) or operator_gain_cap<=0))):
        raise ValueError("Invalid whitening threshold, regularization or operator gain cap")
    c0 = (covs[0]+covs[0].conj().T)/2
    ev,u = linalg.eigh(c0)
    idx = np.flatnonzero(ev>max(ev[-1],1e-30)*epsilon_rel)[-rank:]
    if not len(idx):
        return np.zeros_like(x),{"rank":0},np.zeros((len(x),len(x)),complex)
    ridge=whitening_regularization_rel*max(ev[-1],1e-30)
    b = u[:,idx]/np.sqrt(ev[idx]+ridge)[None,:]
    z = b.conj().T@x
    if ica:
        q,info = complex_fastica_rotation(z,seed=seed,**(ica_options or {}))
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
    operator=(mixing*weights)@demix
    gain=float(linalg.svdvals(operator)[0])
    if operator_gain_cap is not None and gain>operator_gain_cap:
        left,s,right=linalg.svd(operator,full_matrices=False)
        operator=(left*np.minimum(s,operator_gain_cap))@right
    return operator@x,{"rank":len(idx),**info,**score,
        "whitening_eigenvalue_condition":float(ev[idx[-1]]/ev[idx[0]]),
        "whitening_regularization_absolute":float(ridge),
        "whitening_residual":float(np.linalg.norm(b.conj().T@c0@b-np.eye(len(idx)))/np.sqrt(len(idx))),
        "unconstrained_operator_gain_max":gain,"operator_gain_cap":operator_gain_cap,
        "operator_projection":"spectral-norm ball" if operator_gain_cap is not None else "none",
        "operator_gain_max":float(linalg.svdvals(operator)[0])},operator


def complex_fastica_rotation(z,seed=0,max_iter=300,tol=1e-7):
    """Symmetric circular-complex kurtosis ICA in a whitened reduced space."""
    if max_iter<1 or tol<=0 or z.ndim!=2 or not np.isfinite(z).all():
        raise ValueError("Invalid ICA data or stopping settings")
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
    y=w@z
    covariance=z@z.conj().T/z.shape[1]
    pseudo=z@z.T/z.shape[1]
    power=np.mean(abs(y)**2,axis=1)
    cumulant=np.mean(abs(y)**4,axis=1)-2*power**2-abs(np.mean(y*y,axis=1))**2
    return w.conj().T,{"ica_iterations":it+1,"ica_converged":converged,
        "ica_final_change":float(change),"ica_tolerance":float(tol),
        "ica_whitening_error":float(np.linalg.norm(covariance-np.eye(len(z)))/np.sqrt(len(z))),
        "ica_circularity_ratio":float(np.linalg.norm(pseudo)/max(np.linalg.norm(covariance),1e-30)),
        "ica_source_normalized_fourth_cumulants":(cumulant/np.maximum(power**2,1e-30)).tolist()}


def _packet_pipeline(iq,kind,packet=64,centered=False,patch_shape=(8,8),stride=(4,4),
                     tissue_roi=None,rigid=False,rho=0.9,q_min=0.3,rank=8,
                     discard=1,lags=(1,2,4,8,16),use_memory=True,seed=0,
                     ica_options=None,ica_nonconvergence="svd_fallback",
                     phase_compensation=True,motion_phase_error_rad_per_frame=0.,max_operator_bytes=512*1024**2,
                     epsilon_rel=1e-6,whitening_regularization_rel=0.,operator_gain_cap=None):
    if ica_nonconvergence not in ("svd_fallback","retain","error"):
        raise ValueError("Unknown ICA nonconvergence policy")
    if (not 0<=q_min<=1 or rank<1 or int(rank)!=rank or
            not np.isfinite(motion_phase_error_rad_per_frame) or max_operator_bytes<=0 or
            not lags or any(int(l)!=l or l<1 for l in lags)):
        raise ValueError("Invalid packet rank, confidence, lag, phase or memory settings")
    aligned,motion = compensate(iq,tissue_roi,rigid)
    if not phase_compensation:
        motion["phase"]=np.zeros(len(iq))
    motion["phase"]=motion["phase"]+motion_phase_error_rad_per_frame*np.arange(len(iq))
    aligned=apply_motion(iq,motion)
    if packet <= max(lags) or len(iq) % packet and len(iq) % packet <= max(lags):
        # chunks merges only 1-2 frames; avoid silent unequal lag supports.
        if kind in ("sobi","ica"):
            raise ValueError("Packets, including last packet, must exceed the largest lag")
    y = np.zeros_like(aligned)
    den = np.zeros(iq.shape[1:])
    records = []
    plans=[]
    operator_bytes=0
    for patch_id,(sl,w) in enumerate(patches(iq.shape[1:],patch_shape,stride)):
        allx = casorati(aligned[(slice(None),)+sl])
        out = np.empty_like(allx)
        memory = LagMemory((0,)+tuple(lags) if kind != "svd_memory" else (0,),rho)
        packet_plans=[]
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
                operator_bytes+=(ts.stop-ts.start)**2*16
                if operator_bytes>max_operator_bytes:raise MemoryError("Frozen packet operators exceed memory cap")
                packet_plans.append((ts,confidence,True,fb.apply_frozen))
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
                operator=np.eye(len(x))-tissue@tissue.conj().T
                out[:,ts] = operator@x
                info = {"discard":discard}
            else:
                out[:,ts],info,operator = _separate(x,covs,lags,rank,epsilon_rel,ica=kind=="ica",seed=seed,
                    ica_options=ica_options,whitening_regularization_rel=whitening_regularization_rel,
                    operator_gain_cap=operator_gain_cap)
                if kind=="ica" and not info.get("ica_converged",True):
                    if ica_nonconvergence=="error":
                        raise ValueError("ICA did not converge under the frozen stopping rule")
                    if ica_nonconvergence=="svd_fallback":
                        _,u=linalg.eigh((covs[0]+covs[0].conj().T)/2)
                        if not 0<=discard<=len(u):raise ValueError("Invalid fallback rejection rank")
                        tissue=u[:,-discard:] if discard else u[:,:0]
                        operator=np.eye(len(x))-tissue@tissue.conj().T
                        out[:,ts]=operator@x
                        info["ica_candidate_operator_gain_max"]=info["operator_gain_max"]
                        info["operator_gain_max"]=float(linalg.svdvals(operator)[0])
                        info["fallback"]="local_svd_ica_nonconvergence"
            packet_plans.append((ts,confidence,False,operator))
            operator_bytes+=operator.nbytes
            if operator_bytes>max_operator_bytes:raise MemoryError("Frozen packet operators exceed memory cap")
            records.append({"patch":patch_id,"start":ts.start,"q":confidence,"reset":False,**info})
        y[(slice(None),)+sl] += frames(out,(len(iq),)+w.shape)*w
        den[sl] += w
        plans.append((sl,w,packet_plans))
    y,valid = restore(y/den,motion)
    def frozen(other):
        if other.shape!=iq.shape:
            raise ValueError("Frozen packet-operator shape mismatch")
        aligned_other=apply_motion(other,motion)
        output=np.zeros_like(other,dtype=complex)
        for sl,w,packet_plans in plans:
            raw_all=casorati(aligned_other[(slice(None),)+sl])
            filtered=np.empty_like(raw_all)
            mass=0.;mean_sum=np.zeros((len(raw_all),1),complex)
            for ts,confidence,reset,operator in packet_plans:
                raw=raw_all[:,ts]
                if not use_memory or reset:
                    mass=0.;mean_sum[:]=0
                if reset:
                    filtered[:,ts]=casorati(operator(frames(raw,(raw.shape[1],)+w.shape)))
                    continue
                if centered:
                    # Recompute each COMPONENT's linear pooled mean using the
                    # mixture's frozen confidence/reset schedule. Subtracting
                    # the mixture mean from every component would be affine.
                    weight=(1-rho)*confidence
                    mass=rho*mass+weight
                    mean_sum=rho*mean_sum+weight*raw.mean(axis=1,keepdims=True)
                    raw=raw-mean_sum/mass
                filtered[:,ts]=operator@raw
            output[(slice(None),)+sl]+=frames(filtered,(len(other),)+w.shape)*w
        return restore(output/den,motion)[0]
    return FilterResult(y,{"centered":centered,"packets":records,
        "motion_phase_rad":motion["phase"].tolist(),"motion_shifts_pixels":motion["shifts"].tolist(),
        "valid_fraction":float(valid.mean()),"lab_frame_phase_restored":True,
        "registration":"integer rigid 2D" if rigid else "bulk IQ phase only",
        "phase_compensation":phase_compensation,
        "uncentered_axis1_global_phase_invariant":kind=="svd_memory" and not centered,
        "frozen_operator_bytes_estimate":operator_bytes,
        "motion_phase_error_rad_per_frame":float(motion_phase_error_rad_per_frame),
        "component_classifier":"data-only heuristic" if kind != "svd_memory" else None,
        "frozen_parameters":"motion, confidence/reset schedule, covariance bases and component weights",
        "frozen_centering":"component-wise replay of confidence-weighted pooled mean"},frozen)


def motion_packet_filter(iq,**kwargs):
    return _packet_pipeline(iq,"svd_memory",**kwargs)


def multilag_filter(iq,**kwargs):
    return _packet_pipeline(iq,"sobi",**kwargs)


def ica_filter(iq,centered=True,**kwargs):
    """ICA uses zero-mean input by default; uncentered is an explicit ablation."""
    return _packet_pipeline(iq,"ica",use_memory=False,centered=centered,**kwargs)
