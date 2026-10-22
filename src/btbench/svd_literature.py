"""Source-documented SVD/EVD estimators. See docs/SVD_LITERATURE_AUDIT.md.

All rank boundaries are zero-based half-open: reject [0,k), retain [k,h).
These are independent numerical implementations, with explicit fitting choices.
"""
import warnings
import numpy as np
from scipy import linalg, optimize, special
from .core import FilterResult, casorati, frames, prepare, chunks, patches


def spatial_similarity(u):
    """Pearson corr(|U|), including the diagonal (Baranger 2018, Eq. 5).

    A spatially constant envelope has undefined Pearson correlation. Such
    columns are assigned zero correlation and counted, never silently NaN.
    """
    a = np.abs(u).astype(float)
    a -= a.mean(axis=0, keepdims=True)
    norm = np.linalg.norm(a, axis=0)
    valid = norm > 64*np.finfo(float).eps*np.maximum(np.linalg.norm(np.abs(u), axis=0), 1e-30)
    a[:, valid] /= norm[valid]
    a[:, ~valid] = 0
    return np.clip(a.T@a, -1, 1), int((~valid).sum())


def square_fit(ssm, min_tissue=1, min_blood=1, max_tissue=None, max_blood=None):
    """Maximize corr2(SSM, two adjacent all-one squares), Baranger 2023.

    The remaining noise region and the cross-block regions are zero in the
    template. Summed-area tables give O(r^2) candidate search. This objective
    is independently checked against explicitly constructed binary templates.
    """
    c = np.asarray(ssm, float)
    if c.ndim != 2 or c.shape[0] != c.shape[1] or not np.isfinite(c).all():
        raise ValueError("SSM must be a finite square real matrix")
    r = len(c)
    if min_tissue < 0 or min_blood < 1 or min_tissue+min_blood > r:
        raise ValueError("Infeasible tissue/blood rank bounds")
    mt = r-min_blood if max_tissue is None else min(int(max_tissue), r-min_blood)
    mb = r if max_blood is None else int(max_blood)
    if mt < min_tissue or mb < min_blood:
        raise ValueError("Infeasible maximum rank bounds")
    table = np.pad(c, ((1,0),(1,0))).cumsum(0).cumsum(1)
    total = c.sum()
    variance = np.sum((c-c.mean())**2)
    if variance <= np.finfo(float).eps*max(np.sum(c*c), 1e-30):
        return min_tissue, r, {"ssm_informative":False,"fallback":"retain after minimum tissue rank"}
    best = None
    n = r*r
    count = 0
    # Deterministic ties: smaller tissue rank, then smaller upper boundary.
    for k in range(int(min_tissue), mt+1):
        h = np.arange(k+int(min_blood), min(r,k+mb)+1)
        ones = k*k+(h-k)**2
        dot = table[k,k]+table[h,h]-table[k,h]-table[h,k]+table[k,k]
        den = np.sqrt(variance*np.maximum(ones-ones*ones/n, 0))
        score = np.divide(dot-total*ones/n, den, out=np.full(len(h),-np.inf), where=den>0)
        count += len(h)
        if len(h):
            j = int(np.argmax(score))
            if best is None or score[j] > best[0]:
                best = (float(score[j]), k, int(h[j]))
    if best is None or not np.isfinite(best[0]):
        raise ValueError("No nonconstant square template in requested rank range")
    return best[1],best[2],{"ssm_informative":True,"ssm_template_correlation":best[0],
        "ssm_candidates":count,"ssm_modes":r,"ssm_objective":"corr2 with tissue and blood squares",
        "ssm_diagonal_included":True,"ssm_mode_truncation":False}


def baranger_band(u, **options):
    c, undefined = spatial_similarity(u)
    k,h,d = square_fit(c, **options)
    d["constant_envelope_modes"] = undefined
    return k,h,d


def temporal_statistics(vh, prf_hz):
    """Folded mean absolute Doppler frequency; uniform sampling, two-sided PSD.

    The absolute moment avoids cancellation of positive and negative frequency.
    Conjugation of V does not affect this absolute-frequency statistic.
    """
    if prf_hz <= 0:
        raise ValueError("PRF must be positive")
    f = np.fft.fftfreq(vh.shape[1], 1/prf_hz)
    psd = abs(np.fft.fft(vh,axis=1))**2
    return psd@abs(f)/np.maximum(psd.sum(axis=1), 1e-30)


def song_band(s, vh, prf_hz, slope_threshold=0.1, tail_tolerance=0.1,
              noise_rejection=True):
    """Song local spectrum thresholds; numerical conventions fully logged.

    Lower boundary: first log10(sigma/sigma[0]) slope flatter than the supplied
    dB-independent threshold (decades per rank). Upper boundary: sigmoid fit
    of folded mean frequencies, followed by linear singular-value tail fitting.
    The noise transition/tail procedure is described in Ozgun & Byram (2021),
    Section III-A, citing Song (2017). Original author software is unavailable;
    the sigmoid fit domain, deviation rule and fallbacks are benchmark choices.
    """
    r = len(s)
    if r < 3 or slope_threshold <= 0 or tail_tolerance < 0:
        raise ValueError("Need >=3 modes and valid slope/tail parameters")
    curve = np.log10(np.maximum(s/max(s[0],1e-30),1e-15))
    slope = -np.diff(curve)
    flat = np.flatnonzero(slope <= slope_threshold)
    k = int(flat[0]+1) if len(flat) else r-1
    info = {"slope_threshold_decades_per_rank":float(slope_threshold),
        "tail_tolerance_relative":float(tail_tolerance),
        "implementation":"independent Song-style threshold reconstruction",
        "original_author_code_parity":False}
    if not noise_rejection or k>=r-2:
        return k,r,{**info,"noise_rejection":False}
    freq = temporal_statistics(vh, prf_hz)/prf_hz
    n = np.arange(r,dtype=float)
    # Normalize rank to make fitting tolerances independent of packet length.
    z = n/max(r-1,1)
    def sigmoid(q, lo, hi, center, width):
        return lo+(hi-lo)*special.expit((q-center)/width)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", optimize.OptimizeWarning)
            p,_ = optimize.curve_fit(sigmoid,z,freq,p0=[min(freq),max(freq),.5,.1],
                bounds=([0,0,0,.005],[.5,.5,1,1]),maxfev=4000)
        pre = int(np.clip(np.ceil(p[2]*(r-1)), k+1, r-2))
        coef = np.polyfit(n[pre:],s[pre:],1)
        fit = np.polyval(coef,n)
        dev = (s-fit)/np.maximum(abs(fit),s[0]*1e-12)
        above = np.flatnonzero((n>=k)&(n<pre)&(dev>tail_tolerance))
        hi = max(k+1,int(above[-1])+1) if len(above) else max(k+1,pre)
        return k,min(hi,r),{**info,"sigmoid_converged":True,"noise_pre_cutoff":pre,
            "noise_line_slope":float(coef[0]),"noise_line_intercept":float(coef[1])}
    except (RuntimeError,ValueError,optimize.OptimizeWarning):
        return k,r,{**info,"sigmoid_converged":False,"noise_fallback":"no upper truncation"}


def gram_decompose(x):
    """Temporal-Gram EVD: complex Hermitian, descending, nonnegative spectrum."""
    ev,v = linalg.eigh(x.conj().T@x,driver="evd")
    rank = min(x.shape)
    ev,v = np.maximum(ev[::-1][:rank],0),v[:,::-1][:,:rank]
    s = np.sqrt(ev)
    u = x@v
    good = ev > np.finfo(float).eps*max(x.shape)*max(ev[0],1e-30)
    u[:,good] /= s[good]
    u[:,~good] = 0
    return u,s,v.conj().T


def randomized_svd_filter(iq, discard=1, oversampling=5, power_iterations=2,
                          spatial_fraction=1., seed=0, centered=False, window=None):
    """Leading-clutter randomized range finder with optional spatial subsampling.

    Song et al., TUFFC 2017: speed approximations, not new separation criteria.
    Reconstruct by X-X V_k V_k^H, preserving the uncomputed residual spectrum.
    Full-rank range is an exactness control; no oracle rank selection here.
    """
    if not 0<spatial_fraction<=1 or discard<0 or oversampling<0 or power_iterations<0:
        raise ValueError("Invalid randomized SVD settings")
    rng = np.random.default_rng(seed)
    out = np.empty_like(iq,dtype=complex)
    operators,diag = [],[]
    for sl in chunks(len(iq), window or len(iq)):
        x,_ = prepare(casorati(iq[sl]),centered)
        p,t = x.shape
        if discard>min(p,t):
            raise ValueError("Clutter rank exceeds available rank")
        nrows = min(p,max(discard,int(np.ceil(p*spatial_fraction))))
        rows = rng.choice(p,nrows,replace=False)
        a = x[rows]
        qrank = min(min(a.shape),discard+oversampling)
        if discard == 0:
            v = np.zeros((t,0),complex)
        else:
            omega=(rng.standard_normal((t,qrank))+1j*rng.standard_normal((t,qrank)))/np.sqrt(2)
            q,_ = linalg.qr(a@omega,mode="economic")
            for _ in range(int(power_iterations)):
                z,_=linalg.qr(a.conj().T@q,mode="economic")
                q,_=linalg.qr(a@z,mode="economic")
            _,_,vh=linalg.svd(q.conj().T@a,full_matrices=False)
            v=vh[:discard].conj().T
        out[sl]=frames(x-(x@v)@v.conj().T,iq[sl].shape)
        operators.append((sl,v))
        diag.append({"start":sl.start,"stop":sl.stop,"sampled_pixels":nrows,
                     "range_rank":qrank,"discard":int(discard)})
    def apply(other):
        if other.shape!=iq.shape:raise ValueError("Frozen-operator shape mismatch")
        y=np.empty_like(other,dtype=complex)
        for sl,v in operators:
            a,_=prepare(casorati(other[sl]),centered)
            y[sl]=frames(a-(a@v)@v.conj().T,other[sl].shape)
        return y
    return FilterResult(out,{"centered":centered,"seed":int(seed),"windows":diag,
        "spatial_fraction":float(spatial_fraction),"power_iterations":int(power_iterations),
        "oversampling":int(oversampling),"computed_subspace":"leading clutter only"},apply)


def huang_weights(u,s,rule="ncc",beta=1.43,neighbors=5):
    """Huang et al. arXiv:2405.11105 Eqs. 7-9, using lambda=sigma^2.

    NCC is the uncentered normalized inner product of local sums of |U|,
    not Pearson correlation. The sigmoid midpoint uses one-based mode order.
    """
    r=len(s)
    if r<2 or beta<=0 or neighbors<1:raise ValueError("Invalid EVD threshold settings")
    ev=s*s
    d={"rule":rule,"beta":float(beta),"neighbors":int(neighbors)}
    if rule=="median":
        threshold=beta*np.median(ev)
        k=int(np.sum(ev>threshold))
        d["eigenvalue_threshold"]=float(threshold)
    elif rule in ("ncc","sigmoid"):
        a=abs(u)
        cs=np.vstack([np.zeros((1,len(a))),np.cumsum(a.T,axis=0)])
        values=[]
        fn=min(int(neighbors),r//2)
        candidates=range(fn,r-fn+1)
        for k in candidates:
            left=cs[k]-cs[max(0,k-neighbors)]
            right=cs[min(r,k+neighbors)]-cs[k]
            norm=np.linalg.norm(left)*np.linalg.norm(right)
            values.append(float(left@right/norm) if norm>1e-30 else 1.)
        k=list(candidates)[int(np.argmin(values))]
        d["neighbors_actual"]=fn
        d["minimum_image_ncc"]=float(min(values))
    else:
        raise ValueError("Supported Huang rules: median, ncc, sigmoid")
    wt=special.expit(beta*(np.arange(1,r+1)-k)) if rule=="sigmoid" else (np.arange(r)>=k).astype(float)
    return wt,{**d,"discard":k}


def huang_evd_filter(iq,centered=False,rule="ncc",beta=1.43,neighbors=5,
                     window=None,spatial_weight=None):
    """Weighted statistics are fitted to W X, then applied to original X.

    For W=I this is the published global EVD projection. Local IQ combination
    below is a phase-preserving extension; the paper's output is combined PD.
    """
    out=np.empty_like(iq,dtype=complex);ops=[];diag=[]
    w=np.ones(iq.shape[1:]) if spatial_weight is None else np.asarray(spatial_weight,float)
    if w.shape!=iq.shape[1:] or not np.isfinite(w).all() or np.any(w<0) or not w.any():
        raise ValueError("Spatial weight must be a nonnegative nonempty finite image")
    for sl in chunks(len(iq),window or len(iq)):
        x,_=prepare(casorati(iq[sl]),centered)
        u,s,vh=gram_decompose(w.reshape(-1,1)*x)
        wt,d=huang_weights(u,s,rule,beta,neighbors)
        v=vh.conj().T
        out[sl]=frames(((x@v)*wt)@v.conj().T,iq[sl].shape)
        ops.append((sl,v,wt));diag.append(d)
    def apply(other):
        if other.shape!=iq.shape:raise ValueError("Frozen-operator shape mismatch")
        y=np.empty_like(other,dtype=complex)
        for sl,v,wt in ops:
            a,_=prepare(casorati(other[sl]),centered)
            y[sl]=frames(((a@v)*wt)@v.conj().T,other[sl].shape)
        return y
    return FilterResult(out,{"centered":centered,"windows":diag,
        "covariance_weight":"W^2 inside X^H W^2 X","complex_hermitian":True},apply)


def gaussian_local_evd(iq,grid_shape=(2,2),sigma_pixels=None,**kwargs):
    """Huang full-FOV Gaussian covariance localization + explicit IQ extension.

    Published PD: sum_j alpha_j mean_t|W_j X F_j|^2 with alpha_j the reciprocal
    mean subset power. Stored separately from the phase-preserving IQ output.
    """
    z,x=iq.shape[1:]
    grid=tuple(map(int,grid_shape))
    if len(grid)!=2 or min(grid)<1:raise ValueError("grid_shape must contain two positive sizes")
    sigma=np.asarray(sigma_pixels if sigma_pixels is not None else (z/grid[0],x/grid[1]),float)
    if sigma.shape!=(2,) or np.any(sigma<=0):raise ValueError("Need two positive Gaussian widths")
    centers=[np.linspace(0,n-1,g+2)[1:-1] for n,g in zip((z,x),grid)]
    zz,xx=np.mgrid[:z,:x]
    y=np.zeros_like(iq,dtype=complex);den=np.zeros((z,x));pd=np.zeros((z,x));ops=[];diag=[]
    for cz in centers[0]:
        for cx in centers[1]:
            w=np.exp(-.5*((zz-cz)/sigma[0])**2-.5*((xx-cx)/sigma[1])**2)
            r=huang_evd_filter(iq,spatial_weight=w,**kwargs)
            weighted_pd=np.mean(abs(r.iq*w)**2,axis=0)
            alpha=1/max(float(weighted_pd.mean()),1e-30)
            pd+=alpha*weighted_pd
            y+=w*w*r.iq;den+=w*w;ops.append((w*w,r.apply_frozen))
            diag.append({"center_zx_pixels":[float(cz),float(cx)],"pd_alpha":alpha,**r.diagnostics})
    def apply(other):
        return sum(w*op(other) for w,op in ops)/den
    result=FilterResult(y/den,{"subsets":diag,"grid_shape":list(grid),"sigma_pixels":sigma.tolist(),
        "iq_reconstruction":"normalized W^2 overlap-add extension",
        "paper_output":"power_doppler; subset mean-power normalization changes amplitude"},apply)
    result.power_doppler=pd
    return result


def song_local_filter(iq,patch_shape=(8,8),stride=(4,4),**kwargs):
    """Rectangular overlapping blocks and independent local adaptive cutoffs.

    Return equal-weight coherent IQ overlap-add for velocity, with the
    equal-weight incoherent block-power aggregation recorded separately.
    """
    from .filters import svd_filter
    y=np.zeros_like(iq,dtype=complex);den=np.zeros(iq.shape[1:]);pd=den.copy();ops=[];diag=[]
    for sl,_ in patches(iq.shape[1:],patch_shape,stride):
        r=svd_filter(iq[(slice(None),)+sl],selector="song",**kwargs)
        y[(slice(None),)+sl]+=r.iq;den[sl]+=1
        pd[sl]+=np.mean(abs(r.iq)**2,axis=0)
        ops.append((sl,r.apply_frozen));diag.append(r.diagnostics)
    def apply(other):
        out=np.zeros_like(other,dtype=complex)
        for sl,op in ops:out[(slice(None),)+sl]+=op(other[(slice(None),)+sl])
        return out/den
    result=FilterResult(y/den,{"patches":diag,"block_weight":"rectangular",
        "iq_reconstruction":"coherent normalized overlap-add extension",
        "threshold_fidelity":"Song-style independent reconstruction; see audit"},apply)
    result.power_doppler=pd/den
    return result
