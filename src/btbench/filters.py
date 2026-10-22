"""Established algebraic baselines plus explicitly named threshold surrogates."""
import numpy as np
from scipy import linalg, signal

from .core import FilterResult, casorati, frames, prepare, chunks, patches


def decompose(x):
    return linalg.svd(x, full_matrices=False, lapack_driver="gesdd")


def ssm_cutoff(u, max_modes=64):
    """Envelope Pearson SSM + two-block contrast with integral-image sums.

    This is a documented Baranger-inspired surrogate, NOT the published
    three-region sum-table implementation. No use of velocity references.
    """
    r = min(u.shape[1], max_modes)
    if r < 3:
        return min(1, r-1), {"ssm_contrast": 0.0}
    a = np.abs(u[:, :r])
    a -= a.mean(axis=0)
    a /= np.maximum(np.linalg.norm(a, axis=0), 1e-30)
    ssm = a.T @ a
    np.fill_diagonal(ssm, 0)
    ii = np.pad(ssm, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    def area(a, b, c, d):
        return ii[b, d]-ii[a, d]-ii[b, c]+ii[a, c]
    scores = []
    for k in range(1, r):
        within_a = area(0,k,0,k)/max(k*(k-1),1)
        within_b = area(k,r,k,r)/max((r-k)*(r-k-1),1)
        between = area(0,k,k,r)/(k*(r-k))
        scores.append((within_a+within_b)/2-between)
    k = int(np.argmax(scores))+1
    return k, {"ssm_contrast": float(max(scores)), "ssm_search_modes": r}


def choose_rank(u, s, selector="fixed", discard=1, energy_fraction=0.95, **kwargs):
    if selector == "fixed":
        k, diag = int(discard), {}
    elif selector == "energy":
        if not 0 < energy_fraction < 1:
            raise ValueError("energy_fraction must be in (0,1)")
        energy = np.cumsum(s*s)
        k = int(np.searchsorted(energy, energy_fraction*energy[-1]))+1 if energy[-1] else 0
        k, diag = min(k, len(s)-1), {"energy_fraction": energy_fraction}
    elif selector == "ssm_surrogate":
        k, diag = ssm_cutoff(u, **kwargs)
    else:
        raise ValueError("Unknown SVD rank selector: "+selector)
    if not 0 <= k <= len(s):
        raise ValueError("discard outside matrix rank range")
    return k, diag


def svd_filter(iq, centered=False, selector="fixed", discard=1, upper=None,
               window=None, energy_fraction=0.95, selector_options=None, backend="svd"):
    """Retain zero-based singular indices [discard, upper); no mean add-back."""
    y = np.empty_like(iq, dtype=np.complex128)
    info = []
    operators = []
    for sl in chunks(len(iq), window or len(iq)):
        x, _ = prepare(casorati(iq[sl]), centered)
        if backend == "svd":
            u, s, vh = decompose(x)
        elif backend == "gram":
            from .svd_literature import gram_decompose
            u, s, vh = gram_decompose(x)
        else:
            raise ValueError("backend must be svd or gram")
        options = dict(selector_options or {})
        auto_hi = len(s)
        if selector == "baranger":
            from .svd_literature import baranger_band
            k, auto_hi, d = baranger_band(u, **options)
        elif selector == "song":
            from .svd_literature import song_band
            k, auto_hi, d = song_band(s,vh,**options)
        elif selector == "temporal_frequency":
            from .svd_literature import temporal_statistics
            freq = temporal_statistics(vh,options["prf_hz"])
            good = np.flatnonzero(freq > options["cutoff_hz"])
            k = int(good[0]) if len(good) else len(s)
            d = {"cutoff_hz":options["cutoff_hz"],"frequency_statistic":"folded PSD mean absolute frequency"}
        else:
            k, d = choose_rank(u,s,selector,discard,energy_fraction,**options)
        hi = auto_hi if upper is None else int(upper)
        if not k <= hi <= len(s):
            raise ValueError("Require 0 <= discard <= upper <= min(P,T)")
        v = vh[k:hi].conj().T
        y[sl] = frames((x@v)@v.conj().T, iq[sl].shape)
        info.append({"start":sl.start,"stop":sl.stop,"discard":k,"upper":hi,**d})
        operators.append((sl,vh[k:hi].conj().T))
    def frozen(other):
        if other.shape!=iq.shape:raise ValueError("Frozen-operator shape mismatch")
        out=np.empty_like(other,dtype=complex)
        for sl,v in operators:
            x,_=prepare(casorati(other[sl]),centered)
            out[sl]=frames((x@v)@v.conj().T,other[sl].shape)
        return out
    return FilterResult(y,{"centered":centered,"backend":backend,"windows":info},frozen)


class SpatialMoments:
    """Exact streamed second moment / globally centered covariance, complex-safe."""
    def __init__(self, p, max_bytes=512*1024**2):
        if p*p*16 > max_bytes:
            raise MemoryError("Dense spatial covariance exceeds memory cap; use local patches")
        self.n = 0
        self.mean = np.zeros((p,1), complex)
        self.m2 = np.zeros((p,p), complex)

    def update(self, x):
        n = x.shape[1]
        if not n:
            return
        mean = x.mean(axis=1,keepdims=True)
        xc = x-mean
        delta = mean-self.mean
        total = self.n+n
        self.m2 += xc@xc.conj().T + (self.n*n/total)*(delta@delta.conj().T)
        self.mean += delta*(n/total)
        self.n = total

    def covariance(self, centered=False):
        if not self.n:
            raise ValueError("No frames accumulated")
        c = self.m2/self.n
        if not centered:
            c = c+self.mean@self.mean.conj().T
        return (c+c.conj().T)/2


def aer_filter(iq, centered=False, domain="temporal", epsilon_rel=1e-3,
               alpha=0.5, epsilon_abs=None, normalization="samples",
               stream_chunk=64, max_cov_bytes=512*1024**2,epsilon_rule="relative",window=None):
    """Exact regularized spectral remapping. For C=X^H X/P or XX^H/T.

    normalization='samples' matches conventional covariance normalization.
    'common' uses denominator T for both sides to test their mathematical
    equivalence. epsilon_rel is relative to the largest eigenvalue; this is
    an explicit scale-invariant benchmark choice, not Lee's absolute sweep.
    Streaming spatial mode estimates global moments first, applies second pass.
    """
    if alpha < 0 or alpha > 0.5 or epsilon_rel < 0:
        raise ValueError("Require 0<=alpha<=0.5 and nonnegative epsilon_rel")
    if window is not None and window<len(iq):
        out=np.empty_like(iq,dtype=complex)
        ops=[]
        for sl in chunks(len(iq),window):
            r=aer_filter(iq[sl],centered,domain,epsilon_rel,alpha,epsilon_abs,
                         normalization,stream_chunk,max_cov_bytes,epsilon_rule)
            out[sl]=r.iq
            ops.append((sl,r))
        def apply(other):
            y=np.empty_like(other,dtype=complex)
            for sl,r in ops:y[sl]=r.apply_frozen(other[sl])
            return y
        return FilterResult(out,{"window":window,"windows":[r.diagnostics for _,r in ops]},apply)
    x = casorati(iq)
    p, t = x.shape
    if normalization not in ("samples","common"):
        raise ValueError("Unknown covariance normalization")
    if domain == "spatial":
        acc = SpatialMoments(p,max_cov_bytes)
        for sl in chunks(t,stream_chunk):
            acc.update(x[:,sl])
        c = acc.covariance(centered)
        x = x-acc.mean if centered else x
        den = t
    elif domain == "temporal":
        x, _ = prepare(x,centered)
        den = p if normalization == "samples" else t
        if t*t*16 > max_cov_bytes:
            raise MemoryError("Temporal covariance exceeds cap; reduce ensemble or use local spatial AER")
        c = x.conj().T@x/den
    else:
        raise ValueError("domain must be temporal or spatial")
    ev, q = linalg.eigh(c)
    ev = np.maximum(ev,0)
    if epsilon_rule not in ("relative","noise_floor"):
        raise ValueError("Unknown epsilon_rule")
    valid = ev > max(ev[-1],1e-30)*np.finfo(float).eps*max(c.shape)
    if epsilon_abs is not None:
        eps=float(epsilon_abs)
    elif epsilon_rule=="noise_floor":
        nonzero=ev[valid]
        eps=float(np.median(nonzero[:max(len(nonzero)//2,1)])) if len(nonzero) else 1e-30
    else:
        eps=float(epsilon_rel*max(ev[-1],1e-30))
    if eps < 0:
        raise ValueError("epsilon_abs must be nonnegative")
    # Pseudoinverse convention: exact null modes stay zero when epsilon=0.
    g = np.zeros_like(ev)
    if eps > 0:
        g = (ev+eps)**(-alpha)
    else:
        g[valid] = ev[valid]**(-alpha)
    if alpha == 0:
        g[:] = 1
    if domain == "spatial":
        y = q@(g[:,None]*(q.conj().T@x))
    else:
        y = ((x@q)*g)@q.conj().T
    def frozen(other):
        if other.shape!=iq.shape:raise ValueError("Frozen-operator shape mismatch")
        a,_=prepare(casorati(other),centered)
        if domain=="spatial":b=q@(g[:,None]*(q.conj().T@a))
        else:b=((a@q)*g)@q.conj().T
        return frames(b,other.shape)
    return FilterResult(frames(y,iq.shape),{"centered":centered,"domain":domain,
        "epsilon_absolute":eps,"epsilon_rel":epsilon_rel,"epsilon_rule":epsilon_rule,"alpha":alpha,
        "covariance_denominator":den,"numerical_rank":int(valid.sum()),
        "amplitude_calibration":"changed; no truth-based gain normalization"},frozen)


def local_filter(iq, method="svd", patch_shape=(8,8), stride=(4,4), **kwargs):
    y = np.zeros_like(iq,dtype=complex)
    den = np.zeros(iq.shape[1:])
    diagnostics = []
    operators=[]
    fn = {"svd":svd_filter,"aer":aer_filter}[method]
    for sl, w in patches(iq.shape[1:],patch_shape,stride):
        result = fn(iq[(slice(None),)+sl], **kwargs)
        y[(slice(None),)+sl] += result.iq*w
        den[sl] += w
        diagnostics.append(result.diagnostics)
        operators.append((sl,w,result.apply_frozen))
    def frozen(other):
        out=np.zeros_like(other,dtype=complex)
        for sl,w,op in operators:
            out[(slice(None),)+sl]+=op(other[(slice(None),)+sl])*w
        return out/den
    return FilterResult(y/den,{"patches":diagnostics,"overlap_min_weight":float(den.min())},frozen)


def sliding_filter(iq,method="svd",window=256,hop=64,**kwargs):
    """Contiguous long/short ensembles with positive temporal overlap/add.

    Noncausal offline reconstruction. Edge locations use fewer ensembles.
    The same wrapper can give SVD and AER identical temporal support/budgets.
    """
    n=len(iq);window=min(window,n)
    if window<3 or not 1<=hop<=window:
        raise ValueError("Require >=3-frame window and 1<=hop<=window")
    starts=list(range(0,n-window+1,hop))
    if starts[-1]!=n-window:starts.append(n-window)
    w=np.exp(-2*np.linspace(-1,1,window)**2)[:,None,None]
    y=np.zeros_like(iq,dtype=complex);mass=np.zeros((n,1,1))
    ops=[]
    fn={"svd":svd_filter,"aer":aer_filter}[method]
    for start in starts:
        sl=slice(start,start+window)
        r=fn(iq[sl],**kwargs)
        y[sl]+=w*r.iq;mass[sl]+=w
        ops.append((sl,r))
    def frozen(other):
        out=np.zeros_like(other,dtype=complex)
        for sl,r in ops:out[sl]+=w*r.apply_frozen(other[sl])
        return out/mass
    return FilterResult(y/mass,{"window":window,"hop":hop,"offline":True,
        "windows":[{"start":sl.start,"stop":sl.stop,**r.diagnostics} for sl,r in ops]},frozen)


def temporal_highpass(iq, cutoff_hz, prf_hz, order=3, phase_mode="zero_phase"):
    """Pixelwise complex-IQ Butterworth wall filter along slow time.

    ``zero_phase`` uses forward/backward filtering (squared magnitude response).
    ``causal`` uses a one-pass filter initialized at each pixel's first sample.
    Both are linear, including their edge initialization/padding operations.
    """
    iq = np.asarray(iq)
    if (iq.ndim < 2 or not np.iscomplexobj(iq) or not np.isfinite(iq).all()
            or not np.isfinite([cutoff_hz, prf_hz]).all()
            or not 0 < cutoff_hz < prf_hz/2 or int(order) != order or order < 1):
        raise ValueError("Need finite complex IQ, integer order >=1, 0<cutoff<PRF/2")
    if phase_mode not in ("zero_phase", "causal"):
        raise ValueError("phase_mode must be zero_phase or causal")
    sos = signal.butter(int(order), cutoff_hz, fs=prf_hz, btype="highpass", output="sos")
    padlen = 3*(2*len(sos)+1-min((sos[:, 2] == 0).sum(), (sos[:, 5] == 0).sum()))
    if phase_mode == "zero_phase" and len(iq) <= padlen:
        raise ValueError(f"Zero-phase filter needs more than {padlen} frames")
    def apply(other):
        if other.shape != iq.shape:
            raise ValueError("Frozen-operator shape mismatch")
        if phase_mode == "zero_phase":
            return signal.sosfiltfilt(sos, other, axis=0, padlen=int(padlen))
        zi = signal.sosfilt_zi(sos).reshape((len(sos), 2)+(1,)*(other.ndim-1))
        return signal.sosfilt(sos, other, axis=0, zi=zi*other[0])[0]
    return FilterResult(apply(iq), {"cutoff_hz":float(cutoff_hz), "prf_hz":float(prf_hz),
        "order":int(order), "phase_mode":phase_mode, "zero_phase_offline":phase_mode=="zero_phase",
        "padding_frames":int(padlen) if phase_mode=="zero_phase" else 0,
        "initialization":"odd extension" if phase_mode=="zero_phase" else "first-sample steady state",
        "magnitude_response_power":2 if phase_mode=="zero_phase" else 1}, apply)


def rpca_filter(iq, centered=False, sparse_lambda=None, max_iter=300, tol=1e-7):
    """Complex principal-component pursuit, inexact ALM; sparse output is candidate blood."""
    x, _ = prepare(casorati(iq),centered)
    lam = 1/np.sqrt(max(x.shape)) if sparse_lambda is None else sparse_lambda
    if lam <= 0:
        raise ValueError("sparse_lambda must be positive")
    norm = np.linalg.norm(x,"fro")
    if norm == 0:
        return FilterResult(np.zeros_like(iq),{"iterations":0,"converged":True})
    spectral = linalg.svdvals(x)[0]
    dual = x/max(spectral,np.abs(x).max()/lam)
    mu = 1.25/max(spectral,1e-30)
    mu_max = mu*1e7
    sparse = np.zeros_like(x)
    converged = False
    for it in range(max_iter):
        u,s,vh = decompose(x-sparse+dual/mu)
        low = (u*np.maximum(s-1/mu,0))@vh
        a = x-low+dual/mu
        sparse = np.maximum(1-lam/(mu*np.maximum(np.abs(a),1e-30)),0)*a
        residual = x-low-sparse
        rel = np.linalg.norm(residual,"fro")/norm
        if rel < tol:
            converged = True
            break
        dual += mu*residual
        mu = min(mu*1.5,mu_max)
    return FilterResult(frames(sparse,iq.shape),{"iterations":it+1,"converged":converged,
        "relative_residual":float(rel),"lambda":lam,
        "assumption":"blood sparse in pixel-time; dense blood violates this model"})
