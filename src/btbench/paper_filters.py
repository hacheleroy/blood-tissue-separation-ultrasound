"""PDF-audited Song 2017, Maresca 2018 and Lee 2026 algorithm variants.

Unspecified turning-point tolerances are explicit numerical choices, not
claims of author-software parity. See docs/PDF_METHOD_AUDIT.md.
"""
import numpy as np
from .core import FilterResult, chunks, patches


def lag_one_mode_frequency(vh,prf_hz):
    """Magnitude of lag-one autocorrelation frequency, Song 2017 Fig. 2(d)."""
    if prf_hz<=0 or vh.ndim!=2 or vh.shape[1]<2:
        raise ValueError("Need temporal singular vectors and positive revisit PRF")
    return abs(np.angle(np.sum(vh[:,1:]*vh[:,:-1].conj(),axis=1)))*prf_hz/(2*np.pi)


def song_pdf_band(s,vh,prf_hz,tissue_cutoff_hz=10.,slope_threshold=0.1,
                  noise_tail_fraction=.25,tail_tolerance_decades=.05,noise_rejection=True):
    """Two lower cutoffs; lag-one noise plateau; straight fit to log singular values.

    The paper does not give a numerical turning-point/plateau detector. Here
    the slope tolerance, tail fraction and deviation tolerance are frozen
    implementation choices. No sigmoid fit is attributed to Song's PDF.
    """
    s=np.asarray(s,float);r=len(s)
    if (r<3 or vh.shape[0]!=r or not 0<=tissue_cutoff_hz<prf_hz/2 or
            slope_threshold<=0 or not 0<noise_tail_fraction<=1 or tail_tolerance_decades<0):
        raise ValueError("Invalid Song PDF threshold settings")
    curve=np.log10(np.maximum(s/max(s[0],1e-30),1e-15))
    flat=np.flatnonzero(-np.diff(curve)<=slope_threshold)
    lower_slope=int(flat[0]+1) if len(flat) else r-1
    freq=lag_one_mode_frequency(vh,prf_hz)
    moving=np.flatnonzero(freq>tissue_cutoff_hz)
    lower_frequency=int(moving[0]) if len(moving) else r
    k=max(lower_slope,lower_frequency)
    info={"source":"Song 2017 pp.252-253 and Eq.1","original_author_code_parity":False,
          "lower_slope":lower_slope,"lower_frequency":lower_frequency,
          "frequency_statistic":"absolute lag-one autocorrelation frequency",
          "tissue_cutoff_hz":float(tissue_cutoff_hz),"mode_frequency_hz":freq.tolist(),
          "slope_threshold_decades_per_rank":float(slope_threshold),
          "noise_tail_fraction":float(noise_tail_fraction),
          "tail_tolerance_decades":float(tail_tolerance_decades),
          "unspecified_paper_choices":"slope tolerance, plateau support and deviation tolerance"}
    if not noise_rejection or k>=r-2:
        return k,r,{**info,"noise_rejection_applied":False}
    tail=max(k+1,r-max(2,int(np.ceil(noise_tail_fraction*r))))
    noise_frequency=float(freq[tail:].mean())
    hits=np.flatnonzero((np.arange(r)>=k+1)&(freq>=noise_frequency))
    pre=int(np.clip(hits[0] if len(hits) else tail,k+1,r-2))
    ranks=np.arange(r)
    coef=np.polyfit(ranks[pre:],curve[pre:],1)
    deviation=curve-np.polyval(coef,ranks)
    above=np.flatnonzero((ranks>=k)&(ranks<pre)&(deviation>tail_tolerance_decades))
    h=max(k+1,int(above[-1]+1)) if len(above) else max(k+1,pre)
    return k,min(h,r),{**info,"noise_rejection_applied":True,"noise_pre_cutoff":pre,
        "noise_plateau_frequency_hz":noise_frequency,"noise_line_domain":"log10 singular values",
        "noise_line_slope":float(coef[0]),"noise_line_intercept":float(coef[1])}


def maresca_cutoff(s,offset=26,log_floor=1e-12):
    """Maresca 2018 Eqs.2-3: one-based curvature-peak index plus fixed k.

    k=26 was manually calibrated on the authors' first 60-frame ensemble.
    It is an operator choice, not a universal data-only threshold. No upper
    noise threshold is introduced. Interior second differences are explicit.
    """
    s=np.asarray(s,float)
    if len(s)<3 or int(offset)!=offset or not 0<log_floor<1:
        raise ValueError("Need >=3 modes, integer offset and positive log floor")
    c=s/max(s.sum(),1e-30)
    db=20*np.log10(np.maximum(c,log_floor))
    peak=int(np.argmax(np.diff(db,n=2)))+2  # one-based interior mode index
    k=int(np.clip(peak+offset,0,len(s)))
    return k,{"source":"Maresca 2018 Eqs.2-3","curvature_peak_one_based":peak,
              "operator_selected_offset":int(offset),"discard_unclipped":peak+int(offset),
              "normalized_singular_values_db":db.tolist(),"upper_noise_cutoff":False,
              "offset_calibration":"fixed operator choice; original paper used k=26 at 60 frames"}


def song_pdf_local_filter(iq,patch_shape=(8,8),stride=(4,4),window=None,
                          centered=False,selector_options=None,normalization="paper"):
    """Song Eq.1: normalize each block by sum of retained sigma, coherently average.

    normalization='none' is a separately labeled amplitude-preserving ablation.
    Power Doppler is calculated AFTER coherent averaging, as in the paper.
    """
    from .filters import svd_filter
    if normalization not in ("paper","none"):
        raise ValueError("Unknown Song block normalization")
    y=np.zeros_like(iq,dtype=complex);den=np.zeros(iq.shape[1:]);ops=[];diag=[]
    for sl,_ in patches(iq.shape[1:],patch_shape,stride):
        den[sl]+=1
        for ts in chunks(len(iq),window or len(iq)):
            block=iq[(ts,)+sl]
            r=svd_filter(block,centered=centered,selector="song_pdf",selector_options=selector_options)
            mass=r.diagnostics["windows"][0]["retained_sigma_sum"]
            scale=1./mass if normalization=="paper" and mass>0 else 1.
            y[(ts,)+sl]+=r.iq*scale
            ops.append((ts,sl,scale,r.apply_frozen))
            diag.append({"start":ts.start,"stop":ts.stop,"block_start":[s.start for s in sl],
                         "normalization_scale":scale,**r.diagnostics})
    def frozen(other):
        if other.shape!=iq.shape:raise ValueError("Frozen Song shape mismatch")
        out=np.zeros_like(other,dtype=complex)
        for ts,sl,scale,op in ops:out[(ts,)+sl]+=op(other[(ts,)+sl])*scale
        return out/den
    result=FilterResult(y/den,{"blocks":diag,"source":"Song 2017 Eq.1",
        "normalization":normalization,"iq_reconstruction":"coherent equal-weight block average",
        "published_block_normalization":normalization=="paper"},frozen)
    result.power_doppler=np.mean(abs(result.iq)**2,axis=0)
    return result


def lee_aer_filter(iq,domain="temporal",**kwargs):
    """Lee et al. 2026 Algorithms 1-2, centered, regularized ZCA.

    Algorithm 1 divides temporal covariance by P; Algorithm 2 divides globally
    centered spatial covariance by T-1. Absolute lambda retains input units;
    relative regularization and unit-operator scaling are labeled extensions.
    """
    from .filters import aer_filter
    if kwargs.pop("centered",True) is not True or kwargs.pop("alpha",.5)!=.5:
        raise ValueError("Published Lee ZCA uses centered data and alpha=0.5")
    if kwargs.pop("normalization","samples")!="samples":
        raise ValueError("Use aer_filter for a common-denominator ablation")
    result=aer_filter(iq,domain=domain,centered=True,alpha=.5,normalization="samples",
                      covariance_ddof=1 if domain=="spatial" else 0,**kwargs)
    result.diagnostics.update({"source":"Lee et al. 2026 Algorithm 1" if domain=="temporal" else "Lee et al. 2026 Algorithm 2",
        "paper_absolute_regularizer":kwargs.get("epsilon_abs") is not None,
        "relative_regularizer_extension":kwargs.get("epsilon_abs") is None})
    return result
