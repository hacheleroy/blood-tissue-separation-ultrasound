"""Source-blind angle/aperture filters and explicitly source-aware diagnostics.

Layout is (slow_time, angle, spatial_pixel[, receive_subaperture]).  The
three-mode filters below reproduce the *factorizations* in Ozgun--Byram and
Jiang et al.; the deliberately fixed rank/cutoff rules are study adaptations,
not claims to reproduce the authors' trained/manual cutoff choices.
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class Frozen:
    apply: object
    settings: dict


def svd_unfolding(x, discard=1, layout="joint"):
    """SVD of pixels x slow time (compound) or pixels*view x slow time.

    The grouping of angles into spatial samples is Jiang's GA-SVD structure;
    the identical unfolding with receive channels is Ozgun's SVD-a control.
    """
    x = np.asarray(x)
    if x.ndim not in (2, 3) or not 0 <= discard < x.shape[0]:
        raise ValueError("Need (time,pixel[,view]) and a nonnegative discard")
    if layout == "compound" and x.ndim == 3:
        x = x.mean(axis=-1)
    elif layout != "joint" and layout != "compound":
        raise ValueError(layout)
    mat = x.reshape(len(x), -1)
    u, _, _ = np.linalg.svd(mat, full_matrices=False)
    if discard >= u.shape[1]:
        raise ValueError("Rank exceeds available slow-time samples")
    v = u[:,discard:]
    return Frozen(lambda other: (v @ (v.conj().T @ other.reshape(len(other), -1)))
                  .reshape(other.shape), {"unfolding": layout, "discard": discard})


def factor_tensor(x):
    """Complex Tucker/HOSVD on (time,pixel,angle-or-aperture)."""
    x = np.asarray(x)
    if x.ndim != 3 or min(x.shape) < 2:
        raise ValueError("HOSVD expects (time,pixel,view) with each axis >=2")
    # Columns are ordered by descending singular value in each unfolding.
    u = [np.linalg.svd(np.moveaxis(x, mode, 0).reshape(x.shape[mode], -1),
                       full_matrices=False)[0] for mode in range(3)]
    return u


def tensor_gate(x, keep, u=None):
    """Project a tensor through a frozen HOSVD core mask (complex linear)."""
    if u is None:
        u = factor_tensor(x)
    dims = np.shape(x)
    if np.shape(keep) != dims:
        raise ValueError("Core mask has wrong shape")
    def apply(other):
        if np.shape(other) != dims:
            raise ValueError("Frozen HOSVD geometry mismatch")
        core = np.einsum('ti,pj,ak,tpa->ijk',u[0].conj(),u[1].conj(),
                         u[2].conj(),other,optimize=True)
        return np.einsum('ti,pj,ak,ijk->tpa',u[0],u[1],u[2],
                         core*keep,optimize=True)
    return Frozen(apply, {"core_kept":int(np.count_nonzero(keep)),
                          "core_total":int(np.size(keep))})


def paper_tensor(x, temporal_discard=1, family="ozgun", noise_cut=True):
    """Paper-structured HOSVD with frozen, explicit surrogate cutoff rules.

    Ozgun: temporal low-mode rejection + aperture phase plateau + spatial
    Nakagami noise gate. Jiang: temporal low-mode rejection + angular
    singular-value knee; its spatial PCC classifier is omitted in this
    conservative adaptation. These fixed rules are not the original papers'
    full, sometimes manually tuned, classifiers.
    """
    if family not in ("ozgun", "jiang") or temporal_discard not in (0, 1, 2):
        raise ValueError("Invalid paper family/rank")
    u = factor_tensor(x)
    t,p,a = np.shape(x)
    amps2=np.abs(u[1])**2
    m2=amps2.mean(axis=0)
    m4=(amps2**2).mean(axis=0)
    inv=(m2*m2)/np.maximum(m4-m2*m2,1e-30)
    stat=np.sqrt(p)*(inv-1.)
    # Ozgun's spatial classifier: 5-point smoothing, first noise-like mode
    # within three Gaussian standard deviations (sigma=sqrt(4) for m=1).
    smooth=np.convolve(stat,np.ones(5)/5,mode='same')
    hits=np.flatnonzero(np.abs(smooth)<6.)
    spatial_keep=int(hits[0]) if family=='ozgun' and noise_cut and len(hits) else p
    spatial_keep=max(1,spatial_keep)
    rng=np.random.default_rng(7731+a)
    if family=='ozgun':
        # Paper phase-deviation rule uses a Monte Carlo confidence interval
        # around the random-phase plateau pi/sqrt(3); the small-group pilot
        # uses 10k fixed samples and no per-scene tuning.
        phase=np.angle(u[2])
        deviation=np.std(phase,axis=0)
        null=np.std(rng.uniform(-np.pi,np.pi,(10000,a)),axis=1)
        floor=float(np.quantile(null,.025))
        below=np.flatnonzero(deviation<floor)
        angular_keep=int(below[-1]+1) if noise_cut and len(below) else a
    else:
        s=np.linalg.svd(np.moveaxis(x,2,0).reshape(a,-1),compute_uv=False)
        ratios=s[:-1]/np.maximum(s[1:],1e-30)
        angular_keep=int(np.argmax(ratios)+1) if noise_cut and len(ratios) else a
        floor=None
    angular_keep=max(1,min(a,angular_keep))
    keep=np.ones((t,p,a),bool)
    keep[:temporal_discard]=False
    if noise_cut:
        keep[:,spatial_keep:]=False
        keep[:,:,angular_keep:]=False
    op=tensor_gate(x,keep,u)
    op.settings.update({"paper_family":family,"temporal_discard":temporal_discard,
                        "spatial_keep":spatial_keep,"angular_or_aperture_keep":angular_keep,
                        "phase_deviation_null_025":floor,
                        "jiang_spatial_pcc_classification_omitted":family=='jiang',
                        "classifier_is_study_adaptation":True})
    return op


def intersection_filter(x, temporal_rank=1, angular_rank=1):
    """New simple candidate: remove only the joint low-time/low-angle core.

    Unlike the papers' rectangular OR gates, this tests whether a blood
    component with the same temporal factor but different angular factor
    can be retained. Fits observed mixture only; no source/velocity/ROI.
    """
    u=factor_tensor(x)
    t,p,a=x.shape
    keep=np.ones((t,p,a),bool)
    keep[:temporal_rank,:,:angular_rank]=False
    op=tensor_gate(x,keep,u)
    op.settings.update({"temporal_rank":temporal_rank,"angular_rank":angular_rank,
                        "candidate":"intersection of low-time and low-view HOSVD factors"})
    return op


def oracle_joint(tissue, rank=1):
    """Evaluation-only projector fitted to *isolated tissue* in joint coordinates.

    Input shape (time,pixel[,angle[,aperture]]). All non-time axes after pixel
    are retained as one joint snapshot coordinate, rather than collapsed into
    the spatial side of a standard Casorati matrix.
    """
    t,p=tissue.shape[:2]
    mat=np.moveaxis(tissue,1,0).reshape(p,-1)
    _,s,vh=np.linalg.svd(mat,full_matrices=False)
    v=vh[:rank].conj().T
    def removed(x):
        other=np.moveaxis(x,1,0).reshape(p,-1)
        return np.moveaxis(((other@v)@v.conj().T).reshape((p,t,*tissue.shape[2:])),0,1)
    return removed,float(np.sum(s[:rank]**2)/np.sum(s*s))


def minimum_gain(background, removed, tissue, limit=.1):
    """Smallest subtraction strength in [0,1] meeting the background target."""
    a=float(np.vdot(removed,removed).real)
    b=float(np.vdot(background,removed).real)
    c=float(np.vdot(background,background).real-limit*np.vdot(tissue,tissue).real)
    if c <= 0:return 0.
    d=b*b-a*c
    if a<=0 or d<0:return None
    lo=(b-np.sqrt(d))/a
    hi=(b+np.sqrt(d))/a
    gain=max(0.,float(lo))
    return gain if gain<=min(1.,float(hi))+1e-10 else None
