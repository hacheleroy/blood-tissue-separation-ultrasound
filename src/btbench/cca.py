"""Optional downstream fUS denoising, kept separate from IQ clutter filtering."""
import numpy as np
from scipy import linalg


def cca_denoise(functional,noise_reference,min_common_correlation=0.5,rank=None):
    """Faure-style whitening + concatenated SVD + common-mode regression.

    Inputs are real power-Doppler time courses (T,Pf)/(T,Pn), sampled together.
    Noise ROI must be independently identified. Subspace overlap defines the
    selection rule here: s_concat^2=1+rho. This explicit rule is a benchmark
    choice, not a claimed replication of the paper's threshold-selection study.
    Filtering is applied to ORIGINAL centered amplitudes, then mean restored.
    Anatomical/functional signal shared with reference can be removed.
    """
    f,n=np.asarray(functional),np.asarray(noise_reference)
    if f.ndim!=2 or n.ndim!=2 or len(f)!=len(n) or np.iscomplexobj(f) or np.iscomplexobj(n):
        raise ValueError("Need real, time-aligned power-Doppler matrices")
    if not -1<min_common_correlation<1:
        raise ValueError("Invalid common-correlation threshold")
    fm=f.mean(0);nm=n.mean(0)
    def whiten(x):
        u,s,vh=linalg.svd(x,full_matrices=False)
        r=int(np.sum(s>max(s[0],1e-30)*1e-10))
        if rank is not None:r=min(r,rank)
        return u[:,:r]@vh[:r]
    fw,nw=whiten(f-fm),whiten(n-nm)
    u,s,_=linalg.svd(np.concatenate((fw,nw),axis=1),full_matrices=False)
    common=s*s>1+min_common_correlation
    basis=u[:,common]
    clean=f-basis@(basis.T@(f-fm))
    return clean,{"removed_components":int(common.sum()),"common_correlation_threshold":min_common_correlation,
                  "cca_stage":"power-Doppler movie, after IQ clutter filtering",
                  "centered_amplitudes_restored":True,"negative_values_possible":True}
