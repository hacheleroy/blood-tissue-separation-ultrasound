"""Source-blind complex-IQ low-rank/sparse comparators.

All three fits operate on a pixel-by-time Casorati matrix. The full-grid
Cauchy arm implements Sui et al.'s reweighted ADMM core, without their
randomized spatial partition and parallel reassembly. The batch IRLS arm
implements a complex extension of the factorized batch IRLS updates in
Rontogiannis et al. (2020); it is a reconstruction of the family used by
Chen et al. (2025), not verified parity with their static implementation.
"""
import numpy as np

from .core import casorati, frames


def _input(iq):
    a = np.asarray(iq)
    if a.ndim != 3 or a.shape[0] < 3 or not np.iscomplexobj(a) or not np.isfinite(a).all():
        raise ValueError("Expected finite complex IQ (time,z,x), at least three frames")
    x = casorati(a)
    scale = float(np.sqrt(np.mean(np.abs(x)**2)))
    return a, x / max(scale, 1e-30), scale


def _soft(x, threshold):
    return x * np.maximum(1 - threshold / np.maximum(np.abs(x), 1e-30), 0)


def _result(a, s, scale, method, iterations, converged, relative_change, residual, extra):
    output = frames(scale * s, a.shape)
    if not np.isfinite(output).all():
        raise FloatingPointError(f"{method} generated nonfinite IQ")
    return output, dict(method=method, iterations=iterations, converged=bool(converged),
                        relative_change=float(relative_change),
                        relative_residual=float(residual), normalization_rms=scale, **extra)


def stable_rpca(iq, tau, sparse_lambda=None, max_iter=120, tol=1e-5):
    """Minimize 0.5||X-L-S||_F² + tau||L||_* + tau*lambda||S||_1.

    Alternating exact proximal block updates. The residual is the noise term.
    A single normalized amplitude scale is calculated from the observed IQ.
    """
    a, x, scale = _input(iq)
    lam = 1 / np.sqrt(max(x.shape)) if sparse_lambda is None else float(sparse_lambda)
    if tau <= 0 or lam <= 0 or max_iter < 1 or tol <= 0:
        raise ValueError("Invalid stable RPCA parameters")
    l = np.zeros_like(x); s = np.zeros_like(x); change = np.inf
    for k in range(1, max_iter + 1):
        u, singular, vh = np.linalg.svd(x-s, full_matrices=False)
        shrunk = np.maximum(singular-tau, 0)
        next_l = (u*shrunk)@vh
        next_s = _soft(x-next_l, tau*lam)
        change = np.linalg.norm(next_l-l)+np.linalg.norm(next_s-s)
        change /= max(np.linalg.norm(next_l)+np.linalg.norm(next_s), 1e-30)
        l, s = next_l, next_s
        if change < tol: break
    residual = np.linalg.norm(x-l-s)/max(np.linalg.norm(x), 1e-30)
    return _result(a,s,scale,"stable_rpca",k,change<tol,change,residual,
                   dict(tau=float(tau),sparse_lambda=lam,noise_model="quadratic residual"))


def cauchy_rpca_core(iq, sparse_lambda, gamma=0.1, max_iter=120, tol=1e-5):
    """Full-grid reweighted Cauchy RPCA ADMM, complex soft thresholds.

    The Cauchy penalty is gamma²*log(1+|z|²/gamma²), yielding the weights
    2*|z|*gamma²/(gamma²+|z|²) written in Sui et al.'s algorithm. Both the
    singular values and sparse entries are reweighted. The paper's default
    gamma=0.01 is amplitude-unit dependent; IQ is RMS-normalized here and
    gamma is explicitly selected in those units.
    """
    a,x,scale=_input(iq)
    if sparse_lambda<=0 or gamma<=0 or max_iter<1 or tol<=0:
        raise ValueError("Invalid Cauchy RPCA parameters")
    mu=x.size/(4*max(np.sum(np.abs(x)),1e-30))
    l=np.zeros_like(x);s=np.zeros_like(x);dual=np.zeros_like(x)
    weights_s=np.ones(x.shape);weights_l=np.ones(min(x.shape))
    change=np.inf
    for k in range(1,max_iter+1):
        u,singular,vh=np.linalg.svd(x-s+dual/mu,full_matrices=False)
        shrunk=np.maximum(singular-weights_l/mu,0)
        l=(u*shrunk)@vh
        weights_l=2*shrunk*gamma**2/(shrunk**2+gamma**2)
        proposal=x-l+dual/mu
        next_s=_soft(proposal,sparse_lambda*weights_s/mu)
        weights_s=2*np.abs(next_s)*gamma**2/(np.abs(next_s)**2+gamma**2)
        change=np.linalg.norm(next_s-s)/max(np.linalg.norm(next_s),1e-30)
        s=next_s
        residual_matrix=x-l-s
        dual+=mu*residual_matrix
        residual=np.linalg.norm(residual_matrix)/max(np.linalg.norm(x),1e-30)
        if residual<tol and change<tol:break
    return _result(a,s,scale,"cauchy_rpca_core",k,residual<tol and change<tol,
                   change,residual,dict(sparse_lambda=float(sparse_lambda),
                   gamma=float(gamma),admm_mu=float(mu),
                   randomized_downsampling=False,noise_model="absorbed into low-rank/sparse"))


def batch_irls_rpca(iq, low_rank_lambda, sparse_mu, rank=12,
                    epsilon=0.02, max_iter=120, tol=1e-5):
    """Complex batch factorized IRLS with a separate quadratic noise residual.

    U,V and S follow Rontogiannis et al.'s batch closed-form updates with
    Hermitian transposes. The unverified Chen-specific constraints/initialization
    are deliberately absent.
    """
    a,x,scale=_input(iq)
    if min(low_rank_lambda,sparse_mu,epsilon)<=0 or not 1<=rank<=min(x.shape) or max_iter<1 or tol<=0:
        raise ValueError("Invalid IRLS parameters")
    u0,values,vh=np.linalg.svd(x,full_matrices=False)
    root=np.sqrt(values[:rank])
    u=u0[:,:rank]*root
    v=vh[:rank].conj().T*root
    s=np.zeros_like(x)
    eye=np.eye(rank)
    change=np.inf
    for k in range(1,max_iter+1):
        w=1/np.sqrt(np.abs(s)**2+epsilon**2)
        next_s=(x-u@v.conj().T)/(1+sparse_mu*w)
        d=1/np.sqrt(np.sum(np.abs(u)**2,axis=0)+np.sum(np.abs(v)**2,axis=0)+epsilon**2)
        gram_u=u.conj().T@u+low_rank_lambda*np.diag(d)
        next_v=np.linalg.solve(gram_u.T,((x-next_s).conj().T@u).T).T
        gram_v=next_v.conj().T@next_v+low_rank_lambda*np.diag(d)
        next_u=np.linalg.solve(gram_v.T,((x-next_s)@next_v).T).T
        old_l=u@v.conj().T
        new_l=next_u@next_v.conj().T
        change=(np.linalg.norm(new_l-old_l)+np.linalg.norm(next_s-s))
        change/=max(np.linalg.norm(new_l)+np.linalg.norm(next_s),1e-30)
        u,v,s=next_u,next_v,next_s
        if change<tol:break
    residual=np.linalg.norm(x-u@v.conj().T-s)/max(np.linalg.norm(x),1e-30)
    return _result(a,s,scale,"batch_irls_reconstruction",k,change<tol,change,residual,
                   dict(low_rank_lambda=float(low_rank_lambda),sparse_mu=float(sparse_mu),
                        factor_rank=int(rank),epsilon=float(epsilon),
                        chen_paper_parity_verified=False,noise_model="quadratic residual"))
