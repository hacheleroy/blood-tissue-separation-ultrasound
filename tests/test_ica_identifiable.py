import numpy as np
from btbench.axes import complex_fastica_rotation
from btbench.axes import ica_filter
from btbench.synthetic import poiseuille
from numpy.testing import assert_allclose


def test_complex_ica_recovers_nongaussian_sources():
    rng=np.random.default_rng(11)
    s=rng.exponential(size=(3,20000))*np.exp(1j*rng.uniform(0,2*np.pi,size=(3,20000)))
    s-=s.mean(1,keepdims=True)
    mixing=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    x=mixing@s
    ev,u=np.linalg.eigh(x@x.conj().T/x.shape[1])
    b=u/np.sqrt(ev)
    z=b.conj().T@x
    q,info=complex_fastica_rotation(z)
    y=q.conj().T@z
    corr=abs(y@s.conj().T)/np.sqrt(np.sum(abs(y)**2,axis=1)[:,None]*np.sum(abs(s)**2,axis=1)[None,:])
    assert info["ica_converged"]
    assert np.all(np.max(corr,axis=1)>.99)
    assert len(set(np.argmax(corr,axis=1)))==3


def test_ica_nonconvergence_fallback_is_visible_and_linear():
    a=poiseuille(n_frames=64,shape=(5,7))
    r=ica_filter(a.iq,packet=32,rank=4,lags=(1,2,4),patch_shape=(5,7),stride=(5,7),
                 ica_options={"max_iter":1},q_min=0.)
    packets=r.diagnostics["packets"]
    assert all(not p["ica_converged"] for p in packets)
    assert all(p["fallback"]=="local_svd_ica_nonconvergence" for p in packets)
    assert all(p["ica_whitening_error"]<1e-10 for p in packets)
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=1e-11)
