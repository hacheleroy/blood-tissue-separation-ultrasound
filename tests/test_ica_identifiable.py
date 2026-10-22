import numpy as np
from btbench.axes import complex_fastica_rotation


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
