import numpy as np
from numpy.testing import assert_allclose
from btbench.cca import cca_denoise


def test_common_mode_removed_original_gain_preserved():
    rng=np.random.default_rng(42)
    a=rng.normal(size=(200,4));a-=a.mean(0)
    q,_=np.linalg.qr(a)
    noise=q[:,0]
    signal=q[:,1]
    f=np.stack([4*signal+30*noise,2*signal-20*noise],axis=1)+100
    n=np.stack([3*noise,4*q[:,2]],axis=1)
    clean,info=cca_denoise(f,n)
    assert_allclose(clean,np.stack([4*signal,2*signal],axis=1)+100,atol=1e-12)
    assert info["removed_components"]==1
