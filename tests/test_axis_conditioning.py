import numpy as np
from numpy.testing import assert_allclose
from btbench.axes import _separate,multilag_filter,motion_packet_filter


def test_uncentered_axis1_is_invariant_to_global_phase_compensation():
    rng=np.random.default_rng(22);x=rng.normal(size=(96,4,5))+1j*rng.normal(size=(96,4,5))
    x+=20*(rng.normal(size=(1,4,5))+1j*rng.normal(size=(1,4,5)))
    x*=np.exp(1j*(.13*np.arange(96)+.0007*np.arange(96)**2))[:,None,None]
    settings=dict(packet=32,patch_shape=(4,5),stride=(4,5),q_min=0.,rho=.9)
    a=motion_packet_filter(x,phase_compensation=True,**settings)
    b=motion_packet_filter(x,phase_compensation=False,**settings)
    assert_allclose(a.iq,b.iq,atol=2e-12)
    assert a.diagnostics['uncentered_axis1_global_phase_invariant']


def test_sobi_gain_cap_bounds_ill_conditioned_mixture_without_hiding_original_gain():
    n=512;t=np.arange(n)
    sources=np.array([np.exp(.003j*t),np.exp(.43j*t)])
    x=np.array([[1.,1.],[0.,.015]])@sources
    lags=(1,2,4,8)
    covs={lag:x[:,lag:]@x[:,:n-lag].conj().T/(n-lag) for lag in (0,)+lags}
    _,raw,a=_separate(x,covs,lags,rank=2,epsilon_rel=1e-10)
    y,d,b=_separate(x,covs,lags,rank=2,epsilon_rel=1e-10,operator_gain_cap=1.)
    assert raw['operator_gain_max']>2
    assert d['unconstrained_operator_gain_max']>2
    assert np.linalg.norm(b,2)<=1+2e-13
    assert np.linalg.norm(y)<=np.linalg.norm(x)*(1+1e-13)
    assert d['operator_projection']=='spectral-norm ball'


def test_regularized_capped_packet_operator_replays_components_linearly():
    rng=np.random.default_rng(13)
    blood=rng.normal(size=(96,4,5))+1j*rng.normal(size=(96,4,5))
    tissue=50*(rng.normal(size=(1,4,5))+1j*rng.normal(size=(1,4,5)))*np.exp(.01j*np.arange(96))[:,None,None]
    x=blood+tissue
    for centered in (False,True):
        r=multilag_filter(x,packet=32,patch_shape=(4,5),stride=(4,5),lags=(1,2,4),
            q_min=0.,rank=4,epsilon_rel=1e-3,whitening_regularization_rel=1e-3,
            operator_gain_cap=1.,centered=centered)
        assert_allclose(r.apply_frozen(blood)+r.apply_frozen(tissue),r.iq,atol=1e-12)
        assert all(p['operator_gain_max']<=1+1e-12 for p in r.diagnostics['packets'])
