import numpy as np
import pytest
from numpy.testing import assert_allclose
from btbench.filters import aer_filter,svd_filter,sliding_filter
from btbench.stap import lcmv_weights
from btbench.synthetic import poiseuille
from btbench.benchmark import oracle
from btbench.axes import compensate,restore,shift_complex


def test_windowed_aer_and_frozen_operator():
    a=poiseuille(n_frames=64,shape=(4,6))
    r=aer_filter(a.iq,window=32,centered=True)
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=1e-12)
    for k in (0,32):
        ref=aer_filter(a.iq[k:k+32],centered=True).iq
        assert_allclose(r.iq[k:k+32],ref)


def test_data_only_noise_floor_has_finite_weights():
    a=poiseuille(n_frames=64,shape=(4,6))
    r=aer_filter(a.iq,epsilon_rule="noise_floor",centered=True)
    assert np.isfinite(r.iq).all()
    assert r.diagnostics["epsilon_absolute"]>0


def test_overlapping_stap_constraints_rejected():
    s=np.ones((8,1),complex)
    with pytest.raises(ValueError):lcmv_weights(s,s,np.eye(8))


def test_exhaustive_lower_cutoff_grid():
    a=poiseuille(n_frames=20,shape=(3,4))
    _,d=oracle(a,{"method":"svd","params":{}},{"discard":"all"})
    assert len(d["curve"])==12


def test_svd_upper_cutoff_is_real_control():
    rng=np.random.default_rng(3)
    x=rng.normal(size=(30,4,6))+1j*rng.normal(size=(30,4,6))
    y=svd_filter(x,discard=1,upper=3).iq.reshape(30,-1)
    assert np.linalg.matrix_rank(y)==2


def test_integer_rigid_shift_and_boundary_validity():
    rng=np.random.default_rng(22)
    ref=rng.normal(size=(12,14))+1j*rng.normal(size=(12,14))
    iq=np.stack([ref,shift_complex(ref,(1,1)),shift_complex(ref,(1,1))])
    aligned,motion=compensate(iq,np.ones(ref.shape,bool),rigid=True)
    assert_allclose(motion["shifts"][1:],np.array([[-1,-1],[-1,-1]]))
    restored,valid=restore(aligned,motion)
    assert not valid.all()
    assert_allclose(restored[valid],iq[valid],atol=1e-12)


def test_sliding_identity_and_frozen_sum():
    a=poiseuille(n_frames=100,shape=(4,6))
    r=sliding_filter(a.iq,window=40,hop=17,discard=0)
    assert_allclose(r.iq,a.iq,atol=1e-12)
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=1e-12)


def test_sliding_aer_equation_shape():
    a=poiseuille(n_frames=100,shape=(4,6))
    r=sliding_filter(a.iq,method="aer",window=40,hop=20)
    assert np.isfinite(r.iq).all()
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=1e-12)


def test_json_shape_can_be_a_list():
    a=poiseuille(n_frames=20,shape=[4,6])
    assert a.iq.shape==(20,4,6)
