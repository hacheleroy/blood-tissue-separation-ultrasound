import pytest
import numpy as np
from numpy.testing import assert_allclose
from btbench.filters import svd_filter,aer_filter,local_filter
from btbench.measurements import frozen_component_metrics
from btbench.synthetic import poiseuille
from btbench.axes import motion_packet_filter,multilag_filter,ica_filter
from btbench.core import FilterResult


@pytest.mark.parametrize("method",[svd_filter,aer_filter,local_filter])
@pytest.mark.parametrize("centered",[False,True])
def test_frozen_linear_decomposition(method,centered):
    a=poiseuille(n_frames=32,shape=(6,8))
    r=method(a.iq,centered=centered)
    noise=a.iq-a.blood_true-a.tissue_true
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=2e-12)
    assert_allclose(r.apply_frozen(a.blood_true)+r.apply_frozen(a.tissue_true)+r.apply_frozen(noise),r.iq,atol=2e-12)
    m=frozen_component_metrics(r,a)
    assert m["frozen_component_metrics_available"]


@pytest.mark.parametrize("method",[motion_packet_filter,multilag_filter,ica_filter])
@pytest.mark.parametrize("centered",[False,True])
@pytest.mark.parametrize("rigid",[False,True])
def test_packet_components_replay_same_motion_mean_and_weights(method,centered,rigid):
    a=poiseuille(n_frames=64,shape=(5,7),motion_event=True)
    r=method(a.iq,packet=32,lags=(1,2,4),rank=3,centered=centered,
             rigid=rigid,patch_shape=(4,5),stride=(2,3),q_min=0.)
    noise=a.iq-a.blood_true-a.tissue_true
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=2e-11)
    assert_allclose(r.apply_frozen(a.blood_true)+r.apply_frozen(a.tissue_true)+r.apply_frozen(noise),r.iq,atol=2e-11)
    assert_allclose(r.apply_frozen(2j*a.blood_true),2j*r.apply_frozen(a.blood_true),atol=2e-11)
    assert frozen_component_metrics(r,a)["frozen_component_closure_nrmse"]<1e-12


def test_packet_low_confidence_fallback_is_frozen():
    a=poiseuille(n_frames=64,shape=(5,7),motion_event=True)
    r=motion_packet_filter(a.iq,packet=32,centered=True,q_min=1.,
                          patch_shape=(5,7),stride=(5,7))
    assert all(p["reset"] for p in r.diagnostics["packets"])
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=2e-11)


def test_pooled_mean_restarts_after_a_motion_event():
    a=poiseuille(n_frames=96,shape=(5,7),motion_event=True)
    r=motion_packet_filter(a.iq,packet=32,centered=True,q_min=.9,
                          patch_shape=(5,7),stride=(5,7))
    assert any(p["reset"] for p in r.diagnostics["packets"])
    assert any(not p["reset"] for p in r.diagnostics["packets"])
    assert_allclose(r.apply_frozen(a.iq),r.iq,atol=2e-11)
    n=a.iq-a.blood_true-a.tissue_true
    assert_allclose(r.apply_frozen(a.blood_true)+r.apply_frozen(a.tissue_true)+r.apply_frozen(n),r.iq,atol=2e-11)


def test_retained_power_does_not_hide_phase_or_temporal_distortion():
    a=poiseuille(n_frames=32,shape=(4,6))
    identity=FilterResult(a.iq,apply_frozen=lambda x:x)
    rotated=FilterResult(1j*a.iq,apply_frozen=lambda x:1j*x)
    m=frozen_component_metrics(identity,a);r=frozen_component_metrics(rotated,a)
    assert_allclose(m["frozen_blood_power_retention"],r["frozen_blood_power_retention"])
    assert m["frozen_blood_iq_nrmse"]==0
    assert_allclose(r["frozen_blood_iq_nrmse"],np.sqrt(2))
    ramp=np.exp(2j*np.pi*np.arange(32)/4)[:,None,None]
    shifted=FilterResult(a.iq*ramp,apply_frozen=lambda x:x*ramp)
    s=frozen_component_metrics(shifted,a)
    assert s["frozen_blood_spectral_tv"]>.5
    assert s["frozen_blood_lag1_covariance_nrmse"]>1


def test_empty_output_has_no_valid_spectral_shape_or_leakage_ratio():
    a=poiseuille(n_frames=32,shape=(4,6))
    r=FilterResult(np.zeros_like(a.iq),apply_frozen=np.zeros_like)
    m=frozen_component_metrics(r,a)
    assert m["frozen_blood_iq_nrmse"]==1
    assert m["frozen_blood_power_retention"]==0
    assert m["frozen_blood_output_zero"]
    assert not m["frozen_blood_spectral_shape_valid"]
    assert m["frozen_tissue_to_blood_power_ratio"] is None
    assert m["frozen_component_closure_nrmse"]==0
