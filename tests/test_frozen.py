import pytest
import numpy as np
from numpy.testing import assert_allclose
from btbench.filters import svd_filter,aer_filter,local_filter
from btbench.measurements import frozen_component_metrics
from btbench.synthetic import poiseuille


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
