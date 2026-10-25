import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'examples'))
from esdiva_phase7_reference_transport import transported_reference  # noqa: E402


def test_reference_transport_marks_out_of_bounds_without_wrap():
    reference=np.zeros((32,16),complex)
    reference[15,0]=1+1j
    moved,valid=transported_reference(reference,[0.,-.5])
    assert not valid[15,0]
    assert valid[15,15]
    assert moved[15,15]==0


def test_reference_transport_zero_displacement_is_identity():
    rng=np.random.default_rng(42)
    reference=rng.normal(size=(32,16))+1j*rng.normal(size=(32,16))
    moved,valid=transported_reference(reference,[0.,0.])
    assert valid.all()
    np.testing.assert_allclose(moved,reference,atol=0,rtol=0)
