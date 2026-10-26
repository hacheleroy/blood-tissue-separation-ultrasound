"""Exact fixed-projector trade-off calculation, including complex IQ."""
import sys
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"examples"))
from axis1_oracle_limit_study import minimum_gain_for_tissue_limit


def test_minimum_gain_meets_tissue_limit_and_minimizes_blood_error_on_ray():
    background=np.array([1+1j,2-1j])
    removed=background.copy()
    gain=minimum_gain_for_tissue_limit(background,removed,
        float(np.linalg.norm(background)**2),.1)
    assert_allclose(gain,1-np.sqrt(.1),atol=1e-12)
    assert_allclose(np.linalg.norm(background-gain*removed)**2 /
                    np.linalg.norm(background)**2,.1,atol=1e-12)
    assert minimum_gain_for_tissue_limit(background,np.zeros_like(background),
        float(np.linalg.norm(background)**2),.1) is None
    assert minimum_gain_for_tissue_limit(background,background,
        float(np.linalg.norm(background)**2),.5) < gain
