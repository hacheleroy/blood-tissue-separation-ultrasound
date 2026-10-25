"""Numerical bridge from the Example 22 geometry to a measured vessel axis."""

from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))

from esdiva_phase10_support import phase10_support  # noqa: E402
from phantom_phase10_port import structural_support, calibration_support  # noqa: E402


def test_structural_axis_and_calibration_replay_simulated_geometry():
    z=np.arange(14.4,17.6,.1)[:32]
    x=np.arange(-8.,8.)
    axis=np.tan(np.deg2rad(60.))*(z-16.)
    wide,core=structural_support(z,x,axis)
    rng=np.random.default_rng(654)
    tissue=rng.normal(size=(32,16))+1j*rng.normal(size=(32,16))
    cal=np.broadcast_to(tissue,(16,3,32,16)).copy()
    cal+=.005*(rng.normal(size=cal.shape)+1j*rng.normal(size=cal.shape))
    original_wide,original_core,gate=phase10_support(cal,z,x,np.arange(12,20))
    assert np.array_equal(wide,original_wide)
    assert np.array_equal(core,original_core)
    assert calibration_support(cal,wide,core)==gate


def test_structural_support_rejects_unmeasured_or_insufficient_axis():
    z=np.arange(14.4,17.6,.1)[:32]
    x=np.arange(-8.,8.)
    with pytest.raises(ValueError):
        structural_support(z,x,np.full(32,np.nan))
    with pytest.raises(ValueError):
        structural_support(z,np.linspace(-1.,1.,16),np.zeros(32))
