import numpy as np
import pytest
from numpy.testing import assert_allclose

from btbench.motion_first import estimate_bmode_motion
from btbench.native_subtraction import native_tissue_subtraction
from btbench.stable_tissue import bright_reference_subtraction, smooth_motion_path


def test_bright_fit_has_linear_frozen_replay_and_full_field_identity():
    rng = np.random.default_rng(92)
    ref = rng.normal(size=(16, 24)) + 1j*rng.normal(size=(16, 24))
    x = np.stack([ref*np.exp(0.035j*t) + .15*(rng.normal(size=ref.shape) +
                       1j*rng.normal(size=ref.shape)) for t in range(16)])
    motion = smooth_motion_path(estimate_bmode_motion(x))
    assert_allclose(motion.forward[0], np.indices(ref.shape), atol=1e-12)
    fitted = bright_reference_subtraction(x, motion)
    assert fitted.diagnostics["training_pixels"] == 96
    assert_allclose(fitted.iq, fitted.apply_frozen(x), atol=1e-12)
    assert_allclose(fitted.apply_frozen(x[:]*.3)+fitted.apply_frozen(x[:]*.7),
                    fitted.iq, atol=1e-12)
    full = bright_reference_subtraction(x, motion, fraction=1)
    old = native_tissue_subtraction(x, motion)
    assert_allclose(full.iq, old.iq, atol=1e-11)
    with pytest.raises(ValueError):fitted.apply_frozen(x[:-1])
    with pytest.raises(ValueError):bright_reference_subtraction(x, motion, fraction=0)
    with pytest.raises(ValueError):smooth_motion_path(motion, sigma_frames=0)
