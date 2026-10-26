import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import ndimage

from btbench.motion_first import (estimate_bmode_motion, registration_order_filter,
                                  motion_windows)


def test_envelope_motion_and_phase_estimation_uses_no_truth():
    rng = np.random.default_rng(78)
    texture = ndimage.gaussian_filter(rng.normal(size=(32, 36)), .8)
    texture = texture + 1j * ndimage.gaussian_filter(rng.normal(size=(32, 36)), .8)
    shifts = np.array([(0., 0.), (.63, -.42), (-.51, .35)])
    iq = np.stack([ndimage.shift(texture.real, d, order=1, mode="constant")
                   + 1j * ndimage.shift(texture.imag, d, order=1, mode="constant")
                   for d in shifts]) * np.exp(1j*np.array([0., .3, -.2]))[:, None, None]
    fit = estimate_bmode_motion(iq)
    got = fit.forward[:, :, 16, 18] - fit.forward[0, :, 16, 18]
    # Envelope formation after interpolation biases small subpixel shifts.
    assert_allclose(got, shifts, atol=.26)
    assert_allclose(fit.phase[:, 0], [0., .3, -.2], atol=.12)
    assert fit.diagnostics["failed_frames"] == 0


def test_order_filter_replays_frozen_motion_and_linear_source_attribution():
    rng = np.random.default_rng(8)
    tissue = ndimage.gaussian_filter(rng.normal(size=(26, 30)), .7)
    base = np.stack([ndimage.shift(tissue, (.4*np.sin(i/3), .3*np.cos(i/3)),
                      order=1, mode="constant")*np.exp(.05j*i) for i in range(20)])
    blood = np.zeros_like(base, complex)
    blood[:, 9:18, 12:18] = .03*np.exp(1j*np.arange(20)[:, None, None]*.3)
    mixture = base+blood
    motion = estimate_bmode_motion(mixture)
    for order in ("before", "after"):
        for phase in (False, True):
            r = registration_order_filter(mixture, motion, order, phase=phase)
            assert_allclose(r.apply_frozen(mixture), r.iq, atol=1e-12)
            assert_allclose(r.apply_frozen(base)+r.apply_frozen(blood), r.iq, atol=1e-12)
    with pytest.raises(ValueError):motion.align(mixture[:3])


def test_motion_limited_windows_and_fixed_comparator_share_support():
    rng = np.random.default_rng(7)
    iq = rng.normal(size=(48, 16, 20)) + 1j*rng.normal(size=(48, 16, 20))
    motion = estimate_bmode_motion(iq*.01 + 5*np.ones_like(iq))
    # Inject a known estimated path without exposing this motion to a filter.
    motion.forward[:, 0] += np.arange(48)[:, None, None]*.09
    adaptive = motion_windows(iq, motion, 1000, 5e6, .0002, .0002,
                              window=32, min_window=12, hop=6)
    fixed = motion_windows(iq, motion, 1000, 5e6, .0002, .0002,
                           window=32, min_window=12, hop=6, adaptive=False)
    assert min(w["stop"]-w["start"] for w in adaptive.diagnostics["windows"]) < 32
    assert all(w["stop"]-w["start"] == min(32, 48-w["start"])
               for w in fixed.diagnostics["windows"])
    assert_allclose(adaptive.apply_frozen(iq), adaptive.iq, atol=1e-12)
