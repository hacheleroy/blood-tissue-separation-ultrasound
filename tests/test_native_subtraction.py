import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import ndimage

from btbench.motion_first import estimate_bmode_motion, registration_order_filter
from btbench.native_subtraction import native_tissue_subtraction


def test_native_subtraction_is_original_minus_frozen_tissue_estimate():
    rng = np.random.default_rng(59)
    texture = ndimage.gaussian_filter(rng.normal(size=(24, 30)), .9)
    tissue = np.stack([ndimage.shift(texture, (.6*np.sin(t/5), .4*np.cos(t/5)),
                        order=1, mode="constant") * np.exp(.08j*t) for t in range(24)])
    blood = np.zeros_like(tissue)
    blood[:, 7:16, 14:19] = .04*np.exp(1j*np.arange(24)[:, None, None]*.31)
    mixture = tissue+blood
    motion = estimate_bmode_motion(mixture)
    native = native_tissue_subtraction(mixture, motion)
    returned = registration_order_filter(mixture, motion, order="after", phase=True)
    assert_allclose(native.apply_frozen(mixture), native.iq, atol=1e-12)
    assert_allclose(native.apply_frozen(tissue)+native.apply_frozen(blood), native.iq, atol=1e-12)
    assert_allclose(native.iq-returned.iq,
                    mixture-motion.restore(motion.align(mixture)), atol=1e-12)
    assert_allclose(native.apply_frozen(blood)-returned.apply_frozen(blood),
                    blood-motion.restore(motion.align(blood)), atol=1e-12)
    with pytest.raises(ValueError):native.apply_frozen(mixture[:3])
    with pytest.raises(ValueError):native_tissue_subtraction(mixture,motion,discard=0)
