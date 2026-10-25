"""Check the new calibration guard and reversible, open-boundary operator."""

from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"examples"))

from esdiva_phase5_local_phase import tissue_mask  # noqa: E402
from esdiva_phase8_selector import choose_wide_guard, open_unitary_shift  # noqa: E402


def test_guard_selection_uses_training_only_and_checks_validation():
    rng = np.random.default_rng(81)
    z = 14.4+.1*np.arange(32)
    x = np.arange(-8, 8, dtype=float)
    base = rng.normal(size=(32, 16))+1j*rng.normal(size=(32, 16))
    cal = np.broadcast_to(base, (16, 3, 32, 16)).copy()
    bad = tissue_mask(z, x, 5.)
    cal[:, 1, bad] = .6*cal[:, 1, bad]+.8*(
        rng.normal(size=(16, bad.sum()))+1j*rng.normal(size=(16, bad.sum())))
    _, gate = choose_wide_guard(cal, z, x)
    assert gate["selected"]["guard_mm"] == 4.
    assert gate["selection_uses_train_only"]
    cal[8:, 1] = rng.normal(size=cal[8:, 1].shape)+1j*rng.normal(
        size=cal[8:, 1].shape)
    _, failed = choose_wide_guard(cal, z, x)
    assert failed["selected"]["guard_mm"] == gate["selected"]["guard_mm"]
    assert not failed["pass"]


def test_nonperiodic_shift_is_reversible_without_edge_wrap():
    rng = np.random.default_rng(82)
    frame = rng.normal(size=(8, 16))+1j*rng.normal(size=(8, 16))
    displacement = [.7, -.4]
    moved = open_unitary_shift(frame, displacement)
    assert_allclose(open_unitary_shift(moved, -np.asarray(displacement)),
                    frame, rtol=0, atol=2e-14)
    impulse = np.zeros((8, 16), complex)
    impulse[0, 0] = 1
    shifted = open_unitary_shift(impulse, [0., .4])
    assert abs(shifted[0, -1]) < 1e-12
