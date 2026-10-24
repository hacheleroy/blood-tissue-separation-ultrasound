"""Data-only phase confidence and complex-IQ preservation checks."""

from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"examples"))

from btbench.shared import PublicInput, fit_shared, shared_scene  # noqa: E402
from esdiva_phase4_differential_context import (  # noqa: E402
    DifferentialPhaseContext, differential_phase)


def test_disjoint_validation_can_reduce_confidence_without_refitting_phase():
    rng = np.random.default_rng(732)
    ref = rng.normal(size=(3, 8, 16))+1j*rng.normal(size=(3, 8, 16))
    roi = np.ones((8, 16), bool)
    g = {"prf_hz": 1000., "f0_hz": 5e6, "sound_speed_mps": 1540.,
         "phase_sign": 1, "tx_rad": [-.1, 0., .1], "rx_rad": [0., 0., 0.],
         "view_offsets_s": [0., .0001, .0002]}
    known_phase = np.array([0., .4, -.4])
    views = np.broadcast_to(ref[None]*np.exp(1j*known_phase)[None,:,None,None],
                            (16, *ref.shape)).copy()
    cal = np.broadcast_to(ref[None], (16, *ref.shape)).copy()
    original = PublicInput(views, cal, roi, g)
    first, audit = differential_phase(original, np.ones((16,3,1,1)), ref)
    z,x = np.indices(roi.shape)
    corrupted = views.copy()
    corrupted[:,:,((z+2*x)%3 == 0)] *= np.exp(1.1j)
    altered = PublicInput(corrupted, cal, roi, g)
    second, new_audit = differential_phase(altered, np.ones((16,3,1,1)), ref)
    assert audit["fitted_phase_sha256"] == new_audit["fitted_phase_sha256"]
    assert audit["fit_mask_sha256"] != audit["validation_mask_sha256"]
    assert np.mean(abs(second)) < np.mean(abs(first))*.1


def test_differential_phase_roundtrip_preserves_complex_blood():
    scene = shared_scene(seed=93, n_frames=32, calibration_frames=16,
                         shape=(8,16), motion="translation")
    context = DifferentialPhaseContext(scene.public, packet=16, score_start=16)
    result = fit_shared(context, {"method":"roundtrip", "horizon":16})
    assert_allclose(result.views, scene.public.views[16:], atol=1e-13, rtol=1e-13)
    assert_allclose(result.apply_frozen(scene.blood), scene.blood[16:],
                    atol=1e-13, rtol=1e-13)
