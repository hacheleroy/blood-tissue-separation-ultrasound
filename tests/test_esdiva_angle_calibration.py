"""Check the calibration gate with view signals of known coherence."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from esdiva_angle_calibration import angular_diagnostics, calibration_gate  # noqa: E402
from esdiva_gated_stap import fit_gated_loaded_stap  # noqa: E402
from btbench.shared import SharedContext, fit_shared, shared_scene


def test_regression_gain_tracks_coherence_while_rms_gain_tracks_amplitude():
    rng = np.random.default_rng(103)
    ref = rng.standard_normal((8, 4, 8)) + 1j*rng.standard_normal((8, 4, 8))
    independent = rng.standard_normal(ref.shape) + 1j*rng.standard_normal(ref.shape)
    independent -= np.vdot(ref, independent)/np.vdot(ref, ref)*ref
    independent *= np.linalg.norm(ref)/np.linalg.norm(independent)
    views = np.stack([ref, .8*ref + .6*independent, .1*ref + np.sqrt(.99)*independent], axis=1)
    rows = angular_diagnostics(views, views, np.ones(ref.shape[1:], bool), (4, 3, 5))
    assert np.isclose(rows[1]["gain_rms_magnitude"], 1.)
    assert np.isclose(rows[1]["gain_regression_magnitude"], .8)
    assert np.isclose(rows[1]["coherence_validation"], .8)
    assert rows[2]["regression_input_noise_power_gain"] > 90
    assert rows[2]["bounded_input_noise_power_gain"] < 1.01
    assert not calibration_gate(rows, (4, 3, 5))[1]
    assert calibration_gate(rows[:2], (4, 3))[1]
    assert not calibration_gate(rows[:2], (3, 4))[1]


def test_insufficient_local_stap_support_falls_back_to_orthogonal():
    scene = shared_scene(seed=63, shape=(4, 8), n_frames=32, calibration_frames=16)
    context = SharedContext(scene.public, packet=16, score_start=16)
    spec = {"method": "stap", "horizon": 16, "adaptive": True,
            "loading_rel": .01}
    result, gate = fit_gated_loaded_stap(context, spec)
    reference = fit_shared(context, {**spec, "adaptive": False})
    assert gate["selected"] == "orthogonal_fallback"
    assert not gate["support_pass"] and gate["data_only"]
    np.testing.assert_allclose(result.views, reference.views, atol=0, rtol=0)
    np.testing.assert_allclose(result.apply_frozen(scene.blood),
                               reference.apply_frozen(scene.blood), atol=0, rtol=0)
