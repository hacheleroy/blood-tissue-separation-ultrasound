"""Keep spent evaluation seeds out of the development angle survey."""

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"examples"))
from esdiva_phase3_angle_survey import survey  # noqa: E402
from esdiva_phase3_input_gate import candidate_rois, select_diagnostics  # noqa: E402
from esdiva_phase3_bounded_context import Phase3Context  # noqa: E402
from btbench.shared import fit_shared, shared_scene  # noqa: E402


@pytest.mark.parametrize("seed", [2026, 2031, 2053])
def test_rejects_old_development_and_spent_heldout_seeds(tmp_path, seed):
    with pytest.raises(ValueError, match="new development seed"):
        survey(tmp_path, tmp_path/"unused.json", seed)


def test_rejects_nonreference_first_or_repeated_view(tmp_path):
    with pytest.raises(ValueError, match="0-degree reference"):
        survey(tmp_path, tmp_path/"unused.json", 2071, indices=(3, 4, 5))
    with pytest.raises(ValueError, match="0-degree reference"):
        survey(tmp_path, tmp_path/"unused.json", 2071, indices=(4, 3, 3))


def test_training_selection_cannot_use_validation_to_reselect():
    def record(train, validation, noise=1.):
        return {"coherence_train": train, "coherence_validation": validation,
                "bounded_input_noise_power_gain": noise,
                "bounded_residual_validation": 1-validation,
                "regression_residual_validation": 1-validation}
    reference = {"angle_deg": 0., **record(1., 1.)}
    def roi(coherence, validation, pixels=24, noise=1.):
        return {"pixels": pixels, "angles": [reference,
                {"angle_deg": -5., **record(coherence, validation, noise)},
                {"angle_deg": 5., **record(coherence, validation, noise)}]}
    rows = {"outer_both": roi(.84, .83), "outer_left": roi(.91, .64),
            "outer_right": roi(.98, .98, pixels=16), "wide_both": roi(.89, .88, noise=4.1)}
    chosen = select_diagnostics(rows)
    assert chosen["selected_roi"] == "outer_left"
    assert not chosen["input_gate_pass"]
    assert chosen["min_roi_pixels"] == 24


def test_candidate_regions_have_declared_spatial_support():
    import numpy as np
    rois = candidate_rois(np.arange(-8., 8.), 8)
    assert [int(rois[name].sum()) for name in rois] == [40, 24, 16, 56]


def test_bounded_context_retains_identity_and_uses_only_public_input():
    import numpy as np
    scene = shared_scene(seed=19, n_frames=32, calibration_frames=16,
                         shape=(8, 16), motion="none")
    context = Phase3Context(scene.public, packet=16, score_start=16)
    identity = fit_shared(context, {"method":"identity", "horizon":16})
    np.testing.assert_allclose(identity.views, scene.public.views[16:], rtol=1e-13, atol=1e-13)
    assert context.calibration["gain_model"] == "phase_times_rms"
    assert context.calibration["train_frames"] == 8
