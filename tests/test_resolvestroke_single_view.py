"""Numerical checks for the separate Resolve Stroke geometry diagnostic."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from resolvestroke_single_view import power_maps, score_power  # noqa: E402


def test_tube_auc_handles_ties_and_depth_matched_background():
    power = np.ones((3, 5))
    positive = np.zeros_like(power, dtype=bool)
    positive[0, :3] = True
    outside = ~positive
    assert score_power(power, positive, outside)['auc'] == .5
    power[positive] = 5
    score = score_power(power, positive, outside)
    assert score['auc'] == 1
    assert score['tpr_at_outside_99th_percentile'] == 1


def test_complex_filters_suppress_static_clutter_and_keep_oscillating_target():
    t = np.arange(64)
    iq = np.ones((64, 16, 16), dtype=np.complex64) * (10 + 2j)
    iq[:, 6:10, 6:10] += 3 * np.exp(2j * np.pi * .1 * t)[:, None, None]
    power = power_maps(iq, 4000)
    outside = np.ones((16, 16), dtype=bool)
    outside[6:10, 6:10] = False
    for filtered in power[1:]:
        assert filtered[outside].mean() < 1e-3 * power[0, outside].mean()
        assert filtered[6:10, 6:10].mean() > .02 * power[0, 6:10, 6:10].mean()
