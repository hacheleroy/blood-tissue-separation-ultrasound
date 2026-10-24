"""Choose a tissue training region from RF-derived calibration IQ only.

Angular speckle varies with the sampled scatterer cloud. A scalar regression
gain shrinks with decorrelation, unnecessarily boosting receiver noise when
the two views have different but similarly powered tissue speckle. The bounded
phase-plus-RMS gain keeps the measured view amplitude. Its distinct validation
residual is recorded: this is an input-stability candidate, not evidence of
blood preservation or tissue rejection by a filter.
"""

import numpy as np

from esdiva_angle_calibration import angular_diagnostics


ROI_ORDER = ("outer_both", "outer_left", "outer_right", "wide_both")


def candidate_rois(x_mm, n_depth):
    """Geometric lateral tissue guards in millimetres around the vessel."""
    x_mm = np.asarray(x_mm)
    masks = (abs(x_mm) >= 6., x_mm <= -6., x_mm >= 6., abs(x_mm) >= 5.)
    return {name: np.broadcast_to(mask[None, :], (n_depth,len(x_mm))).copy()
            for name,mask in zip(ROI_ORDER,masks)}


def select_diagnostics(rows_by_roi, coherence_min=.7, max_noise_power_gain=4.,
                       min_roi_pixels=24):
    """Select from first-half calibration statistics and then audit the second."""
    if tuple(rows_by_roi) != ROI_ORDER:
        raise ValueError("Expected the four ordered geometric ROI candidates")
    candidates = {}
    for name, roi in rows_by_roi.items():
        rows = roi["angles"]
        if len(rows) != 3 or [r["angle_deg"] for r in rows] != [0., -5., 5.]:
            raise ValueError("Expected three physical 0, -5, +5 degree views")
        nonreference = rows[1:]
        candidates[name] = {"roi_pixels": int(roi["pixels"]),
                            "min_coherence_train": min(r["coherence_train"] for r in nonreference),
                            "min_coherence_validation": min(r["coherence_validation"] for r in nonreference),
                            "max_bounded_noise_power_gain_train": max(
                                r["bounded_input_noise_power_gain"] for r in nonreference),
                            "bounded_residual_validation_max": max(
                                r["bounded_residual_validation"] for r in nonreference),
                            "regression_residual_validation_max": max(
                                r["regression_residual_validation"] for r in nonreference)}
    eligible = [name for name in ROI_ORDER if
                candidates[name]["roi_pixels"] >= min_roi_pixels and
                candidates[name]["min_coherence_train"] >= coherence_min and
                candidates[name]["max_bounded_noise_power_gain_train"] <= max_noise_power_gain]
    selected = max(eligible, key=lambda name: candidates[name]["min_coherence_train"]) if eligible else None
    passed = bool(selected is not None and
                  candidates[selected]["min_coherence_validation"] >= coherence_min)
    return {"selected_roi": selected, "input_gate_pass": passed,
            "selection_uses_train_only": True,
            "validation_shares_stationary_tissue_scatterers": True,
            "coherence_min": coherence_min,
            "max_bounded_input_noise_power_gain": max_noise_power_gain,
            "min_roi_pixels": min_roi_pixels,
            "candidates": candidates}


def choose_input(calibration, x_mm, indices=(4, 3, 5), coherence_min=.7,
                 max_noise_power_gain=4., min_roi_pixels=24):
    """Fit on first-half IQ; gate the chosen ROI on independent-time IQ.

    The second half never selects or reselects a region. No blood, source truth,
    displacement or evaluation-frame IQ enters this function. Temporal halves
    still image the same stationary tissue scatterers, so the second half is
    an independent-noise check and not an independent tissue sample.
    """
    indices = tuple(indices)
    if (indices != (4, 3, 5) or calibration.ndim != 4 or
            calibration.shape[1] != 3 or len(calibration) != 16 or
            calibration.shape[2:] != (8, len(x_mm)) or
            not np.iscomplexobj(calibration) or not np.isfinite(calibration).all()):
        raise ValueError("Require the phase-3 0, -5, +5 degree 16-frame calibration tile")
    rois = candidate_rois(x_mm, calibration.shape[2])
    rows_by_roi = {}
    for name, roi in rois.items():
        if not roi.any():
            raise ValueError("Empty geometric tissue guard")
        rows = angular_diagnostics(calibration[:8], calibration[8:], roi, indices)
        rows_by_roi[name] = {"pixels": int(roi.sum()), "angles": rows}
    audit = select_diagnostics(rows_by_roi, coherence_min, max_noise_power_gain,
                               min_roi_pixels)
    selected = audit["selected_roi"]
    return rois[selected] if selected else None, audit
