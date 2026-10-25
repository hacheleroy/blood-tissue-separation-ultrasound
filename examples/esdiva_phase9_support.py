"""Guarded source-independent calibration eligibility and fixed phase-9 motion.

The calibration phase map is a gate-only stationary-tissue diagnostic. It is
fitted on the first eight prior tissue calibration frames and checked on the
last eight. The regular filter still receives unmodified mixture IQ and its
prior calibration. Blood truth and prescribed motion are never inputs here.
"""

import numpy as np
from scipy.ndimage import gaussian_filter

from btbench.shared import array_hash
from esdiva_angle_calibration import angular_diagnostics
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase5_local_phase import (
    LocalPhaseOptions, fit_local_phase, tissue_mask,
)
from esdiva_phase6_unit_shift import ShiftOptions, fit_shifts
from esdiva_phase7_reference_transport import ReferenceTransportContext
from esdiva_phase8_selector import open_unitary_shift, select_mode


def phase_coherent_support(calibration, z_mm, x_mm, core_rows):
    """Certify the fixed 5 mm tissue region after train-only angular phase fit."""
    roi = tissue_mask(z_mm, x_mm, 5.)
    core_rows = np.asarray(core_rows)
    if (calibration.shape != (16, 3, 32, 16) or
            not np.array_equal(core_rows, np.arange(12,20)) or
            roi.sum() != 192 or roi[core_rows].sum() != 48):
        raise ValueError("Changed phase-9 calibration geometry")
    registered = calibration.copy()
    phase = np.zeros((3,32,16), dtype=float)
    smooth_mask = gaussian_filter(roi.astype(float), .4, mode="constant")
    for j in (1,2):
        cross = np.mean(calibration[:8,j]*calibration[:8,0].conj(), axis=0)
        local = gaussian_filter(cross*roi, .4, mode="constant")/np.maximum(
            smooth_mask, 1e-15)
        phase[j] = np.angle(local)
        registered[:,j] *= np.exp(-1j*phase[j])[None]
    result = {"model":"train-only 0.4-pixel angular phase map; gate only",
              "wide_pixels":int(roi.sum()),"core_pixels":int(roi[core_rows].sum()),
              "phase_sha256":array_hash(phase),
              "calibration_sha256":array_hash(calibration),
              "fit_frames":8,"validation_frames":8,
              "uses_source_truth_or_evaluation_IQ":False}
    for name, array, mask in (
            ("wide", registered, roi),
            ("core", registered[:,:,core_rows], roi[core_rows])):
        rows = angular_diagnostics(array[:8], array[8:], mask, (4,3,5))[1:]
        train = min(row["coherence_train"] for row in rows)
        validation = min(row["coherence_validation"] for row in rows)
        noise = max(row["bounded_input_noise_power_gain"] for row in rows)
        result[name] = {
            "minimum_coherence_train":train,
            "minimum_coherence_validation":validation,
            "maximum_bounded_noise_power_gain_train":noise,
            "pass":bool(
                train >= (.7 if name=="wide" else .5) and
                validation >= (.7 if name=="wide" else .5) and
                noise <= 4. and
                (name=="wide" or validation >= .85*train)),
        }
    result["pass"] = result["wide"]["pass"] and result["core"]["pass"]
    return roi, roi[core_rows].copy(), result


class Phase9HybridContext(Phase3Context):
    """Unitary open-boundary shift at 0.5 plus accepted local phase at 0.6."""

    def __init__(self, core, wide, core_rows, z_mm, x_mm, score_start=64):
        expected = tissue_mask(z_mm, x_mm, 5.)
        if (not np.array_equal(wide.training_roi, expected) or
                not np.array_equal(core.training_roi, expected[core_rows]) or
                not np.array_equal(core.calibration_views,
                                   wide.calibration_views[:,:,core_rows]) or
                not np.array_equal(core.views, wide.views[:,:,core_rows]) or
                core.geometry != wide.geometry):
            raise ValueError("Phase-9 core and wider mixture must match exactly")
        super().__init__(core, packet=16, score_start=score_start, motion="none")
        shift, shift_audit = fit_shifts(wide, z_mm, x_mm, ShiftOptions())
        phase, phase_audit = fit_local_phase(
            wide, core_rows, z_mm, x_mm,
            LocalPhaseOptions(vessel_guard_mm=5., phase_fraction=.6,
                              spatial_sigma_pixels=4.,
                              minimum_validation_phase_gain=.01))
        self.shifts = .5*shift
        self.modulation = self.modulation*np.exp(1j*phase)
        self.motion = "phase9_hybrid"
        self.motion_summary = {
            "model":"accepted local phase 0.6 plus half nonperiodic unitary IQ shift",
            "shift_fit":{**shift_audit,
                         "application":"nonperiodic skew derivative",
                         "inherited_periodic_application_not_used":True},
            "shift_sha256":array_hash(self.shifts),
            "phase_fit":phase_audit,
        }
        self.motion_diagnostics = self.motion_summary
        self.calibration = {**self.calibration, "phase9_hybrid":self.motion_summary}
        self.aligned = self.align(core.views)
        self.x = self.aligned.reshape(len(core.views), -1).T
        self.cache = {}

    def align(self, other):
        value = other/self.modulation
        if not hasattr(self, "shifts"):
            return value
        return np.stack([np.stack([open_unitary_shift(frame, delta)
                                   for frame, delta in zip(row, deltas)])
                         for row, deltas in zip(value, self.shifts)])

    def restore_block(self, aligned, start):
        if not hasattr(self, "shifts"):
            return super().restore_block(aligned, start)
        native = np.stack([np.stack([open_unitary_shift(frame, -delta)
                                     for frame, delta in zip(row, deltas)])
                           for row, deltas in
                           zip(aligned, self.shifts[start:start+len(aligned)])])
        return native*self.modulation[start:start+len(aligned)]


def choose_phase9_mode(core, wide, core_rows, z_mm, x_mm):
    """Retain the frozen phase-8 mixture-only phase-consistency switch."""
    selected, audit = select_mode(core, wide, core_rows, z_mm, x_mm,
                                  score_start=64)
    if selected not in ("reference", "hybrid"):
        raise ValueError("The fixed 5 mm guard needs reference or hybrid mode")
    return selected, audit


def build_phase9_context(mode, core, wide, core_rows, z_mm, x_mm):
    if mode == "none":
        return Phase3Context(core, packet=16, score_start=64, motion="none")
    if mode == "hybrid":
        return Phase9HybridContext(core, wide, core_rows, z_mm, x_mm)
    if mode == "reference":
        return ReferenceTransportContext(
            core, wide, core_rows, z_mm, x_mm, packet=16, score_start=64,
            transport_options={"transport_strength":1.1,
                               "reference_amplitude_fraction":.2,
                               "minimum_validation_phase_gain":.01,
                               "minimum_valid_check_pixels":30})
    raise ValueError("Unknown phase-9 mode")
