"""Predeclared phase-8 calibration support and measured phase selection.

The selector sees mixture IQ, its prior calibration, and the fixed vessel
geometry. Only frames before the scored packet enter the mode decision.
Neither source components nor prescribed simulator motion enter this module.
"""

import numpy as np
from scipy.linalg import expm

from btbench.shared import array_hash
from esdiva_angle_calibration import angular_diagnostics
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate
from esdiva_phase5_local_phase import (
    LocalPhaseContext, LocalPhaseOptions, fit_local_phase, tissue_mask,
)
from esdiva_phase6_unit_shift import ShiftOptions, fit_shifts
from esdiva_phase7_reference_transport import ReferenceTransportContext


def choose_wide_guard(calibration, z_mm, x_mm):
    """Pick guard 4 or 5 by train coherence, then gate its validation half."""
    rows = []
    for guard in (4., 5.):
        roi = tissue_mask(z_mm, x_mm, guard)
        stats = angular_diagnostics(calibration[:8], calibration[8:], roi, [4, 3, 5])[1:]
        rows.append({
            "guard_mm": guard,
            "pixels": int(roi.sum()),
            "minimum_coherence_train": min(r["coherence_train"] for r in stats),
            "minimum_coherence_validation": min(r["coherence_validation"] for r in stats),
            "maximum_bounded_noise_power_gain_train": max(
                r["bounded_input_noise_power_gain"] for r in stats),
        })
    selected = max(rows, key=lambda row: row["minimum_coherence_train"])
    return tissue_mask(z_mm, x_mm, selected["guard_mm"]), {
        "candidates": rows,
        "selected": selected,
        "pass": bool(selected["pixels"] >= 128 and
                     min(selected["minimum_coherence_train"],
                         selected["minimum_coherence_validation"]) >= .700 and
                     selected["maximum_bounded_noise_power_gain_train"] <= 4.),
        "selection_uses_train_only": True,
    }


def open_unitary_shift(frame, displacement):
    """Exact inverse at negative displacement; open, nonperiodic core edges."""
    operators = []
    for size, amount in zip(frame.shape, displacement):
        derivative = (np.diag(np.full(size-1, .5), -1) -
                      np.diag(np.full(size-1, .5), 1))
        operators.append(expm(float(amount)*derivative))
    return operators[0] @ frame @ operators[1].T


class HybridContext(Phase3Context):
    """Half-strength nonperiodic rigid IQ transport and accepted local phase."""

    def __init__(self, core, wide, core_rows, z_mm, x_mm, packet=16, score_start=64):
        if (not np.array_equal(core.calibration_views,
                               wide.calibration_views[:, :, core_rows]) or
                not np.array_equal(core.views, wide.views[:, :, core_rows]) or
                core.geometry != wide.geometry or
                not np.array_equal(wide.training_roi, tissue_mask(z_mm, x_mm, 5.))):
            raise ValueError("Hybrid needs the exact guarded wider replay")
        super().__init__(core, packet=packet, score_start=score_start, motion="none")
        shift, audit = fit_shifts(wide, z_mm, x_mm, ShiftOptions())
        phase, phase_audit = fit_local_phase(
            wide, core_rows, z_mm, x_mm,
            LocalPhaseOptions(vessel_guard_mm=5., phase_fraction=.5,
                              spatial_sigma_pixels=4.,
                              minimum_validation_phase_gain=.01))
        self.shifts = .5*shift
        self.modulation = self.modulation*np.exp(1j*phase)
        self.motion = "phase8_hybrid"
        self.motion_summary = {
            "model": "validated local phase plus half nonperiodic unitary IQ shift",
            "shift_fit": {**audit, "application": "nonperiodic skew derivative",
                          "inherited_periodic_application_not_used": True},
            "shift_sha256": array_hash(self.shifts),
            "phase_fit": phase_audit,
        }
        self.motion_diagnostics = self.motion_summary
        self.calibration = {**self.calibration, "phase8_hybrid": self.motion_summary}
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


def build_context(mode, core, wide, core_rows, z_mm, x_mm, score_start=64):
    if mode == "none":
        return Phase3Context(core, packet=16, score_start=score_start, motion="none")
    guard = 4. if np.array_equal(wide.training_roi, tissue_mask(z_mm, x_mm, 4.)) else 5.
    if mode == "local":
        return LocalPhaseContext(
            core, wide, core_rows, z_mm, x_mm, packet=16, score_start=score_start,
            phase_options={"vessel_guard_mm":guard, "phase_fraction":.5,
                           "spatial_sigma_pixels":4.,
                           "minimum_validation_phase_gain":.01})
    if guard != 5.:
        raise ValueError("Reference and hybrid modes require guard 5")
    if mode == "reference":
        return ReferenceTransportContext(
            core, wide, core_rows, z_mm, x_mm, packet=16, score_start=score_start,
            transport_options={"transport_strength":1.1,
                               "reference_amplitude_fraction":.2,
                               "minimum_validation_phase_gain":.01,
                               "minimum_valid_check_pixels":30})
    if mode == "hybrid":
        return HybridContext(core, wide, core_rows, z_mm, x_mm,
                             packet=16, score_start=score_start)
    raise ValueError("Unknown phase-8 mode")


def phase_loss(context, reference, roi, stop=64):
    """Weighted 1-cos phase discordance on the pre-score check pixels."""
    mixed = context.aligned[:stop, :, roi]
    ref = np.broadcast_to(reference[:, roi], mixed.shape)
    cross = mixed*ref.conj()
    return float(np.sum(abs(cross)-cross.real)/max(np.sum(abs(cross)), 1e-30))


def select_mode(core, wide, core_rows, z_mm, x_mm, score_start=64):
    """Rule fixed on seen clouds: reference only in high, limited excess loss."""
    if score_start != 64:
        raise ValueError("Phase-8 mode must use the 64 pre-score frames")
    guard = 4. if np.array_equal(wide.training_roi, tissue_mask(z_mm, x_mm, 4.)) else 5.
    if guard == 4.:
        return "local", {"guard_mm":guard, "rule":"guard-4 local only",
                         "pre_score_frames":64}
    z, x = np.indices(wide.training_roi.shape)
    check = (wide.training_roi & ((z+2*x)%3 == 0))[core_rows]
    local = build_context("local", core, wide, core_rows, z_mm, x_mm, score_start)
    reference = build_context("reference", core, wide, core_rows, z_mm, x_mm, score_start)
    _, _, _, refs = bounded_calibrate(core)
    local_loss = phase_loss(local, refs, check, score_start)
    reference_loss = phase_loss(reference, refs, check, score_start)
    excess = reference_loss-local_loss
    mode = ("reference" if local_loss > .28 and excess < .18
            else "hybrid")
    return mode, {
        "guard_mm":guard, "rule":"reference if local_loss > 0.28 and reference excess < 0.18; otherwise hybrid",
        "pre_score_frames":64, "local_phase_loss":local_loss,
        "reference_phase_loss":reference_loss, "reference_excess_loss":excess,
        "check_pixels":int(check.sum()), "source_truth_used":False,
    }
