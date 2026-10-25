"""Boundary-aware complex tissue reference transport for phase-7 development.

Fit displacement on observed wide-tile tissue magnitude. Interpolate only the
prior tissue calibration reference at in-bounds coordinates. Its phase change
defines a unit-magnitude, exactly reversible correction on the original core.
No current IQ, isolated blood, or known motion is spatially resampled.
"""

from dataclasses import asdict, dataclass

import numpy as np
from scipy.ndimage import map_coordinates

from btbench.shared import array_hash
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate
from esdiva_phase5_local_phase import tissue_mask
from esdiva_phase6_unit_shift import ShiftOptions, fit_shifts


@dataclass(frozen=True)
class ReferenceTransportOptions:
    transport_strength: float = 1.1
    reference_amplitude_fraction: float = 0.2
    minimum_validation_phase_gain: float = 0.01
    minimum_valid_check_pixels: int = 30


def transported_reference(reference, displacement):
    """Return interpolated complex reference and its exact measured support."""
    z, x = np.indices(reference.shape, dtype=float)
    coords = (z+displacement[0], x+displacement[1])
    valid = ((coords[0] >= 0) & (coords[0] <= reference.shape[0]-1) &
             (coords[1] >= 0) & (coords[1] <= reference.shape[1]-1))
    # Nearest supplies finite values outside the tile, but the validity mask
    # prevents those values from entering a correction or validation score.
    moved = (map_coordinates(reference.real, coords, order=1, mode="nearest",
                             prefilter=False) +
             1j*map_coordinates(reference.imag, coords, order=1, mode="nearest",
                                 prefilter=False))
    return moved, valid


def fit_reference_phase(wide, core_rows, z_mm, x_mm, options=None):
    options = options or ReferenceTransportOptions()
    if (not 0 < options.transport_strength <= 1.2 or
            not 0 < options.reference_amplitude_fraction < 1 or
            options.minimum_validation_phase_gain < 0 or
            options.minimum_valid_check_pixels < 1):
        raise ValueError("Invalid complex reference transport options")
    roi = tissue_mask(z_mm, x_mm)
    core_rows = np.asarray(core_rows)
    if (not np.array_equal(roi, wide.training_roi) or
            not np.array_equal(core_rows, np.arange(12,20))):
        raise ValueError("Changed tissue mask or central score tile")
    zz, xx = np.indices(roi.shape)
    fit = roi & ((zz+2*xx)%3 != 0)
    check = roi & ~fit
    if fit.sum() < 80 or check.sum() < 40:
        raise ValueError("Insufficient disjoint tissue support")
    displacement, shift_audit = fit_shifts(wide, z_mm, x_mm, ShiftOptions())
    # The inherited fitter also documents a cyclic IQ operation used by its
    # separate phase-6 context. That operation is never called here.
    shift_audit = {key: value for key, value in shift_audit.items()
                   if key not in ("model", "edge_assumption")}
    shift_audit["model"] = "disjoint-pixel observed-magnitude translation only"
    modulation, _, calibration, reference_views = bounded_calibrate(wide)
    normalized = wide.views/modulation
    phase = np.zeros(wide.views.shape, dtype=float)
    gains = np.zeros(wide.views.shape[:2], dtype=float)
    accepted = np.zeros_like(gains, dtype=bool)
    inbounds = np.zeros_like(gains)
    reliable = np.zeros_like(gains)
    for t in range(len(normalized)):
        for j in range(normalized.shape[1]):
            reference = reference_views[j]
            moved, valid = transported_reference(reference, displacement[t,j])
            threshold = options.reference_amplitude_fraction*np.median(abs(reference[roi]))
            support = valid & (abs(reference)>=threshold) & (abs(moved)>=threshold)
            inbounds[t,j] = np.mean(valid[core_rows])
            reliable[t,j] = np.mean(support[core_rows])
            if (np.sum(fit & support) < 80 or
                    np.sum(check & support) < options.minimum_valid_check_pixels):
                continue
            proposal = np.where(support,
                options.transport_strength*np.angle(moved*reference.conj()),0.)
            observed = np.angle(normalized[t,j]*reference.conj())
            weights = abs(normalized[t,j]*reference.conj())*(check & support)
            gain = np.sum(weights*(np.cos(observed-proposal)-np.cos(observed)))/max(
                np.sum(weights),1e-30)
            gains[t,j] = gain
            if gain > options.minimum_validation_phase_gain:
                phase[t,j] = proposal
                accepted[t,j] = True
    applied = phase[:,:,core_rows]
    audit = {"model":"in-bounds complex calibration-reference transport; phase-only core",
             "options":asdict(options),"fit_pixels":int(fit.sum()),
             "validation_pixels":int(check.sum()),
             "wide_calibration_sha256":calibration["calibration_sha256"],
             "magnitude_translation_audit":shift_audit,
             "phase_sha256":array_hash(applied),
             "accepted_view_frames":int(accepted.sum()),
             "view_frames":int(accepted.size),
             "median_validation_gain":float(np.median(gains)),
             "mean_core_inbounds_fraction":float(np.mean(inbounds)),
             "mean_core_reliable_reference_fraction":float(np.mean(reliable)),
             "boundary_rule":"no cyclic wrap; interpolate reference only; zero phase on invalid or weak support",
             "uses_observed_IQ_and_prior_calibration_only":True}
    return applied,audit


class ReferenceTransportContext(Phase3Context):
    def __init__(self,core,wide,core_rows,z_mm,x_mm,packet=16,score_start=64,
                 transport_options=None):
        if (not np.array_equal(core.calibration_views,wide.calibration_views[:,:,core_rows]) or
                not np.array_equal(core.views,wide.views[:,:,core_rows]) or
                core.geometry!=wide.geometry):
            raise ValueError("Wider observation must reproduce the core exactly")
        super().__init__(core,packet=packet,score_start=score_start,motion="none")
        phase,audit=fit_reference_phase(
            wide,core_rows,z_mm,x_mm,
            ReferenceTransportOptions(**transport_options) if transport_options else None)
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion="reference_transport"
        self.motion_summary=audit
        self.motion_diagnostics=audit
        self.calibration={**self.calibration,"reference_transport":audit}
        self.aligned=self.align(core.views)
        self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}
