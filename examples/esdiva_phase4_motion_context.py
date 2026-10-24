"""Fit motion from the tissue envelope with disjoint spatial validation.

The same phase-and-RMS angular calibration and geometric tissue ROI feed the
phase-3 filters. For envelope registration only, phase changes from a small
axial scatterer displacement do not count as image mismatch. The fitted map
still passes the existing independent-pixel motion gate, and exactly zero
motion is applied when it fails. No isolated source or prescribed displacement
enters the fit or its gate.
"""

import hashlib
import json

import numpy as np

from btbench.shared import MotionGate, estimate_displacement, inverse_coordinates
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate


class EnvelopeMotionContext(Phase3Context):
    """Use calibration-envelope geometry with the shared linear filters."""

    def __init__(self, public, packet=16, score_start=64, motion="smooth",
                 motion_options=None, gain_model="phase_times_rms"):
        if motion not in ("rigid", "smooth"):
            raise ValueError("Envelope motion requires rigid or smooth registration")
        super().__init__(public, packet=packet, score_start=score_start, motion="none",
                         gain_model=gain_model)
        modulation, reference, _, _ = bounded_calibrate(public, gain_model)
        if not np.array_equal(modulation, self.modulation):
            raise AssertionError("Calibration changed between envelope and filtering")
        displacement, diagnostics = estimate_displacement(
            public.views/self.modulation, reference, public.training_roi, motion,
            reference_views=None,
            options=MotionGate(**motion_options) if motion_options else None)
        self.displacement = displacement
        self.motion = motion
        self.motion_diagnostics = diagnostics
        if "audit_rows" in diagnostics:
            diagnostics["audit_sha256"] = hashlib.sha256(json.dumps(
                diagnostics["audit_rows"], separators=(",", ":"),
                allow_nan=False).encode()).hexdigest()
        self.motion_summary = {key:value for key,value in diagnostics.items()
                               if key != "audit_rows"}
        self.has_displacement = bool(np.any(displacement))
        grid = np.indices(public.views.shape[2:], dtype=float)
        self.alignment_coordinates = grid[None]+displacement
        self.restoration_coordinates = (np.stack([inverse_coordinates(d) for d in displacement])
                                        if self.has_displacement else self.alignment_coordinates)
        self.aligned = self.align(public.views)
        self.x = self.aligned.reshape(len(public.views),-1).T
        self.cache = {}
