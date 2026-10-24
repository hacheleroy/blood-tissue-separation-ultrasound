"""Reversible local phase registration fitted on the wider Example 22 tissue tile.

The phase field is inferred from geometric tissue pixels of the current IQ
observation and a disjoint-time tissue calibration reference. Its fit and
validation pixels do not overlap. It multiplies the original eight-row scoring
tile by a unit-magnitude factor, so no interpolation of blood is performed.
Isolated sources and prescribed motion are never supplied to this module.
"""

from dataclasses import asdict, dataclass

import numpy as np
from scipy.ndimage import gaussian_filter

from btbench.shared import array_hash
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate


@dataclass(frozen=True)
class LocalPhaseOptions:
    vessel_guard_mm: float = 5.
    spatial_sigma_pixels: float = 4.
    phase_fraction: float = .5
    minimum_validation_phase_gain: float = .01


def tissue_mask(z_mm, x_mm, guard_mm=5.):
    """Example 22's known 60-degree vessel axis, with a fixed tissue guard."""
    z_mm, x_mm = np.asarray(z_mm), np.asarray(x_mm)
    if (z_mm.shape != (32,) or x_mm.shape != (16,) or
            not np.allclose(np.diff(z_mm), .1, atol=1e-10) or
            not np.allclose(z_mm[[0,-1]], [14.4,17.5], atol=1e-10) or
            not np.allclose(x_mm, np.arange(-8,8), atol=1e-10) or guard_mm <= 0):
        raise ValueError("Expected the fixed Example 22 wider beamforming grid")
    return abs(x_mm[None,:]-np.tan(np.deg2rad(60.))*(z_mm[:,None]-16.))>=guard_mm


def fit_local_phase(wide, core_rows, z_mm, x_mm, options=None):
    """Fit/check phase separately for each observed frame and physical view."""
    options = options or LocalPhaseOptions()
    if not 0 < options.phase_fraction <= 1 or options.spatial_sigma_pixels <= 0:
        raise ValueError("Invalid local phase options")
    roi=tissue_mask(z_mm,x_mm,options.vessel_guard_mm)
    core_rows=np.asarray(core_rows)
    if (not np.array_equal(core_rows,np.arange(12,20)) or
            not np.array_equal(roi,wide.training_roi)):
        raise ValueError("Wider and scoring tiles or tissue ROI differ")
    zz,xx=np.indices(roi.shape)
    fit=roi&((zz+2*xx)%3!=0)
    check=roi&~fit
    if fit.sum()<80 or check.sum()<40:
        raise ValueError("Too few disjoint tissue pixels")
    modulation,_,calibration,reference=bounded_calibrate(wide)
    normalized=wide.views/modulation
    phase=np.empty(wide.views.shape,float)
    sigma=options.spatial_sigma_pixels
    for t in range(len(normalized)):
        for j in range(normalized.shape[1]):
            cross=normalized[t,j]*reference[j].conj()
            real=gaussian_filter((cross*fit).real,sigma,mode="nearest")
            imag=gaussian_filter((cross*fit).imag,sigma,mode="nearest")
            phase[t,j]=np.angle(real+1j*imag)
    # Half-angle correction needs a continuous phase across the depth axis.
    phase=np.unwrap(phase,axis=2)
    observed=np.angle(normalized*reference.conj()[None])[:,:,check]
    proposal=phase[:,:,check]*options.phase_fraction
    weights=abs(normalized*reference.conj()[None])[:,:,check]
    validation_gain=np.sum(weights*(np.cos(observed-proposal)-np.cos(observed)),axis=2)/np.maximum(
        weights.sum(axis=2),1e-30)
    applied=np.where(validation_gain>options.minimum_validation_phase_gain,
                     options.phase_fraction,0.)[:,:,None,None]*phase[:,:,core_rows]
    audit={"model":"masked, spatially smoothed local tissue phase; no IQ resampling",
           "options":asdict(options),"fit_pixels":int(fit.sum()),"validation_pixels":int(check.sum()),
           "fit_mask_sha256":array_hash(fit),"validation_mask_sha256":array_hash(check),
           "wide_calibration_sha256":calibration["calibration_sha256"],
           "phase_sha256":array_hash(applied),"accepted_view_frames":int(np.sum(
               validation_gain>options.minimum_validation_phase_gain)),
           "view_frames":int(validation_gain.size),
           "median_validation_phase_gain":float(np.median(validation_gain)),
           "minimum_validation_phase_gain_observed":float(validation_gain.min()),
           "mean_applied_phase_abs_rad":float(np.mean(abs(applied))),
           "uses_current_frame_and_prior_calibration_only":True}
    return applied,audit


class LocalPhaseContext(Phase3Context):
    """Apply the validated, unitary phase field on the original score tile."""

    def __init__(self, core, wide, core_rows, z_mm, x_mm, packet=16,
                 score_start=64, phase_options=None):
        if (not np.array_equal(core.calibration_views,wide.calibration_views[:,:,core_rows]) or
                not np.array_equal(core.views,wide.views[:,:,core_rows]) or
                core.geometry!=wide.geometry):
            raise ValueError("Wider observations must reproduce the original core exactly")
        super().__init__(core,packet=packet,score_start=score_start,motion="none")
        phase,audit=fit_local_phase(wide,core_rows,z_mm,x_mm,
                                    LocalPhaseOptions(**phase_options) if phase_options else None)
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion="local_phase"
        self.motion_summary=audit
        self.motion_diagnostics=audit
        self.calibration={**self.calibration,"local_phase":audit}
        self.aligned=self.align(core.views)
        self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}
