"""Track differential tissue phase across nearby ultrasound views.

Rigid and smooth sub-wavelength scatterer motion changes the IQ carrier
phase before the envelope moves a full pixel. A common phase across views
does not alter a temporal SVD subspace, so this candidate removes only the
view-specific phase around their circular mean. The fit uses geometric tissue
pixels; disjoint pixels attenuate a poorly supported correction continuously.
It applies a unit-modulus factor and never spatially resamples blood IQ.
"""

from dataclasses import asdict, dataclass

import numpy as np

from btbench.shared import array_hash
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate


@dataclass(frozen=True)
class DifferentialPhaseGate:
    phase_deadband_rad: float = .02
    phase_ramp_rad: float = .25
    validation_coherence_low: float = .65
    validation_coherence_ramp: float = .15
    phase_disagreement_low_rad: float = .1
    phase_disagreement_ramp_rad: float = .1


def differential_phase(public, modulation, reference_views, gate=None):
    """Fit per-view phase on one spatial subset; weight it on disjoint pixels."""
    gate=gate or DifferentialPhaseGate()
    roi=np.asarray(public.training_roi,bool)
    z,x=np.indices(roi.shape)
    check=roi&((z+2*x)%3==0)
    fit=roi&~check
    if fit.sum()<12 or check.sum()<8:
        raise ValueError("Too few tissue pixels for independent phase validation")
    normalized=public.views/modulation

    def cross(mask):
        ref=reference_views[:,mask]
        obs=normalized[:,:,mask]
        v=np.einsum("jp,tjp->tj",ref.conj(),obs)
        c=abs(v)/np.maximum(np.sqrt(np.sum(abs(ref)**2,axis=1)[None]*
                                    np.sum(abs(obs)**2,axis=2)),1e-30)
        return v,c

    trained,_=cross(fit)
    checked,coherence=cross(check)
    fitted_phase=np.angle(trained)
    validated_phase=np.angle(checked)
    disagreement=abs(np.angle(np.exp(1j*(fitted_phase-validated_phase))))
    common=np.angle(np.sum(np.exp(1j*fitted_phase),axis=1))
    differential=np.angle(np.exp(1j*(fitted_phase-common[:,None])))
    amplitude=np.clip((abs(differential)-gate.phase_deadband_rad)/gate.phase_ramp_rad,0.,1.)
    confidence=np.clip((coherence-gate.validation_coherence_low)/
                       gate.validation_coherence_ramp,0.,1.)
    agreement=np.clip((gate.phase_disagreement_low_rad+
                       gate.phase_disagreement_ramp_rad-disagreement)/
                      gate.phase_disagreement_ramp_rad,0.,1.)
    weight=amplitude*confidence*agreement
    applied=differential*weight
    audit={"model":"differential view phase; no spatial displacement or resampling",
           "options":asdict(gate),"fit_pixels":int(fit.sum()),
           "validation_pixels":int(check.sum()),
           "fit_mask_sha256":array_hash(fit),
           "validation_mask_sha256":array_hash(check),
           "fitted_phase_sha256":array_hash(fitted_phase),
           "applied_phase_sha256":array_hash(applied),
           "mean_applied_phase_abs_rad":float(np.mean(abs(applied))),
           "median_validation_coherence":float(np.median(coherence)),
           "median_phase_disagreement_rad":float(np.median(disagreement)),
           "view_frames_nonzero":int(np.sum(weight>0)),
           "view_frames_full_weight":int(np.sum(weight>=1.)),
           "total_view_frames":int(weight.size),
           "uses_current_frame_and_prior_calibration_only":True}
    return applied,audit


class DifferentialPhaseContext(Phase3Context):
    """Shared filters with data-only, confidence-weighted angular phase."""

    def __init__(self, public, packet=16, score_start=64, gain_model="phase_times_rms",
                 phase_options=None):
        super().__init__(public,packet=packet,score_start=score_start,motion="none",
                         gain_model=gain_model)
        modulation,_,_,reference=bounded_calibrate(public,gain_model)
        if not np.array_equal(modulation,self.modulation):
            raise AssertionError("Differential phase changed calibration")
        phase,audit=differential_phase(public,modulation,reference,
                                      DifferentialPhaseGate(**phase_options) if phase_options else None)
        self.modulation=self.modulation*np.exp(1j*phase)[:,:,None,None]
        self.motion="differential_phase"
        self.motion_summary=audit
        self.motion_diagnostics=audit
        self.calibration={**self.calibration,"differential_phase":audit}
        self.aligned=self.align(public.views)
        self.x=self.aligned.reshape(len(public.views),-1).T
        self.cache={}
