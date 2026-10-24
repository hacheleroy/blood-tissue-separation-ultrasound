"""Gate bulk tissue phase tracking without spatially resampling blood IQ.

Each uncompounded view gets a per-frame unit-modulus tissue phase measured on
geometric training pixels. Disjoint lateral/depth pixels check phase agreement,
coherence and residual improvement. Rejected estimates become exactly zero
phase. The calibration/reference and all decisions use observed IQ only; the
isolated blood and prescribed scatterer displacement stay with the evaluator.
"""

from dataclasses import asdict, dataclass

import numpy as np

from btbench.shared import array_hash
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate


@dataclass(frozen=True)
class PhaseGate:
    minimum_phase_rad: float = .15
    maximum_phase_disagreement_rad: float = .1
    minimum_validation_coherence: float = .8
    minimum_validation_nrmse_improvement: float = .01


def fit_phase(public, modulation, reference_views, options=None):
    """Fit on one spatial subset and accept on its disjoint complement."""
    options=options or PhaseGate()
    roi=np.asarray(public.training_roi,bool)
    z,x=np.indices(roi.shape)
    validation=roi&((z+2*x)%3==0)
    training=roi&~validation
    if training.sum()<12 or validation.sum()<8:
        raise ValueError("Insufficient independent tissue phase pixels")
    observed=public.views/modulation

    def summary(mask):
        ref=reference_views[:,mask]
        obs=observed[:,:,mask]
        cross=np.einsum("jp,tjp->tj",ref.conj(),obs)
        denominator=np.sqrt(np.sum(abs(ref)**2,axis=1)[None]*
                            np.sum(abs(obs)**2,axis=2))
        return cross,abs(cross)/np.maximum(denominator,1e-30)

    fitted,training_coherence=summary(training)
    checked,validation_coherence=summary(validation)
    phase=np.angle(fitted)
    disagreement=abs(np.angle(np.exp(1j*(phase-np.angle(checked)))))
    ref=reference_views[:,validation]
    obs=observed[:,:,validation]
    denominator=np.maximum(np.linalg.norm(obs,axis=2),1e-30)
    before=np.linalg.norm(obs-ref[None],axis=2)/denominator
    after=np.linalg.norm(obs-ref[None]*np.exp(1j*phase)[:,:,None],axis=2)/denominator
    accepted=(abs(phase)>=options.minimum_phase_rad)&(
        disagreement<=options.maximum_phase_disagreement_rad)&(
        validation_coherence>=options.minimum_validation_coherence)&(
        before-after>=options.minimum_validation_nrmse_improvement)
    applied=np.where(accepted,phase,0.)
    audit={"model":"per-frame, per-view tissue phase; no spatial warp",
           "options":asdict(options), "fit_pixels":int(training.sum()),
           "validation_pixels":int(validation.sum()),
           "fit_mask_sha256":array_hash(training),
           "validation_mask_sha256":array_hash(validation),
           "accepted_view_frames":int(accepted.sum()),
           "total_view_frames":int(accepted.size),
           "accepted_frames_any_view":int(accepted.any(axis=1).sum()),
           "accepted_frames_all_views":int(accepted.all(axis=1).sum()),
           "phase_sha256":array_hash(applied),
           "median_phase_disagreement_rad":float(np.median(disagreement)),
           "median_validation_coherence":float(np.median(validation_coherence)),
           "median_validation_improvement":float(np.median(before-after)),
           "fit_uses_calibration_and_current_frame_only":True}
    return applied,audit


class PhaseOnlyContext(Phase3Context):
    """Use exact unit-modulus phase demodulation in the shared filter API."""

    def __init__(self, public, packet=16, score_start=64, gain_model="phase_times_rms",
                 phase_options=None):
        super().__init__(public,packet=packet,score_start=score_start,motion="none",
                         gain_model=gain_model)
        modulation,_,_,reference_views=bounded_calibrate(public,gain_model)
        if not np.array_equal(modulation,self.modulation):
            raise AssertionError("The phase reference changed after calibration")
        phase,audit=fit_phase(public,modulation,reference_views,
                              PhaseGate(**phase_options) if phase_options else None)
        self.modulation=self.modulation*np.exp(1j*phase)[:,:,None,None]
        self.motion="phase_only"
        self.motion_summary=audit
        self.motion_diagnostics=audit
        self.calibration={**self.calibration,"phase_gate":audit}
        self.aligned=self.align(public.views)
        self.x=self.aligned.reshape(len(public.views),-1).T
        self.cache={}
