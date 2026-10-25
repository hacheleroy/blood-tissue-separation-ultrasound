"""Development-only data-fitted reversible rigid shift on the phase-5 core.

Estimate translation from held-out tissue pixels on the wider observed IQ,
then apply an exactly invertible Fourier shift on the eight-row score tile.
The Fourier operation assumes periodic boundaries on that small tile; this
limitation must be checked in source-resolved scores before any freeze.
"""

from dataclasses import asdict, dataclass

import numpy as np
from scipy.ndimage import fourier_shift, shift
from scipy.optimize import minimize

from btbench.shared import array_hash
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate
from esdiva_phase5_local_phase import LocalPhaseOptions, fit_local_phase, tissue_mask


@dataclass(frozen=True)
class ShiftOptions:
    maximum_axial_pixels: float = 1.25
    maximum_lateral_pixels: float = .75
    minimum_validation_gain: float = .01


def unitary_shift(frame, displacement):
    """Invertible cyclic translation of a complex 2-D IQ frame."""
    return np.fft.ifftn(fourier_shift(np.fft.fftn(frame),displacement))


def fit_shifts(wide, z_mm, x_mm, options=None):
    options=options or ShiftOptions()
    roi=tissue_mask(z_mm,x_mm)
    if not np.array_equal(roi,wide.training_roi):
        raise ValueError('Wider tissue mask differs from fixed geometry')
    zz,xx=np.indices(roi.shape)
    fit=roi&((zz+2*xx)%3!=0)
    check=roi&~fit
    if fit.sum()<80 or check.sum()<40:
        raise ValueError('Insufficient disjoint translation support')
    modulation,_,calibration,reference_views=bounded_calibrate(wide)
    reference=abs(reference_views)
    movements=np.zeros((len(wide.views),wide.views.shape[1],2),float)
    improvements=np.zeros(movements.shape[:2])
    candidates=np.zeros_like(movements)
    bounds=[(-options.maximum_axial_pixels,options.maximum_axial_pixels),
            (-options.maximum_lateral_pixels,options.maximum_lateral_pixels)]
    for t in range(len(wide.views)):
        for j in range(wide.views.shape[1]):
            current=abs(wide.views[t,j]/modulation[t,j])
            target=reference[j]
            def loss(d,mask):
                moved=shift(current,d,order=1,mode='nearest',prefilter=False)
                observed=moved[mask];expected=target[mask]
                scale=np.vdot(expected,observed).real/max(np.vdot(observed,observed).real,1e-30)
                return float(np.mean((scale*observed-expected)**2)/
                             max(np.mean(expected**2),1e-30))
            candidate=minimize(lambda d:loss(d,fit),[0.,0.],method='Powell',bounds=bounds,
                               options={'maxiter':30,'xtol':.03,'ftol':.001})
            e0=loss([0.,0.],check)
            e1=loss(candidate.x,check)
            gain=(e0-e1)/max(e0,1e-30)
            candidates[t,j]=candidate.x
            improvements[t,j]=gain
            if candidate.success and gain>options.minimum_validation_gain:
                movements[t,j]=candidate.x
    audit={'model':'disjoint-pixel magnitude translation; cyclic unitary IQ shift',
           'options':asdict(options),'fit_pixels':int(fit.sum()),
           'validation_pixels':int(check.sum()),
           'wide_calibration_sha256':calibration['calibration_sha256'],
           'shift_sha256':array_hash(movements),'candidate_sha256':array_hash(candidates),
           'accepted_view_frames':int(np.any(movements!=0,axis=-1).sum()),
           'view_frames':int(improvements.size),
           'median_validation_gain':float(np.median(improvements)),
           'mean_applied_shift_pixels':float(np.mean(np.linalg.norm(movements,axis=-1))),
           'uses_observed_IQ_and_prior_calibration_only':True,
           'edge_assumption':'periodic core, not an anatomical boundary condition'}
    return movements,audit


class UnitShiftContext(Phase3Context):
    def __init__(self,core,wide,core_rows,z_mm,x_mm,packet=16,score_start=64,
                 shift_options=None,phase_options=None):
        if (not np.array_equal(core.calibration_views,wide.calibration_views[:,:,core_rows]) or
                not np.array_equal(core.views,wide.views[:,:,core_rows]) or
                core.geometry!=wide.geometry):
            raise ValueError('Wider observations must replay the original core')
        super().__init__(core,packet=packet,score_start=score_start,motion='none')
        movements,audit=fit_shifts(wide,z_mm,x_mm,
                                   ShiftOptions(**shift_options) if shift_options else None)
        if phase_options is not None:
            phase,phase_audit=fit_local_phase(wide,core_rows,z_mm,x_mm,
                                             LocalPhaseOptions(**phase_options))
            self.modulation=self.modulation*np.exp(1j*phase)
            audit['local_phase']=phase_audit
        self.shifts=movements
        self.motion='unit_shift'
        self.motion_summary=audit
        self.motion_diagnostics=audit
        self.calibration={**self.calibration,'unit_shift':audit}
        self.aligned=self.align(core.views)
        self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}

    def align(self,other):
        value=other/self.modulation
        if not hasattr(self,'shifts'):
            return value
        return np.stack([np.stack([unitary_shift(frame,d) for frame,d in zip(row,ds)])
                         for row,ds in zip(value,self.shifts)])

    def restore_block(self,aligned,start):
        if not hasattr(self,'shifts'):
            return super().restore_block(aligned,start)
        native=np.stack([np.stack([unitary_shift(frame,-d) for frame,d in zip(row,ds)])
                         for row,ds in zip(aligned,self.shifts[start:start+len(aligned)])])
        return native*self.modulation[start:start+len(aligned)]
