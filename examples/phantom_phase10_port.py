"""Geometry-aware phase-10 candidate for a future independently referenced phantom.

This port is separate from the frozen eSDIva scorer. The vessel centerline is
measured from a preacquisition structural image and never fitted to flow IQ.
No physical capture is bundled; equivalence on Example 22 is a prerequisite.
"""
from dataclasses import asdict

import numpy as np
from scipy.ndimage import gaussian_filter, shift
from scipy.optimize import minimize

from btbench.shared import array_hash
from esdiva_angle_calibration import angular_diagnostics
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate
from esdiva_phase5_local_phase import LocalPhaseOptions
from esdiva_phase6_unit_shift import ShiftOptions
from esdiva_phase7_reference_transport import ReferenceTransportOptions, transported_reference
from esdiva_phase8_selector import open_unitary_shift, phase_loss


CORE_ROWS = np.arange(12,20)


def structural_support(z_mm,x_mm,axis_x_mm):
    """Fixed lateral-distance guards around a structurally measured vessel axis."""
    z,x,axis=map(np.asarray,(z_mm,x_mm,axis_x_mm))
    if (z.shape!=(32,) or x.shape!=(16,) or axis.shape!=(32,) or
            not np.isfinite(z).all() or not np.isfinite(x).all() or
            not np.isfinite(axis).all() or not np.all(np.diff(z)>0) or
            not np.all(np.diff(x)>0)):
        raise ValueError('Physical grid and independent vessel axis are required')
    distance=abs(x[None,:]-axis[:,None])
    wide=distance>=5.
    core=(distance>=5.5)[CORE_ROWS]
    if wide.sum()<192 or core.sum()<40 or np.any(core&~wide[CORE_ROWS]):
        raise ValueError('Insufficient physical phase-10 tissue support')
    return wide,core


def calibration_support(calibration,wide_roi,core_roi):
    """Apply the unchanged held-half angular coherence gate, without IQ replacement."""
    if (calibration.shape!=(16,3,32,16) or wide_roi.shape!=(32,16) or
            core_roi.shape!=(8,16)):
        raise ValueError('Expected 16 prior frames, three views, and the 32x16 tile')
    registered=calibration.copy()
    phase=np.zeros((3,32,16),float)
    denominator=gaussian_filter(wide_roi.astype(float),.4,mode='constant')
    for j in (1,2):
        cross=np.mean(calibration[:8,j]*calibration[:8,0].conj(),axis=0)
        phase[j]=np.angle(gaussian_filter(cross*wide_roi,.4,mode='constant')/
                          np.maximum(denominator,1e-15))
        registered[:,j]*=np.exp(-1j*phase[j])[None]
    result={'model':'train-only angular phase map; gate only',
            'wide_guard_mm':5.,'core_guard_mm':5.5,
            'wide_pixels':int(wide_roi.sum()),'core_pixels':int(core_roi.sum()),
            'phase_sha256':array_hash(phase),
            'calibration_sha256':array_hash(calibration),
            'fit_frames':8,'validation_frames':8,
            'uses_source_truth_or_evaluation_IQ':False}
    for name,array,roi in (('wide',registered,wide_roi),
                           ('core',registered[:,:,CORE_ROWS],core_roi)):
        stats=angular_diagnostics(array[:8],array[8:],roi,(4,3,5))[1:]
        train=min(row['coherence_train'] for row in stats)
        validation=min(row['coherence_validation'] for row in stats)
        noise=max(row['bounded_input_noise_power_gain'] for row in stats)
        floor=.7 if name=='wide' else .5
        result[name]={'minimum_coherence_train':train,
                      'minimum_coherence_validation':validation,
                      'maximum_bounded_noise_power_gain_train':noise,
                      'pass':bool(train>=floor and validation>=floor and
                                  noise<=4. and
                                  (name=='wide' or validation>=.85*train))}
    result['pass']=result['wide']['pass'] and result['core']['pass']
    return result


def _fit_check(roi):
    zz,xx=np.indices(roi.shape)
    fit=roi&((zz+2*xx)%3!=0)
    check=roi&~fit
    if fit.sum()<80 or check.sum()<40:
        raise ValueError('Insufficient disjoint tissue support')
    return fit,check


def fit_shifts(wide,wide_roi,options=None):
    options=options or ShiftOptions()
    if not np.array_equal(wide.training_roi,wide_roi):
        raise ValueError('Wide structural mask mismatch')
    fit,check=_fit_check(wide_roi)
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


def fit_local_phase(wide,wide_roi,options=None):
    options=options or LocalPhaseOptions(vessel_guard_mm=5.)
    if (not 0<options.phase_fraction<=1 or options.spatial_sigma_pixels<=0 or
            not np.array_equal(wide.training_roi,wide_roi)):
        raise ValueError('Invalid local phase support/options')
    fit,check=_fit_check(wide_roi)
    modulation,_,calibration,reference=bounded_calibrate(wide)
    normalized=wide.views/modulation
    phase=np.empty(wide.views.shape,float)
    sigma=options.spatial_sigma_pixels
    for t in range(len(normalized)):
        for j in range(normalized.shape[1]):
            cross=normalized[t,j]*reference[j].conj()
            real=gaussian_filter((cross*fit).real,sigma,mode='nearest')
            imag=gaussian_filter((cross*fit).imag,sigma,mode='nearest')
            phase[t,j]=np.angle(real+1j*imag)
    phase=np.unwrap(phase,axis=2)
    observed=np.angle(normalized*reference.conj()[None])[:,:,check]
    proposal=phase[:,:,check]*options.phase_fraction
    weights=abs(normalized*reference.conj()[None])[:,:,check]
    validation_gain=np.sum(weights*(np.cos(observed-proposal)-np.cos(observed)),axis=2)/np.maximum(
        weights.sum(axis=2),1e-30)
    applied=np.where(validation_gain>options.minimum_validation_phase_gain,
                     options.phase_fraction,0.)[:,:,None,None]*phase[:,:,CORE_ROWS]
    audit={'model':'masked, spatially smoothed local tissue phase; no IQ resampling',
           'options':asdict(options),'fit_pixels':int(fit.sum()),'validation_pixels':int(check.sum()),
           'fit_mask_sha256':array_hash(fit),'validation_mask_sha256':array_hash(check),
           'wide_calibration_sha256':calibration['calibration_sha256'],
           'phase_sha256':array_hash(applied),'accepted_view_frames':int(np.sum(
               validation_gain>options.minimum_validation_phase_gain)),
           'view_frames':int(validation_gain.size),
           'median_validation_phase_gain':float(np.median(validation_gain)),
           'minimum_validation_phase_gain_observed':float(validation_gain.min()),
           'mean_applied_phase_abs_rad':float(np.mean(abs(applied))),
           'uses_current_frame_and_prior_calibration_only':True}
    return applied,audit


def fit_reference_phase(wide,wide_roi,options=None):
    options=options or ReferenceTransportOptions()
    if (not 0<options.transport_strength<=1.2 or
            not 0<options.reference_amplitude_fraction<1 or
            options.minimum_validation_phase_gain<0 or
            options.minimum_valid_check_pixels<1 or
            not np.array_equal(wide.training_roi,wide_roi)):
        raise ValueError('Invalid reference transport support/options')
    fit,check=_fit_check(wide_roi)
    displacement,shift_audit=fit_shifts(wide,wide_roi,ShiftOptions())
    shift_audit={key:value for key,value in shift_audit.items()
                 if key not in ('model','edge_assumption')}
    shift_audit['model']='disjoint-pixel observed-magnitude translation only'
    modulation,_,calibration,reference_views=bounded_calibrate(wide)
    normalized=wide.views/modulation
    phase=np.zeros(wide.views.shape,float)
    gains=np.zeros(wide.views.shape[:2],float)
    accepted=np.zeros_like(gains,dtype=bool)
    inbounds=np.zeros_like(gains)
    reliable=np.zeros_like(gains)
    for t in range(len(normalized)):
        for j in range(normalized.shape[1]):
            reference=reference_views[j]
            moved,valid=transported_reference(reference,displacement[t,j])
            threshold=options.reference_amplitude_fraction*np.median(abs(reference[wide_roi]))
            support=valid&(abs(reference)>=threshold)&(abs(moved)>=threshold)
            inbounds[t,j]=np.mean(valid[CORE_ROWS])
            reliable[t,j]=np.mean(support[CORE_ROWS])
            if (np.sum(fit&support)<80 or
                    np.sum(check&support)<options.minimum_valid_check_pixels):
                continue
            proposal=np.where(support,
                options.transport_strength*np.angle(moved*reference.conj()),0.)
            observed=np.angle(normalized[t,j]*reference.conj())
            weights=abs(normalized[t,j]*reference.conj())*(check&support)
            gain=np.sum(weights*(np.cos(observed-proposal)-np.cos(observed)))/max(
                np.sum(weights),1e-30)
            gains[t,j]=gain
            if gain>options.minimum_validation_phase_gain:
                phase[t,j]=proposal
                accepted[t,j]=True
    applied=phase[:,:,CORE_ROWS]
    audit={'model':'in-bounds complex calibration-reference transport; phase-only core',
           'options':asdict(options),'fit_pixels':int(fit.sum()),
           'validation_pixels':int(check.sum()),
           'wide_calibration_sha256':calibration['calibration_sha256'],
           'magnitude_translation_audit':shift_audit,
           'phase_sha256':array_hash(applied),
           'accepted_view_frames':int(accepted.sum()),
           'view_frames':int(accepted.size),
           'median_validation_gain':float(np.median(gains)),
           'mean_core_inbounds_fraction':float(np.mean(inbounds)),
           'mean_core_reliable_reference_fraction':float(np.mean(reliable)),
           'boundary_rule':'no cyclic wrap; interpolate reference only; zero phase on invalid or weak support',
           'uses_observed_IQ_and_prior_calibration_only':True}
    return applied,audit


def _check_inputs(core,wide,wide_roi,core_roi):
    if (not np.array_equal(wide.training_roi,wide_roi) or
            not np.array_equal(core.training_roi,core_roi) or
            not np.array_equal(core.calibration_views,wide.calibration_views[:,:,CORE_ROWS]) or
            not np.array_equal(core.views,wide.views[:,:,CORE_ROWS]) or
            core.geometry!=wide.geometry or
            np.any(core_roi&~wide_roi[CORE_ROWS])):
        raise ValueError('Physical core and wide observations/masks differ')


class LocalContext(Phase3Context):
    def __init__(self,core,wide,wide_roi,core_roi):
        _check_inputs(core,wide,wide_roi,core_roi)
        super().__init__(core,packet=16,score_start=64,motion='none')
        phase,audit=fit_local_phase(wide,wide_roi,
            LocalPhaseOptions(vessel_guard_mm=5.,phase_fraction=.5,
                              spatial_sigma_pixels=4.,minimum_validation_phase_gain=.01))
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion='local_phase';self.motion_summary=audit;self.motion_diagnostics=audit
        self.calibration={**self.calibration,'local_phase':audit}
        self.aligned=self.align(core.views);self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}


class ReferenceContext(Phase3Context):
    def __init__(self,core,wide,wide_roi,core_roi):
        _check_inputs(core,wide,wide_roi,core_roi)
        super().__init__(core,packet=16,score_start=64,motion='none')
        phase,audit=fit_reference_phase(wide,wide_roi,
            ReferenceTransportOptions(transport_strength=1.1,
                reference_amplitude_fraction=.2,
                minimum_validation_phase_gain=.01,minimum_valid_check_pixels=30))
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion='reference_transport';self.motion_summary=audit;self.motion_diagnostics=audit
        self.calibration={**self.calibration,'reference_transport':audit}
        self.aligned=self.align(core.views);self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}


class HybridContext(Phase3Context):
    """Selected local phase and open-boundary unitary shift on the original core."""
    def __init__(self,core,wide,wide_roi,core_roi):
        _check_inputs(core,wide,wide_roi,core_roi)
        super().__init__(core,packet=16,score_start=64,motion='none')
        shifts,shift_audit=fit_shifts(wide,wide_roi,ShiftOptions())
        phase,phase_audit=fit_local_phase(wide,wide_roi,
            LocalPhaseOptions(vessel_guard_mm=5.,phase_fraction=.75,
                              spatial_sigma_pixels=4.,minimum_validation_phase_gain=.01))
        self.shifts=.5*shifts
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion='phase10_hybrid'
        self.motion_summary={
            'model':'accepted local phase plus nonperiodic unitary IQ shift',
            'phase_fraction':.75,'shift_scale':.5,
            'shift_fit':{**shift_audit,'application':'nonperiodic skew derivative',
                         'inherited_periodic_application_not_used':True},
            'shift_sha256':array_hash(self.shifts),'phase_fit':phase_audit}
        self.motion_diagnostics=self.motion_summary
        self.calibration={**self.calibration,'phase10_hybrid':self.motion_summary}
        self.aligned=self.align(core.views);self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}

    def align(self,other):
        value=other/self.modulation
        if not hasattr(self,'shifts'):
            return value
        return np.stack([np.stack([open_unitary_shift(frame,d)
                                   for frame,d in zip(row,ds)])
                         for row,ds in zip(value,self.shifts)])

    def restore_block(self,aligned,start):
        if not hasattr(self,'shifts'):
            return super().restore_block(aligned,start)
        native=np.stack([np.stack([open_unitary_shift(frame,-d)
                                   for frame,d in zip(row,ds)])
                         for row,ds in zip(aligned,self.shifts[start:start+len(aligned)])])
        return native*self.modulation[start:start+len(aligned)]


def choose_mode(core,wide,wide_roi,core_roi):
    """Apply the same pre-score phase loss and prior-reference safeguard."""
    _check_inputs(core,wide,wide_roi,core_roi)
    _,_,_,refs=bounded_calibrate(core)
    check=(wide_roi&((np.indices(wide_roi.shape)[0]+2*np.indices(wide_roi.shape)[1])%3==0))[CORE_ROWS]
    local=LocalContext(core,wide,wide_roi,core_roi)
    reference=ReferenceContext(core,wide,wide_roi,core_roi)
    local_loss=phase_loss(local,refs,check,64)
    reference_loss=phase_loss(reference,refs,check,64)
    excess=reference_loss-local_loss
    requested='reference' if local_loss>.28 and excess<.18 else 'hybrid'
    _,_,_,references=bounded_calibrate(wide)
    fractions=[]
    for ref in references:
        threshold=.2*np.median(abs(ref[wide_roi]))
        fractions.append(float(np.mean((abs(ref[CORE_ROWS])>=threshold)[core_roi])))
    minimum=min(fractions)
    selected='reference' if requested=='reference' and minimum>=.9 else 'hybrid'
    return selected,{'guard_mm':5.,'rule':'reference if local_loss > 0.28 and reference excess < 0.18; otherwise hybrid',
                     'pre_score_frames':64,'local_phase_loss':local_loss,
                     'reference_phase_loss':reference_loss,'reference_excess_loss':excess,
                     'check_pixels':int(check.sum()),'source_truth_used':False,
                     'requested_mode':requested,
                     'minimum_core_reference_amplitude_support':minimum,
                     'reference_support_by_view':fractions,
                     'minimum_reference_support_for_reference_mode':.9,
                     'selection_uses_prior_calibration_and_pre_score_mixture_only':True,
                     'selected_mode':selected}


def build_context(mode,core,wide,wide_roi,core_roi):
    _check_inputs(core,wide,wide_roi,core_roi)
    if mode=='none':return Phase3Context(core,packet=16,score_start=64,motion='none')
    if mode=='hybrid':return HybridContext(core,wide,wide_roi,core_roi)
    if mode=='reference':return ReferenceContext(core,wide,wide_roi,core_roi)
    raise ValueError('Unknown phase-10 mode')
