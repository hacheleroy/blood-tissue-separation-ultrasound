"""Development-only core support and measured reference-mode eligibility.

The wide motion fit retains the 5 mm anatomical guard. The scored eight-row
core uses a 5.5 mm guard, chosen on now-seen clouds without isolated sources
as an input. Calibration and the 64 pre-score mixture frames alone determine
eligibility and mode. Blood truth and known motion are evaluator-only.
"""

import numpy as np
from scipy.ndimage import gaussian_filter

from btbench.shared import array_hash
from esdiva_angle_calibration import angular_diagnostics
from esdiva_phase3_bounded_context import Phase3Context, bounded_calibrate
from esdiva_phase5_local_phase import LocalPhaseOptions, fit_local_phase, tissue_mask
from esdiva_phase6_unit_shift import ShiftOptions, fit_shifts
from esdiva_phase7_reference_transport import ReferenceTransportContext
from esdiva_phase8_selector import select_mode
from esdiva_phase9_support import Phase9HybridContext


def phase10_support(calibration, z_mm, x_mm, core_rows):
    """Train-only angular phase gate on wide-5 and core-5.5 tissue support."""
    wide=tissue_mask(z_mm,x_mm,5.)
    core=tissue_mask(z_mm,x_mm,5.5)[core_rows]
    if (calibration.shape!=(16,3,32,16) or
            not np.array_equal(core_rows,np.arange(12,20)) or
            wide.sum()!=192 or core.sum()!=40):
        raise ValueError('Changed phase-10 grid or anatomical support')
    registered=calibration.copy()
    phase=np.zeros((3,32,16),float)
    denominator=gaussian_filter(wide.astype(float),.4,mode='constant')
    for j in (1,2):
        cross=np.mean(calibration[:8,j]*calibration[:8,0].conj(),axis=0)
        phase[j]=np.angle(gaussian_filter(cross*wide,.4,mode='constant')/
                          np.maximum(denominator,1e-15))
        registered[:,j]*=np.exp(-1j*phase[j])[None]
    result={'model':'train-only angular phase map; gate only',
            'wide_guard_mm':5.,'core_guard_mm':5.5,
            'wide_pixels':int(wide.sum()),'core_pixels':int(core.sum()),
            'phase_sha256':array_hash(phase),
            'calibration_sha256':array_hash(calibration),
            'fit_frames':8,'validation_frames':8,
            'uses_source_truth_or_evaluation_IQ':False}
    for name,array,roi in (('wide',registered,wide),
                           ('core',registered[:,:,core_rows],core)):
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
    return wide,core,result


class Phase10HybridContext(Phase9HybridContext):
    """Same reversible hybrid on a source-independent narrower core mask."""

    def __init__(self,core,wide,core_rows,z_mm,x_mm,score_start=64,
                 phase_fraction=.75,shift_scale=.5):
        if not (0<phase_fraction<=1 and 0<shift_scale<=1):
            raise ValueError('Invalid phase-10 hybrid strength')
        if (not np.array_equal(wide.training_roi,tissue_mask(z_mm,x_mm,5.)) or
                not np.array_equal(core.training_roi,
                                   tissue_mask(z_mm,x_mm,5.5)[core_rows]) or
                not np.array_equal(core.calibration_views,
                                   wide.calibration_views[:,:,core_rows]) or
                not np.array_equal(core.views,wide.views[:,:,core_rows]) or
                core.geometry!=wide.geometry):
            raise ValueError('Phase-10 core and wide input mismatch')
        Phase3Context.__init__(self,core,packet=16,score_start=score_start,
                               motion='none')
        shifts,shift_audit=fit_shifts(wide,z_mm,x_mm,ShiftOptions())
        phase,phase_audit=fit_local_phase(
            wide,core_rows,z_mm,x_mm,
            LocalPhaseOptions(vessel_guard_mm=5.,phase_fraction=phase_fraction,
                              spatial_sigma_pixels=4.,
                              minimum_validation_phase_gain=.01))
        self.shifts=shift_scale*shifts
        self.modulation=self.modulation*np.exp(1j*phase)
        self.motion='phase10_hybrid'
        self.motion_summary={
            'model':'accepted local phase plus nonperiodic unitary IQ shift',
            'phase_fraction':phase_fraction,'shift_scale':shift_scale,
            'shift_fit':{**shift_audit,'application':'nonperiodic skew derivative',
                         'inherited_periodic_application_not_used':True},
            'shift_sha256':array_hash(self.shifts),'phase_fit':phase_audit}
        self.motion_diagnostics=self.motion_summary
        self.calibration={**self.calibration,'phase10_hybrid':self.motion_summary}
        self.aligned=self.align(core.views)
        self.x=self.aligned.reshape(len(core.views),-1).T
        self.cache={}


def choose_phase10_mode(core,wide,core_rows,z_mm,x_mm):
    """Require a strong calibration reference before allowing reference mode."""
    requested,audit=select_mode(core,wide,core_rows,z_mm,x_mm,score_start=64)
    if requested not in ('hybrid','reference'):
        raise ValueError('Phase-10 fixed wide guard requires hybrid/reference')
    _,_,_,references=bounded_calibrate(wide)
    core_roi=tissue_mask(z_mm,x_mm,5.5)[core_rows]
    fractions=[]
    for reference in references:
        threshold=.2*np.median(abs(reference[wide.training_roi]))
        fractions.append(float(np.mean(
            (abs(reference[core_rows])>=threshold)[core_roi])))
    minimum=min(fractions)
    selected='reference' if requested=='reference' and minimum>=.9 else 'hybrid'
    return selected,{**audit,'requested_mode':requested,
                     'minimum_core_reference_amplitude_support':minimum,
                     'reference_support_by_view':fractions,
                     'minimum_reference_support_for_reference_mode':.9,
                     'selection_uses_prior_calibration_and_pre_score_mixture_only':True,
                     'selected_mode':selected}


def build_phase10_context(mode,core,wide,core_rows,z_mm,x_mm):
    if mode=='none':
        return Phase3Context(core,packet=16,score_start=64,motion='none')
    if mode=='hybrid':
        return Phase10HybridContext(core,wide,core_rows,z_mm,x_mm)
    if mode=='reference':
        return ReferenceTransportContext(
            core,wide,core_rows,z_mm,x_mm,packet=16,score_start=64,
            transport_options={'transport_strength':1.1,
                               'reference_amplitude_fraction':.2,
                               'minimum_validation_phase_gain':.01,
                               'minimum_valid_check_pixels':30})
    raise ValueError('Unknown phase-10 mode')
