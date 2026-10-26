"""Candidate phase-11 hybrid strengths chosen from pre-score magnitude validation.

This revision changes only the hybrid phase/shift strengths. Phase-10 source
support, calibration, reference selection, reference transport, and numerical
filter settings remain as frozen. Strengths are selected once from the first
64 mixture frames and shared with the evaluator's separate background fit.
"""

import numpy as np

from btbench.shared import PublicInput
from esdiva_phase6_unit_shift import fit_shifts
from esdiva_phase10_support import Phase10HybridContext, build_phase10_context


def choose_hybrid_strength(wide, z_mm, x_mm):
    if len(wide.views) != 80 or len(wide.calibration_views) != 16:
        raise ValueError('Need 64 pre-score frames and a disjoint 16-frame calibration')
    pre = PublicInput(wide.views[:64], wide.calibration_views,
                      wide.training_roi, wide.geometry)
    _, audit = fit_shifts(pre, z_mm, x_mm)
    gain = audit['median_validation_gain']
    if not np.isfinite(gain):
        raise ValueError('Nonfinite pre-score magnitude validation gain')
    weight = float(np.clip((gain - .25)/.20, 0., 1.))
    return {'pre_score_frames':64, 'pre_score_shift_median_validation_gain':gain,
            'pre_score_shift_audit':audit,
            'weight':weight,'phase_fraction':1.0-.5*weight,
            'shift_scale':.25+.5*weight,
            'rule':'clip((median first-64 magnitude validation gain - 0.25)/0.20,0,1); phase=1-0.5w, shift=0.25+0.5w',
            'uses_source_truth_or_scored_frames':False}


def build_phase11_context(mode, core, wide, core_rows, z_mm, x_mm, strength=None):
    if mode == 'hybrid':
        if strength is None:
            raise ValueError('Hybrid needs frozen mixture-only strength')
        return Phase10HybridContext(core, wide, core_rows, z_mm, x_mm,
                phase_fraction=strength['phase_fraction'],shift_scale=strength['shift_scale'])
    if mode in ('none','reference'):
        return build_phase10_context(mode, core, wide, core_rows, z_mm, x_mm)
    raise ValueError('Unknown phase-11 mode')
