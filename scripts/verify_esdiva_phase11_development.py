"""Verify the phase-11 seen-cloud development grid and source fingerprints."""

import hashlib
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    folder=ROOT/'runs/esdiva_phase11_development'
    report=json.loads((folder/'report.json').read_text())
    manifest=json.loads((folder/'manifest.json').read_text())
    if (manifest['files_sha256'] != {key:sha(ROOT/key) for key in manifest['files_sha256']}
            or set(manifest['files_sha256']) != {
                'examples/esdiva_phase11_adaptive.py',
                'examples/esdiva_phase11_development_screen.py',
                'examples/esdiva_phase10_support.py',
                'configs/esdiva_phase10_stress.json',
                'runs/esdiva_phase10_stress/decision.json',
                'runs/esdiva_phase10_stress/manifest.json',
                'runs/esdiva_phase11_development/report.json'}
            or report['source_sha256'] != sha(ROOT/'examples/esdiva_phase11_development_screen.py')
            or report['adaptive_source_sha256'] != sha(ROOT/'examples/esdiva_phase11_adaptive.py')
            or report['stress_config_sha256'] != sha(ROOT/'configs/esdiva_phase10_stress.json')
            or report['parent_support_sha256'] != sha(ROOT/'examples/esdiva_phase10_support.py')
            or report['seen_seeds'] != [3607,3619]
            or report['scenes'] != ['rigid_fast','rigid_slow','smooth_slow']
            or report['phase_fraction_grid'] != [.5,.6,.75,1.]
            or report['shift_scale_grid'] != [.25,.5,.75]
            or report['selected_development_setting'] is not None
            or not report['no_new_scatterer_seed_opened']):
        raise ValueError('Changed phase-11 development source or bank')
    grid={(r['phase_fraction'],r['shift_scale'],r['seed'],r['scene'])
          for r in report['all_scene_rows']}
    if (len(grid)!=72 or len(report['candidate_rows'])!=12
            or any(c['development_gate_pass'] for c in report['candidate_rows'])):
        raise ValueError('Incomplete or relabeled fixed-strength grid')
    adaptive={(r['seed'],r['scene']) for r in report['adaptive_rows']}
    if adaptive!={(s,c) for s in (3607,3619) for c in report['scenes']}:
        raise ValueError('Incomplete adaptive development matrix')
    for row in report['adaptive_rows']:
        strength=row['strength']
        gain=strength['pre_score_shift_median_validation_gain']
        weight=max(0.,min(1.,(gain-.25)/.20))
        if (strength['pre_score_frames']!=64 or strength['uses_source_truth_or_scored_frames']
                or abs(strength['weight']-weight)>1e-12
                or abs(strength['phase_fraction']-(1.-.5*weight))>1e-12
                or abs(strength['shift_scale']-(.25+.5*weight))>1e-12):
            raise ValueError('Development strength changed')
    summary=report['adaptive_summary']
    blood=max(r['methods'][m]['blood_nrmse_increase'] for r in report['adaptive_rows']
              for m in ('svd_cut1','axis1_cut1'))
    tissue=max(r['methods'][m]['tissue_leakage_ratio'] for r in report['adaptive_rows']
               for m in ('svd_cut1','axis1_cut1'))
    roundtrip=max(r['methods']['roundtrip']['blood_nrmse'] for r in report['adaptive_rows'])
    if (not summary['development_gate_pass'] or blood!=summary['worst_blood_nrmse_increase']
            or tissue!=summary['worst_tissue_leakage_ratio']
            or roundtrip!=summary['worst_roundtrip_nrmse']
            or blood>.02 or tissue>.9 or roundtrip>1e-10):
        raise ValueError('Adaptive development gate does not reproduce')
    return {'status':'seen-cloud development archive verified','grid_settings':12,
            'adaptive_scenes':6,'worst_tissue_ratio':tissue,'worst_blood_delta':blood}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
