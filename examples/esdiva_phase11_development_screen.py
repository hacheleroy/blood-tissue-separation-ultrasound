"""Screen hybrid strengths on the two now-seen phase-10 stress clouds.

The six source-resolved scenes are development data. A candidate selected from
this grid must be frozen and evaluated on a subsequently declared RF bank.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, relative_error
from esdiva_phase10_support import Phase10HybridContext, phase10_support
from esdiva_phase11_adaptive import choose_hybrid_strength, build_phase11_context


ROOT = Path(__file__).resolve().parents[1]
PHASE_FRACTIONS = (.5, .6, .75, 1.)
SHIFT_SCALES = (.25, .5, .75)
SEEDS = (3607, 3619)
SCENES = ('rigid_fast', 'rigid_slow', 'smooth_slow')
METHODS = ('roundtrip', 'svd_cut1', 'axis1_cut1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(a, b):
    return float(np.mean(abs(a)**2) / max(np.mean(abs(b)**2), 1e-30))


def screen(acquisitions, output):
    acquisitions, output = Path(acquisitions), Path(output)
    if output.exists():
        raise FileExistsError(output)
    config = json.loads((ROOT/'configs/esdiva_phase10_stress.json').read_text())
    if (list(config['independent_test_seeds']) != list(SEEDS)
            or list(config['scenes']) != list(SCENES)
            or sha(ROOT/'examples/esdiva_phase10_support.py') != config['support_source_sha256']):
        raise ValueError('Development bank or fixed support differs')
    cases = []
    for seed in SEEDS:
        for scene in SCENES:
            folder = acquisitions/str(seed)/scene
            report = json.loads((folder/'report.json').read_text())
            audit = json.loads((folder/'broad_iq.json').read_text())
            score = json.loads((folder/'score.json').read_text())
            setting = config['scenes'][scene]
            if (report['seed'] != seed or report['motion'] != setting['motion']
                    or report['flow_scale'] != setting['flow_scale']
                    or report['script_sha256'] != config['acquisition_source_sha256']
                    or audit['source_report_sha256'] != sha(folder/'report.json')
                    or audit['script_sha256'] != config['rebeam_source_sha256']
                    or score['acquisition_sha256'] != sha(folder/'report.json')
                    or score['broad_audit_sha256'] != sha(folder/'broad_iq.json')
                    or score['config_sha256'] != sha(ROOT/'configs/esdiva_phase10_stress.json')
                    or score['selected_mode'] != 'hybrid'
                    or not score['calibration_support']['pass']):
                raise ValueError(f'Changed development input: {seed}/{scene}')
            with np.load(folder/'iq_sources.npz') as old, np.load(folder/'broad_iq.npz') as wide:
                core = {key:old[key] for key in ('views','calibration_views','tissue','blood',
                                                 'noise','z_mm','x_mm','training_roi')}
                large = {key:wide[key] for key in ('views','calibration_views','tissue','blood',
                                                   'noise','z_mm','x_mm','core_rows')}
            for key, label in (('views','mixture_sha256'),('calibration_views','calibration_sha256'),
                               ('tissue','tissue_sha256'),('blood','blood_sha256'),
                               ('noise','noise_sha256')):
                if digest(core[key]) != report[label]:
                    raise ValueError(f'Changed core IQ {seed}/{scene}/{key}')
            for key in ('views','calibration_views','tissue','blood','noise'):
                if (digest(large[key]) != audit['broad_iq_sha256'][key]
                        or not np.array_equal(large[key][:,:,large['core_rows']], core[key])):
                    raise ValueError(f'Changed wide IQ {seed}/{scene}/{key}')
            wide_roi,core_roi,gate=phase10_support(
                large['calibration_views'],large['z_mm'],large['x_mm'],large['core_rows'])
            if (not gate['pass'] or not np.array_equal(core_roi,core['training_roi'])
                    or ratio(core['blood'][:,:,core_roi],core['tissue'][:,:,core_roi])>.1
                    or ratio(large['blood'][:,:,wide_roi],large['tissue'][:,:,wide_roi])>.1):
                raise ValueError(f'Ineligible development input {seed}/{scene}')
            geo=report['geometry']
            inputs={label:(PublicInput(c,core['calibration_views'],core_roi,geo),
                           PublicInput(w,large['calibration_views'],wide_roi,geo))
                    for label,(c,w) in {
                        'mixture':(core['views'],large['views']),
                        'background':(core['tissue']+core['noise'],
                                      large['tissue']+large['noise'])}.items()}
            cases.append((seed,scene,folder,core,large,inputs,score))
    rows=[]
    for phase_fraction in PHASE_FRACTIONS:
        for shift_scale in SHIFT_SCALES:
            for seed,scene,folder,core,large,inputs,prior in cases:
                contexts={label:Phase10HybridContext(ci,wi,large['core_rows'],
                                  large['z_mm'],large['x_mm'],phase_fraction=phase_fraction,
                                  shift_scale=shift_scale)
                          for label,(ci,wi) in inputs.items()}
                if (contexts['mixture'].calibration['calibration_sha256'] !=
                        contexts['background'].calibration['calibration_sha256']):
                    raise ValueError('Paired calibration differs')
                methods={}
                for name in METHODS:
                    spec=config['methods'][name]
                    mixed,bg=(fit_shared(contexts[label],spec)
                              for label in ('mixture','background'))
                    blood=relative_error(mixed.views-bg.views,core['blood'][64:])
                    tissue=ratio(bg.views,core['tissue'][64:])
                    baseline=prior['variants']['none']['methods'][name]
                    methods[name]={'blood_nrmse':blood,
                                   'blood_nrmse_increase':
                                       blood-baseline['paired_added_blood_iq_nrmse'],
                                   'tissue_leakage_ratio':
                                       tissue/baseline['refitted_background_output_to_input_tissue_power']}
                rows.append({'seed':seed,'scene':scene,'phase_fraction':phase_fraction,
                             'shift_scale':shift_scale,'methods':methods,
                             'phase_accepted_view_frames':
                                 contexts['mixture'].motion_summary['phase_fit']['accepted_view_frames'],
                             'shift_accepted_view_frames':
                                 contexts['mixture'].motion_summary['shift_fit']['accepted_view_frames']})
            print(f'completed phase={phase_fraction} shift={shift_scale}',flush=True)
    limits=config['comparison']
    candidates=[]
    for phase_fraction in PHASE_FRACTIONS:
        for shift_scale in SHIFT_SCALES:
            sub=[row for row in rows if row['phase_fraction']==phase_fraction
                 and row['shift_scale']==shift_scale]
            worst_blood=max(row['methods'][m]['blood_nrmse_increase'] for row in sub
                            for m in ('svd_cut1','axis1_cut1'))
            worst_tissue=max(row['methods'][m]['tissue_leakage_ratio'] for row in sub
                             for m in ('svd_cut1','axis1_cut1'))
            worst_roundtrip=max(row['methods']['roundtrip']['blood_nrmse']
                                for row in sub)
            candidates.append({'phase_fraction':phase_fraction,'shift_scale':shift_scale,
                               'worst_blood_nrmse_increase':worst_blood,
                               'worst_tissue_leakage_ratio':worst_tissue,
                               'worst_roundtrip_increase':worst_roundtrip,
                               'development_gate_pass':bool(
                                   worst_blood<=limits['blood_nrmse_increase_max_motion']
                                   and worst_tissue<=limits['moving_tissue_leakage_ratio_max']
                                   and worst_roundtrip<=limits['roundtrip_nrmse_max'])})
    admitted=[x for x in candidates if x['development_gate_pass']]
    selected=min(admitted,key=lambda x:(x['worst_tissue_leakage_ratio'],
                                        x['worst_blood_nrmse_increase'],
                                        x['phase_fraction'],x['shift_scale'])) if admitted else None
    adaptive_rows=[]
    for seed,scene,folder,core,large,inputs,prior in cases:
        strength=choose_hybrid_strength(inputs['mixture'][1],large['z_mm'],large['x_mm'])
        contexts={label:build_phase11_context('hybrid',ci,wi,large['core_rows'],
                         large['z_mm'],large['x_mm'],strength)
                  for label,(ci,wi) in inputs.items()}
        methods={}
        for name in METHODS:
            mixed,bg=(fit_shared(contexts[label],config['methods'][name])
                      for label in ('mixture','background'))
            blood=relative_error(mixed.views-bg.views,core['blood'][64:])
            tissue=ratio(bg.views,core['tissue'][64:])
            baseline=prior['variants']['none']['methods'][name]
            methods[name]={'blood_nrmse':blood,
                           'blood_nrmse_increase':blood-baseline['paired_added_blood_iq_nrmse'],
                           'tissue_leakage_ratio':
                               tissue/baseline['refitted_background_output_to_input_tissue_power']}
        adaptive_rows.append({'seed':seed,'scene':scene,'strength':strength,'methods':methods})
    adaptive_worst_blood=max(row['methods'][m]['blood_nrmse_increase'] for row in adaptive_rows
                             for m in ('svd_cut1','axis1_cut1'))
    adaptive_worst_tissue=max(row['methods'][m]['tissue_leakage_ratio'] for row in adaptive_rows
                              for m in ('svd_cut1','axis1_cut1'))
    adaptive_worst_roundtrip=max(row['methods']['roundtrip']['blood_nrmse']
                                 for row in adaptive_rows)
    adaptive_summary={'worst_blood_nrmse_increase':adaptive_worst_blood,
                      'worst_tissue_leakage_ratio':adaptive_worst_tissue,
                      'worst_roundtrip_nrmse':adaptive_worst_roundtrip,
                      'development_gate_pass':bool(
                          adaptive_worst_blood<=limits['blood_nrmse_increase_max_motion']
                          and adaptive_worst_tissue<=limits['moving_tissue_leakage_ratio_max']
                          and adaptive_worst_roundtrip<=limits['roundtrip_nrmse_max'])}
    result={'status':'development grid complete; new holdout unopened',
            'source_sha256':sha(__file__),
            'adaptive_source_sha256':sha(ROOT/'examples/esdiva_phase11_adaptive.py'),
            'stress_config_sha256':sha(ROOT/'configs/esdiva_phase10_stress.json'),
            'parent_support_sha256':config['support_source_sha256'],
            'seen_seeds':list(SEEDS),'scenes':list(SCENES),
            'phase_fraction_grid':list(PHASE_FRACTIONS),'shift_scale_grid':list(SHIFT_SCALES),
            'selection_rule':'Among gate-passing settings, minimize worst tissue ratio, then worst blood increase, phase fraction, shift scale. No selection if none passes.',
            'candidate_rows':candidates,'all_scene_rows':rows,
            'selected_development_setting':selected,
            'adaptive_rows':adaptive_rows,'adaptive_summary':adaptive_summary,
            'no_new_scatterer_seed_opened':True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--acquisitions',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    result=screen(a.acquisitions,a.output)
    print(json.dumps({'status':result['status'],
                      'selected_development_setting':result['selected_development_setting'],
                      'adaptive_summary':result['adaptive_summary']}))
