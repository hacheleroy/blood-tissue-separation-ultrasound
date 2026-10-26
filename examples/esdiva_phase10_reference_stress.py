"""Pre-score counterfactual of the frozen phase-10 reference context on rigid RF.

This never affects the primary mode selector or independent stress decision.
Source truth enters only paired evaluation, as in the frozen phase-10 scorer.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, relative_error
from esdiva_phase10_support import build_phase10_context, phase10_support


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(a, b):
    return float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2), 1e-30))


def evaluate(acquisition, broad, primary_score, config_path, output):
    acquisition, broad, primary_score, config_path, output = map(Path,
        (acquisition, broad, primary_score, config_path, output))
    if output.exists():
        raise FileExistsError(output)
    config = json.loads(config_path.read_text())
    stress = json.loads((ROOT/'configs/esdiva_phase10_stress.json').read_text())
    report = json.loads((acquisition/'report.json').read_text())
    audit = json.loads(broad.with_suffix('.json').read_text())
    score = json.loads(primary_score.read_text())
    scene = next((name for name, setting in stress['scenes'].items()
                  if setting['motion'] == report['motion']
                  and setting['flow_scale'] == report['flow_scale']), None)
    if (sha(ROOT/'configs/esdiva_phase10_stress.json') != config['primary_config_sha256']
            or config['source_sha256'] != sha(__file__)
            or config['rigid_scenes'] != ['rigid_fast', 'rigid_slow']
            or report['seed'] not in stress['independent_test_seeds']
            or scene not in config['rigid_scenes']
            or report['script_sha256'] != stress['acquisition_source_sha256']
            or audit['source_report_sha256'] != sha(acquisition/'report.json')
            or audit['script_sha256'] != stress['rebeam_source_sha256']
            or score['status'] != 'phase-10 paired source-resolved independent stress diagnostic'
            or score['acquisition_sha256'] != sha(acquisition/'report.json')
            or score['broad_audit_sha256'] != sha(broad.with_suffix('.json'))
            or score['config_sha256'] != sha(ROOT/'configs/esdiva_phase10_stress.json')
            or score['script_sha256'] != stress['score_script_sha256']
            or not score['calibration_support']['pass']
            or not report['source_eligibility_pass']):
        raise ValueError('Counterfactual input or frozen parent identity changed')
    with np.load(acquisition/'iq_sources.npz') as old, np.load(broad) as wide:
        core = {key:old[key] for key in ('views','calibration_views','tissue','blood','noise',
                                        'z_mm','x_mm','training_roi')}
        large = {key:wide[key] for key in ('views','calibration_views','tissue','blood',
                                          'noise','z_mm','x_mm','core_rows')}
    if (any(digest(core[key]) != report[label] for key,label in
            (('views','mixture_sha256'),('calibration_views','calibration_sha256'),
             ('tissue','tissue_sha256'),('blood','blood_sha256'),('noise','noise_sha256')))
            or any(digest(large[key]) != audit['broad_iq_sha256'][key]
                   or not np.array_equal(large[key][:,:,large['core_rows']],core[key])
                   for key in ('views','calibration_views','tissue','blood','noise'))):
        raise ValueError('RF-replayed core or source IQ mismatch')
    wide_roi, core_roi, gate = phase10_support(large['calibration_views'],
                            large['z_mm'],large['x_mm'],large['core_rows'])
    if (not gate['pass'] or not np.array_equal(core_roi,core['training_roi'])
            or ratio(core['blood'][:,:,core_roi],core['tissue'][:,:,core_roi]) > .1
            or ratio(large['blood'][:,:,wide_roi],large['tissue'][:,:,wide_roi]) > .1):
        raise ValueError('Counterfactual scene is ineligible')
    geo = report['geometry']
    inputs = {}
    for label, c, w in (
        ('mixture', core['views'], large['views']),
        ('background', core['tissue']+core['noise'], large['tissue']+large['noise'])):
        inputs[label] = (PublicInput(c,core['calibration_views'],core_roi,geo),
                         PublicInput(w,large['calibration_views'],wide_roi,geo))
    contexts = {label:build_phase10_context('reference',ci,wi,large['core_rows'],
                                             large['z_mm'],large['x_mm'])
                for label,(ci,wi) in inputs.items()}
    if contexts['mixture'].calibration['calibration_sha256'] != \
            contexts['background'].calibration['calibration_sha256']:
        raise ValueError('Paired calibration differs')
    methods = {}
    for method in ('roundtrip','svd_cut1','axis1_cut1'):
        mixed = fit_shared(contexts['mixture'],stress['methods'][method])
        background = fit_shared(contexts['background'],stress['methods'][method])
        blood_nrmse = relative_error(mixed.views-background.views,core['blood'][64:])
        tissue_leakage = ratio(background.views,core['tissue'][64:])
        baseline = score['variants']['none']['methods'][method]
        methods[method] = {
            'paired_added_blood_iq_nrmse':blood_nrmse,
            'refitted_background_output_to_input_tissue_power':tissue_leakage,
            'blood_nrmse_increase_over_none':
                blood_nrmse-baseline['paired_added_blood_iq_nrmse'],
            'tissue_leakage_ratio_to_none':
                tissue_leakage/baseline['refitted_background_output_to_input_tissue_power']}
    result = {'status':'reference counterfactual; not the selected-mode decision',
              'seed':report['seed'],'scene':scene,
              'naturally_selected_reference':score['selected_mode']=='reference',
              'requested_mode':score['selection']['requested_mode'],
              'minimum_core_reference_amplitude_support':
                  score['selection']['minimum_core_reference_amplitude_support'],
              'primary_score_sha256':sha(primary_score),'broad_audit_sha256':sha(broad.with_suffix('.json')),
              'config_sha256':sha(config_path),'source_sha256':sha(__file__),
              'fit':contexts['mixture'].motion_summary,
              'background_fit':contexts['background'].motion_summary,
              'methods':methods,'uses_source_truth_for_fit_or_selection':False,
              'interpretation':'Forced reference is a branch sensitivity diagnostic; it cannot repair missing natural coverage.'}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--acquisition',type=Path,required=True)
    p.add_argument('--broad',type=Path,required=True)
    p.add_argument('--primary-score',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    r=evaluate(a.acquisition,a.broad,a.primary_score,a.config,a.output)
    print(json.dumps({k:r[k] for k in ('seed','scene','naturally_selected_reference')}))
