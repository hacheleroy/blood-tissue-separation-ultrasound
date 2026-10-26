"""Reproduce four pre-score reference counterfactual diagnostics."""

import argparse
import json
from pathlib import Path

from review_esdiva_phase10_stress import ROOT, sha


def review(folder, output):
    folder, output = Path(folder), Path(output)
    config_path = ROOT/'configs/esdiva_phase10_reference_ablation.json'
    config = json.loads(config_path.read_text())
    if (config['source_sha256'] != sha(ROOT/'examples/esdiva_phase10_reference_stress.py')
            or config['primary_config_sha256'] != sha(ROOT/'configs/esdiva_phase10_stress.json')
            or config['independent_test_seeds'] != [3607,3619]
            or config['rigid_scenes'] != ['rigid_fast','rigid_slow']):
        raise ValueError('Reference sensitivity source or scene contract changed')
    missing, rows = [], []
    for seed in config['independent_test_seeds']:
        for scene in config['rigid_scenes']:
            base = folder/str(seed)/scene
            path = base/'reference_ablation.json'
            primary_path = base/'score.json'
            if primary_path.is_file():
                primary = json.loads(primary_path.read_text())
                if primary['status'] == 'input gate failed; no filter scores':
                    if path.exists():
                        raise ValueError(f'Ineligible scene has forced reference scores: {seed}/{scene}')
                    rows.append({'seed':seed,'scene':scene,'status':'input ineligible; reference not fitted',
                                 'primary_score_sha256':sha(primary_path)})
                    continue
            if not path.is_file():
                missing.append(f'{seed}/{scene}')
                continue
            row = json.loads(path.read_text())
            primary = json.loads(primary_path.read_text())
            if (row['seed'] != seed or row['scene'] != scene
                    or row['config_sha256'] != sha(config_path)
                    or row['source_sha256'] != config['source_sha256']
                    or row['primary_score_sha256'] != sha(primary_path)
                    or row['broad_audit_sha256'] != sha(base/'wide_audit.json')
                    or row['naturally_selected_reference'] != (primary['selected_mode']=='reference')
                    or row['uses_source_truth_for_fit_or_selection']
                    or row['fit']['model'] != 'in-bounds complex calibration-reference transport; phase-only core'
                    or not row['fit']['uses_observed_IQ_and_prior_calibration_only']):
                raise ValueError(f'Changed reference branch input/source in {seed}/{scene}')
            for method in ('roundtrip','svd_cut1','axis1_cut1'):
                actual = row['methods'][method]
                baseline = primary['variants']['none']['methods'][method]
                if (abs(actual['blood_nrmse_increase_over_none'] -
                        (actual['paired_added_blood_iq_nrmse']-
                         baseline['paired_added_blood_iq_nrmse'])) > 1e-12
                        or abs(actual['tissue_leakage_ratio_to_none']-
                               actual['refitted_background_output_to_input_tissue_power']/
                               baseline['refitted_background_output_to_input_tissue_power']) > 1e-12):
                    raise ValueError(f'Reference metric arithmetic changed in {seed}/{scene}')
                if row['naturally_selected_reference']:
                    selected = primary['variants']['selected']['methods'][method]
                    if (abs(actual['paired_added_blood_iq_nrmse']-
                            selected['paired_added_blood_iq_nrmse']) > 1e-12
                            or abs(actual['refitted_background_output_to_input_tissue_power']-
                                   selected['refitted_background_output_to_input_tissue_power']) > 1e-12):
                        raise ValueError('Natural reference and counterfactual disagree')
            rows.append({'seed':seed,'scene':scene,
                         'naturally_selected_reference':row['naturally_selected_reference'],
                         'minimum_core_reference_amplitude_support':
                             row['minimum_core_reference_amplitude_support'],
                         'accepted_view_frames':row['fit']['accepted_view_frames'],
                         'view_frames':row['fit']['view_frames'],
                         'mean_core_reliable_reference_fraction':
                             row['fit']['mean_core_reliable_reference_fraction'],
                         'methods':{key:row['methods'][key] for key in ('roundtrip','svd_cut1','axis1_cut1')}})
    result = {'status':('incomplete counterfactual reference diagnostic' if missing else
                         'four-scene reference branch eligibility and diagnostics recorded'),
              'config_sha256':sha(config_path),'missing':missing,'rows':rows,
              'interpretation':'Forced reference rows are engineering sensitivity only; natural selection is assessed by the separate primary decision.'}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(review(a.folder,a.output)['status']))
