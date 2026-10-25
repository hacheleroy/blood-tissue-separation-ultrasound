"""Audit every declared phase-6 design trial; never use reserved seeds here."""

import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(root, acquisitions, configurations, scores, output):
    root,acquisitions,configurations,scores,output=map(lambda p:Path(p).resolve(),
        (root,acquisitions,configurations,scores,output))
    if output.exists():raise FileExistsError(output)
    split=json.loads((root/'configs/esdiva_phase6_split.json').read_text())
    assert split['design_development_seeds']==[2309,2311]
    rules=json.loads((root/'configs/esdiva_phase6_trial_rules.json').read_text())
    expected={f"{trial['gate']}_{trial['motion']}_f{round(trial['fraction']*100):03d}_s{round(trial['sigma'])}.json"
              for trial in rules['trials']}
    actual={path.name for path in configurations.glob('*.json')}
    if len(expected)!=len(rules['trials']) or actual!=expected:
        raise ValueError('Incomplete or changed phase-6 trial set')
    freeze=json.loads((root/'configs/esdiva_phase5_freeze.json').read_text())
    records=[]
    for path in sorted(configurations.glob('*.json')):
        config=json.loads(path.read_text())
        if (config['admitted_development_seeds']!=[2309,2311] or
                config['reserved_heldout_seeds']!=[2401,2411] or
                config['phase6_split_sha256']!=sha(root/'configs/esdiva_phase6_split.json') or
                config['phase6_trial_rules_sha256']!=sha(root/'configs/esdiva_phase6_trial_rules.json')):
            raise ValueError('Trial configuration changes the declared split')
        failures=[];rows=[]
        for seed in split['design_development_seeds']:
            for scene in split['scenes']:
                folder=acquisitions/str(seed)/scene
                score_path=scores/str(seed)/scene/path.name
                if not score_path.is_file():
                    failures.append(f'missing {seed}/{scene}')
                    continue
                report_path=folder/'report.json'
                audit_path=folder/'broad_iq.json'
                item=json.loads(score_path.read_text())
                report=json.loads(report_path.read_text())
                audit=json.loads(audit_path.read_text())
                if (report['seed']!=seed or report['script_sha256']!=
                        freeze['files_sha256']['examples/esdiva_phase2_acquire.py'] or
                        report['benchmark_source_sha256']!=freeze['numerical_source_sha256'] or
                        audit['source_report_sha256']!=sha(report_path) or
                        not audit['rf_sources_and_noise_verified'] or
                        item['seed']!=seed or item['scene']!=scene or
                        item['acquisition_sha256']!=sha(report_path) or
                        item['broad_audit_sha256']!=sha(audit_path) or
                        item['config_sha256']!=sha(path) or
                        item['script_sha256']!=sha(root/'examples/esdiva_phase6_score.py') or
                        item['local_phase_context_sha256']!=sha(root/'examples/esdiva_phase5_local_phase.py') or
                        item['unit_shift_context_sha256']!=sha(root/'examples/esdiva_phase6_unit_shift.py') or
                        item['rebeam_sha256']!=sha(root/'examples/esdiva_phase5_rebeam.py') or
                        item['bounded_context_sha256']!=sha(root/'examples/esdiva_phase3_bounded_context.py') or
                        item['input_gate_sha256']!=sha(root/'examples/esdiva_phase3_input_gate.py') or
                        item['numerical_source_sha256']!=freeze['numerical_source_sha256']):
                    raise ValueError(f'Changed source or score: {path.name} {seed}/{scene}')
                if (item['status']!='phase-6 paired source-resolved development comparison' or
                        not item['input_gate_core']['input_gate_pass'] or
                        not item['input_gate_wide']['pass'] or
                        item['blood_to_tissue_wide_roi_power_evaluator_only']>
                           config['evaluator_only_blood_to_tissue_roi_power_max']):
                    failures.append(f'input/source {seed}/{scene}')
                    continue
                none=item['variants']['none']['methods']
                phase=item['variants'][config['motion_variant']]['methods']
                for mode,methods in (('none',none),(config['motion_variant'],phase)):
                    if (methods['roundtrip']['paired_added_blood_iq_nrmse']>
                            config['comparison']['roundtrip_nrmse_max']):
                        failures.append(f'roundtrip {mode} {seed}/{scene}')
                for method in ('svd_cut1','axis1_cut1'):
                    blood=(phase[method]['paired_added_blood_iq_nrmse']-
                           none[method]['paired_added_blood_iq_nrmse'])
                    leakage=(phase[method]['refitted_background_output_to_input_tissue_power']/
                             none[method]['refitted_background_output_to_input_tissue_power'])
                    rows.append({'seed':seed,'scene':scene,'method':method,
                                 'blood_nrmse_increase':blood,'tissue_leakage_ratio':leakage})
                    if blood>config['comparison']['blood_nrmse_increase_max_motion']:
                        failures.append(f'blood {method} {seed}/{scene}')
                    if (scene in ('rigid_fast','smooth_slow') and leakage>
                            config['comparison']['moving_tissue_leakage_ratio_max']):
                        failures.append(f'leakage {method} {seed}/{scene}')
                if scene=='smooth_slow':
                    slow=none['axis2_slow_r4_cap040']
                    if (slow['paired_added_blood_iq_nrmse']>
                            config['comparison']['slow_axis2_blood_nrmse_max'] or
                            slow['refitted_background_output_to_input_tissue_power']>
                            config['comparison']['slow_axis2_tissue_leakage_max']):
                        failures.append(f'Axis 2 slow {seed}/{scene}')
        moving=[row['tissue_leakage_ratio'] for row in rows
                if row['scene'] in ('rigid_fast','smooth_slow')]
        records.append({'trial':path.name,'config_sha256':sha(path),
                        'motion_variant':config['motion_variant'],
                        'status':'pass' if not failures else 'fail',
                        'failures':failures,'rows':rows,
                        'worst_moving_leakage_ratio':max(moving) if moving else None})
    passing=[r for r in records if r['status']=='pass']
    passing.sort(key=lambda r:(r['worst_moving_leakage_ratio'],
        json.loads((configurations/r['trial']).read_text())['local_phase']['phase_fraction'],
        -json.loads((configurations/r['trial']).read_text())['local_phase']['spatial_sigma_pixels']))
    result={'status':'design candidate found' if passing else 'no design candidate passes',
            'phase6_split_sha256':sha(root/'configs/esdiva_phase6_split.json'),
            'source_freeze_sha256':sha(root/'configs/esdiva_phase5_freeze.json'),
            'selected_trial':passing[0]['trial'] if passing else None,
            'prospective_seed_unopened':2333,'reserved_seeds_unopened':[2401,2411],
            'trials':records}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return {'status':result['status'],'selected_trial':result['selected_trial'],
            'trial_count':len(records),'failures_per_trial':{r['trial']:len(r['failures']) for r in records}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--acquisitions',required=True,type=Path)
    p.add_argument('--configurations',required=True,type=Path)
    p.add_argument('--scores',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    print(json.dumps(review(a.root,a.acquisitions,a.configurations,a.scores,a.output)))
