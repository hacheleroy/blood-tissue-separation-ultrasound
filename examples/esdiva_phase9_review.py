"""Review retrospective phase-9 screening and each fresh development scene."""

import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(root, acquisitions, scores, output):
    root,acquisitions,scores,output=map(lambda p:Path(p).resolve(),
                                       (root,acquisitions,scores,output))
    if output.exists():raise FileExistsError(output)
    config_path=root/'configs/esdiva_phase9_candidate.json'
    config=json.loads(config_path.read_text())
    split=json.loads((root/'configs/esdiva_phase9_split.json').read_text())
    freeze=json.loads((root/'configs/esdiva_phase5_freeze.json').read_text())
    screen_path=root/'runs/esdiva_phase9_design/retrospective.json'
    retrospective=json.loads(screen_path.read_text())
    old=split['previously_seen_development_seeds']
    fresh=split['new_development_confirmation_seeds']
    if (old!=[2309,2311,2503,2521,2707,2729] or fresh!=[3001,3019] or
            config['admitted_development_seeds']!=fresh or
            config['phase9_split_sha256']!=sha(root/'configs/esdiva_phase9_split.json') or
            config['retrospective_screen_sha256']!=sha(screen_path) or
            config['screen_script_sha256']!=sha(root/'examples/esdiva_phase9_screen_seen.py') or
            config['replay_source_sha256']!=sha(root/'examples/esdiva_phase9_replay_seen.py') or
            config['support_source_sha256']!=sha(root/'examples/esdiva_phase9_support.py') or
            config['acquisition_source_sha256']!=sha(root/'examples/esdiva_phase9_acquire.py') or
            config['selector_source_sha256']!=sha(root/'examples/esdiva_phase8_selector.py') or
            config['transport_source_sha256']!=sha(root/'examples/esdiva_phase7_reference_transport.py') or
            config['score_script_sha256']!=sha(root/'examples/esdiva_phase9_score.py') or
            config['review_script_sha256']!=sha(__file__) or
            config['magnitude_shift_fit_sha256']!=sha(root/'examples/esdiva_phase6_unit_shift.py') or
            config['numerical_source_sha256']!=freeze['numerical_source_sha256']):
        raise ValueError('Phase-9 candidate or split is not the fixed source')
    if (retrospective['seen_seeds']!=old or len(retrospective['rows'])!=24 or
            retrospective['source_sha256']!=config['screen_script_sha256'] or
            retrospective['support_sha256']!=config['support_source_sha256'] or
            retrospective['replay_sha256']!=config['replay_source_sha256']):
        raise ValueError('Retrospective screening identity or scope differs')
    rows=[]; failures=list(retrospective['failures']); missing=[]
    for seed in fresh:
        for scene,setting in split['scenes'].items():
            folder=acquisitions/str(seed)/scene
            score_path=scores/str(seed)/scene/'score.json'
            if not score_path.is_file():
                missing.append(f'{seed}/{scene}')
                continue
            report_path=folder/'report.json'
            audit_path=folder/'broad_iq.json'
            report=json.loads(report_path.read_text())
            audit=json.loads(audit_path.read_text())
            score=json.loads(score_path.read_text())
            if (report['seed']!=seed or report['motion']!=setting['motion'] or
                    report['flow_scale']!=setting['flow_scale'] or
                    report['script_sha256']!=config['acquisition_source_sha256'] or
                    report['benchmark_source_sha256']!=freeze['numerical_source_sha256'] or
                    audit['source_report_sha256']!=sha(report_path) or
                    audit['script_sha256']!=sha(root/'examples/esdiva_phase5_rebeam.py') or
                    not audit['rf_sources_and_noise_verified'] or
                    max(audit['core_replay_relative_errors'].values())>3e-5 or
                    score['seed']!=seed or score['scene']!=scene or
                    score['config_sha256']!=sha(config_path) or
                    score['script_sha256']!=sha(root/'examples/esdiva_phase9_score.py') or
                    score['selector_context_sha256']!=config['selector_source_sha256'] or
                    score['support_context_sha256']!=config['support_source_sha256'] or
                    score['acquisition_script_sha256']!=config['acquisition_source_sha256'] or
                    score['transport_context_sha256']!=config['transport_source_sha256'] or
                    score['magnitude_shift_fit_sha256']!=config['magnitude_shift_fit_sha256'] or
                    score['acquisition_sha256']!=sha(report_path) or
                    score['broad_audit_sha256']!=sha(audit_path) or
                    score['numerical_source_sha256']!=freeze['numerical_source_sha256']):
                raise ValueError(f'Source mismatch in {seed}/{scene}')
            if (score['status']!='phase-9 paired source-resolved development diagnostic' or
                    not score['calibration_support']['pass'] or
                    not report['source_eligibility_pass'] or
                    score['blood_to_tissue_core_roi_power_evaluator_only']>
                        config['evaluator_only_blood_to_tissue_roi_power_max'] or
                    score['blood_to_tissue_wide_roi_power_evaluator_only']>
                        config['evaluator_only_blood_to_tissue_roi_power_max']):
                failures.append(f'input/source {seed}/{scene}')
                continue
            base=score['variants']['none']['methods']
            moved=score['variants']['selected']['methods']
            bound=config['comparison']
            for name,methods in (('none',base),('selected',moved)):
                if methods['roundtrip']['paired_added_blood_iq_nrmse']>bound['roundtrip_nrmse_max']:
                    failures.append(f'roundtrip {name} {seed}/{scene}')
            for method in ('svd_cut1','axis1_cut1'):
                delta=(moved[method]['paired_added_blood_iq_nrmse']-
                       base[method]['paired_added_blood_iq_nrmse'])
                leakage=(moved[method]['refitted_background_output_to_input_tissue_power']/
                         base[method]['refitted_background_output_to_input_tissue_power'])
                rows.append({'seed':seed,'scene':scene,'method':method,
                             'blood_nrmse_increase':delta,'tissue_leakage_ratio':leakage})
                if delta>bound['blood_nrmse_increase_max_motion']:
                    failures.append(f'blood {method} {seed}/{scene}')
                if (scene in ('rigid_fast','smooth_slow') and
                        leakage>bound['moving_tissue_leakage_ratio_max']):
                    failures.append(f'leakage {method} {seed}/{scene}')
            if scene=='smooth_slow':
                axis2=base['axis2_slow_r4_cap040']
                if (axis2['paired_added_blood_iq_nrmse']>
                        bound['slow_axis2_blood_nrmse_max'] or
                        axis2['refitted_background_output_to_input_tissue_power']>
                        bound['slow_axis2_tissue_leakage_max']):
                    failures.append(f'Axis 2 slow {seed}/{scene}')
    decision={
        'status':('incomplete development comparison' if missing else
                  'qualified for prospective development 2333' if not failures else
                  'no qualified candidate; prospective and reserve unopened'),
        'candidate_config_sha256':sha(config_path),
        'phase9_split_sha256':sha(root/'configs/esdiva_phase9_split.json'),
        'retrospective_screen_sha256':sha(screen_path),
        'review_script_sha256':sha(__file__),
        'phase5_source_freeze_sha256':sha(root/'configs/esdiva_phase5_freeze.json'),
        'seen_design_seeds':old,'new_development_seeds':fresh,
        'prospective_seed_unopened':2333,'reserved_seeds_unopened':[2401,2411],
        'missing':missing,'failures':failures,'rows':rows}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(decision,indent=2,allow_nan=False)+'\n')
    return {key:decision[key] for key in ('status','missing','failures')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--acquisitions',type=Path,required=True)
    p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(review(a.root,a.acquisitions,a.scores,a.output)))
