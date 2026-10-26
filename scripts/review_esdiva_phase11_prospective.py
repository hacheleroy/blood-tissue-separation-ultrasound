"""Apply the pre-RF phase-11 gates to every declared prospective scene."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(folder, output):
    folder, output = Path(folder), Path(output)
    config_path = ROOT / 'configs/esdiva_phase11_prospective.json'
    config = json.loads(config_path.read_text())
    parent = json.loads((ROOT / 'configs/esdiva_phase10_reserve_candidate.json').read_text())
    development = json.loads((ROOT / 'runs/esdiva_phase11_development/report.json').read_text())
    if (config['independent_test_seeds'] != [3907, 3919]
            or list(config['scenes']) != ['stationary_fast','rigid_fast','rigid_slow',
                                          'smooth_slow','smooth_zero']
            or config['methods'] != parent['methods']
            or config['comparison'] != parent['comparison']
            or sha(ROOT / 'configs/esdiva_phase10_reserve_candidate.json') != config['parent_candidate_sha256']
            or sha(ROOT / 'runs/esdiva_phase10_stress/decision.json') != config['prior_stress_decision_sha256']
            or sha(ROOT / 'runs/esdiva_phase10_stress/manifest.json') != config['prior_stress_manifest_sha256']
            or sha(ROOT / 'runs/esdiva_phase11_development/report.json') != config['development_report_sha256']
            or not development['adaptive_summary']['development_gate_pass']
            or development['seen_seeds'] != [3607,3619]
            or development['source_sha256'] != config['development_screen_sha256']
            or development['adaptive_source_sha256'] != config['adaptive_source_sha256']
            or any(sha(ROOT / name) != config[key] for key, name in (
                ('acquisition_source_sha256', 'examples/esdiva_phase10_acquire.py'),
                ('pilot_adapter_sha256','examples/esdiva_example22_pilot.py'),
                ('support_source_sha256', 'examples/esdiva_phase10_support.py'),
                ('adaptive_source_sha256', 'examples/esdiva_phase11_adaptive.py'),
                ('development_screen_sha256','examples/esdiva_phase11_development_screen.py'),
                ('rebeam_source_sha256', 'examples/esdiva_phase5_rebeam.py'),
                ('selector_source_sha256', 'examples/esdiva_phase8_selector.py'),
                ('transport_source_sha256', 'examples/esdiva_phase7_reference_transport.py'),
                ('magnitude_shift_fit_sha256', 'examples/esdiva_phase6_unit_shift.py'),
                ('score_script_sha256', 'examples/esdiva_phase11_prospective_score.py'),
                ('review_script_sha256','scripts/review_esdiva_phase11_prospective.py')))):
        raise ValueError('Frozen source, parent settings, or scene matrix changed')
    missing, failures, rows, modes = [], [], [], []
    bound = config['comparison']
    for seed in config['independent_test_seeds']:
        for scene, setting in config['scenes'].items():
            base = folder / str(seed) / scene
            paths = [base / name for name in ('report.json', 'wide_audit.json', 'score.json')]
            if any(not path.is_file() for path in paths):
                missing.append(f'{seed}/{scene}')
                continue
            report, audit, score = (json.loads(path.read_text()) for path in paths)
            if (report['seed'] != seed or report['motion'] != setting['motion']
                    or report['flow_scale'] != setting['flow_scale']
                    or report['script_sha256'] != config['acquisition_source_sha256']
                    or report['pilot_adapter_sha256'] != config['pilot_adapter_sha256']
                    or report['benchmark_source_sha256'] != config['numerical_source_sha256']
                    or report['upstream_commit'] != config['upstream_commit']
                    or report['scatterer_fraction'] != config['scatterer_fraction']
                    or report['calibration_frames'] != config['calibration_frames']
                    or report['evaluation_frames'] != config['evaluation_frames']
                    or report['view_indices'] != config['view_indices']
                    or audit['source_report_sha256'] != sha(paths[0])
                    or audit['script_sha256'] != config['rebeam_source_sha256']
                    or not audit['rf_sources_and_noise_verified']
                    or max(audit['core_replay_relative_errors'].values()) > 3e-5
                    or score['acquisition_sha256'] != sha(paths[0])
                    or score['broad_audit_sha256'] != sha(paths[1])
                    or score['config_sha256'] != sha(config_path)
                    or score['script_sha256'] != config['score_script_sha256']
                    or score['support_context_sha256'] != config['support_source_sha256']
                    or score['adaptive_context_sha256'] != config['adaptive_source_sha256']
                    or score['selector_context_sha256'] != config['selector_source_sha256']
                    or score['transport_context_sha256'] != config['transport_source_sha256']
                    or score['numerical_source_sha256'] != config['numerical_source_sha256']
                    or score['score_frames'] != 16):
                raise ValueError(f'Changed input or source identity in {seed}/{scene}')
            if (score['status'] != 'phase-11 paired source-resolved prospective diagnostic'
                    or not report['source_eligibility_pass']
                    or not score['calibration_support']['pass']
                    or score['blood_to_tissue_core_roi_power_evaluator_only'] > .1
                    or score['blood_to_tissue_wide_roi_power_evaluator_only'] > .1):
                failures.append(f'input/source {seed}/{scene}')
                continue
            selection = score['selection']
            mode = score['selected_mode']
            expected_request = ('reference' if selection['local_phase_loss'] > .28
                                and selection['reference_excess_loss'] < .18 else 'hybrid')
            expected_mode = ('reference' if expected_request == 'reference'
                             and selection['minimum_core_reference_amplitude_support'] >= .9
                             else 'hybrid')
            if (selection['selected_mode'] != mode
                    or selection['requested_mode'] != expected_request
                    or mode != expected_mode
                    or selection['source_truth_used']
                    or not selection['selection_uses_prior_calibration_and_pre_score_mixture_only']
                    or selection['pre_score_frames'] != 64
                    or (mode == 'reference' and not (
                        selection['requested_mode'] == 'reference'
                        and selection['minimum_core_reference_amplitude_support'] >= .9))):
                raise ValueError(f'Invalid frozen mode selection in {seed}/{scene}')
            strength = score['strength_selection']
            if mode == 'hybrid':
                if (strength is None or strength['pre_score_frames'] != 64
                        or strength['uses_source_truth_or_scored_frames']
                        or not 0 <= strength['weight'] <= 1
                        or abs(strength['weight']-max(0.,min(1.,
                            (strength['pre_score_shift_median_validation_gain']-.25)/.20)))>1e-12
                        or abs(strength['phase_fraction']-(1.-.5*strength['weight']))>1e-12
                        or abs(strength['shift_scale']-(.25+.5*strength['weight']))>1e-12):
                    raise ValueError(f'Invalid pre-score strength in {seed}/{scene}')
            elif strength is not None:
                raise ValueError(f'Reference selected but hybrid strength fitted in {seed}/{scene}')
            modes.append({'seed': seed, 'scene': scene, 'mode': mode,
                          'requested_mode': selection['requested_mode'],
                          'local_phase_loss': selection['local_phase_loss'],
                          'reference_excess_loss': selection['reference_excess_loss'],
                          'minimum_core_reference_amplitude_support':
                              selection['minimum_core_reference_amplitude_support'],
                          'pre_score_shift_gain':
                              None if strength is None else strength['pre_score_shift_median_validation_gain'],
                          'phase_fraction':None if strength is None else strength['phase_fraction'],
                          'shift_scale':None if strength is None else strength['shift_scale']})
            normal = score['variants']['none']['methods']
            selected = score['variants']['selected']['methods']
            for variant, methods in (('none', normal), ('selected', selected)):
                if methods['roundtrip']['paired_added_blood_iq_nrmse'] > bound['roundtrip_nrmse_max']:
                    failures.append(f'roundtrip {variant} {seed}/{scene}')
            for method in ('svd_cut1', 'axis1_cut1'):
                delta = (selected[method]['paired_added_blood_iq_nrmse']
                         - normal[method]['paired_added_blood_iq_nrmse'])
                leakage = (selected[method]['refitted_background_output_to_input_tissue_power']
                           / normal[method]['refitted_background_output_to_input_tissue_power'])
                rows.append({'seed': seed, 'scene': scene, 'mode': mode, 'method': method,
                             'blood_nrmse_increase': delta, 'tissue_leakage_ratio': leakage,
                             'selected_blood_nrmse': selected[method]['paired_added_blood_iq_nrmse'],
                             'none_blood_nrmse': normal[method]['paired_added_blood_iq_nrmse']})
                if delta > bound['blood_nrmse_increase_max_motion']:
                    failures.append(f'blood {method} {seed}/{scene}')
                if scene in ('rigid_fast','rigid_slow','smooth_slow') and \
                        leakage > bound['moving_tissue_leakage_ratio_max']:
                    failures.append(f'leakage {method} {seed}/{scene}')
            if scene.endswith('_slow'):
                axis2 = normal['axis2_slow_r4_cap040']
                if (axis2['paired_added_blood_iq_nrmse'] > bound['slow_axis2_blood_nrmse_max']
                        or axis2['refitted_background_output_to_input_tissue_power']
                        > bound['slow_axis2_tissue_leakage_max']):
                    failures.append(f'Axis 2 slow {seed}/{scene}')
    coverage = any(item['mode'] == 'reference' for item in modes)
    status = ('incomplete phase-11 prospective comparison' if missing else
              'phase-11 prospective input or numerical gate failed' if failures else
              'phase-11 prospective numerical gates passed')
    result = {'status': status,
              'config_sha256': sha(config_path),
              'prior_stress_decision_sha256': config['prior_stress_decision_sha256'],
              'independent_test_seeds': config['independent_test_seeds'],
              'missing': missing, 'failures': failures,
              'reference_coverage_achieved': coverage, 'mode_rows': modes,
              'paired_rows': rows,
              'slow_axis2': [
                  {'seed': seed, 'scene': scene, 'blood_nrmse':
                   json.loads((folder / str(seed) / scene / 'score.json').read_text())
                   ['variants']['none']['methods']['axis2_slow_r4_cap040']
                   ['paired_added_blood_iq_nrmse'], 'tissue_leakage':
                   json.loads((folder / str(seed) / scene / 'score.json').read_text())
                   ['variants']['none']['methods']['axis2_slow_r4_cap040']
                   ['refitted_background_output_to_input_tissue_power']}
                  for seed in config['independent_test_seeds']
                  for scene in ('rigid_slow', 'smooth_slow')
                  if (folder / str(seed) / scene / 'score.json').is_file()
                  and 'variants' in json.loads((folder / str(seed) / scene / 'score.json').read_text())]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(review(args.folder, args.output)['status']))
