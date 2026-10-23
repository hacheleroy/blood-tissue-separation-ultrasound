"""Check support, calibration, split and selection invariants in the shared run."""
import hashlib
import json
import lzma
from pathlib import Path
from btbench.shared_experiment import validate_shared_config, select_global_operating_points


def verify(root):
    root=Path(root);manifest=json.loads((root/'runs/shared_v05_manifest.json').read_text());item=manifest['archives'][0]
    raw=lzma.decompress(b''.join((root/p['path']).read_bytes() for p in item['parts']));report=json.loads(raw)
    config=validate_shared_config(json.loads((root/item['config_path']).read_text()))
    expected={(c['id'],s,split) for c in config['conditions'] for split,key in
              [('training','training_seeds'),('evaluation','evaluation_seeds')] for s in config[key]}
    records=report['records']
    if len(records)!=len(expected) or {(c['condition'],c['seed'],c['split']) for c in records}!=expected:
        raise ValueError('Case bank or seed split differs from frozen configuration')
    if report['evaluation_seed_truth_used_for_selection']:raise ValueError('Held-out truth reached selection')
    specs={m['id']:m for m in config['methods']};axis2_count=0;nonconverged=0;loaded_counts=[]
    for case in records:
        if len(case['results'])!=len(specs) or {r['id'] for r in case['results']}!=set(specs):raise ValueError('Incomplete case grid')
        if case['config_sha256']!=report['config_sha256'] or case['code_sha256']!=report['code_sha256']:raise ValueError('Mixed freeze')
        if case['source_metadata'].get('generator')!='specified angular IQ with imposed image-coordinate motion; not acoustic truth':
            raise ValueError('Unexpected data provenance')
        for mode,contexts in case['contexts'].items():
            calibration=[ctx['calibration'] for ctx in contexts.values()]
            if any(c!=calibration[0] for c in calibration):raise ValueError('Counterfactual calibration differs')
            if any(c['calibration_sha256']!=case['calibration_sha256'] for c in calibration):raise ValueError('Calibration identity differs')
            if contexts['actual']['input_sha256']!=case['input_sha256']:raise ValueError('Actual observation identity differs')
            if mode=='known' and any(ctx['evaluation_only_displacement_rmse_pixels']!=0 for ctx in contexts.values()):
                raise ValueError('Known-warp control differs from prescribed coordinates')
        for row in case['results']:
            spec=specs[row['id']];d=row['diagnostics']
            if row['status']!='ok' or row['evaluation_only']!=spec['evaluation_only']:raise ValueError('Unsuccessful/mislabelled output')
            if d['input_views']!=3 or d['spatial_tile']!=config['generator']['shape']:raise ValueError('Information support changed')
            if d['horizon_frames']!=spec['horizon'] or d['latency_max_frames']!=config['packet_frames']:raise ValueError('Temporal support changed')
            if d['calibration']['calibration_sha256']!=case['calibration_sha256']:raise ValueError('Method used a different calibration')
            block_sets=[d['blocks'],*row['counterfactual_blocks'].values()]
            expected_stops=list(range(config['score_start']+config['packet_frames'],config['generator']['n_frames']+1,config['packet_frames']))
            for blocks in block_sets:
                if [b['stop'] for b in blocks]!=expected_stops:raise ValueError('Scoring interval changed')
                for b in blocks:
                    if b['start']!=b['stop']-config['packet_frames'] or b['history_start']!=max(0,b['stop']-spec['horizon']):raise ValueError('Packet/history support changed')
                    if spec['method']=='axis2':
                        axis2_count+=1;nonconverged+=not b['ajd_converged']
                        if b['operator_gain_max']>1+1e-12:raise ValueError('Axis 2 gain cap failed')
                    if spec['method']=='stap' and spec.get('adaptive'):
                        loaded_counts.append(b['training_snapshots_min']/b['snapshot_dimension'])
            if spec['method']=='axis2':
                diag=row['metrics']['axis2_basis_diagnostic']
                if not diag['evaluation_only'] or diag['selection_use']:raise ValueError('Truth diagnostic not isolated')
                if not diag['solver_success'] or diag['best_diagonal_uncapped_nrmse']>diag['heuristic_uncapped_nrmse']+1e-7:
                    raise ValueError('Diagonal-weight diagnostic failed')
    selected=select_global_operating_points(records,config['methods'],config['background_targets'])
    if selected!=report['operating_points']:raise ValueError('Training-only global selection is not reproducible')
    if hashlib.sha256((root/item['config_path']).read_bytes()).hexdigest()!=report['config_sha256']:raise ValueError('Config hash differs')
    return dict(cases=len(records),outputs=sum(len(c['results']) for c in records),axis2_packet_fits=axis2_count,
                axis2_iteration_limit_fits=nonconverged,loaded_stap_min_K_over_d=min(loaded_counts),
                loaded_stap_max_of_min_K_over_d=max(loaded_counts),selection_reproduced=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('root',nargs='?',default='.')
    print(json.dumps(verify(p.parse_args().root),indent=2))
