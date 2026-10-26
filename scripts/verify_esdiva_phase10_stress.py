"""Verify the independent phase-10 stress archive without large RF checkpoints."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from review_esdiva_phase10_stress import ROOT, review, sha
from review_esdiva_phase10_reference_ablation import review as review_ablation


def verify():
    folder = ROOT / 'runs/esdiva_phase10_stress'
    manifest = json.loads((folder / 'manifest.json').read_text())
    required = {
        'configs/esdiva_phase10_stress.json',
        'configs/esdiva_phase10_reference_ablation.json',
        'configs/esdiva_phase10_reserve_candidate.json',
        'docs/ESDIVA_PHASE10_STRESS_FREEZE.md',
        'docs/ESDIVA_PHASE10_STRESS_SCHEMA_REPAIR.md',
        'docs/ESDIVA_PHASE10_REFERENCE_ABLATION_FREEZE.md',
        'examples/esdiva_phase10_stress_score.py',
        'examples/esdiva_phase10_reference_stress.py',
        'examples/esdiva_phase10_acquire.py',
        'examples/esdiva_phase5_rebeam.py',
        'examples/esdiva_phase10_support.py',
        'examples/esdiva_phase8_selector.py',
        'examples/esdiva_phase7_reference_transport.py',
        'examples/esdiva_phase6_unit_shift.py',
        'scripts/review_esdiva_phase10_stress.py',
        'scripts/review_esdiva_phase10_reference_ablation.py',
        'scripts/verify_esdiva_phase10_stress.py',
        'runs/esdiva_phase10_reserve/decision.json',
        'runs/esdiva_phase10_reserve/manifest.json',
    }
    required.update(f'runs/esdiva_phase10_stress/fresh/{seed}/{scene}/{name}'
                    for seed in (3607, 3619)
                    for scene in ('rigid_fast', 'rigid_slow', 'smooth_slow')
                    for name in ('report.json', 'wide_audit.json', 'score.json'))
    for seed in (3607, 3619):
        for scene in ('rigid_fast', 'rigid_slow'):
            score_path = folder / 'fresh' / str(seed) / scene / 'score.json'
            if json.loads(score_path.read_text())['status'] != 'input gate failed; no filter scores':
                required.add(f'runs/esdiva_phase10_stress/fresh/{seed}/{scene}/reference_ablation.json')
    if (manifest['pre_rf_freeze_commit'] != 'f12dc0d4195030bb4f14eadcb9588a7a2068190c'
            or set(manifest['files_sha256']) != required):
        raise ValueError('Incomplete independent stress manifest')
    for name, expected in manifest['files_sha256'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha(ROOT / name) != expected:
            raise ValueError(f'Stress archive source or result changed: {name}')
    with TemporaryDirectory() as temp:
        path = Path(temp) / 'decision.json'
        ablation_path = Path(temp) / 'ablation_decision.json'
        reconstructed = review(folder / 'fresh', path)
        archived = json.loads((folder / 'decision.json').read_text())
        if path.read_bytes() != (folder / 'decision.json').read_bytes() or reconstructed != archived:
            raise ValueError('Independent stress decision does not reproduce')
        reproduced_ablation = review_ablation(folder / 'fresh', ablation_path)
        ablation = json.loads((folder / 'ablation_decision.json').read_text())
        if (ablation_path.read_bytes() != (folder / 'ablation_decision.json').read_bytes()
                or reproduced_ablation != ablation):
            raise ValueError('Reference branch diagnostic does not reproduce')
    if manifest['decision_sha256'] != sha(folder / 'decision.json'):
        raise ValueError('Independent stress decision hash changed')
    if manifest['ablation_decision_sha256'] != sha(folder / 'ablation_decision.json'):
        raise ValueError('Reference branch decision hash changed')
    return {'status': archived['status'], 'files': len(manifest['files_sha256']),
            'reference_coverage_achieved': archived['reference_coverage_achieved'],
            'paired_rows': len(archived['paired_rows']),
            'counterfactual_status': ablation['status']}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
