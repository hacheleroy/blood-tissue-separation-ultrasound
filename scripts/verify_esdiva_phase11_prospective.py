"""Verify the source-bound phase-11 prospective archive without raw RF files."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from review_esdiva_phase11_prospective import ROOT, review, sha


def verify():
    folder = ROOT / 'runs/esdiva_phase11_prospective'
    manifest = json.loads((folder / 'manifest.json').read_text())
    freeze = json.loads((folder / 'local_freeze_manifest.json').read_text())
    config = json.loads((ROOT / 'configs/esdiva_phase11_prospective.json').read_text())
    if (manifest['pre_rf_local_freeze_commit'] !=
            '5447696d12f62b7b1022997b8cf2e4b85cddbfe2'
            or freeze['upstream_esdiva_commit'] != config['upstream_commit']
            or freeze['files_sha256']['configs/esdiva_phase11_prospective.json'] !=
            sha(ROOT / 'configs/esdiva_phase11_prospective.json')
            or set(freeze['files_sha256']) != set(manifest['frozen_files_sha256'])
            or freeze['files_sha256'] != manifest['frozen_files_sha256']):
        raise ValueError('Local pre-RF freeze identity changed')
    for name, expected in freeze['files_sha256'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha(ROOT / name) != expected:
            raise ValueError(f'Frozen phase-11 source changed: {name}')
    required = {
        'configs/esdiva_phase11_prospective.json',
        'docs/ESDIVA_PHASE11_DEVELOPMENT.md',
        'docs/ESDIVA_PHASE11_PROSPECTIVE_FREEZE.md',
        'docs/ESDIVA_PHASE11_PROSPECTIVE_RESULT.md',
        'examples/esdiva_phase11_adaptive.py',
        'examples/esdiva_phase11_development_screen.py',
        'examples/esdiva_phase11_prospective_score.py',
        'scripts/review_esdiva_phase11_prospective.py',
        'scripts/verify_esdiva_phase11_development.py',
        'scripts/verify_esdiva_phase11_prospective.py',
        'runs/esdiva_phase11_development/report.json',
        'runs/esdiva_phase11_development/manifest.json',
        'runs/esdiva_phase10_stress/decision.json',
        'runs/esdiva_phase10_stress/manifest.json',
        'runs/esdiva_phase11_prospective/local_freeze_manifest.json',
    }
    required.update(
        f'runs/esdiva_phase11_prospective/fresh/{seed}/{scene}/{name}'
        for seed in config['independent_test_seeds']
        for scene in config['scenes']
        for name in ('report.json', 'wide_audit.json', 'score.json')
    )
    if set(manifest['files_sha256']) != required:
        raise ValueError('Incomplete phase-11 prospective archive manifest')
    for name, expected in manifest['files_sha256'].items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha(ROOT / name) != expected:
            raise ValueError(f'Prospective source or result changed: {name}')
    with TemporaryDirectory() as temp:
        path = Path(temp) / 'decision.json'
        reproduced = review(folder / 'fresh', path)
        archived = json.loads((folder / 'decision.json').read_text())
        if path.read_bytes() != (folder / 'decision.json').read_bytes() or reproduced != archived:
            raise ValueError('Prospective decision does not reproduce')
    if manifest['decision_sha256'] != sha(folder / 'decision.json'):
        raise ValueError('Prospective decision hash changed')
    if archived['missing'] or len(archived['mode_rows']) + sum(
            failure.startswith('input/source') for failure in archived['failures']) != 10:
        raise ValueError('Incomplete prospective scene matrix')
    return {'status': archived['status'], 'files': len(manifest['files_sha256']),
            'scenes': 10, 'paired_rows': len(archived['paired_rows']),
            'failures': len(archived['failures']),
            'reference_coverage_achieved': archived['reference_coverage_achieved']}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
