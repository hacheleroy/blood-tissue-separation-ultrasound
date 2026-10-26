"""Verify the archived synthetic phase-10 handoff, never physical data."""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / 'runs/phantom_phase10_dryrun'


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def main():
    manifest = json.loads((FOLDER / 'manifest.json').read_text())
    for name, expected in manifest['artifact_sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError(f'Changed handoff artifact {name}')
    report = json.loads((FOLDER / 'report.json').read_text())
    if (report['status'] != 'SIMULATED SOFTWARE FIXTURE PASSED; NO PHYSICAL VALIDATION' or
            report['physical_eligible'] is not False or
            report['plan']['capture_slots'] != 180 or
            report['plan']['physical_eligible'] is not False or
            report['window_starts'] != list(range(0, 177, 16)) or
            report['assembled_score_frames'] != [64, 256] or
            report['slow_band_fft_bins'] < 6 or
            report['pump_response_identifiable_on_fixture'] is not True):
        raise ValueError('Fixture scope or packet contract changed')
    if len(report['negative_checks']) != 8 or not all(x['rejected'] for x in report['negative_checks']):
        raise ValueError('Incomplete negative controls')
    for kind in ('tissue_only', 'flow_off', 'flow_on'):
        row = report['window_results'][kind]
        if (len(row['selected_modes']) != 12 or
                set(row['selected_modes']) != {'hybrid'} or
                row['maximum_roundtrip_relative_error'] > 1e-10):
            raise ValueError('Changed selected-mode IQ/roundtrip result')
    for variant in ('selected', 'none'):
        metrics = report['paired_metric_plumbing'][variant]
        if (not np.isfinite(metrics['pump_response_retention']) or
                not np.isfinite(metrics['tissue_only_leakage'])):
            raise ValueError('Invalid synthetic metric arithmetic')
    for name, expected in report['source_sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError(f'Changed fixture or port source {name}')
    print(json.dumps({'status': 'synthetic handoff archive verified',
                      'physical_eligible': False, 'windows_per_capture': 12,
                      'negative_checks': 8}))


if __name__ == '__main__':
    main()
