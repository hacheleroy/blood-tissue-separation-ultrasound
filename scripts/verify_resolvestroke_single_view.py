"""Verify the released Resolve Stroke aggregate and optional local checkpoints."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from resolvestroke_single_view import score_power, tube_versus_background  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'runs/resolvestroke_single_view/manifest.json'


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def close(actual, expected):
    if not np.isclose(actual, expected, rtol=1e-7, atol=1e-9):
        raise ValueError(f'Archive score mismatch: {actual} != {expected}')


def verify(checkpoints=None):
    manifest = json.loads(MANIFEST.read_text())
    for relative, expected in manifest['artifact_sha256'].items():
        if digest(ROOT / relative) != expected:
            raise ValueError(f'Changed released artifact: {relative}')
    cfg = json.loads((ROOT / 'configs/resolvestroke_single_view.json').read_text())
    report = json.loads((ROOT / 'runs/resolvestroke_single_view/report.json').read_text())
    for field, relative in [('script_sha256', 'scripts/resolvestroke_single_view.py'),
                            ('pipeline_sha256', 'configs/resolvestroke_complex.yaml')]:
        if cfg[field] != digest(ROOT / relative) or report[field] != cfg[field]:
            raise ValueError(f'Changed frozen {field}')
    if (report['protocol_sha256'] != digest(ROOT / 'configs/resolvestroke_single_view.json') or
            report['dataset_sha256'] != manifest['source_sha256'] or
            report['dataset_sha256'] != cfg['source_sha256'] or
            report['power_map_sha256'] != manifest['artifact_sha256'][
                'runs/resolvestroke_single_view/clip_power_maps.npz'] or
            digest(ROOT / 'src/btbench/filters.py') != manifest['filter_source_sha256']):
        raise ValueError('Changed source or archive lineage')
    if (report['methods'] != ['raw', 'svd_cut1', 'highpass_100hz'] or
            len(report['clip_scores']) != 5 or len(report['window_scores']) != 75 or
            report['n_windows_per_clip'] != 15):
        raise ValueError('Incomplete result')
    with np.load(ROOT / 'runs/resolvestroke_single_view/clip_power_maps.npz') as d:
        power, grid, mask = d['power'], d['coordinates'], d['tube_mask']
    if (power.shape != (5, 3, 57, 77) or grid.shape != (57, 77, 3) or
            mask.shape != (57, 77) or np.sum(mask == 1) != 35 or
            np.sum(mask == 2) != 9 or not np.isfinite(power).all() or
            np.any(power < 0) or not np.isclose(grid[0, 0, 0], -0.0048)):
        raise ValueError('Invalid released tube plane or power maps')
    rois = {str(k): tube_versus_background(mask, grid, cfg, k) for k in (1, 2)}
    for i, clip in enumerate(cfg['clips']):
        entry = report['clip_scores'][i]
        if entry['clip'] != clip:
            raise ValueError('Changed clip order')
        for j, name in enumerate(report['methods']):
            for label, roi in rois.items():
                recalculated = score_power(power[i, j], *roi)
                scored = entry['scores'][name][label]
                for key, value in recalculated.items():
                    close(scored[key], value)
                baseline = report['clip_scores'][0]['scores'][name][label]['mean_tube_power']
                close(scored['tube_power_vs_baseline_db'],
                      10 * np.log10(scored['mean_tube_power'] / baseline))
                if name != 'raw':
                    raw_outside = entry['scores']['raw'][label]['mean_outside_power']
                    close(scored['outside_power_relative_to_raw_db'],
                          10 * np.log10(scored['mean_outside_power'] / raw_outside))
        for k, offset in enumerate(cfg['window_starts_relative']):
            row = report['window_scores'][i * 15 + k]
            expected = [i * 4000 + offset, i * 4000 + offset + 256]
            if row['clip'] != clip or row['window'] != k or row['absolute_frames'] != expected:
                raise ValueError(f'Changed window schedule: {clip}, {k}')
    if checkpoints:
        folder = Path(checkpoints)
        for stem, pair in manifest['local_checkpoint_sha256'].items():
            for ext, expected in pair.items():
                if digest(folder / f'{stem}.{ext}') != expected:
                    raise ValueError(f'Changed local checkpoint {stem}.{ext}')
            row = json.loads((folder / f'{stem}.json').read_text())
            i, k = (int(v) for v in stem.split('_'))
            if (row['scores'] != report['window_scores'][i * 15 + k]['scores'] or
                    row['power_sha256'] != pair['npz']):
                raise ValueError(f'Changed checkpoint score {stem}')
        for name, expected in manifest['local_auxiliary_sha256'].items():
            if digest(folder / name) != expected:
                raise ValueError(f'Changed local auxiliary file {name}')
    print(json.dumps({'status': 'verified', 'clips': 5, 'windows': 75,
                      'local_checkpoints': bool(checkpoints)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-checkpoints', type=Path)
    args = parser.parse_args()
    verify(args.check_checkpoints)
