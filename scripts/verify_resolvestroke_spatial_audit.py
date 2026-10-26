"""Verify the portable five-plane Resolve Stroke archive and its source links."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from resolvestroke_single_view import tube_versus_background  # noqa: E402
from resolvestroke_spatial_audit import score_variants  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / 'runs/resolvestroke_spatial_audit'


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def same(a, b):
    if not np.isclose(a, b, rtol=1e-7, atol=1e-9):
        raise ValueError(f'Archive score differs: {a} vs {b}')


def verify(checkpoints=None):
    manifest = json.loads((ARCHIVE / 'manifest.json').read_text())
    for relative, expected in manifest['artifact_sha256'].items():
        if digest(ROOT / relative) != expected:
            raise ValueError(f'Changed artifact {relative}')
    cfg = json.loads((ROOT / 'configs/resolvestroke_spatial_audit.json').read_text())
    old = json.loads((ROOT / 'runs/resolvestroke_single_view/report.json').read_text())
    report = json.loads((ARCHIVE / 'report.json').read_text())
    ids = report['identity']
    expected_ids = {'config_sha256': digest(ROOT / 'configs/resolvestroke_spatial_audit.json'),
                    'script_sha256': digest(ROOT / 'scripts/resolvestroke_spatial_audit.py'),
                    'pipeline_sha256': digest(ROOT / 'configs/resolvestroke_complex.yaml'),
                    'prior_report_sha256': digest(ROOT / 'runs/resolvestroke_single_view/report.json'),
                    'dataset_sha256': cfg['source_sha256']}
    if ids != expected_ids or cfg['script_sha256'] != ids['script_sha256'] or \
       cfg['pipeline_sha256'] != ids['pipeline_sha256'] or \
       cfg['prior_report_sha256'] != ids['prior_report_sha256'] or \
       cfg['prior_script_sha256'] != digest(ROOT / 'scripts/resolvestroke_single_view.py') or \
       old['dataset_sha256'] != ids['dataset_sha256'] or \
       report['clip_plane_power_sha256'] != digest(ARCHIVE / 'clip_plane_power_maps.npz'):
        raise ValueError('Frozen source lineage mismatch')
    with np.load(ARCHIVE / 'clip_plane_power_maps.npz') as d:
        power, grid, masks = d['power'], d['coordinates'], d['tube_mask']
    if (power.shape != (5, 5, 3, 57, 77) or grid.shape != (5, 57, 77, 3) or
            masks.shape != (5, 57, 77) or not np.isfinite(power).all() or np.any(power < 0) or
            report['x_indices'] != [25, 52, 77, 105, 132] or
            report['window_indices'] != [0, 7, 14] or
            report['methods'] != ['raw', 'svd_cut1', 'highpass_100hz']):
        raise ValueError('Changed five-plane result shape or schedule')
    rois = [{str(k): tube_versus_background(masks[i], grid[i], cfg, k)
             for k in (1, 2)} for i in range(5)]
    for ci, clip in enumerate(cfg['clips']):
        row = report['spatial_scores'][ci]
        if row['clip'] != clip or len(row['planes']) != 5:
            raise ValueError('Missing clip or structural plane')
        for pi in range(5):
            plane = row['planes'][pi]
            if (plane['x_index'] != cfg['x_indices'][pi] or
                    not np.isclose(plane['x_m'], grid[pi, 0, 0, 0])):
                raise ValueError('Changed structural x plane')
            for k in (1, 2):
                if plane['tube_pixels'][str(k)] != int(np.sum(masks[pi] == k)):
                    raise ValueError('Changed tube mask geometry')
            for mi, method in enumerate(report['methods']):
                for label, roi in rois[pi].items():
                    saved = plane['scores'][method][label]
                    recalculated = score_variants(power[ci, pi, mi], *roi)
                    if saved.keys() != recalculated.keys():
                        raise ValueError('Missing mask perturbation')
                    for name, result in recalculated.items():
                        for key, value in result.items():
                            same(saved[name][key], value)
                    current = saved['original']
                    baseline = report['spatial_scores'][0]['planes'][pi][
                        'scores'][method][label]['original']['mean_tube_power']
                    same(current['tube_power_vs_baseline_db'],
                         10 * np.log10(current['mean_tube_power'] / baseline))
                    if method != 'raw':
                        raw = plane['scores']['raw'][label]['original']['mean_outside_power']
                        same(current['outside_power_relative_to_raw_db'],
                             10 * np.log10(current['mean_outside_power'] / raw))
    if len(report['central_plane_temporal_power']) != 75:
        raise ValueError('Missing fifteen windows in any clip')
    for new, prior in zip(report['central_plane_temporal_power'], old['window_scores']):
        if new['clip'] != prior['clip'] or new['window'] != prior['window'] or \
           new['absolute_frames'] != prior['absolute_frames']:
            raise ValueError('Changed central-plane temporal schedule')
        for method in report['methods']:
            for label in ('1', '2'):
                same(new['tube_power'][method][label],
                     prior['scores'][method][label]['mean_tube_power'])
    replay = json.loads((ARCHIVE / 'central_plane_replay.json').read_text())
    old_manifest = json.loads((ROOT / 'runs/resolvestroke_single_view/manifest.json').read_text())
    if (replay['prior_report_sha256'] != ids['prior_report_sha256'] or
            replay['new_report_sha256'] != digest(ARCHIVE / 'report.json') or
            len(replay['comparisons']) != 15):
        raise ValueError('Incomplete central-plane replay')
    for item in replay['comparisons']:
        stem = item['stem']
        if (item['previous_power_sha256'] != old_manifest['local_checkpoint_sha256'][stem]['npz'] or
                item['new_five_plane_checkpoint_sha256'] != manifest['local_checkpoint_sha256'][stem]['npz'] or
                item['relative_l2_error'] != 0 or item['maximum_absolute_power_difference'] != 0):
            raise ValueError('Changed exact central-plane replay')
    if checkpoints:
        for stem, pair in manifest['local_checkpoint_sha256'].items():
            for suffix, expected in pair.items():
                if digest(Path(checkpoints) / f'{stem}.{suffix}') != expected:
                    raise ValueError(f'Changed local checkpoint {stem}.{suffix}')
    print(json.dumps({'status': 'verified', 'clips': 5, 'planes': 5,
                      'sampled_windows': 15, 'central_replays': 15,
                      'local_checkpoints': bool(checkpoints)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-checkpoints', type=Path)
    args = parser.parse_args()
    verify(args.check_checkpoints)
