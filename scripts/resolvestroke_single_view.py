"""Exploratory single-transmit CEUS tube-mask diagnostic on Resolve Stroke IQ.

This deliberately does not run or score the frozen three-view phase-10 method.
The known tube mask fixes the anatomical crop and is otherwise evaluator-only.
Requires the optional zea==0.1.6 CPU/JAX environment for beamforming.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.stats import rankdata

from btbench.filters import svd_filter, temporal_highpass


ROOT = Path(__file__).resolve().parents[1]
METHODS = ('raw', 'svd_cut1', 'highpass_100hz')


def digest(path):
    with Path(path).open('rb') as stream:
        h = sha256()
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def protocol(path):
    path = Path(path)
    cfg = json.loads(path.read_text())
    if (cfg['window_frames'] % cfg['beamform_block_frames'] != 0 or
            len(cfg['clips']) != 5 or len(cfg['window_starts_relative']) != 15 or
            cfg['window_starts_relative'][-1] + cfg['window_frames'] > cfg['source_frames_per_clip']):
        raise ValueError('Changed or invalid ensemble schedule')
    return cfg


def geometry(h5, cfg):
    coords = h5['custom/computed_references/coordinates']
    x = coords[0, 0, :, 0]
    y = coords[0, :, 0, 1]
    z = coords[:, 0, 0, 2]
    ix = int(np.argmin(abs(x - cfg['grid']['x_target_m'])))
    iy = np.flatnonzero(abs(y) <= cfg['grid']['y_absolute_max_m'])
    iz = np.flatnonzero((z >= cfg['grid']['z_min_m']) & (z <= cfg['grid']['z_max_m']))
    if not len(iy) or not len(iz) or np.any(np.diff(iy) != 1) or np.any(np.diff(iz) != 1):
        raise ValueError('Empty or noncontiguous anatomical grid')
    zs, ys = slice(int(iz[0]), int(iz[-1]) + 1), slice(int(iy[0]), int(iy[-1]) + 1)
    grid = np.asarray(coords[zs, ys, ix, :])
    mask = np.asarray(h5['custom/computed_references/tube_mask'][zs, ys, ix])
    if grid.shape[:2] != mask.shape or set(np.unique(mask)) != {0, 1, 2}:
        raise ValueError('Expected both known tubes in the anatomical plane')
    return grid, mask, {'x_index': ix, 'x_m': float(x[ix]), 'y_indices': [int(iy[0]), int(iy[-1])],
                        'z_indices': [int(iz[0]), int(iz[-1])],
                        'tube_pixels': {str(k): int(np.sum(mask == k)) for k in (1, 2)}}


def tube_versus_background(mask, grid, cfg, label):
    tube = mask == label
    z = grid[..., 2]
    margin = cfg['mask']['depth_match_margin_m']
    depth = (z >= np.min(z[tube]) - margin) & (z <= np.max(z[tube]) + margin)
    outside = (mask == 0) & (distance_transform_edt(mask == 0) > cfg['mask']['outside_guard_pixels']) & depth
    if tube.sum() < 5 or outside.sum() < 50:
        raise ValueError('Insufficient tube or depth-matched outside pixels')
    return tube, outside


def score_power(power, tube, outside):
    pos = np.asarray(power[tube], dtype=np.float64)
    neg = np.asarray(power[outside], dtype=np.float64)
    if not (np.isfinite(pos).all() and np.isfinite(neg).all() and np.min(pos) >= 0 and np.min(neg) >= 0):
        raise ValueError('Nonfinite or negative Doppler power')
    rank = rankdata(np.concatenate([pos, neg]), method='average')
    n_pos, n_neg = len(pos), len(neg)
    auc = (rank[:n_pos].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    threshold = np.quantile(neg, .99)
    return {'auc': float(auc), 'tpr_at_outside_99th_percentile': float(np.mean(pos > threshold)),
            'mean_tube_power': float(np.mean(pos)), 'mean_outside_power': float(np.mean(neg)),
            'mean_power_contrast_db': float(10 * np.log10(max(np.mean(pos), 1e-30) / max(np.mean(neg), 1e-30))),
            'n_tube_pixels': n_pos, 'n_outside_pixels': n_neg}


def power_maps(iq, prf):
    if iq.ndim != 3 or not np.iscomplexobj(iq) or not np.isfinite(iq).all():
        raise ValueError('Need finite complex beamformed [time,z,y] IQ')
    svd = svd_filter(iq, centered=False, selector='fixed', discard=1).iq
    hp = temporal_highpass(iq, cutoff_hz=100.0, prf_hz=prf,
                           order=3, phase_mode='zero_phase').iq
    return np.stack([np.mean(abs(a)**2, axis=0) for a in (iq, svd, hp)]).astype(np.float32)


def check_source(h5, cfg):
    raw = h5['tracks/track_0/data/raw_data']
    scan = h5['tracks/track_0/scan']
    if raw.shape != (20000, 1, 320, 256, 2):
        raise ValueError('Changed raw IQ shape')
    if float(scan['center_frequency'][()]) != 2031250.0 or float(scan['sampling_frequency'][()]) != 2031250.0:
        raise ValueError('Changed channel IQ frequencies')
    if scan['tgc_gain_curve'].shape != (320,) or np.min(scan['tgc_gain_curve'][:]) <= 0:
        raise ValueError('Missing or invalid hardware TGC')
    if abs(float(scan['sound_speed'][()]) - 1540) > .01:
        raise ValueError('Changed beamformer sound speed')
    labels = h5['metadata/annotations/label']
    dt = scan['time_to_next_transmit'][:, 0]
    for i, name in enumerate(cfg['clips']):
        start = i * cfg['source_frames_per_clip']
        if labels[start] != name.encode() or labels[start + 3999] != name.encode():
            raise ValueError('Changed clip label or ordering')
        # The near-uniform 4 kHz segments must not be treated as one contiguous ensemble.
        segment = dt[start:start + 3999]
        if not np.allclose(segment, 1 / cfg['frame_rate_hz'], rtol=.02, atol=1e-6):
            raise ValueError('Unexpected within-clip transmit timing')
    return raw


def process(source, output, cfg_path, pipeline_path):
    import h5py
    import hdf5plugin  # noqa: F401; registers published HDF5 compression filters
    import zea
    from zea import Config, File, Pipeline

    cfg = protocol(cfg_path)
    source, output = Path(source), Path(output)
    if zea.__version__ != cfg['zea_version']:
        raise ValueError('Changed Zea reconstruction version')
    if source.stat().st_size != cfg['source_size_bytes'] or digest(source) != cfg['source_sha256']:
        raise ValueError('Dataset bytes differ from the declared release')
    output.mkdir(parents=True, exist_ok=True)
    protocol_sha = digest(cfg_path)
    script_sha = digest(__file__)
    pipeline_sha = digest(pipeline_path)
    if script_sha != cfg['script_sha256'] or pipeline_sha != cfg['pipeline_sha256']:
        raise ValueError('Analysis or beamformer source changed after the protocol freeze')
    config = Config.from_path(str(pipeline_path))
    zea.init_device()
    with h5py.File(source, 'r') as h5, File(str(source)) as zf:
        raw = check_source(h5, cfg)
        grid, mask, grid_audit = geometry(h5, cfg)
        parameters = zf.load_parameters(**config.parameters)
        pipeline = Pipeline.from_config(config)
        prepared = pipeline.prepare_parameters(parameters, grid=grid, flatgrid=grid.reshape(-1, 3))
        tgc = np.asarray(parameters.tgc_gain_curve, dtype=np.float32).reshape(1, 1, -1, 1, 1)
        rois = {str(k): tube_versus_background(mask, grid, cfg, k) for k in (1, 2)}
        np.savez_compressed(output / 'grid.npz', coordinates=grid, tube_mask=mask)
        (output / 'source_audit.json').write_text(json.dumps({
            'dataset_sha256': cfg['source_sha256'], 'protocol_sha256': protocol_sha,
            'script_sha256': script_sha, 'pipeline_sha256': pipeline_sha,
            'grid': grid_audit, 'methods': list(METHODS),
            'zea_version': zea.__version__}, indent=2) + '\n')
        block_size = cfg['beamform_block_frames']
        for clip_idx, clip in enumerate(cfg['clips']):
            for window_idx, offset in enumerate(cfg['window_starts_relative']):
                stem = f'{clip_idx:02d}_{window_idx:02d}'
                target, audit = output / f'{stem}.npz', output / f'{stem}.json'
                if target.exists() and audit.exists():
                    record = json.loads(audit.read_text())
                    if (record['protocol_sha256'] != protocol_sha or
                            record['script_sha256'] != script_sha or
                            record['power_sha256'] != digest(target)):
                        raise ValueError(f'Changed checkpoint {stem}')
                    print(clip, window_idx + 1, 'checkpoint verified', flush=True)
                    continue
                start = clip_idx * cfg['source_frames_per_clip'] + offset
                end = start + cfg['window_frames']
                if end > (clip_idx + 1) * cfg['source_frames_per_clip']:
                    raise ValueError('Window crosses the clip gap')
                pieces = []
                for pos in range(start, end, block_size):
                    channel = np.asarray(raw[pos:pos + block_size], dtype=np.float32)
                    channel /= tgc
                    result = pipeline(**{pipeline.key: channel}, **prepared, return_numpy=True)
                    array = np.asarray(result[pipeline.output_key])
                    if array.shape != (block_size, *mask.shape, 2):
                        raise ValueError('Unexpected complex beamformer output shape')
                    pieces.append(np.asarray(array[..., 0] + 1j * array[..., 1], dtype=np.complex64))
                iq = np.concatenate(pieces)
                power = power_maps(iq, cfg['frame_rate_hz'])
                scores = {name: {label: score_power(power[j], *roi)
                                 for label, roi in rois.items()}
                          for j, name in enumerate(METHODS)}
                for name in ('svd_cut1', 'highpass_100hz'):
                    for label in rois:
                        scores[name][label]['outside_power_relative_to_raw_db'] = float(
                            10 * np.log10(max(scores[name][label]['mean_outside_power'], 1e-30) /
                                          max(scores['raw'][label]['mean_outside_power'], 1e-30)))
                np.savez_compressed(target, power=power)
                record = {'clip': clip, 'clip_index': clip_idx, 'window_index': window_idx,
                          'absolute_frames': [start, end], 'protocol_sha256': protocol_sha,
                          'script_sha256': script_sha, 'pipeline_sha256': pipeline_sha,
                          'power_sha256': digest(target), 'scores': scores}
                audit.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
                print(clip, window_idx + 1, '/', len(cfg['window_starts_relative']),
                      'tube AUC', round(scores['svd_cut1']['1']['auc'], 3),
                      round(scores['svd_cut1']['2']['auc'], 3), flush=True)
    return {'status': 'all windows acquired', 'windows': len(cfg['clips']) * len(cfg['window_starts_relative'])}


def finalize(output, cfg_path):
    cfg = protocol(cfg_path)
    output = Path(output)
    audit = json.loads((output / 'source_audit.json').read_text())
    if audit['protocol_sha256'] != digest(cfg_path) or audit['script_sha256'] != digest(__file__):
        raise ValueError('Source or analysis changed after acquisition')
    if audit['pipeline_sha256'] != cfg['pipeline_sha256']:
        raise ValueError('Beamformer configuration changed after acquisition')
    with np.load(output / 'grid.npz') as g:
        grid, mask = g['coordinates'], g['tube_mask']
    rois = {str(k): tube_versus_background(mask, grid, cfg, k) for k in (1, 2)}
    maps, window_rows, clip_scores = [], [], []
    for clip_idx, clip in enumerate(cfg['clips']):
        powers = []
        for window_idx, offset in enumerate(cfg['window_starts_relative']):
            stem = f'{clip_idx:02d}_{window_idx:02d}'
            target, p = output / f'{stem}.npz', output / f'{stem}.json'
            row = json.loads(p.read_text())
            if (row['clip'] != clip or row['absolute_frames'] !=
                    [clip_idx * 4000 + offset, clip_idx * 4000 + offset + 256] or
                    row['protocol_sha256'] != audit['protocol_sha256'] or
                    row['script_sha256'] != audit['script_sha256'] or
                    row['power_sha256'] != digest(target)):
                raise ValueError(f'Changed or missing checkpoint {stem}')
            with np.load(target) as d:
                powers.append(d['power'])
            window_rows.append({'clip': clip, 'window': window_idx,
                                'absolute_frames': row['absolute_frames'], 'scores': row['scores']})
        mean_power = np.mean(powers, axis=0)
        maps.append(mean_power)
        clip_scores.append({'clip': clip, 'scores': {name: {
            label: score_power(mean_power[j], *roi) for label, roi in rois.items()}
            for j, name in enumerate(METHODS)}})
    for entry in clip_scores:
        for label in rois:
            for name in ('svd_cut1', 'highpass_100hz'):
                entry['scores'][name][label]['outside_power_relative_to_raw_db'] = float(
                    10 * np.log10(max(entry['scores'][name][label]['mean_outside_power'], 1e-30) /
                                  max(entry['scores']['raw'][label]['mean_outside_power'], 1e-30)))
    for label in rois:
        for j, name in enumerate(METHODS):
            baseline = clip_scores[0]['scores'][name][label]['mean_tube_power']
            for entry in clip_scores:
                entry['scores'][name][label]['tube_power_vs_baseline_db'] = float(
                    10 * np.log10(max(entry['scores'][name][label]['mean_tube_power'], 1e-30) /
                                  max(baseline, 1e-30)))
    np.savez_compressed(output / 'clip_power_maps.npz', power=np.asarray(maps, dtype=np.float32),
                        coordinates=grid, tube_mask=mask)
    report = {'status': 'exploratory single-view CEUS tube-mask diagnostic completed',
              'interpretation': cfg['limitations'], 'dataset_sha256': cfg['source_sha256'],
              'protocol_sha256': audit['protocol_sha256'], 'script_sha256': audit['script_sha256'],
              'pipeline_sha256': audit['pipeline_sha256'], 'power_map_sha256': digest(output / 'clip_power_maps.npz'),
              'grid': audit['grid'], 'n_windows_per_clip': len(cfg['window_starts_relative']),
              'methods': list(METHODS), 'clip_scores': clip_scores, 'window_scores': window_rows}
    (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--config', type=Path, default=ROOT / 'configs/resolvestroke_single_view.json')
    p.add_argument('--pipeline', type=Path, default=ROOT / 'configs/resolvestroke_complex.yaml')
    p.add_argument('--finalize-only', action='store_true')
    args = p.parse_args()
    if args.finalize_only:
        report = finalize(args.output, args.config)
        print(json.dumps({'status': report['status'], 'clips': len(report['clip_scores'])}))
    else:
        if args.input is None:
            p.error('--input required for beamforming')
        print(json.dumps(process(args.input, args.output, args.config, args.pipeline)))
        report = finalize(args.output, args.config)
        print(json.dumps({'status': report['status'], 'clips': len(report['clip_scores'])}))


if __name__ == '__main__':
    main()
