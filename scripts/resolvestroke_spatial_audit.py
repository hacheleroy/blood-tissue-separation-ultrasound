"""Frozen exploratory longitudinal/geometry audit of Resolve Stroke CEUS IQ.

Five structural x planes and three fixed windows per clip are analyzed with
the already frozen single-view reconstruction and filters. The tube labels
are used only for fixed anatomical selection and evaluation.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation

from resolvestroke_single_view import (
    check_source, digest, power_maps, score_power, tube_versus_background,
)


ROOT = Path(__file__).resolve().parents[1]
METHODS = ('raw', 'svd_cut1', 'highpass_100hz')


def read_config(path):
    cfg = json.loads(Path(path).read_text())
    if cfg['x_indices'] != [25, 52, 77, 105, 132] or cfg['window_indices'] != [0, 7, 14]:
        raise ValueError('Changed predeclared plane/window schedule')
    return cfg


def file_sha(path):
    return digest(path)


def source_grid(h5, cfg):
    coords = h5['custom/computed_references/coordinates']
    x = coords[0, 0, :, 0]
    y = coords[0, :, 0, 1]
    z = coords[:, 0, 0, 2]
    subset = cfg['grid']
    iy = np.flatnonzero(abs(y) <= subset['y_absolute_max_m'])
    iz = np.flatnonzero((z >= subset['z_min_m']) & (z <= subset['z_max_m']))
    zs, ys = slice(int(iz[0]), int(iz[-1]) + 1), slice(int(iy[0]), int(iy[-1]) + 1)
    grids = np.stack([coords[zs, ys, ix, :] for ix in cfg['x_indices']])
    tube = h5['custom/computed_references/tube_mask']
    masks = np.stack([tube[zs, ys, ix] for ix in cfg['x_indices']])
    if grids.shape != (5, 57, 77, 3) or masks.shape != (5, 57, 77):
        raise ValueError('Changed structural grid')
    counts = []
    for plane, mask in enumerate(masks):
        if set(np.unique(mask)) != {0, 1, 2}:
            raise ValueError(f'Missing tube in structural plane {plane}')
        counts.append({str(k): int(np.sum(mask == k)) for k in (1, 2)})
    return grids, masks, {'x_m': [float(v) for v in x[cfg['x_indices']]],
                          'tube_pixels': counts, 'shape': list(masks.shape)}


def shifted(mask, dz, dy):
    out = np.zeros_like(mask, dtype=bool)
    z0, z1 = max(0, dz), mask.shape[0] + min(0, dz)
    y0, y1 = max(0, dy), mask.shape[1] + min(0, dy)
    out[z0:z1, y0:y1] = mask[z0-dz:z1-dz, y0-dy:y1-dy]
    return out


def score_variants(power, tube, outside):
    struct = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    variants = {'original': tube, 'dilate_cross_1px': binary_dilation(tube, structure=struct),
                'shift_z_minus_1px': shifted(tube, -1, 0),
                'shift_z_plus_1px': shifted(tube, 1, 0),
                'shift_y_minus_1px': shifted(tube, 0, -1),
                'shift_y_plus_1px': shifted(tube, 0, 1)}
    if any(np.any(pos & outside) for pos in variants.values()):
        raise ValueError('Perturbed tube intersects the fixed negative region')
    return {name: score_power(power, pos, outside) for name, pos in variants.items()}


def identities(cfg, cfg_path, pipeline_path, prior_path):
    if file_sha(__file__) != cfg['script_sha256']:
        raise ValueError('Changed frozen analysis source')
    if file_sha(pipeline_path) != cfg['pipeline_sha256']:
        raise ValueError('Changed frozen reconstruction configuration')
    if file_sha(ROOT / 'scripts/resolvestroke_single_view.py') != cfg['prior_script_sha256']:
        raise ValueError('Changed predecessor method source')
    if file_sha(prior_path) != cfg['prior_report_sha256']:
        raise ValueError('Changed previous 75-window report')
    prior = json.loads(Path(prior_path).read_text())
    if prior['dataset_sha256'] != cfg['source_sha256'] or prior['methods'] != list(METHODS):
        raise ValueError('Prior report/source mismatch')
    return {'config_sha256': file_sha(cfg_path), 'script_sha256': file_sha(__file__),
            'pipeline_sha256': cfg['pipeline_sha256'], 'prior_report_sha256': cfg['prior_report_sha256'],
            'dataset_sha256': cfg['source_sha256']}


def process(source, output, cfg_path, pipeline_path, prior_path):
    import h5py
    import hdf5plugin  # noqa: F401
    import zea
    from zea import Config, File, Pipeline

    cfg = read_config(cfg_path)
    identity = identities(cfg, cfg_path, pipeline_path, prior_path)
    source, output = Path(source), Path(output)
    if source.stat().st_size != cfg['source_size_bytes'] or file_sha(source) != cfg['source_sha256']:
        raise ValueError('Changed 5.25 GB source HDF5')
    if zea.__version__ != cfg['zea_version']:
        raise ValueError('Changed beamformer version')
    output.mkdir(parents=True, exist_ok=True)
    config = Config.from_path(str(pipeline_path))
    zea.init_device()
    with h5py.File(source, 'r', locking=False) as h5, File(str(source)) as zf:
        raw = check_source(h5, cfg)
        grids, masks, grid_audit = source_grid(h5, cfg)
        params = zf.load_parameters(**config.parameters)
        pipeline = Pipeline.from_config(config)
        prepared = pipeline.prepare_parameters(params, grid=grids, flatgrid=grids.reshape(-1, 3))
        tgc = np.asarray(params.tgc_gain_curve, dtype=np.float32).reshape(1, 1, -1, 1, 1)
        np.savez_compressed(output / 'structural_planes.npz', coordinates=grids, tube_mask=masks)
        (output / 'source_audit.json').write_text(json.dumps({**identity, 'grid': grid_audit,
                                                             'zea_version': zea.__version__}, indent=2) + '\n')
        for clip_index, clip in enumerate(cfg['clips']):
            for window_index in cfg['window_indices']:
                stem = f'{clip_index:02d}_{window_index:02d}'
                target = output / f'{stem}.npz'
                audit = output / f'{stem}.json'
                if target.exists() and audit.exists():
                    record = json.loads(audit.read_text())
                    if record['identity'] != identity or record['power_sha256'] != file_sha(target):
                        raise ValueError(f'Changed checkpoint {stem}')
                    print(clip, window_index, 'checkpoint verified', flush=True)
                    continue
                offset = cfg['window_starts_relative'][window_index]
                start = clip_index * 4000 + offset
                pieces = []
                for pos in range(start, start + cfg['window_frames'], cfg['beamform_block_frames']):
                    channel = np.asarray(raw[pos:pos + cfg['beamform_block_frames']], dtype=np.float32)
                    channel /= tgc
                    result = pipeline(**{pipeline.key: channel}, **prepared, return_numpy=True)
                    array = np.asarray(result[pipeline.output_key])
                    if array.shape != (cfg['beamform_block_frames'], 5, 57, 77, 2):
                        raise ValueError('Changed complex-IQ output shape')
                    pieces.append(np.asarray(array[..., 0] + 1j * array[..., 1], dtype=np.complex64))
                iq = np.concatenate(pieces)
                power = np.stack([power_maps(iq[:, i], cfg['frame_rate_hz']) for i in range(5)])
                if power.shape != (5, 3, 57, 77) or not np.isfinite(power).all():
                    raise ValueError('Invalid power maps')
                np.savez_compressed(target, power=power)
                audit.write_text(json.dumps({'clip': clip, 'window_index': window_index,
                                             'absolute_frames': [start, start + cfg['window_frames']],
                                             'power_sha256': file_sha(target), 'identity': identity},
                                            indent=2) + '\n')
                print(clip, window_index, 'acquired', flush=True)
    return {'status': 'all fixed plane/windows acquired', 'windows': 15}


def finalize(output, cfg_path, pipeline_path, prior_path):
    cfg = read_config(cfg_path)
    identity = identities(cfg, cfg_path, pipeline_path, prior_path)
    output = Path(output)
    audit = json.loads((output / 'source_audit.json').read_text())
    if any(audit[k] != v for k, v in identity.items()):
        raise ValueError('Source audit differs from freeze')
    with np.load(output / 'structural_planes.npz') as d:
        grids, masks = d['coordinates'], d['tube_mask']
    if grids.shape != (5, 57, 77, 3) or masks.shape != (5, 57, 77):
        raise ValueError('Changed structural plane archive')
    rois = [{str(k): tube_versus_background(masks[i], grids[i], cfg, k)
             for k in (1, 2)} for i in range(5)]
    clip_rows = []
    for clip_index, clip in enumerate(cfg['clips']):
        windows = []
        for window_index in cfg['window_indices']:
            stem = f'{clip_index:02d}_{window_index:02d}'
            record = json.loads((output / f'{stem}.json').read_text())
            start = clip_index * 4000 + cfg['window_starts_relative'][window_index]
            if (record['identity'] != identity or record['clip'] != clip or
                    record['absolute_frames'] != [start, start + cfg['window_frames']] or
                    record['power_sha256'] != file_sha(output / f'{stem}.npz')):
                raise ValueError(f'Changed fixed window {stem}')
            with np.load(output / f'{stem}.npz') as d:
                windows.append(d['power'])
        mean_power = np.mean(windows, axis=0).astype(np.float32)
        clip_rows.append(mean_power)
    stack = np.stack(clip_rows)
    np.savez_compressed(output / 'clip_plane_power_maps.npz', power=stack,
                        coordinates=grids, tube_mask=masks)
    spatial = []
    for clip_index, clip in enumerate(cfg['clips']):
        plane_rows = []
        for i in range(5):
            scores = {}
            for j, method in enumerate(METHODS):
                scores[method] = {label: score_variants(stack[clip_index, i, j], *roi)
                                  for label, roi in rois[i].items()}
            for label in ('1', '2'):
                raw = scores['raw'][label]['original']['mean_outside_power']
                for method in METHODS[1:]:
                    score = scores[method][label]['original']
                    score['outside_power_relative_to_raw_db'] = float(
                        10 * np.log10(max(score['mean_outside_power'], 1e-30) / max(raw, 1e-30)))
            plane_rows.append({'x_index': cfg['x_indices'][i], 'x_m': float(grids[i, 0, 0, 0]),
                               'tube_pixels': {str(k): int(np.sum(masks[i] == k)) for k in (1, 2)},
                               'scores': scores})
        spatial.append({'clip': clip, 'planes': plane_rows})
    for i in range(5):
        for method in METHODS:
            for label in ('1', '2'):
                baseline = spatial[0]['planes'][i]['scores'][method][label]['original']['mean_tube_power']
                for row in spatial:
                    score = row['planes'][i]['scores'][method][label]['original']
                    score['tube_power_vs_baseline_db'] = float(
                        10 * np.log10(max(score['mean_tube_power'], 1e-30) / max(baseline, 1e-30)))
    prior = json.loads(Path(prior_path).read_text())
    temporal = []
    for row in prior['window_scores']:
        temporal.append({'clip': row['clip'], 'window': row['window'],
                         'absolute_frames': row['absolute_frames'],
                         'tube_power': {method: {label: row['scores'][method][label]['mean_tube_power']
                                                for label in ('1', '2')} for method in METHODS}})
    report = {'status': 'exploratory structural-plane robustness completed',
              'identity': identity, 'grid': audit['grid'], 'methods': list(METHODS),
              'x_indices': cfg['x_indices'], 'window_indices': cfg['window_indices'],
              'clip_plane_power_sha256': file_sha(output / 'clip_plane_power_maps.npz'),
              'spatial_scores': spatial, 'central_plane_temporal_power': temporal,
              'interpretation': 'Known tube anatomy tests localization only; one CEUS acquisition, single transmit, no independent blood or moving-tissue truth.'}
    (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/resolvestroke_spatial_audit.json')
    parser.add_argument('--pipeline', type=Path, default=ROOT / 'configs/resolvestroke_complex.yaml')
    parser.add_argument('--prior-report', type=Path, default=ROOT / 'runs/resolvestroke_single_view/report.json')
    parser.add_argument('--finalize-only', action='store_true')
    args = parser.parse_args()
    if not args.finalize_only:
        if args.input is None:
            parser.error('--input required')
        print(json.dumps(process(args.input, args.output, args.config, args.pipeline, args.prior_report)))
    result = finalize(args.output, args.config, args.pipeline, args.prior_report)
    print(json.dumps({'status': result['status'], 'clips': len(result['spatial_scores'])}))
