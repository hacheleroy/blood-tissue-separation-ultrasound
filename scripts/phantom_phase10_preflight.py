"""Input-only checks for the phase-10 physical inventory and a labeled fixture.

This never calls a filter or computes a physical result. An incomplete card is
always rejected by `preflight`; `fixture` validates one explicitly synthetic
triplet for software plumbing only.
"""

import argparse
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import binary_dilation


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'examples')]
from phantom_phase10_port import (  # noqa: E402
    _fit_check, calibration_support, structural_support,
)

MOTIONS = ('static', 'rigid', 'nonrigid')
SPEEDS = {'slow': .003, 'nominal': .020}
CAPTURES = ('tissue_only', 'flow_off', 'flow_on')
CORE_ROWS = slice(12, 20)


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def safe_file(base, relative, expected=None):
    path = (Path(base) / relative).resolve()
    if not path.is_relative_to(Path(base).resolve()) or not path.is_file():
        raise ValueError(f'Missing or escaping file: {relative}')
    if expected is not None and digest(path) != expected:
        raise ValueError(f'SHA-256 mismatch: {relative}')
    return path


def hexhash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def inventory_shape(inventory, card):
    bridge = json.loads((ROOT / 'configs/phantom_phase10_bridge.json').read_text())
    sources = {'reserve_decision_sha256': 'runs/esdiva_phase10_reserve/decision.json',
               'reserve_manifest_sha256': 'runs/esdiva_phase10_reserve/manifest.json',
               'phase10_reserve_split_sha256': 'configs/esdiva_phase10_reserve_split.json',
               'phase10_reserve_candidate_sha256': 'configs/esdiva_phase10_reserve_candidate.json',
               'phantom_v07_freeze_sha256': 'docs/PHANTOM_FREEZE_V07.json',
               'phase10_bridge_sha256': 'docs/PHANTOM_PHASE10_BRIDGE.md',
               'phase10_bridge_config_sha256': 'configs/phantom_phase10_bridge.json',
               'phase10_port_source_sha256': 'examples/phantom_phase10_port.py',
               'phase10_equivalence_script_sha256': 'scripts/phantom_phase10_equivalence.py',
               'phase10_equivalence_report_sha256': 'runs/phantom_phase10_port/equivalence.json'}
    for key, path in sources.items():
        expected = digest(ROOT / path)
        if card.get(key) != expected or inventory.get(key) != expected:
            raise ValueError(f'Bridge provenance mismatch: {key}')
    if (bridge['reserve_decision_sha256'] != card['reserve_decision_sha256'] or
            inventory.get('required_wide_iq_shape') != [256, 3, 32, 16] or
            inventory.get('required_calibration_iq_shape') != [64, 3, 32, 16] or
            inventory.get('required_scored_windows_start') != list(range(0, 177, 16))):
        raise ValueError('Changed phase-10 bridge shapes or packet schedule')
    sessions = inventory.get('sessions', [])
    if len(sessions) != 2 or {s.get('split') for s in sessions} != {'development', 'heldout'}:
        raise ValueError('Two distinct development/held-out sessions required')
    expected = {(split, motion, speed, repeat) for split in ('development', 'heldout')
                for motion in MOTIONS for speed in SPEEDS for repeat in range(1, 6)}
    keys, paths = set(), set()
    for session in sessions:
        if len(session.get('triplets', [])) != 30:
            raise ValueError('Incomplete 3×2×5 session')
        for item in session['triplets']:
            key = (session['split'], item.get('motion'), item.get('speed'), item.get('repeat'))
            if key in keys or key not in expected or set(item.get('captures', {})) != set(CAPTURES):
                raise ValueError('Duplicate cell or missing matched control')
            keys.add(key)
            for kind in CAPTURES:
                capture = item['captures'][kind]
                for name in ('raw_rf_path', 'wide_iq_path', 'independent_reference_path'):
                    path = capture.get(name)
                    if not isinstance(path, str) or not path or path in paths:
                        raise ValueError('Missing or reused capture/reference path')
                    paths.add(path)
    if keys != expected or len(paths) != 540:
        raise ValueError('Missing factorial cell or capture path')
    return {'sessions': 2, 'triplets': 60, 'capture_slots': 180}


def plan(card_path, inventory_path):
    card = json.loads(Path(card_path).read_text())
    inventory = json.loads(Path(inventory_path).read_text())
    counts = inventory_shape(inventory, card)
    if (not card.get('status', '').startswith('INCOMPLETE') or
            not inventory.get('status', '').startswith('INCOMPLETE') or
            card.get('scanner') is not None or inventory.get('acquisition_card_sha256') is not None):
        raise ValueError('Plan check expects the original incomplete scaffold')
    return {**counts, 'status': 'scaffold structure verified; physical inputs missing',
            'physical_eligible': False}


def parse_utc(value):
    try:
        stamp = datetime.fromisoformat(value.replace('Z', '+00:00'))
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError('Missing ISO acquisition timestamp') from exc
    if stamp.tzinfo is None:
        raise ValueError('Acquisition time must include timezone')
    return stamp


def card_geometry(card, base, *, fixture=False):
    if card.get('status') != ('SIMULATED_FIXTURE_NOT_PHYSICAL' if fixture else 'MEASURED_CARD'):
        raise ValueError('Incomplete or simulated card is ineligible for physical preflight')
    if not parse_utc(card.get('preregistered_utc')) < parse_utc(card.get('first_capture_utc')):
        raise ValueError('Card must predate the first capture')
    required = ('scanner', 'probe', 'transmit_waveform', 'receive_aperture',
                'measured_psf', 'optical_registration', 'stage_program',
                'phantom_medium_and_temperature', 'pump_and_flowmeter_calibration',
                'piv_calibration', 'clock_synchronization', 'capture_order_randomization',
                'fresh_remount_plan')
    if any(not card.get(k) for k in required) or not hexhash(card.get('beamformer_sha256')):
        raise ValueError('Missing hardware, independent reference or beamformer field')
    area = card.get('channel_area_m2')
    if not isinstance(area, (int, float)) or not np.isfinite(area) or area <= 0:
        raise ValueError('Missing independently measured channel area')
    g = card['measured_geometry']
    for key in ('same_view_prf_hz', 'f0_hz', 'sound_speed_mps'):
        if not isinstance(g.get(key), (int, float)) or not np.isfinite(g[key]) or g[key] <= 0:
            raise ValueError('Missing measured acquisition constant: ' + key)
    if g.get('phase_sign') not in (-1, 1):
        raise ValueError('Missing measured phase sign')
    nominal = np.deg2rad([0., -5., 5.])
    for key in ('tx_rad', 'rx_rad', 'view_offsets_s'):
        a = np.asarray(g.get(key), dtype=float)
        if a.shape != (3,) or not np.isfinite(a).all():
            raise ValueError('Exactly three true view geometries required')
    if (np.max(abs(np.asarray(g['tx_rad']) - nominal)) > np.deg2rad(.5) or
            not np.all(np.diff(g['view_offsets_s']) > 0) or
            min(g['view_offsets_s']) < 0 or max(g['view_offsets_s']) >= 1/g['same_view_prf_hz']):
        raise ValueError('Measured angles/order or transmit offsets violate the bridge')
    if (card['wide_crop_and_core_rows'].get('wide_shape_zx') != [32, 16] or
            card['wide_crop_and_core_rows'].get('core_rows_zero_based') != list(range(12, 20)) or
            not card['wide_crop_and_core_rows'].get('crop_origin_mm')):
        raise ValueError('Changed wide/core crop')
    z, x = (np.asarray(g.get(key), dtype=float) for key in ('z_mm', 'x_mm'))
    source = safe_file(base, card['structural_source_path'],
                       card.get('structural_centerline_and_wall_source_sha256'))
    masks_file = safe_file(base, card['structural_masks_path'], card.get('structural_masks_sha256'))
    with np.load(source, allow_pickle=False) as d:
        if set(d.files) != {'axis_x_mm', 'lumen_mask_wide'}:
            raise ValueError('Independent structural centerline and lumen mask required')
        axis, lumen = d['axis_x_mm'], d['lumen_mask_wide']
    with np.load(masks_file, allow_pickle=False) as d:
        required_masks = {'wide_tissue_roi', 'core_tissue_roi', 'flow_roi', 'tissue_score_roi'}
        if set(d.files) != required_masks:
            raise ValueError('Fixed structural and score masks required')
        masks = {key: d[key] for key in required_masks}
    wide, core = structural_support(z, x, axis)
    fit, check = _fit_check(wide)
    flow, tissue = masks['flow_roi'], masks['tissue_score_roi']
    if (lumen.shape != (32, 16) or lumen.dtype != bool or
            masks['wide_tissue_roi'].dtype != bool or masks['core_tissue_roi'].dtype != bool or
            not np.array_equal(wide, masks['wide_tissue_roi']) or
            not np.array_equal(core, masks['core_tissue_roi']) or
            np.any(wide & lumen) or flow.shape != (8, 16) or tissue.shape != (8, 16) or
            flow.dtype != bool or tissue.dtype != bool or flow.sum() < 8 or tissue.sum() < 8 or
            np.any(flow & ~lumen[CORE_ROWS]) or np.any(tissue & ~core) or
            np.any(tissue & fit[CORE_ROWS]) or
            np.any(binary_dilation(flow, iterations=1) & tissue)):
        raise ValueError('Lumen/tissue anatomy, independent score pixels or flow guard invalid')
    return {'geometry': {'prf_hz': g['same_view_prf_hz'], 'f0_hz': g['f0_hz'],
                         'sound_speed_mps': g['sound_speed_mps'], 'phase_sign': g['phase_sign'],
                         'tx_rad': np.asarray(g['tx_rad']), 'rx_rad': np.asarray(g['rx_rad']),
                         'view_offsets_s': np.asarray(g['view_offsets_s'])},
            'wide_roi': wide, 'core_roi': core, 'flow_roi': flow, 'tissue_score_roi': tissue,
            'support_counts': {'wide': int(wide.sum()), 'core': int(core.sum()),
                               'fit': int(fit.sum()), 'check': int(check.sum()),
                               'flow': int(flow.sum()), 'tissue_score': int(tissue.sum())}}


def capture_input(base, item, card, support, motion, speed, kind):
    for key in ('raw_rf_sha256', 'wide_iq_sha256', 'independent_reference_sha256'):
        if not hexhash(item.get(key)):
            raise ValueError('Missing physical RF/IQ/reference checksum: ' + key)
    safe_file(base, item['raw_rf_path'], item['raw_rf_sha256'])
    iq_file = safe_file(base, item['wide_iq_path'], item['wide_iq_sha256'])
    reference_file = safe_file(base, item['independent_reference_path'], item['independent_reference_sha256'])
    source = item.get('reference_source')
    if (source not in ('optical_DIC', 'optical_fiducials', 'stage_encoder') or
            (motion == 'nonrigid' and source != 'optical_DIC') or
            item.get('calibration_pump_off') is not True or
            item.get('calibration_motion_rest') is not True):
        raise ValueError('Independent motion source or prior calibration state missing')
    with np.load(iq_file, allow_pickle=False) as d:
        keys = {'views', 'calibration_views', 'timestamps_s', 'calibration_timestamps_s'}
        if set(d.files) != keys:
            raise ValueError('Wide IQ capture contract changed')
        iq = {key: d[key] for key in keys}
    views, cal = iq['views'], iq['calibration_views']
    if (views.shape != (256, 3, 32, 16) or cal.shape != (64, 3, 32, 16) or
            not np.iscomplexobj(views) or not np.iscomplexobj(cal) or
            not np.isfinite(views).all() or not np.isfinite(cal).all()):
        raise ValueError('Missing finite unfiltered complex wide IQ or calibration')
    rate = support['geometry']['prf_hz']
    t, ct = iq['timestamps_s'], iq['calibration_timestamps_s']
    if (t.shape != (256,) or ct.shape != (64,) or not np.isfinite(t).all() or
            not np.isfinite(ct).all() or not ct[-1] < t[0] or
            not np.allclose(np.diff(t), 1/rate, rtol=1e-5, atol=1e-7/rate) or
            not np.allclose(np.diff(ct), 1/rate, rtol=1e-5, atol=1e-7/rate)):
        raise ValueError('Invalid same-view clocks or prior calibration')
    cal_gate = calibration_support(cal[:16], support['wide_roi'], support['core_roi'])
    if not cal_gate['pass']:
        raise ValueError('Prior angular calibration support failed')
    with np.load(reference_file, allow_pickle=False) as d:
        if set(d.files) != {'displacement', 'valid', 'sigma_pixels', 'timestamps_s'}:
            raise ValueError('Independent motion reference contract changed')
        ref = {key: d[key] for key in d.files}
    d, valid = ref['displacement'], ref['valid']
    if (d.shape != (256, 2, 8, 16) or valid.shape != (256, 8, 16) or valid.dtype != bool or
            not np.isfinite(d).all() or not valid[:, 1:-1, 1:-1].all() or
            not valid[64:, support['flow_roi'] | support['tissue_score_roi']].all() or
            not np.array_equal(ref['timestamps_s'], t) or
            not np.isfinite(float(ref['sigma_pixels'])) or not 0 <= float(ref['sigma_pixels']) <= .1):
        raise ValueError('Independent optical/stage field, uncertainty or synchronization invalid')
    interior = d[64:, :, 1:-1, 1:-1]
    rms = float(np.sqrt(np.mean(interior**2)))
    if motion == 'static' and rms > .1:
        raise ValueError('Static reference moved')
    if motion != 'static' and (rms < .2 or np.max(abs(interior)) > 1.25):
        raise ValueError('Moving reference outside preregistered range')
    if motion == 'nonrigid' and np.sqrt(np.mean((interior-interior.mean(axis=(-2,-1), keepdims=True))**2)) < .05:
        raise ValueError('Missing dense nonrigid variation')
    if kind == 'flow_on':
        measured = item.get('flowmeter_mps')
        if not isinstance(measured, (int, float)) or not .8*SPEEDS[speed] <= measured <= 1.2*SPEEDS[speed]:
            raise ValueError('Independent flowmeter outside target')
    elif item.get('flowmeter_mps') != 'not_applicable':
        raise ValueError('Pump-off/tissue-only flowmeter state changed')
    return {'input_iq_sha256': item['wide_iq_sha256'], 'reference_sha256': item['independent_reference_sha256'],
            'calibration_min_wide_validation': cal_gate['wide']['minimum_coherence_validation'],
            'calibration_min_core_validation': cal_gate['core']['minimum_coherence_validation'],
            'reference_rms_pixels': rms, 'iq': iq}


def preflight(card_path, inventory_path):
    card_path, inventory_path = Path(card_path).resolve(), Path(inventory_path).resolve()
    base = inventory_path.parent
    inventory, card = json.loads(inventory_path.read_text()), json.loads(card_path.read_text())
    counts = inventory_shape(inventory, card)
    if (inventory.get('status') != 'MEASURED_INVENTORY' or
            inventory.get('acquisition_card_path') != str(card_path.relative_to(base)) or
            inventory.get('acquisition_card_sha256') != digest(card_path)):
        raise ValueError('Unsealed, incomplete or simulated acquisition inventory')
    support = card_geometry(card, base)
    seen_iq, seen_rf = set(), set()
    ids = set()
    for session in inventory['sessions']:
        if (not session.get('session_id') or session['session_id'] in ids or
                session.get('geometry_matches_card') is not True or
                session.get('structural_masks_sha256') != card['structural_masks_sha256']):
            raise ValueError('Independent remount/session geometry not documented')
        ids.add(session['session_id'])
        for triplet in session['triplets']:
            for kind in CAPTURES:
                item = triplet['captures'][kind]
                for key, seen in [('wide_iq_sha256', seen_iq), ('raw_rf_sha256', seen_rf)]:
                    if item.get(key) in seen:
                        raise ValueError('Repeated capture/reference across independent controls')
                    seen.add(item.get(key))
                capture_input(base, item, card, support, triplet['motion'], triplet['speed'], kind)
    return {**counts, 'status': 'input schema and independent-control preflight passed',
            'analysis_input_eligible': True, 'experimental_origin_verified': False,
            'support_counts': support['support_counts']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'preflight'))
    parser.add_argument('--card', type=Path, required=True)
    parser.add_argument('--inventory', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(plan(args.card, args.inventory) if args.command == 'plan'
                     else preflight(args.card, args.inventory), indent=2))
