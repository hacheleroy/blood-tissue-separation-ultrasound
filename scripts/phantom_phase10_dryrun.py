"""Exercise the phase-10 physical port on one explicitly simulated triplet.

The 180-slot acquisition card remains incomplete. This fixture is a software
contract check, never an eligible physical capture or blood-retention result.
"""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import gaussian_filter


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'examples'), str(ROOT / 'scripts')]
from btbench.shared import PublicInput, array_hash, fit_shared  # noqa: E402
from phantom_phase10_port import (  # noqa: E402
    _fit_check, build_context, choose_mode, structural_support,
)
from phantom_phase10_preflight import (  # noqa: E402
    card_geometry, capture_input, digest, inventory_shape, plan,
)


def expect_reject(label, action):
    try:
        action()
    except (ValueError, FileNotFoundError, KeyError, TypeError) as exc:
        return {'case': label, 'rejected': True, 'reason': str(exc)}
    raise AssertionError('Unsafe fixture was accepted: ' + label)


def make_fixture(folder, scaffold_card):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(7612)
    z = 14.4 + .1*np.arange(32)
    x = np.arange(-8., 8.)
    axis = np.tan(np.deg2rad(60.))*(z-16.)
    distance = abs(x[None, :] - axis[:, None])
    lumen = distance < 1.6
    wide_roi, core_roi = structural_support(z, x, axis)
    fit, check = _fit_check(wide_roi)
    tissue_score = np.zeros((8, 16), bool)
    eligible_tissue = np.argwhere(core_roi & check[12:20])
    tissue_score[tuple(eligible_tissue[:8].T)] = True
    flow = np.zeros((8, 16), bool)
    flow[tuple(np.argwhere(lumen[12:20])[:12].T)] = True
    np.savez_compressed(folder / 'structural_source.npz', axis_x_mm=axis,
                        lumen_mask_wide=lumen)
    np.savez_compressed(folder / 'structural_masks.npz', wide_tissue_roi=wide_roi,
                        core_tissue_roi=core_roi, flow_roi=flow, tissue_score_roi=tissue_score)
    rate = 3700.
    card = deepcopy(scaffold_card)
    card.update({'status': 'SIMULATED_FIXTURE_NOT_PHYSICAL',
                 'preregistered_utc': '2026-10-25T09:00:00Z',
                 'first_capture_utc': '2026-10-25T10:00:00Z',
                 'scanner': 'SIMULATED', 'probe': 'SIMULATED',
                 'transmit_waveform': 'SIMULATED three measured event placeholders',
                 'receive_aperture': 'SIMULATED', 'beamformer_sha256': '1'*64,
                 'structural_source_path': 'structural_source.npz',
                 'structural_centerline_and_wall_source_sha256': digest(folder / 'structural_source.npz'),
                 'structural_masks_path': 'structural_masks.npz',
                 'structural_masks_sha256': digest(folder / 'structural_masks.npz'),
                 'measured_psf': 'SIMULATED', 'optical_registration': 'SIMULATED',
                 'stage_program': 'SIMULATED static',
                 'phantom_medium_and_temperature': 'SIMULATED',
                 'channel_area_m2': float(np.pi*(.002)**2),
                 'pump_and_flowmeter_calibration': 'SIMULATED',
                 'piv_calibration': 'SIMULATED', 'clock_synchronization': 'SIMULATED',
                 'capture_order_randomization': 'SIMULATED order only',
                 'fresh_remount_plan': 'SIMULATED; no second physical session'})
    card['measured_geometry'] = {'tx_rad': np.deg2rad([0., -5., 5.]).tolist(),
                                 'rx_rad': [0., 0., 0.],
                                 'view_offsets_s': [0., 1/(3*rate), 2/(3*rate)],
                                 'same_view_prf_hz': rate, 'f0_hz': 5e6,
                                 'sound_speed_mps': 1540., 'phase_sign': 1,
                                 'z_mm': z.tolist(), 'x_mm': x.tolist()}
    card['wide_crop_and_core_rows']['crop_origin_mm'] = [float(z[0]), float(x[0])]
    (folder / 'fixture_card.json').write_text(json.dumps(card, indent=2) + '\n')
    base = gaussian_filter(rng.normal(size=(32,16)) + 1j*rng.normal(size=(32,16)), .7)
    base = base/max(np.sqrt(np.mean(abs(base[wide_roi])**2)), 1e-30)
    base[lumen] *= .1
    blood = gaussian_filter(rng.normal(size=(32,16)) + 1j*rng.normal(size=(32,16)), .5)*lumen
    blood /= max(np.sqrt(np.mean(abs(blood[lumen])**2)), 1e-30)
    gains = np.exp(1j*np.array([0., .08, -.07]))
    times = 1 + np.arange(256)/rate
    cal_times = np.arange(64)/rate
    captures = {}
    for index, kind in enumerate(('tissue_only', 'flow_off', 'flow_on')):
        static = (base + (.12*blood if kind != 'tissue_only' else 0.))
        cal = static[None, None]*gains[None, :, None, None] + .003*(
            rng.normal(size=(64,3,32,16)) + 1j*rng.normal(size=(64,3,32,16)))
        signal = static[None, None]*gains[None, :, None, None] + .005*(
            rng.normal(size=(256,3,32,16)) + 1j*rng.normal(size=(256,3,32,16)))
        if kind == 'flow_on':
            phase = np.exp(2j*np.pi*.02*np.arange(256))
            signal += .8*phase[:, None, None, None]*blood[None, None]*gains[None, :, None, None]
        iq_path = folder / f'{kind}_wide_iq.npz'
        np.savez_compressed(iq_path, views=signal.astype(np.complex64),
                            calibration_views=cal.astype(np.complex64),
                            timestamps_s=times, calibration_timestamps_s=cal_times)
        reference_path = folder / f'{kind}_independent_reference.npz'
        np.savez_compressed(reference_path, displacement=np.zeros((256,2,8,16),np.float32),
                            valid=np.ones((256,8,16),bool), sigma_pixels=.02+index*.001,
                            timestamps_s=times)
        raw = folder / f'{kind}_SIMULATED_raw_rf.bin'
        raw.write_bytes(('SIMULATED RAW RF PLACEHOLDER ' + kind).encode())
        captures[kind] = {'raw_rf_path': raw.name, 'raw_rf_sha256': digest(raw),
                          'wide_iq_path': iq_path.name, 'wide_iq_sha256': digest(iq_path),
                          'independent_reference_path': reference_path.name,
                          'independent_reference_sha256': digest(reference_path),
                          'reference_source': 'optical_fiducials',
                          'calibration_pump_off': True, 'calibration_motion_rest': True,
                          'flowmeter_mps': .003 if kind == 'flow_on' else 'not_applicable'}
    return card, captures


def run(plan_directory, fixture_directory, output):
    plan_directory = Path(plan_directory)
    card_path = plan_directory / 'acquisition_card.json'
    inventory_path = plan_directory / 'capture_inventory.json'
    original = json.loads(card_path.read_text())
    planned = plan(card_path, inventory_path)
    card, captures = make_fixture(fixture_directory, original)
    support = card_geometry(card, fixture_directory, fixture=True)
    checked = {kind: capture_input(fixture_directory, entry, card, support,
                                   'static', 'slow', kind)
               for kind, entry in captures.items()}
    negative = []
    missing = json.loads(inventory_path.read_text())
    del missing['sessions'][0]['triplets'][0]['captures']['flow_off']
    negative.append(expect_reject('missing matched control', lambda: inventory_shape(missing, original)))
    wrong_angle = deepcopy(card)
    wrong_angle['measured_geometry']['tx_rad'][1] += np.deg2rad(1.)
    negative.append(expect_reject('transmit angle > 0.5 degree',
                                  lambda: card_geometry(wrong_angle, fixture_directory, fixture=True)))
    missing_reference = deepcopy(captures['flow_on'])
    missing_reference['independent_reference_path'] = 'absent.npz'
    negative.append(expect_reject('missing independent optical reference',
                                  lambda: capture_input(fixture_directory, missing_reference, card,
                                                        support, 'static', 'slow', 'flow_on')))
    missing_calibration = deepcopy(captures['flow_on'])
    missing_calibration['calibration_motion_rest'] = None
    negative.append(expect_reject('unverified prior calibration state',
                                  lambda: capture_input(fixture_directory, missing_calibration, card,
                                                        support, 'static', 'slow', 'flow_on')))
    bad_meter = deepcopy(captures['flow_on'])
    bad_meter['flowmeter_mps'] = None
    negative.append(expect_reject('missing independent flowmeter',
                                  lambda: capture_input(fixture_directory, bad_meter, card,
                                                        support, 'static', 'slow', 'flow_on')))
    bad_clock = deepcopy(captures['flow_on'])
    with np.load(Path(fixture_directory) / bad_clock['wide_iq_path']) as d:
        shifted_iq = {k:d[k] for k in d.files}
    shifted_iq['calibration_timestamps_s'] = shifted_iq['calibration_timestamps_s'] + 2
    clock_file = Path(fixture_directory) / 'SIMULATED_bad_clock.npz'
    np.savez_compressed(clock_file, **shifted_iq)
    bad_clock['wide_iq_path'] = clock_file.name
    bad_clock['wide_iq_sha256'] = digest(clock_file)
    negative.append(expect_reject('calibration occurs after evaluation',
                                  lambda: capture_input(fixture_directory, bad_clock, card,
                                                        support, 'static', 'slow', 'flow_on')))
    bad_masks = deepcopy(card)
    with np.load(Path(fixture_directory) / card['structural_masks_path']) as d:
        altered_masks = {k:d[k].copy() for k in d.files}
    altered_masks['tissue_score_roi'] = altered_masks['flow_roi'].copy()
    mask_file = Path(fixture_directory) / 'SIMULATED_bad_masks.npz'
    np.savez_compressed(mask_file, **altered_masks)
    bad_masks['structural_masks_path'] = mask_file.name
    bad_masks['structural_masks_sha256'] = digest(mask_file)
    negative.append(expect_reject('tissue score overlaps the flow lumen',
                                  lambda: card_geometry(bad_masks, fixture_directory, fixture=True)))
    negative.append(expect_reject('simulated card cannot pass physical preflight',
                                  lambda: card_geometry(card, fixture_directory, fixture=False)))
    from phantom_phase10_port import build_context, choose_mode
    rows, outputs, raw_core = {}, {}, {}
    starts = list(range(0, 177, 16))
    spec = {'method':'axis1', 'horizon':64, 'discard':1, 'rho':.75}
    for kind in ('tissue_only', 'flow_off', 'flow_on'):
        iq = checked[kind]['iq']
        selected, none, choices, modes, errors = [], [], [], [], []
        for start in starts:
            core = PublicInput(iq['views'][start:start+80,:,12:20],
                               iq['calibration_views'][:16,:,12:20],
                               support['core_roi'], support['geometry'])
            wide = PublicInput(iq['views'][start:start+80], iq['calibration_views'][:16],
                               support['wide_roi'], support['geometry'])
            mode, choice = choose_mode(core, wide, support['wide_roi'], support['core_roi'])
            ctx = build_context(mode, core, wide, support['wide_roi'], support['core_roi'])
            roundtrip = fit_shared(ctx, {'method':'roundtrip','horizon':64}).views
            axis1 = fit_shared(ctx, spec).views
            none_axis1 = fit_shared(build_context('none', core, wide, support['wide_roi'],
                                                 support['core_roi']), spec).views
            true_core = core.views[64:80]
            error = float(np.linalg.norm(roundtrip-true_core)/max(np.linalg.norm(true_core),1e-30))
            if (roundtrip.shape != (16,3,8,16) or axis1.shape != roundtrip.shape or
                    none_axis1.shape != roundtrip.shape or not np.isfinite(axis1).all() or
                    not np.isfinite(none_axis1).all() or error > 1e-10):
                raise ValueError('Physical port failed simulated selected-mode roundtrip')
            selected.append(axis1)
            none.append(none_axis1)
            choices.append(choice)
            modes.append(mode)
            errors.append(error)
        outputs[kind] = {'selected': np.concatenate(selected), 'none': np.concatenate(none)}
        raw_core[kind] = iq['views'][64:256,:,12:20]
        if outputs[kind]['selected'].shape != raw_core[kind].shape:
            raise ValueError('Packet assembly omitted or repeated scored frames')
        rows[kind] = {'selected_modes': modes, 'mode_choice_first_packet': choices[0],
                      'maximum_roundtrip_relative_error': max(errors),
                      'axis1_output_sha256': array_hash(outputs[kind]['selected']),
                      'axis1_none_output_sha256': array_hash(outputs[kind]['none']),
                      'calibration_min_wide_validation':checked[kind]['calibration_min_wide_validation'],
                      'calibration_min_core_validation':checked[kind]['calibration_min_core_validation']}
    from phantom_v07 import band_power
    band = (.005, .05)
    bins = int(np.count_nonzero((abs(np.fft.fftfreq(192)) >= band[0]) &
                                (abs(np.fft.fftfreq(192)) <= band[1])))
    if bins < 6:
        raise ValueError('Insufficient slow-band bins after packet assembly')
    flow = support['flow_roi']
    input_off = band_power(raw_core['flow_off'], flow, band)
    input_excess = band_power(raw_core['flow_on'], flow, band)-input_off
    identifiable = input_excess > .5*input_off
    response = {}
    for variant in ('selected', 'none'):
        off = band_power(outputs['flow_off'][variant], flow, band)
        on = band_power(outputs['flow_on'][variant], flow, band)
        response[variant] = {'pump_response_retention':(on-off)/input_excess if identifiable else None,
                             'tissue_only_leakage':float(np.mean(abs(outputs['tissue_only'][variant][:,:,support['tissue_score_roi']])**2)/
                                                         np.mean(abs(raw_core['tissue_only'][:,:,support['tissue_score_roi']])**2))}
    report = {'status':'SIMULATED SOFTWARE FIXTURE PASSED; NO PHYSICAL VALIDATION',
              'physical_eligible':False, 'plan':planned,
              'fixture_label':card['status'], 'triplet':{'motion':'static','speed':'slow',
                                                       'captures':list(rows)},
              'source_sha256':{name:digest(ROOT/name) for name in
                               ('examples/phantom_phase10_port.py','scripts/phantom_phase10_preflight.py',
                                'scripts/phantom_phase10_dryrun.py','configs/phantom_phase10_bridge.json')},
              'support_counts':support['support_counts'],
              'window_starts':starts, 'assembled_score_frames':[64,256],
              'window_results':rows, 'slow_band_cycles_per_same_view_frame':list(band),
              'slow_band_fft_bins':bins, 'input_pump_excess_over_off':input_excess/input_off,
              'pump_response_identifiable_on_fixture':bool(identifiable),
              'paired_metric_plumbing':response, 'negative_checks':negative,
              'limitations':'One synthetic static/slow triplet over twelve windows. The physical 180-capture card remains incomplete; fixture metrics cannot validate flow preservation or moving-tissue rejection.'}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-directory', type=Path, required=True)
    parser.add_argument('--fixture-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.plan_directory, args.fixture_directory, args.output)
    print(json.dumps({'status':result['status'], 'mode_counts':{
        k:{m:v['selected_modes'].count(m) for m in set(v['selected_modes'])}
        for k,v in result['window_results'].items()},
        'windows_per_capture':len(result['window_starts']),
        'negative_checks':len(result['negative_checks'])}))
