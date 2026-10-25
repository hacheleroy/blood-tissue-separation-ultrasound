"""Write an intentionally incomplete phase-10 physical acquisition inventory.

This is a lab planning aid, not a physical scorer or evidence of acquisition.
The phase-10 geometry-aware filter port must be frozen separately before capture.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MOTIONS = ('static', 'rigid', 'nonrigid')
SPEEDS = ('slow', 'nominal')
CAPTURES = ('tissue_only', 'flow_off', 'flow_on')


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def scaffold(output):
    decision_path = ROOT / 'runs/esdiva_phase10_reserve/decision.json'
    manifest_path = ROOT / 'runs/esdiva_phase10_reserve/manifest.json'
    if not decision_path.is_file() or not manifest_path.is_file():
        raise ValueError('Reserve decision and archive are required before physical preparation')
    decision = json.loads(decision_path.read_text())
    if (decision['status'] !=
            'reserved bank passed; physical bridge eligible for acquisition preparation'
            or decision['missing'] or decision['failures'] or
            len(decision['rows']) != 16):
        raise ValueError('The complete frozen reserve bank must pass before physical preparation')
    reserve_manifest = json.loads(manifest_path.read_text())
    for relative,expected in reserve_manifest['files_sha256'].items():
        if digest(ROOT / relative) != expected:
            raise ValueError('Reserve archive source changed: '+relative)
    if decision['candidate_config_sha256'] != digest(
            ROOT / 'configs/esdiva_phase10_reserve_candidate.json'):
        raise ValueError('Reserve decision is not bound to the candidate')
    port=ROOT / 'examples/phantom_phase10_port.py'
    equivalence_script=ROOT / 'scripts/phantom_phase10_equivalence.py'
    equivalence_path=ROOT / 'runs/phantom_phase10_port/equivalence.json'
    equivalence=json.loads(equivalence_path.read_text())
    if (equivalence['port_sha256'] != digest(port) or
            equivalence['script_sha256'] != digest(equivalence_script) or
            len(equivalence['cases']) < 3 or
            {case['selected_mode'] for case in equivalence['cases']} != {'hybrid','reference'}):
        raise ValueError('The phase-10 physical port equivalence is not fixed')
    bridge_config_path=ROOT / 'configs/phantom_phase10_bridge.json'
    bridge_config=json.loads(bridge_config_path.read_text())
    if (bridge_config['reserve_decision_sha256']!=digest(decision_path) or
            bridge_config['reserve_manifest_sha256']!=digest(manifest_path) or
            bridge_config['reserve_candidate_sha256']!=digest(
                ROOT / 'configs/esdiva_phase10_reserve_candidate.json') or
            bridge_config['physical_port_source_sha256']!=digest(port) or
            bridge_config['equivalence_script_sha256']!=digest(equivalence_script) or
            bridge_config['equivalence_report_sha256']!=digest(equivalence_path)):
        raise ValueError('The physical bridge configuration is not fixed')
    provenance = {
        'reserve_decision_sha256': digest(decision_path),
        'reserve_manifest_sha256': digest(manifest_path),
        'phase10_reserve_split_sha256': digest(ROOT / 'configs/esdiva_phase10_reserve_split.json'),
        'phase10_reserve_candidate_sha256': digest(ROOT / 'configs/esdiva_phase10_reserve_candidate.json'),
        'phantom_v07_freeze_sha256': digest(ROOT / 'docs/PHANTOM_FREEZE_V07.json'),
        'phase10_bridge_sha256': digest(ROOT / 'docs/PHANTOM_PHASE10_BRIDGE.md'),
        'phase10_bridge_config_sha256': digest(bridge_config_path),
        'phase10_port_source_sha256': digest(port),
        'phase10_equivalence_script_sha256': digest(equivalence_script),
        'phase10_equivalence_report_sha256': digest(equivalence_path),
    }
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    card = {
        'status': 'INCOMPLETE: measured card and geometry-aware port required before acquisition',
        **provenance,
        'preregistered_utc': None,
        'first_capture_utc': None,
        'scanner': None,
        'probe': None,
        'transmit_waveform': None,
        'receive_aperture': None,
        'measured_geometry': {
            'tx_rad': [None]*3, 'rx_rad': [None]*3,
            'view_offsets_s': [None]*3, 'same_view_prf_hz': None,
            'f0_hz': None, 'sound_speed_mps': None, 'phase_sign': None,
            'z_mm': None, 'x_mm': None,
        },
        'beamformer_sha256': None,
        'structural_centerline_and_wall_source_sha256': None,
        'structural_masks_sha256': None,
        'wide_crop_and_core_rows': {'wide_shape_zx': [32,16],
                                    'core_rows_zero_based': list(range(12,20)),
                                    'crop_origin_mm': None},
        'measured_psf': None,
        'optical_registration': None,
        'stage_program': None,
        'phantom_medium_and_temperature': None,
        'channel_area_m2': None,
        'pump_and_flowmeter_calibration': None,
        'piv_calibration': None,
        'clock_synchronization': None,
        'capture_order_randomization': None,
        'fresh_remount_plan': None,
    }
    sessions = []
    for split in ('development', 'heldout'):
        triplets = []
        for motion in MOTIONS:
            for speed in SPEEDS:
                for repeat in range(1,6):
                    prefix = f'{split}_{motion}_{speed}_{repeat:02d}'
                    captures = {}
                    for kind in CAPTURES:
                        entry = {
                            'raw_rf_path': f'raw/{prefix}_{kind}.bin',
                            'raw_rf_sha256': None,
                            'wide_iq_path': f'iq/{prefix}_{kind}.npz',
                            'wide_iq_sha256': None,
                            'independent_reference_path': f'reference/{prefix}_{kind}.npz',
                            'independent_reference_sha256': None,
                            'reference_source': None,
                            'calibration_pump_off': None,
                            'calibration_motion_rest': None,
                            'flowmeter_mps': None if kind == 'flow_on' else 'not_applicable',
                        }
                        captures[kind] = entry
                    triplets.append({'motion': motion, 'speed': speed,
                                     'repeat': repeat, 'captures': captures})
        sessions.append({'split': split, 'session_id': None,
                         'geometry_matches_card': None, 'structural_masks_sha256': None,
                         'triplets': triplets})
    inventory = {
        'status': 'INCOMPLETE: no measured capture is represented by a filled placeholder',
        **provenance,
        'acquisition_card_path': 'acquisition_card.json',
        'acquisition_card_sha256': None,
        'required_wide_iq_shape': [256,3,32,16],
        'required_calibration_iq_shape': [64,3,32,16],
        'required_scored_windows_start': list(range(0,177,16)),
        'sessions': sessions,
    }
    (output / 'acquisition_card.json').write_text(json.dumps(card,indent=2)+'\n')
    (output / 'capture_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    return {'output': str(output), 'sessions': len(sessions),
            'triplets': sum(len(s['triplets']) for s in sessions),
            'capture_slots': sum(len(t['captures']) for s in sessions for t in s['triplets']),
            'ready_for_acquisition': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(scaffold(parser.parse_args().output)))


if __name__ == '__main__':
    main()
