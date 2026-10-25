"""Verify the conditional phase-5 physical bridge and its pre-capture card.

This preflight checks metadata and structural masks, not phantom data or a
filter result. A physical release also needs the acquired triplets and the
full analysis specified in docs/PHANTOM_PHASE5_BRIDGE.md.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_dilation


ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_contract(root=ROOT):
    root=Path(root).resolve()
    freeze=json.loads((root/'configs/phantom_phase5_bridge_freeze.json').read_text())
    config=json.loads((root/'configs/phantom_phase5_bridge.json').read_text())
    for name,expected in freeze['files_sha256'].items():
        if digest(root/name)!=expected:
            raise ValueError(f'Bridge source changed: {name}')
    if (config['phase6_split_sha256']!=digest(root/'configs/esdiva_phase6_split.json') or
            config['phantom_v07_freeze_sha256']!=digest(root/'docs/PHANTOM_FREEZE_V07.json') or
            config['core_rows_zero_based']!=list(range(12,20)) or
            config['wide_shape_zx']!=[32,16] or
            config['calibration_used_frames_zero_based']!=list(range(16)) or
            config['calibration_acquired_frames']!=64 or config['evaluation_frames']!=256 or
            config['history_frames']!=64 or config['score_start_frame']!=64 or
            config['packet_frames']!=16):
        raise ValueError('Bridge mapping or earlier frozen protocol differs')
    return config


def verified_file(base, item):
    if not isinstance(item,dict) or not item.get('path') or not item.get('sha256'):
        raise ValueError('Missing file path or SHA256 in acquisition card')
    path=(base/item['path']).resolve()
    if not path.is_relative_to(base.resolve()) or not path.is_file() or digest(path)!=item['sha256']:
        raise ValueError(f'Missing, external or changed acquisition-card file: {path}')
    return path


def preflight(manifest_path,root=ROOT):
    config=check_contract(root)
    manifest_path=Path(manifest_path).resolve()
    base=manifest_path.parent
    manifest=json.loads(manifest_path.read_text())
    card=json.loads(verified_file(base,manifest['acquisition_card']).read_text())
    if (card.get('bridge_freeze_sha256')!=digest(Path(root)/'configs/phantom_phase5_bridge_freeze.json') or
            not all(card.get(key) for key in ('preregistered_at_utc','first_capture_after_utc',
                'scanner','probe','waveform','receive_aperture','beamformer',
                'optical_registration','stage_program','pump_medium','flowmeter_calibration'))):
        raise ValueError('Acquisition card is incomplete or bound to another bridge')
    decision=json.loads(verified_file(base,card['phase6_reserved_decision']).read_text())
    if (decision.get('status')!='pass' or decision.get('reserved_heldout_seeds')!=[2401,2411]):
        raise ValueError('Fresh phase-6 simulation reserve has not passed')
    verified_file(base,card['candidate_config'])
    verified_file(base,card['physical_port_source'])
    g=card['geometry']
    tx=np.asarray(g['tx_angles_deg'],float)
    rx=np.asarray(g['rx_angles_deg'],float)
    z=np.asarray(g['z_mm'],float)
    x=np.asarray(g['x_mm'],float)
    offsets=np.asarray(g['view_offsets_s'],float)
    if (tx.shape!=(3,) or rx.shape!=(3,) or offsets.shape!=(3,) or
            z.shape!=(32,) or x.shape!=(16,) or not all(np.isfinite(a).all()
            for a in (tx,rx,offsets,z,x)) or
            np.max(abs(tx-np.asarray(config['target_tx_angles_deg'])))>
                config['maximum_tx_angle_error_deg'] or
            not (np.all(np.diff(z)>0) and np.all(np.diff(x)>0) and
                 np.all(np.diff(offsets)>0)) or
            any(float(g[k])<=0 for k in ('prf_hz','carrier_hz','sound_speed_mps')) or
            int(g['phase_sign']) not in (-1,1)):
        raise ValueError('Three-view physical geometry does not meet the bridge')
    masks=verified_file(base,card['structural_masks'])
    with np.load(masks,allow_pickle=False) as data:
        arrays={key:data[key] for key in ('registration_roi','tissue_score_roi','flow_roi')}
    if any(a.shape!=(32,16) or a.dtype!=np.bool_ for a in arrays.values()):
        raise ValueError('Structural masks must be Boolean on the measured 32x16 grid')
    reg,tissue,flow=(arrays[k] for k in ('registration_roi','tissue_score_roi','flow_roi'))
    zz,xx=np.indices(reg.shape)
    fit=reg&((zz+2*xx)%3!=0)
    check=reg&~fit
    if (reg.sum()<128 or reg[12:20].sum()<40 or fit.sum()<80 or check.sum()<40 or
            tissue[12:20].sum()<8 or flow[12:20].sum()<8 or
            np.any(reg&tissue) or np.any(reg&binary_dilation(flow,np.ones((3,3),bool))) or
            np.any(tissue&flow)):
        raise ValueError('Structural ROI support, flow guard or disjointness failed')
    return {'status':'card structurally eligible; physical captures still required',
            'bridge_sha256':digest(Path(root)/'configs/phantom_phase5_bridge_freeze.json'),
            'card_sha256':digest(verified_file(base,manifest['acquisition_card'])),
            'registration_pixels':int(reg.sum()),'fit_pixels':int(fit.sum()),
            'check_pixels':int(check.sum()),'core_tissue_pixels':int(reg[12:20].sum())}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('check-contract','preflight'))
    p.add_argument('--manifest',type=Path)
    args=p.parse_args()
    if args.command=='check-contract':result={'status':'conditional bridge frozen',
                                             'config':check_contract()['study']}
    else:
        if args.manifest is None:p.error('--manifest is required for preflight')
        result=preflight(args.manifest)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
