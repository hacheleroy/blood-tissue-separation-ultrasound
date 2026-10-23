"""Physical-control scoring checks; these arrays are fixtures, not phantom data."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('phantom_v07',ROOT/'scripts/phantom_v07.py')
phantom=importlib.util.module_from_spec(spec);spec.loader.exec_module(phantom)


def fixture(tmp_path):
    config=json.loads((ROOT/'configs/phantom_v07.json').read_text())
    config.update(n_frames=32,calibration_frames=16,score_start=16,repeats_per_cell=1,
                  methods=[dict(id='identity_none_h16',method='identity',motion='none',horizon=16)])
    config['metrics']['minimum_band_bins']=1
    reg=np.zeros((8,16),bool);reg[1:7,:6]=True
    tissue=np.zeros((8,16),bool);tissue[1:7,7:9]=True
    flow=np.zeros((8,16),bool);flow[2:6,11:14]=True
    masks=dict(registration_roi=reg,tissue_score_roi=tissue,flow_roi=flow)
    np.savez(tmp_path/'masks.npz',**masks)
    rng=np.random.default_rng(31);base=rng.normal(size=(3,8,16))+1j*rng.normal(size=(3,8,16))
    noise=.03*(rng.normal(size=(32,3,8,16))+1j*rng.normal(size=(32,3,8,16)))
    files={};reference=np.zeros((32,2,8,16))
    for kind in ('tissue_only','flow_off','flow_on'):
        views=base[None]+noise.copy()
        if kind=='flow_on':
            for t in range(32):views[t,:,flow]+=4*np.exp(2j*np.pi*t/8)
        path=tmp_path/f'{kind}.npz'
        np.savez(path,views=views,calibration_views=np.broadcast_to(base,(16,3,8,16)),
                 timestamps_s=np.arange(32)/1000.,calibration_timestamps_s=np.arange(-16,0)/1000.,
                 reference_displacement=reference,reference_valid=np.ones((32,8,16),bool),
                 reference_sigma_pixels=np.array(.01))
        files[kind]=dict(path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    session=dict(id='fixture',split='heldout',geometry=dict(prf_hz=1000.,f0_hz=5e6,sound_speed_mps=1540.,
        phase_sign=1,tx_rad=[-.15,0.,.15],rx_rad=[-.08,0.,.08],view_offsets_s=[0.,1e-6,2e-6]))
    item=dict(motion='static',speed='nominal',replicate=1,files=files)
    return config,masks,session,item


def test_pump_on_off_proxy_and_tissue_negative_control(tmp_path):
    config,masks,session,item=fixture(tmp_path)
    result=phantom.evaluate_triplet(item,session,tmp_path,config,masks)
    row=result['rows'][0]
    assert row['id']=='identity_none_h16' and row['identifiable']
    assert row['pump_response_retention']==pytest.approx(1.,abs=1e-12)
    assert row['tissue_leakage_ratio']==pytest.approx(1.,abs=1e-12)
    assert row['geometry_rmse_pixels']==0
    # Evaluation-only optical coordinates can change the score, but not the
    # data-only filtered output or its hash.
    source=item['files']['flow_on'];path=tmp_path/source['path']
    with np.load(path,allow_pickle=False) as f:values={k:f[k] for k in f.files}
    values['reference_displacement']=np.full((32,2,8,16),.04)
    np.savez(path,**values);source['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    changed=phantom.evaluate_triplet(item,session,tmp_path,config,masks)['rows'][0]
    assert changed['output_sha256']==row['output_sha256']


def test_preflight_rejects_clock_and_optical_support_corruption(tmp_path):
    config,_,_,item=fixture(tmp_path);path=tmp_path/item['files']['tissue_only']['path']
    with np.load(path,allow_pickle=False) as f:values={k:f[k] for k in f.files}
    values['timestamps_s'][18]+=.0001
    np.savez(path,**values)
    with pytest.raises(ValueError,match='Nonuniform'):
        phantom._capture(path,config,1000.,'static')
    values['timestamps_s']=np.arange(32)/1000.
    values['reference_sigma_pixels']=np.array(.2)
    np.savez(path,**values)
    with pytest.raises(ValueError,match='uncertainty'):
        phantom._capture(path,config,1000.,'static')


def test_heldout_gate_does_not_promote_unidentifiable_pump_response():
    config=json.loads((ROOT/'configs/phantom_v07.json').read_text());config['repeats_per_cell']=1
    records=[]
    for motion in config['motions']:
        for speed in config['flow_speeds_mps']:
            rows=[dict(id='axis1_smooth_h64_k1',identifiable=False,pump_response_retention=None,
                       tissue_leakage_ratio=0.,geometry_rmse_pixels=0.,accepted_fraction=0.),
                  dict(id='axis1_none_h64_k1',geometry_rmse_pixels=.5)]
            records.append(dict(split='heldout',motion=motion,speed=speed,rows=rows))
    decision=phantom.gate(dict(records=records),config)
    assert not decision['passed'] and all(not cell['passed'] for cell in decision['cells'])


def test_complete_small_manifest_requires_independent_controls_and_reference(tmp_path):
    config,masks,session,example=fixture(tmp_path)
    geometry=session['geometry'];card=dict(geometry=geometry,beamformer_sha256='b'*64,
        preregistered_utc='2026-01-01T10:00:00Z',first_capture_utc='2026-01-02T10:00:00Z',
        scanner='test scanner',probe='test probe',transmit_waveform='test wave',
        receive_aperture='test aperture',optical_registration='independent calibration',
        stage_program='preregistered waveform',pump_medium='seeded fluid',
        flowmeter_calibration='independent')
    card_path=tmp_path/'card.json';card_path.write_text(json.dumps(card))
    sessions=[];capture_number=0
    for split in config['splits']:
        triplets=[]
        for motion in config['motions']:
            for speed in config['flow_speeds_mps']:
                files={}
                for kind in config['matched_captures']:
                    with np.load(tmp_path/example['files'][kind]['path'],allow_pickle=False) as f:
                        values={k:f[k] for k in f.files}
                    field=values['reference_displacement'].copy()
                    if motion!='static':
                        field[:,0]+=.3;field[:,1]-=.2
                    if motion=='nonrigid':
                        z,x=np.meshgrid(np.linspace(-1,1,8),np.linspace(-1,1,16),indexing='ij')
                        field[:,0]+=.4*z*x
                    values['reference_displacement']=field
                    capture_number+=1
                    values['views']=values['views']+1e-5*capture_number
                    path=tmp_path/f'{split}_{motion}_{speed}_{kind}.npz';np.savez(path,**values)
                    files[kind]=dict(path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        rf_sha256=hashlib.sha256(path.name.encode()).hexdigest(),beamformer_sha256='b'*64,
                        reference_source='optical_DIC' if motion=='nonrigid' else 'optical_fiducials',
                        calibration_pump_off=True,calibration_motion_rest=True)
                    if kind=='flow_on':files[kind]['flowmeter_mps']=config['flow_speeds_mps'][speed]
                triplets.append(dict(motion=motion,speed=speed,replicate=1,files=files))
        mask_path=tmp_path/'masks.npz'
        sessions.append(dict(id=split,split=split,geometry=geometry,
            masks=dict(path='masks.npz',sha256=hashlib.sha256(mask_path.read_bytes()).hexdigest()),triplets=triplets))
    manifest=dict(protocol_sha256=hashlib.sha256((ROOT/'configs/phantom_v07.json').read_bytes()).hexdigest(),
        acquisition_card=dict(path='card.json',sha256=hashlib.sha256(card_path.read_bytes()).hexdigest()),sessions=sessions)
    path=tmp_path/'manifest.json';path.write_text(json.dumps(manifest))
    _,inventory=phantom.preflight(path,config)
    assert len(inventory)==36
    sessions[1]['triplets'][-1]['files']['flow_on']['reference_source']='ultrasound_tracker'
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='Ultrasound-derived'):
        phantom.preflight(path,config)


def test_scaffold_declares_every_capture_but_never_fabricates_data(tmp_path):
    config=json.loads((ROOT/'configs/phantom_v07.json').read_text())
    directory=tmp_path/'planned'
    assert phantom.scaffold(directory,config)==180
    manifest=json.loads((directory/'manifest.json').read_text())
    assert len(manifest['sessions'])==2
    assert all(len(session['triplets'])==30 for session in manifest['sessions'])
    assert manifest['sessions'][0]['triplets'][0]['files']['flow_on']['sha256'] is None
