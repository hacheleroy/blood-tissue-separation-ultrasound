import json
import numpy as np
import pytest
from numpy.testing import assert_allclose
from btbench.stap import (steering,null_projector,local_stap_views,estimate_tissue_velocity,
    lcmv_velocity_spectrum,lcmv_weights)
from btbench.stap_experiment import angle_scene,view_svd,default_stap_config,run_stap_study


def random(shape,seed=0):
    rng=np.random.default_rng(seed)
    return rng.normal(size=shape)+1j*rng.normal(size=shape)


@pytest.mark.parametrize("sign",[-1,1])
def test_stap_moving_tissue_absolute_snapshot_times(sign):
    tx=np.array([-.25,0,.25]);rx=np.zeros(3);time=np.arange(48)/1000
    gains=np.array([1.,.8*np.exp(.3j),1.2*np.exp(-.2j)])
    vel=np.array([.003,.012]);offsets=np.array([0.,1e-4,2e-4])
    st=steering(vel,tx,rx,time,5e6,phase_sign=sign,view_offsets_s=offsets,gains=gains).reshape(48,3)
    views=st[:,:,None,None]*random((1,1,4,5))
    r=local_stap_views(views,1000,5e6,tx,rx,[vel],phase_sign=sign,
        view_offsets_s=offsets,gains=gains,taps=8,hop=3)
    assert np.linalg.norm(r.views)/np.linalg.norm(views)<1e-14


def test_projector_single_column_and_unidentifiable_blood():
    s=random((12,))
    p=null_projector(s)
    assert p.shape==(12,12)
    assert_allclose(p@s,0,atol=1e-14)
    with pytest.raises(ValueError):lcmv_weights(s,s,np.eye(12))


def test_adaptive_stap_frozen_component_accounting():
    a=angle_scene(n_frames=32,shape=(4,13));g=a.geometry
    views=a.acquisition.iq_views
    r=local_stap_views(views,1000,5e6,g['tx_rad'],g['rx_rad'],[[.001,0]],
        taps=4,hop=2,adaptive=True,training_roi=a.training_roi,training_radius=20)
    assert_allclose(r.apply_frozen(a.blood_views)+r.apply_frozen(a.tissue_views)+r.apply_frozen(a.noise_views),r.views,atol=1e-12)
    assert r.diagnostics['training_samples_min']>0


def test_guard_prevents_cut_from_training_its_own_operator():
    x=random((16,3,5,7));changed=x.copy();changed[:,:,2,3]*=100
    params=dict(prf_hz=1000,f0_hz=5e6,tx_rad=[-.2,0,.2],rx_rad=[0,0,0],
        tissue_velocities_zx_mps=[[.001,.002]],taps=4,hop=4,adaptive=True,
        training_roi=np.ones((5,7),bool),guard=1,training_radius=3)
    a=local_stap_views(x,**params);b=local_stap_views(changed,**params)
    probe=random(x.shape,8)
    assert_allclose(a.apply_frozen(probe)[:,:,2,3],b.apply_frozen(probe)[:,:,2,3],atol=1e-13)


def test_calibration_recovers_velocity_from_independent_packet():
    vel=[.003,.008];tx=np.array([-.3,0,.3]);rx=np.zeros(3)
    s=steering(vel,tx,rx,np.arange(48)/1000,5e6).reshape(48,3)
    views=s[:,:,None,None]*random((1,1,3,4))
    fit=estimate_tissue_velocity(views,np.ones((3,4),bool),tx,rx,1000,5e6)
    assert_allclose(fit,vel,atol=1e-14)


def test_starved_training_is_explicit_and_matches_orthogonal_fallback():
    x=random((16,3,3,4))
    kw=dict(prf_hz=1000,f0_hz=5e6,tx_rad=[-.2,0,.2],rx_rad=[0,0,0],
        tissue_velocities_zx_mps=[[.001,0]],taps=4,hop=4)
    a=local_stap_views(x,**kw)
    b=local_stap_views(x,**kw,adaptive=True,training_roi=np.zeros((3,4),bool))
    assert_allclose(a.views,b.views,atol=1e-14)
    assert b.diagnostics['training_fallbacks']==4*3*4


def test_lcmv_known_response_and_correct_grid_detection():
    tx=np.array([-.3,0,.3]);rx=np.zeros(3);time=np.arange(48)/1000
    st=steering([.001,0],tx,rx,time,5e6).reshape(48,3)
    sb=steering([.025,0],tx,rx,time,5e6).reshape(48,3)
    texture=random((1,1,3,5));flow=np.zeros((3,5));flow[1,2]=1
    data=50*st[:,:,None,None]*texture+sb[:,:,None,None]*flow
    result=lcmv_velocity_spectrum(data,[[.001,0],[.015,0],[.025,0],[.04,0]],[[.001,0]],
        tx,rx,1000,5e6,~flow.astype(bool),(1,2),taps=12,hop=12,guard=0)
    assert result['selected_velocity_zx_mps']==[.025,0.]
    assert result['rejected_snapshot_constraints']>=4
    assert result['constraint_response_error_max']<1e-8


def test_joint_view_svd_identity_and_frozen_shape():
    x=random((32,3,3,4))
    r=view_svd(x,window=16,hop=8,discard=0)
    assert_allclose(r.views,x,atol=1e-13)
    assert_allclose(r.apply_frozen(2*x),2*x,atol=2e-13)


def test_stap_study_writes_frozen_config_and_complete_records(tmp_path):
    cfg=default_stap_config();cfg['seeds']=[9]
    cfg['generator'].update(n_frames=32,shape=[4,13])
    cfg['factorial']={'umax_relative_mps':[.02]};cfg['stress_cases']=[]
    cfg['filter'].update(taps=4,hop=4);cfg['oracle_svd_grid']=[0,1,2]
    cfg['lcmv_candidates_zx_mps']=[[.001,0],[.02,0],[.03,0]]
    r=run_stap_study(cfg,tmp_path/'study')
    assert len(r['records'])==1
    assert all(x['status']=='ok' for x in r['records'][0]['results'])
    assert (tmp_path/'study/frozen_config.json').is_file()
    json.loads((tmp_path/'study/results.json').read_text())
    assert r['records'][0]['lcmv_detector']['grid_selection_uses_test_truth'] is False
    with pytest.raises(FileExistsError):run_stap_study(cfg,tmp_path/'study')
