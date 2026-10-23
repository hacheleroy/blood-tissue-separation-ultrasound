import copy
from dataclasses import fields, replace
import numpy as np
import pytest
from numpy.testing import assert_allclose
from btbench.shared import (PublicInput, SharedContext, bounded_moments, fit_shared,
                            shared_scene, calibrate, estimate_displacement)
from btbench.shared_experiment import (default_shared_config, validate_shared_config,
    select_global_operating_points, run_shared_study)
from btbench.paper_filters import lee_aer_filter


def scene(**kwargs):
    return shared_scene(seed=991,n_frames=64,calibration_frames=32,shape=(4,8),**kwargs)


@pytest.mark.parametrize('centered',[False,True])
@pytest.mark.parametrize('horizon,rho',[(16,.75),(32,.75),(32,1.)])
def test_expiring_moments_equal_direct_weighted_svd(centered,horizon,rho):
    a=scene();c=SharedContext(a.public,packet=16,score_start=32)
    params=dict(horizon=horizon,rho=rho,centered=centered,discard=1)
    axis=fit_shared(c,dict(method='axis1',**params))
    direct=fit_shared(c,dict(method='weighted_svd',**params))
    assert_allclose(axis.views,direct.views,atol=3e-11,rtol=1e-10)
    assert_allclose(axis.apply_frozen(a.blood),direct.apply_frozen(a.blood),atol=3e-11)


def test_bounded_lag_moments_expire_and_use_matching_endpoint_means():
    rng=np.random.default_rng(2);x=rng.normal(size=(5,48))+1j*rng.normal(size=(5,48))
    packet=8;horizon=24;rho=.7;lags=(0,1,4)
    snapshots=bounded_moments(x,packet,horizon,rho,lags)
    for index, state in enumerate(snapshots):
        stop=(index+1)*packet;left=max(0,stop-horizon)
        for lag in lags:
            values=[];weights=[];plus=[];minus=[]
            for start in range(left,stop,packet):
                a=x[:,start+lag:start+packet];b=x[:,start:start+packet-lag]
                values.append(a@b.conj().T/(packet-lag));plus.append(a.mean(1,keepdims=True));minus.append(b.mean(1,keepdims=True))
                weights.append(rho**((stop-start)//packet-1))
            mass,c,p,m=state[lag]
            assert_allclose(mass,sum(weights),atol=1e-14)
            for got,terms in [(c,values),(p,plus),(m,minus)]:
                assert_allclose(got,sum(w*v for w,v in zip(weights,terms)),atol=4e-15)


@pytest.mark.parametrize('kwargs',[dict(packet=0,horizon=16),dict(packet=8,horizon=0),
    dict(packet=8,horizon=16,lags=[]),dict(packet=8,horizon=16,lags=[8]),dict(packet=8,horizon=16,rho=0)])
def test_invalid_moment_contract(kwargs):
    with pytest.raises(ValueError):bounded_moments(np.ones((2,32)),**kwargs)


@pytest.mark.parametrize('method,params',[
    ('svd',dict(discard=1,centered=True)),('axis1',dict(discard=1,centered=True)),
    ('axis2',dict(rank=4,phase_width=.1)),('stap',dict(adaptive=True)),
    ('lee_aer',dict(domain='temporal',epsilon_abs=26.))])
def test_frozen_replay_all_views_closure_and_packet_causality(method,params):
    a=scene(motion='smooth');c=SharedContext(a.public,packet=16,score_start=32,motion='smooth')
    spec=dict(method=method,horizon=32,**params);r=fit_shared(c,spec)
    assert r.views.shape==(32,3,4,8)
    assert_allclose(r.apply_frozen(a.public.views),r.views,atol=5e-11)
    assert_allclose(sum(r.apply_frozen(v) for v in (a.blood,a.tissue,a.noise)),r.views,atol=5e-11)
    changed=a.public.views.copy();changed[48:]*=3*np.exp(.3j)
    other=SharedContext(replace(a.public,views=changed),16,32,'smooth')
    assert_allclose(fit_shared(other,spec).views[:16],r.views[:16],atol=2e-11)


@pytest.mark.parametrize('domain',['temporal','spatial'])
def test_common_view_lee_matches_full_history_published_covariance_output(domain):
    a=scene();c=SharedContext(a.public,16,32)
    r=fit_shared(c,dict(method='lee_aer',horizon=32,domain=domain,epsilon_abs=26.))
    expected=[]
    for stop in c.stops:
        hist=c.aligned[stop-32:stop].reshape(32,3*4,8)
        reference=lee_aer_filter(hist,domain=domain,epsilon_abs=26.).iq[-16:].reshape(16,3,4,8)
        expected.append(c.restore_block(reference,stop-16))
    assert_allclose(r.views,np.concatenate(expected),atol=5e-12)


def test_truth_is_separate_and_calibration_is_identical_for_counterfactual_refits():
    assert {f.name for f in fields(PublicInput)}=={'views','calibration_views','training_roi','geometry'}
    a=scene(motion='smooth');background=replace(a.public,views=a.tissue+a.noise)
    assert_allclose(calibrate(a.public)[0],calibrate(background)[0])
    c=SharedContext(a.public,16,32,'smooth');other=SharedContext(background,16,32,'smooth')
    assert c.calibration==other.calibration
    # Independent fits may both correctly reject motion and return exact zero.
    assert c.motion_diagnostics['input_sha256']!=other.motion_diagnostics['input_sha256']
    assert c.displacement is not other.displacement
    with pytest.raises(ValueError):SharedContext(a.public,16,32,'smooth',a.known_displacement)
    with pytest.raises(ValueError):SharedContext(a.public,16,32,'known')
    with pytest.raises(ValueError):fit_shared(c,dict(method='identity',horizon=16))


def test_known_motion_roundtrip_keeps_boundary_and_interpolation_loss():
    a=scene(motion='translation');c=SharedContext(a.public,16,32,'known',a.known_displacement)
    r=fit_shared(c,dict(method='roundtrip',horizon=16))
    assert np.linalg.norm(r.views-a.public.views[32:])>0
    assert_allclose(r.views,r.apply_frozen(a.public.views),atol=1e-12)
    zero=SharedContext(a.public,16,32,'known',np.zeros_like(a.known_displacement))
    assert_allclose(fit_shared(zero,dict(method='roundtrip',horizon=16)).views,a.public.views[32:],atol=2e-13)


def test_registration_recovers_interior_subpixel_translation_without_truth():
    from scipy.ndimage import shift
    rng=np.random.default_rng(3);reference=rng.uniform(.1,1.,(20,24))
    from scipy.ndimage import gaussian_filter
    reference=gaussian_filter(reference,1)
    frame=shift(reference,(.3,-.25),order=1,mode='constant')
    roi=np.zeros(reference.shape,bool);roi[3:-3,3:-3]=True
    d,_=estimate_displacement(np.broadcast_to(frame,(2,3,20,24)),reference,roi,'rigid',iterations=12)
    assert_allclose(d.mean(axis=(0,2,3)),[.3,-.25],atol=.04)


def test_axis2_truth_weights_are_diagnostic_and_gain_is_capped():
    a=scene();c=SharedContext(a.public,16,32)
    r=fit_shared(c,dict(method='axis2',horizon=32,rank=4,phase_width=.1,gain_cap=.1))
    assert max(b['unconstrained_operator_gain_max'] for b in r.diagnostics['blocks'])>.1
    assert max(b['operator_gain_max'] for b in r.diagnostics['blocks'])<=.1
    # Verify the actual action in the calibrated frame, rather than trusting
    # a reported gain. With no warp/centering each scored block is a contraction.
    assert np.linalg.norm(r.views/c.modulation[32:])<=.1*np.linalg.norm(c.aligned[32:])+1e-10
    diag=r.source_diagnostic(a.blood)
    assert diag['evaluation_only'] and not diag['selection_use'] and diag['solver_success']
    assert diag['best_diagonal_uncapped_nrmse']<=diag['heuristic_uncapped_nrmse']+1e-8
    assert all('ajd_converged' in b for b in r.diagnostics['blocks'])


def test_global_selection_ignores_held_out_truth_and_known_motion():
    def row(mid,e,b):return dict(id=mid,status='ok',metrics=dict(blood_iq_nrmse=e,refitted_background_to_tissue_power=b))
    methods=[dict(id=mid,family='f',method='svd',motion='none',horizon=16,evaluation_only=False) for mid in ['a','b']]
    methods.append(dict(id='oracle',family='f',method='svd',motion='known',horizon=16,evaluation_only=True))
    records=[dict(condition=c,seed=seed,split=split,results=[row('a',ea,.001),row('b',eb,.1),row('oracle',0,0)])
             for c,seed,split,ea,eb in [('c1',1,'training',.4,.1),('c2',1,'training',.8,.2),('c1',2,'evaluation',.9,0),('c2',2,'evaluation',.9,0)]]
    s=select_global_operating_points(records,methods,[.01,.00001])
    assert len(s)==3 and s[0]['selected_id']=='b' and s[1]['selected_id']=='a'
    assert s[2]['status']=='infeasible_on_training'
    other=copy.deepcopy(records)
    for case in other[2:]:case['results']=[row('a',0,0),row('b',100,100),row('oracle',0,0)]
    assert [x['selected_id'] for x in select_global_operating_points(other,methods,[.01,.00001])]==[x['selected_id'] for x in s]


def test_default_support_and_frozen_small_study(tmp_path):
    c=validate_shared_config(default_shared_config())
    assert len(c['methods'])==73 and len(c['conditions'])==8
    assert not set(c['training_seeds'])&set(c['evaluation_seeds'])
    c['training_seeds']=[1];c['evaluation_seeds']=[2];c['conditions']=c['conditions'][:1]
    c['case_workers']=1
    c['generator'].update(n_frames=64,calibration_frames=32,shape=[4,8]);c['score_start']=32
    c['methods']=[m for m in c['methods'] if m['id'] in ['axis1_none_h64_k1','weighted_svd_none_h64_k1','axis2_none_h64_r4_p0.1','stap_none_h64_loaded0.01']]
    for m in c['methods']:m['horizon']=32
    r=run_shared_study(c,tmp_path/'shared')
    assert len(r['records'])==2 and not r['evaluation_seed_truth_used_for_selection']
    assert (tmp_path/'shared'/'frozen_config.json').exists()
    for case in r['records']:
        for row in case['results']:
            assert row['status']=='ok' and row['metrics']['frozen_component_closure_nrmse']<1e-10
            assert row['metrics']['refitted_blood_response_nrmse_alpha0_1']>=0
