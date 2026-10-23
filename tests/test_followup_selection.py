import copy
from btbench.followup_experiment import calibrate_operating_points,default_followup_config


def row(mid,blood,tissue):
    return {'id':mid,'status':'ok','metrics':{'blood_iq_nrmse':blood,'refitted_background_to_tissue_power':tissue}}


def test_matched_rejection_selects_on_calibration_only_and_reports_test_failure():
    methods=[{'id':'a','family':'f'},{'id':'b','family':'f'}]
    records=[{'condition':'c','split':'calibration','seed':11,'results':[row('a',.2,1e-5),row('b',.1,1e-2)]},
             {'condition':'c','split':'evaluation','seed':61,'results':[row('a',.9,.1),row('b',0,0)]}]
    selected=calibrate_operating_points(records,methods,[1e-4,1e-8])
    assert selected[0]['selected_id']=='a'
    assert not selected[0]['evaluation'][0]['meets_target']
    assert selected[1]['status']=='infeasible_on_calibration'
    changed=copy.deepcopy(records);changed[1]['results']=[row('a',0,0),row('b',1,1)]
    assert calibrate_operating_points(changed,methods,[1e-4])[0]['selected_id']=='a'


def test_followup_has_disjoint_seeds_explicit_lee_arms_and_unique_settings():
    c=default_followup_config()
    assert not set(c['calibration_seeds'])&set(c['evaluation_seeds'])
    assert len({m['id'] for m in c['methods']})==len(c['methods'])
    for domain in ('temporal','spatial'):
        assert any(m['method']=='lee_aer' and m['params']['domain']==domain and 'epsilon_abs' in m['params'] for m in c['methods'])
    assert any(m['method']=='song_pdf_local' and m['params']['normalization']=='paper' for m in c['methods'])


def test_followup_execution_records_refitted_primary_and_frozen_secondary_metrics(tmp_path):
    from btbench.followup_experiment import run_followup_study
    c=default_followup_config();c['calibration_seeds']=[11];c['evaluation_seeds']=[61]
    c['conditions']=c['conditions'][:1];c['generator']['n_frames']=64
    ids=['svd_short_c0_k1','axis1_memory_c0_k1','axis2_ridge_cap_c0_e0.0001',
         'lee_temporal_w32_relative0.001','lee_spatial_w32_relative0.001','song_pdf_paper']
    c['methods']=[m for m in c['methods'] if m['id'] in ids]
    assert len(c['methods'])==len(ids)
    r=run_followup_study(c,tmp_path/'run')
    assert len(r['records'])==2
    for case in r['records']:
        for row in case['results']:
            assert row['status']=='ok',row
            assert row['metrics']['refitted_background_to_tissue_power']>=0
            assert row['metrics']['frozen_component_closure_nrmse']<1e-11
    assert not r['evaluation_seed_truth_used_for_selection']


def test_experimental_followup_retains_clock_and_does_not_promote_truth(tmp_path):
    from dataclasses import replace
    from btbench.stap_experiment import angle_scene
    from btbench.followup_experiment import run_openhrf_followup,default_openhrf_followup_config
    c=default_openhrf_followup_config();c['acquisitions']=c['acquisitions'][:1]
    c['calibration_frames']=16;c['evaluation_frames']=32;c['evaluation_spatial_stride']=1
    c['methods']=[{'id':'unfiltered','method':'identity'}]
    a=angle_scene(n_frames=48,shape=(4,13)).acquisition
    a=replace(a,metadata={**a.metadata,'subject_id':'1'})
    report=run_openhrf_followup([a],c,tmp_path/'real')
    record=report['records'][0]
    assert record['acquisition_metadata']['evaluation_time_origin_s']==16/a.prf_hz
    assert not record['acquisition_metadata']['independent_truth_available']
    assert not record['results'][0]['metrics']['frozen_component_metrics_available']
    assert 'blood_iq_nrmse' not in record['results'][0]['metrics']
