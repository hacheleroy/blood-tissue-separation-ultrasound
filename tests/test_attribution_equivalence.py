import numpy as np
from numpy.testing import assert_allclose
from btbench.filters import svd_filter,aer_filter
from btbench.core import Acquisition
from btbench.followup_experiment import refitted_response_metrics


def scene():
    rng=np.random.default_rng(8)
    b=rng.normal(size=(64,4,5))+1j*rng.normal(size=(64,4,5))
    t=30*(rng.normal(size=(1,4,5))+1j*rng.normal(size=(1,4,5)))*np.exp(.01j*np.arange(64))[:,None,None]
    return Acquisition(b+t,1000,5e6,blood_true=b,tissue_true=t)


def test_equal_svd_outputs_can_have_different_frozen_blood_attribution():
    a=scene();right=svd_filter(a.iq,discard=1);left=svd_filter(a.iq,discard=1,projection_side='spatial')
    assert_allclose(left.iq,right.iq,atol=5e-13)
    assert np.linalg.norm(left.apply_frozen(a.blood_true)-right.apply_frozen(a.blood_true))>.1*np.linalg.norm(a.blood_true)
    d=[]
    for side,result in [('temporal',right),('spatial',left)]:
        spec={'method':'svd','params':{'discard':1,'projection_side':side}}
        d.append(refitted_response_metrics(a,spec,result))
    for key in d[0]:assert_allclose(d[0][key],d[1][key],atol=2e-11)


def test_equivalent_aer_outputs_and_refitted_responses_agree_across_sides():
    a=scene();results=[];diagnostics=[]
    for domain in ['temporal','spatial']:
        params={'domain':domain,'centered':True,'normalization':'common','epsilon_rel':1e-3}
        r=aer_filter(a.iq,**params);results.append(r)
        diagnostics.append(refitted_response_metrics(a,{'method':'aer','params':params},r))
    assert_allclose(results[0].iq,results[1].iq,atol=2e-11)
    assert np.linalg.norm(results[0].apply_frozen(a.blood_true)-results[1].apply_frozen(a.blood_true))>1
    for key in diagnostics[0]:assert_allclose(diagnostics[0][key],diagnostics[1][key],atol=2e-11)


def test_axis1_without_memory_equals_local_svd_at_identical_support():
    from btbench.axes import motion_packet_filter
    from btbench.filters import local_filter
    a=scene();kw={'patch_shape':(4,5),'stride':(4,5),'discard':1}
    left=motion_packet_filter(a.iq,packet=32,rho=0.,q_min=0.,centered=False,**kw)
    right=local_filter(a.iq,method='svd',window=32,centered=False,**kw)
    assert_allclose(left.iq,right.iq,atol=1e-12)
    assert np.linalg.norm(left.apply_frozen(a.blood_true)-right.apply_frozen(a.blood_true))>1
