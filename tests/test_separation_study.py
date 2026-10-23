import copy
import json
import numpy as np
from numpy.testing import assert_allclose
from btbench.synthetic import advected_scatterers
from btbench.measurements import injection_response,doppler_velocity
from btbench.separation_experiment import default_separation_config,run_separation_study


def test_transport_phase_has_declared_sign_and_speed():
    for sign in (-1,1):
        a=advected_scatterers(seed=3,n_frames=12,shape=(5,7),umax_mps=.002,
            n_blood=1,n_tissue=10,uniform_flow=True,phase_sign=sign)
        v,_=doppler_velocity(a.blood_true,a.prf_hz,a.f0_hz,phase_sign=sign)
        power=np.mean(abs(a.blood_true)**2,axis=0)
        observable=power>power.max()*1e-10
        assert observable.any()
        assert_allclose(v[observable],.002,atol=1e-14)
        assert a.metadata["fluid_structure_coupling"] is False


def test_transport_static_control_and_reproducibility():
    kw=dict(seed=8,n_frames=8,shape=(4,7),umax_mps=0.,tissue_velocity_mps=0.,n_blood=30,n_tissue=30)
    a=advected_scatterers(**kw);b=advected_scatterers(**kw)
    assert_allclose(a.iq,b.iq)
    assert_allclose(a.blood_true,np.broadcast_to(a.blood_true[0],a.blood_true.shape))
    assert_allclose(a.tissue_true,np.broadcast_to(a.tissue_true[0],a.tissue_true.shape))


def test_counterfactual_injection_is_separately_labeled():
    blood=np.ones((4,2,2),complex);background=np.ones_like(blood)*3
    result=injection_response(lambda x:x,background+blood,blood)
    assert all(r["incremental_blood_iq_nrmse"]<1e-14 for r in result["responses"])
    nonlinear=injection_response(lambda x:x*abs(x),background+blood,blood)
    assert nonlinear["responses"][0]["incremental_blood_iq_nrmse"]!=nonlinear["responses"][1]["incremental_blood_iq_nrmse"]
    assert "no additive" in nonlinear["interpretation"]


def test_small_study_has_frozen_config_and_both_model_families(tmp_path):
    cfg=copy.deepcopy(default_separation_config())
    cfg.update(seeds=[17],methods=[{"id":"identity","method":"identity"}],
        statistical_conditions=cfg["statistical_conditions"][:1],
        scatterer_conditions=[{"label":"static","umax_mps":0.,"tissue_velocity_mps":0.,"n_blood":20,"n_tissue":20}])
    cfg["generator"].update(n_frames=32,shape=[4,7])
    r=run_separation_study(cfg,tmp_path/"run")
    assert len(r["records"])==2
    assert {c["family"] for c in r["records"]}=={"statistical","psf_transport"}
    assert all(c["results"][0]["metrics"]["frozen_blood_iq_nrmse"]==0 for c in r["records"])
    assert json.loads((tmp_path/"run/frozen_config.json").read_text())["seeds"]==[17]
