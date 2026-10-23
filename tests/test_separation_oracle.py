import json
import numpy as np
import pytest
from btbench.core import Acquisition
from btbench.benchmark import oracle,OracleInfeasible,pareto_front,run_benchmark


def orthogonal_scene():
    time=np.arange(32)
    b=np.exp(2j*np.pi*time/8)[:,None,None]*np.array([[[1.,-1.]]])
    tissue=np.ones_like(b)*20
    return Acquisition(b+tissue,1000.,5e6,blood_true=b,tissue_true=tissue)


def test_preservation_constrained_oracle_rejects_silence():
    a=orthogonal_scene()
    r,d=oracle(a,{"method":"svd"},{"discard":[0,1,2]},
        objective="frozen_tissue_power_retention",
        constraints={"frozen_blood_iq_nrmse":{"max":.1}})
    assert d["selected"]["params"]["discard"]==1
    assert not d["curve"][-1]["feasible"]
    assert np.linalg.norm(r.iq-a.blood_true)<1e-11


def test_infeasible_search_preserves_candidate_records(tmp_path):
    spec={"id":"impossible","method":"svd","oracle_grid":{"discard":[2]},
          "objective":"frozen_tissue_power_retention",
          "constraints":{"frozen_blood_iq_nrmse":{"max":.1}}}
    with pytest.raises(OracleInfeasible) as exc:
        oracle(orthogonal_scene(),spec,spec["oracle_grid"],spec["objective"],constraints=spec["constraints"])
    assert exc.value.diagnostics["feasible_candidates"]==0
    report=run_benchmark(orthogonal_scene(),{"methods":[spec]},tmp_path/"run")
    row=report["results"][0]
    assert row["status"]=="infeasible" and row["oracle"]
    assert len(json.loads((tmp_path/"run/results.json").read_text())["results"][0]["curve"])==1


def test_pareto_excludes_dominated_and_unavailable_points():
    rows=[{"metrics":{"frozen_blood_iq_nrmse":x,"frozen_tissue_power_retention":y}}
          for x,y in [(0,1),(.5,.2),(1,.5),(1,0),(None,None)]]
    assert pareto_front(rows)==[0,1,3]
