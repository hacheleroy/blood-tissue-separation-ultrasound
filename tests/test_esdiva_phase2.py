"""Check particle trajectories and disjoint calibration before acoustic RF."""

from pathlib import Path
import json
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"examples"))
from esdiva_phase2_acquire import make_trajectories, tissue_shift_mm  # noqa: E402
from esdiva_phase2_freeze import sha256, verify_freeze  # noqa: E402
from esdiva_example22_pilot import numerical_source_hash  # noqa: E402


class SmallVessel:
    SEED=2026
    VESSEL_CENTER_MM=np.array([0.,0.,16.])
    FLOW_DIR=np.array([1.,0.,0.])
    VESSEL_LENGTH_MM=24.

    @classmethod
    def build_phantom(cls,n_events,prf_emission):
        rng=np.random.default_rng(cls.SEED)
        tissue=np.array([[0.,0.,15.],[2.,0.,17.]])
        blood=np.array([[1.,0.,16.],[-1.,0.,16.]])
        velocities=np.array([.01,.02])
        positions=[]
        for event in range(n_events):
            shift=1e3*velocities*event/prf_emission
            positions.append(np.vstack([tissue,blood+np.column_stack((shift,np.zeros((2,2))))]))
        amps=rng.normal(size=4).astype(np.float32)
        return np.asarray(positions,dtype=np.float32),amps,{
            "n_tissue":2,"n_blood":2,"v_blood":velocities}


def test_phase2_keeps_disjoint_stationary_calibration_and_controls_flow():
    original=SmallVessel.SEED
    kwargs=dict(n_events=18,emission_prf=3000.,seed=37,fraction=1.,
                calibration_frames=2,evaluation_frames=4)
    moving,amps,n_tissue,_=make_trajectories(SmallVessel,motion="smooth",flow_scale=.02,**kwargs)
    zero,amps0,_,_=make_trajectories(SmallVessel,motion="none",flow_scale=0.,**kwargs)
    assert SmallVessel.SEED==original and n_tissue==2
    np.testing.assert_array_equal(amps,amps0)
    np.testing.assert_array_equal(moving[:6,:n_tissue],
                                  np.broadcast_to(moving[0,:n_tissue],moving[:6,:n_tissue].shape))
    assert np.linalg.norm(moving[9,:n_tissue]-moving[0,:n_tissue])>0
    np.testing.assert_allclose(zero[:,n_tissue:],
                               np.broadcast_to(zero[0:1,n_tissue:],zero[:,n_tissue:].shape),atol=1e-6)
    assert np.linalg.norm(moving[-1,n_tissue:]-moving[0,n_tissue:])>0


def test_rigid_and_smooth_tissue_motion_have_declared_mm_geometry():
    p=np.array([[-5.,0.,15.],[5.,0.,17.]])
    np.testing.assert_array_equal(tissue_shift_mm(p,-1,"smooth"),np.zeros((2,2)))
    rigid=tissue_shift_mm(p,.125,"rigid")
    smooth=tissue_shift_mm(p,.125,"smooth")
    np.testing.assert_allclose(rigid[0],rigid[1])
    assert rigid[0,0]>.29 and rigid[0,1]>.04
    assert smooth[0,0]!=smooth[1,0] and smooth[0,1]!=smooth[1,1]


def test_freeze_rejects_mutated_settings_and_unreserved_seed(tmp_path):
    root=tmp_path
    (root/"src/btbench").mkdir(parents=True)
    (root/"src/btbench/x.py").write_text("def f(): return 1\n")
    (root/"configs").mkdir()
    scenarios={"reserved_heldout_seeds":[2031],"scenes":[{"id":"zero","motion":"smooth","flow_scale":0.}],
               "upstream_commit":"pinned","view_indices":[4,3,5],"particle_fraction":.1,
               "calibration_frames":16,"evaluation_frames":80}
    scenario_path=root/"configs/esdiva_phase2_scenarios.json"
    scenario_path.write_text(json.dumps(scenarios))
    manifest={"upstream_commit":"pinned",
              "numerical_source_sha256":numerical_source_hash(root),
              "files_sha256":{"configs/esdiva_phase2_scenarios.json":sha256(scenario_path)}}
    freeze=root/"freeze.json"
    freeze.write_text(json.dumps(manifest))
    assert verify_freeze(root,freeze,2031,"smooth",0.)==manifest
    with pytest.raises(ValueError,match="Seed is not"):
        verify_freeze(root,freeze,2053,"smooth",0.)
    scenarios["scenes"][0]["flow_scale"]=.02
    scenario_path.write_text(json.dumps(scenarios))
    with pytest.raises(ValueError,match="Frozen file changed"):
        verify_freeze(root,freeze,2031,"smooth",0.)
