import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import linalg

from btbench.core import Acquisition,casorati,frames
from btbench.filters import svd_filter,aer_filter,SpatialMoments,local_filter,ssm_cutoff,rpca_filter
from btbench.axes import LagMemory,joint_diagonalize,compensate,restore,multilag_filter,motion_packet_filter
from btbench.measurements import doppler_velocity,vector_doppler,wall_shear_scalar,relative_error
from btbench.stap import steering,null_projector,lcmv_weights,local_stap
from btbench.synthetic import poiseuille
from btbench.io import save_acquisition,load_acquisition,inspect_openhrf
from btbench.benchmark import oracle,run_benchmark


def rand(shape,seed=4):
    rng=np.random.default_rng(seed)
    return rng.normal(size=shape)+1j*rng.normal(size=shape)


def test_casorati_roundtrip():
    a=rand((31,4,5))
    assert_allclose(frames(casorati(a),a.shape),a)


@pytest.mark.parametrize("centered",[False,True])
def test_svd_identity_and_centering(centered):
    a=rand((31,4,5))
    y=svd_filter(a,centered=centered,discard=0).iq
    ref=a-a.mean(0) if centered else a
    assert_allclose(y,ref,atol=1e-12)


def test_mean_subtraction_changes_rank_one():
    a=rand((1,4,5))*np.ones((31,1,1))
    assert np.linalg.norm(svd_filter(a,centered=True,discard=0).iq)<1e-12
    assert np.linalg.norm(svd_filter(a,centered=False,discard=0).iq)>1


@pytest.mark.parametrize("centered",[False,True])
def test_streamed_spatial_moments(centered):
    x=rand((13,41))+3
    m=SpatialMoments(13)
    for a in np.array_split(x,5,axis=1):m.update(a)
    xc=x-x.mean(1,keepdims=True) if centered else x
    assert_allclose(m.covariance(centered),xc@xc.conj().T/41,atol=1e-13)


@pytest.mark.parametrize("centered",[False,True])
def test_temporal_spatial_aer_equivalence(centered):
    # Same denominator and same regularization yield the same filtered data.
    a=rand((31,4,5))
    kw={"centered":centered,"normalization":"common","epsilon_rel":1e-3}
    y=aer_filter(a,domain="temporal",**kw).iq
    z=aer_filter(a,domain="spatial",stream_chunk=7,**kw).iq
    assert_allclose(y,z,atol=3e-12)


def test_zca_spectrum_and_covariance_rank():
    a=rand((31,4,5));x=casorati(a)
    y=casorati(aer_filter(a,domain="spatial",epsilon_rel=0).iq)
    assert_allclose(y@y.conj().T/31,np.eye(20),atol=2e-13)
    s=linalg.svdvals(y)
    assert_allclose(s,np.sqrt(31),atol=3e-13)


def test_regularized_spectrum():
    a=rand((31,4,5));x=casorati(a)
    r=aer_filter(a,domain="spatial",epsilon_abs=0.3)
    ev=linalg.eigvalsh(x@x.conj().T/31)
    y=casorati(r.iq)
    actual=linalg.eigvalsh(y@y.conj().T/31)
    assert_allclose(actual,ev/(ev+0.3),atol=2e-13)


def test_no_whitening_identity():
    a=rand((31,4,5))
    assert_allclose(aer_filter(a,alpha=0).iq,a,atol=2e-13)


def test_memory_cap():
    with pytest.raises(MemoryError):SpatialMoments(100,max_bytes=100)


def test_local_overlap_add_identity():
    a=rand((31,7,9))
    r=local_filter(a,patch_shape=(4,4),stride=(2,2),discard=0)
    assert_allclose(r.iq,a,atol=2e-13)
    assert r.diagnostics["overlap_min_weight"]>0


def test_confidence_memory_no_scale_drift():
    x=rand((3,40))
    memory=LagMemory((0,1),rho=0.9)
    for q in [0.2,0.9,0.3]:c=memory.update(x,q)
    assert_allclose(c[0],x@x.conj().T/40,atol=1e-13)
    assert_allclose(c[1],x[:,1:]@x[:,:-1].conj().T/39,atol=1e-13)


def test_complex_lag_centering():
    x=rand((3,40))+2+3j
    memory=LagMemory((0,2),rho=0)
    c=memory.update(x,1,True)
    a=x[:,2:];b=x[:,:-2]
    assert_allclose(c[2],(a-a.mean(1,keepdims=True))@(b-b.mean(1,keepdims=True)).conj().T/38,atol=1e-13)


def test_complex_joint_diagonalization():
    rng=np.random.default_rng(2)
    u,_=np.linalg.qr(rand((4,4)))
    mats=np.array([u@np.diag(rng.normal(size=4))@u.conj().T for _ in range(8)])
    q,info=joint_diagonalize(mats,max_sweeps=80)
    assert info["ajd_relative_offdiag"]<1e-7
    assert_allclose(q.conj().T@q,np.eye(4),atol=1e-12)


def test_phase_compensation_and_inverse():
    phi=np.cumsum(np.linspace(.01,.15,64))
    a=rand((1,4,5))*np.exp(1j*phi[:,None,None])
    demod,m=compensate(a,np.ones((4,5),bool))
    assert_allclose(demod,demod[:1]*np.ones((64,1,1)),atol=1e-12)
    recovered,valid=restore(demod,m)
    assert_allclose(recovered,a,atol=1e-12)
    assert valid.all()


def test_known_doppler_sign_and_scale():
    speed=0.023;c=1540;f0=5e6;prf=1000
    for sign in (-1,1):
        phase=sign*4*np.pi*f0*speed/c*np.arange(64)/prf
        a=np.exp(1j*phase)[:,None,None]
        v,_=doppler_velocity(a,prf,f0,c,sign)
        assert_allclose(v,speed,atol=1e-14)


def test_vector_doppler_geometry():
    tx=np.array([-.3,.3]);rx=tx
    truth=np.array([.01,.025]);prf=2000;f0=5e6
    a=steering(truth,tx,rx,np.arange(64)/prf,f0).reshape(64,2,1,1)
    v,_=vector_doppler(a,tx,rx,prf,f0)
    assert_allclose(v[0,0],truth,atol=1e-13)


@pytest.mark.parametrize("fit",["free","no_slip"])
def test_poiseuille_wss_exact_gradient(fit):
    z=np.linspace(0,.003,8);x=np.linspace(-.001,.001,101)
    v=np.broadcast_to(.02*(1-(x/.001)**2),(len(z),len(x)))
    wall={"positions_zx_m":[[.001,-.001],[.001,.001]],"inward_normals_zx":[[0,1],[0,-1]],
          "distances_m":[.0001,.0002,.0003,.0004],"viscosity_pa_s":.0035}
    tau=wall_shear_scalar(v,z,x,wall,fit)
    assert_allclose(tau,2*.0035*.02/.001,atol=1e-12)


def test_stap_projectors_and_lcmv():
    st=rand((12,2));sb=rand((12,1))
    c=rand((12,20));c=c@c.conj().T
    for cov in (None,c):
        p=null_projector(st,cov)
        assert_allclose(p@st,0,atol=2e-13)
        assert_allclose(p@p,p,atol=2e-13)
    w=lcmv_weights(sb,st,c)
    assert_allclose(w.conj()@sb[:,0],1,atol=2e-13)
    assert_allclose(w.conj()@st,0,atol=2e-13)


def test_stap_stationary_tissue_null():
    views=np.broadcast_to(rand((1,2,3,4)),(16,2,3,4)).copy()
    result=local_stap(views,1000,5e6,[0,.2],[0,.2],[[0,0]],taps=3)
    # Equal per-view gains assumed by steering; input unequal views leave a
    # view-differential residual, whose coherent mean cancels exactly.
    assert_allclose(result.iq,0,atol=1e-12)


def test_npz_roundtrip_no_pickle(tmp_path):
    a=poiseuille(n_frames=64,shape=(6,8))
    save_acquisition(tmp_path/"a.npz",a)
    b=load_acquisition(tmp_path/"a.npz")
    assert_allclose(a.iq,b.iq)
    assert b.metadata==a.metadata


def test_invalid_predictions_not_removed():
    assert np.isinf(relative_error(np.array([1,np.nan]),np.array([1,1])))


def test_oracle_objective_and_log():
    a=poiseuille(n_frames=64,shape=(6,8))
    _,diag=oracle(a,{"method":"svd","params":{}},{"discard":[0,1,2]},"velocity_nrmse")
    assert diag["selection_uses_test_truth"]
    assert len(diag["curve"])==3
    assert diag["selected"]["params"]["discard"] in (0,1,2)


def test_oracle_rejects_no_truth():
    a=Acquisition(rand((20,3,4)),1000,5e6)
    with pytest.raises(ValueError):oracle(a,{"method":"svd"},{"discard":[1]})


@pytest.mark.parametrize("fn",[motion_packet_filter,multilag_filter])
def test_axes_finite_and_lab_phase_restored(fn):
    a=poiseuille(n_frames=64,shape=(6,8))
    r=fn(a.iq,packet=32,patch_shape=(4,4),stride=(2,2),tissue_roi=a.tissue_roi,
         lags=(1,2,4,8),rank=4)
    assert r.iq.shape==a.iq.shape and np.isfinite(r.iq).all()
    assert r.diagnostics["lab_frame_phase_restored"]


def test_rpca_reconstruction_residual():
    a=rand((20,4,5))*.01
    a+=rand((1,4,5))
    r=rpca_filter(a,max_iter=400)
    assert r.diagnostics["converged"]


def test_reproducible_synthetic():
    assert_allclose(poiseuille(seed=8,n_frames=32).iq,poiseuille(seed=8,n_frames=32).iq)


def test_results_strict_json(tmp_path):
    import json
    a=poiseuille(n_frames=32,shape=(4,6))
    r=run_benchmark(a,{"methods":[{"id":"svd","method":"svd","params":{}}]},tmp_path)
    assert r["results"][0]["status"]=="ok"
    data=json.loads((tmp_path/"results.json").read_text())
    assert "input_sha256" in data
    with pytest.raises(FileExistsError):run_benchmark(a,{"methods":[]},tmp_path)


def test_h5_inspection_without_raw_reads(tmp_path):
    h5py=pytest.importorskip("h5py")
    p=tmp_path/"a.hdf5"
    with h5py.File(p,"w") as f:
        f.create_dataset("data/raw_data",shape=(30000,2,2048,128,1),dtype="f4",chunks=(1,1,2048,128,1))
        f.create_dataset("scan/center_frequency",data=5e6)
    r=inspect_openhrf(p)
    assert r["datasets"]["data/raw_data"]["decoded_bytes"]>1e9
    assert r["scalars"]["scan/center_frequency"]==5e6
