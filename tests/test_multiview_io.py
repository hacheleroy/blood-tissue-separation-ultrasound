import json
import numpy as np
import pytest
from numpy.testing import assert_allclose
from btbench.io import export_openhrf_rf_window,load_acquisition,selected_view_timing
from btbench.benchmark import run_benchmark
from btbench.cli import main
from btbench.measurements import doppler_velocity,frozen_component_metrics
from btbench.stap_experiment import angle_scene
from btbench.stap import local_stap
from btbench.experimental import run_openhrf_pilot


def point_rf_file(path):
    h5py=pytest.importorskip("h5py")
    fs=40e6;f0=5e6;c=1540.;depth=.008;prf=3000.;speed=.02
    tx=np.array([-.2,.2]);t=np.arange(1024)/fs
    raw=np.empty((16,2,1024,1,1))
    for i in range(16):
        for j in range(2):
            d=speed*(i/prf+j/(2*prf))
            delay=(1+np.cos(tx[j]))*(depth+d)/c
            raw[i,j,:,0,0]=np.exp(-((t-delay)/.4e-6)**2)*np.cos(2*np.pi*f0*(t-delay))
    with h5py.File(path,"w") as f:
        p="tracks/track_0/"
        f[p+"data/raw_data"]=raw
        f["probe/probe_geometry"]=np.zeros((1,3))
        for k,v in {"sampling_frequency":fs,"center_frequency":f0,"sound_speed":c,
                    "polar_angles":tx,"initial_times":np.zeros(2),
                    "time_to_next_transmit":np.full((16,2),1/(2*prf))}.items():
            f[p+"scan/"+k]=v
    return {"time_reference_verified":True,"full_receive_aperture_verified":True,
        "frames_start_stop":[0,16],"transmit_indices":[0,1],"root":"tracks/track_0",
        "tx_angles_rad":tx.tolist(),"initial_times_s":[0,0],"tx_reference_times_s":[0,0],
        "z_m":[depth],"x_m":[0],"phase_sign":-1,"revisit_prf_hz":prf}


def test_rf_npz_standard_stap_runner_and_phase(tmp_path):
    source=tmp_path/"rf.h5";g=point_rf_file(source)
    a=export_openhrf_rf_window(source,g,tmp_path/"iq.npz")
    loaded=load_acquisition(tmp_path/"iq.npz")
    assert_allclose(loaded.iq_views,a.iq_views)
    assert_allclose(loaded.iq,a.iq_views.mean(axis=1))
    assert_allclose(loaded.view_offsets_s,[0,1/6000])
    for j,angle in enumerate(a.tx_angles_rad):
        v,_=doppler_velocity(a.iq_views[:,j],a.prf_hz,a.f0_hz,phase_sign=-1,
                             geometry_factor=1+np.cos(angle))
        assert_allclose(v,.02,rtol=.025)
    report=run_benchmark(loaded,{"methods":[{"id":"stap","method":"stap",
        "params":{"tissue_velocities_zx_mps":[[0,0]],"taps":4}}]},tmp_path/"run")
    assert report["results"][0]["status"]=="ok"
    assert not report["results"][0]["metrics"]["frozen_component_metrics_available"]


def test_timing_gaps_and_declared_wrong_prf_rejected(tmp_path):
    dt=np.full((12,2),1/6000);dt[5,1]=.1
    with pytest.raises(ValueError,match="Gapped"):
        selected_view_timing(dt,[0,1])
    p=tmp_path/"rf.h5";g=point_rf_file(p);g["revisit_prf_hz"]=1000
    with pytest.raises(ValueError,match="PRF"):
        export_openhrf_rf_window(p,g,tmp_path/"unused.npz")


def test_stap_component_metrics_use_view_truth():
    s=angle_scene(n_frames=32,shape=(4,13))
    r=local_stap(s.acquisition.iq_views,1000,5e6,s.geometry["tx_rad"],s.geometry["rx_rad"],
        [[.001,0]],view_offsets_s=s.geometry["view_offsets_s"],taps=4)
    m=frozen_component_metrics(r,s.acquisition)
    assert m["frozen_component_input_domain"]=="uncompounded_views_to_coherent_mean"
    assert m["frozen_component_closure_nrmse"]<1e-12


def test_config_cli_uses_actual_revisit_rate(tmp_path):
    path=tmp_path/"cfg.json"
    main(["write-config",str(path),"--prf-hz","3000","--packet","32","--long-window","96"])
    methods=json.loads(path.read_text())["methods"]
    hp=next(m for m in methods if m["method"]=="highpass")
    assert hp["params"]["prf_hz"]==3000
    assert hp["params"]["cutoff_hz"]==30


def test_pilot_retains_clock_origin_and_never_promotes_component_truth(tmp_path):
    a=angle_scene(n_frames=64,shape=(4,13)).acquisition
    assert a.timestamps_s is None
    r=run_openhrf_pilot(a,tmp_path/"pilot",calibration_frames=32,spatial_stride=1)
    assert_allclose(r["acquisition"]["evaluation_time_origin_s"],.032)
    assert all(not row["oracle"] for row in r["results"])
    assert all(not row.get("metrics",{}).get("frozen_component_metrics_available",False) for row in r["results"])
    assert all(row["status"]=="ok" for row in r["results"])
