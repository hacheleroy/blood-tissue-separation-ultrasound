"""Bounded experimental processing pilot, with no promoted reference truth."""
from dataclasses import replace
from pathlib import Path
import hashlib
import json
import numpy as np
from .benchmark import run_benchmark,standard_methods
from .stap import estimate_tissue_velocity


def run_openhrf_pilot(a,output,calibration_frames=32,spatial_stride=2):
    """Disjoint calibration/test packets and matched angular-information arms.

    The brightest calibration quartile is an automatic nuisance ROI, not a
    verified tissue-only annotation. No blood/tissue component truth is created.
    The short duration is intended for pipeline verification, not physiology.
    """
    if (a.iq_views is None or any(v is None for v in (a.tx_angles_rad,a.rx_angles_rad,a.view_offsets_s))
            or calibration_frames<3 or len(a.iq)-calibration_frames<32 or spatial_stride<1):
        raise ValueError("Pilot needs views/geometry, >=3 calibration and >=32 evaluation frames")
    sl=(slice(None),slice(None),slice(None,None,spatial_stride),slice(None,None,spatial_stride))
    views=a.iq_views[sl]
    cal=views[:calibration_frames];test=views[calibration_frames:]
    power=np.mean(abs(cal)**2,axis=(0,1))
    roi=power>=np.quantile(power,.75)
    v=estimate_tissue_velocity(cal,roi,a.tx_angles_rad,a.rx_angles_rad,a.prf_hz,a.f0_hz,
        a.sound_speed_mps,a.phase_sign)
    # The implicit uniform clock starts at zero for the INPUT packet. Slicing
    # must retain its origin for sequential-view STAP/calibrated demodulation.
    clock=np.arange(len(a.iq))/a.prf_hz if a.timestamps_s is None else np.asarray(a.timestamps_s)
    times=clock[calibration_frames:]
    meta={**a.metadata,"pilot":"disjoint calibration and evaluation packets",
        "independent_truth_available":False,"calibration_frames":calibration_frames,
        "evaluation_frames":len(test),"calibration_roi":"brightest quartile, unverified tissue-only assumption",
        "evaluation_time_origin_s":float(times[0]),
        "estimated_tissue_velocity_zx_mps":v.tolist(),"spatial_stride":spatial_stride,
        "calibration_sha256":hashlib.sha256(cal.tobytes()).hexdigest()}
    evaluation=replace(a,iq=test.mean(axis=1),iq_views=test,timestamps_s=times,
        z_m=None if a.z_m is None else a.z_m[::spatial_stride],
        x_m=None if a.x_m is None else a.x_m[::spatial_stride],
        tissue_roi=roi,flow_mask=None,near_wall_mask=None,blood_true=None,tissue_true=None,
        blood_views_true=None,tissue_views_true=None,velocity_true_mps=None,wall=None,metadata=meta)
    methods=[{"id":"unfiltered","method":"identity"}]+standard_methods(packet=32,long_window=len(test),prf_hz=a.prf_hz)
    for name,params in [
        ("view_svd_matched",{"window":4,"hop":2}),
        ("view_svd_long",{"window":len(test),"hop":len(test)}),
        ("calibrated_view_svd_matched",{"window":4,"hop":2,"calibration_velocity_zx_mps":v.tolist()})]:
        methods.append({"id":name,"method":"view_svd","params":params})
    for adaptive in (False,True):
        params={"tissue_velocities_zx_mps":[v.tolist()],"taps":4,"hop":2,
                "adaptive":adaptive,"training_radius":3,"guard":1}
        if adaptive:params["training_roi"]=roi.tolist()
        methods.append({"id":"stap_covariance" if adaptive else "stap_orthogonal","method":"stap","params":params})
    config={"blas_threads":1,"measurement_window":32,"measurement_hop":32,"methods":methods,
        "protocol":"short experimental engineering pilot; no oracle or independent source truth",
        "calibration_frames":calibration_frames,"spatial_stride":spatial_stride}
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError("Choose a new pilot output directory")
    out.mkdir(parents=True,exist_ok=True)
    (out/"frozen_config.json").write_text(json.dumps(config,indent=2,allow_nan=False))
    report=run_benchmark(evaluation,config,out)
    report["config_sha256"]=hashlib.sha256((out/"frozen_config.json").read_bytes()).hexdigest()
    (out/"pilot_manifest.json").write_text(json.dumps({
        "config_sha256":report["config_sha256"],"code_sha256":report["code_sha256"],
        "input_sha256":report["input_sha256"],"acquisition":meta,
        "output_claim":"processing path and numerical diagnostics; no separation accuracy claim"},indent=2))
    return report
