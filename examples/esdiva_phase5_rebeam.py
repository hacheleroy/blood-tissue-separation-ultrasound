"""Reconstruct a wider tissue tile from verified phase-4 RF checkpoints.

The central published eight-row tile is retained as the scoring interval.
Additional rows only inform the motion fit. Isolated components and known
scatterer motion remain evaluator-only, outside the regular PublicInput.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, digest, relative_error, view_rf


def load_rf(folder, n, channels, samples):
    data = np.zeros((n, channels, samples), dtype=np.float32)
    t0 = np.empty(n)
    dt = None
    for event in range(n):
        with np.load(folder/f"rf_event_{event:04d}.npz") as packet:
            rf = packet["rf"]
            if rf.shape[0] != channels or rf.shape[1] > samples:
                raise ValueError("RF checkpoint shape mismatch")
            data[event,:,:rf.shape[1]] = rf
            t0[event] = float(packet["t0"])
            step = float(packet["dt"])
            if dt is not None and dt != step:
                raise ValueError("RF checkpoint sample spacing changed")
            dt = step
    return data, {"dt":dt,"t0_per_event":t0}


def rebeam(acquisition, upstream_root, output, z_start_mm=14.4, n_depth=32):
    acquisition,upstream_root,output = (Path(p).resolve() for p in
                                        (acquisition,upstream_root,output))
    if output.exists():
        raise FileExistsError(output)
    report = json.loads((acquisition/"report.json").read_text())
    head = subprocess.check_output(["git","-C",str(upstream_root),"rev-parse","HEAD"],text=True).strip()
    if (head != UPSTREAM_COMMIT or report["upstream_commit"] != head or
        report["view_indices"] != [4,3,5] or report["rf_shape"][0] != 288 or
        n_depth < 16 or z_start_mm < GRID_MM["z_extent"][0] or
        z_start_mm+(n_depth-1)*GRID_MM["dz"] > GRID_MM["z_extent"][1]):
        raise ValueError("Unsupported RF source or wider grid")
    n,ch,length = report["rf_shape"]
    tissue,coords = load_rf(acquisition/"rf_tissue",n,ch,length)
    blood,bcoords = load_rf(acquisition/"rf_blood",n,ch,length)
    if (digest(tissue) != report["rf_tissue_sha256"] or
        digest(blood) != report["rf_blood_sha256"] or
        digest(coords["t0_per_event"]) != report["rf_t0_sha256"] or
        not np.array_equal(coords["t0_per_event"],bcoords["t0_per_event"]) or
        coords["dt"] != bcoords["dt"]):
        raise ValueError("RF checkpoints disagree with the acquisition report")
    sys.path.insert(0,str(upstream_root/"examples/example22_flow_doppler"))
    import step1_define_flow_phantom as example
    from doppler_tools import rf2iq
    from esdiva.utilities import add_noise, rf_rms
    sim,tx,_ = example.build_simulator()
    del sim
    from step1_define_flow_phantom import compound_events
    _,nine_events = compound_events(tx,1)
    events = [nine_events[i] for i in report["view_indices"]]
    # IQ demodulation uses a Hilbert transform along the full depth grid.
    # Cropping before it changes even the interior eight rows.
    grid = GRID_MM
    clean = tissue+blood
    mixed = add_noise(clean,6.,reference=rf_rms(clean),
                      rng=np.random.default_rng(report["seed"]+101))
    noise = mixed-clean
    if digest(noise) != report["rf_noise_sha256"]:
        raise ValueError("RF noise replay disagrees with the acquisition")
    ncal = report["calibration_frames"]*3
    cal_beam,z_mm,x_mm = view_rf(tissue[:ncal]+noise[:ncal],
                                {"dt":coords["dt"],"t0_per_event":coords["t0_per_event"][:ncal]},
                                events,tx,grid,example.C)
    eval_coords = {"dt":coords["dt"],"t0_per_event":coords["t0_per_event"][ncal:]}
    mix_beam,_,_ = view_rf(mixed[ncal:],eval_coords,events,tx,grid,example.C)
    tissue_beam,_,_ = view_rf(tissue[ncal:],eval_coords,events,tx,grid,example.C)
    blood_beam,_,_ = view_rf(blood[ncal:],eval_coords,events,tx,grid,example.C)
    fs = example.C/(2*GRID_MM["dz"]*1e-3)
    carrier = report["carrier_hz"]
    broad = {"calibration_views":rf2iq(cal_beam,fs=fs,fc=carrier,axis=-2),
             "views":rf2iq(mix_beam,fs=fs,fc=carrier,axis=-2),
             "tissue":rf2iq(tissue_beam,fs=fs,fc=carrier,axis=-2),
             "blood":rf2iq(blood_beam,fs=fs,fc=carrier,axis=-2)}
    start = int(np.argmin(abs(z_mm-z_start_mm)))
    if abs(z_mm[start]-z_start_mm)>1e-10 or start+n_depth>len(z_mm):
        raise ValueError("Wider tile is off the original reconstruction grid")
    broad = {k:v[:,:,start:start+n_depth].copy() for k,v in broad.items()}
    z_mm = z_mm[start:start+n_depth]
    with np.load(acquisition/"iq_sources.npz") as old:
        old_iq = {k:old[k] for k in ("views","calibration_views","tissue","blood","noise","z_mm","x_mm")}
    core_start = int(np.argmin(abs(z_mm-old_iq["z_mm"][0])))
    core = slice(core_start,core_start+len(old_iq["z_mm"]))
    if (len(z_mm) != n_depth or core.stop > n_depth or
        np.max(abs(z_mm[core]-old_iq["z_mm"])) > 1e-10 or
        not np.array_equal(x_mm,old_iq["x_mm"])):
        raise ValueError("Wider tile does not contain the published score tile")
    core_errors = {k:relative_error(v[:,:,core],old_iq[k]) for k,v in broad.items()}
    if max(core_errors.values()) > 3e-5:
        raise ValueError(f"Central IQ replay mismatch: {core_errors}")
    broad["noise"] = broad["views"]-broad["tissue"]-broad["blood"]
    # Preserve the published IQ bit-for-bit for all core processing.
    for k in ("views","calibration_views","tissue","blood","noise"):
        broad[k][:,:,core] = old_iq[k]
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,**broad,z_mm=z_mm,x_mm=x_mm,
                        core_rows=np.arange(core.start,core.stop))
    audit={"seed":report["seed"],"motion":report["motion"],"flow_scale":report["flow_scale"],
           "source_report_sha256":hashlib.sha256((acquisition/"report.json").read_bytes()).hexdigest(),
           "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "broad_shape":list(broad["views"].shape),"core_rows":list(range(core.start,core.stop)),
           "core_replay_relative_errors":core_errors,
           "rf_sources_and_noise_verified":True,"source_truth_evaluator_only":True,
           "broad_iq_sha256":{k:digest(v) for k,v in broad.items()}}
    output.with_suffix(".json").write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n")
    return audit


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",required=True,type=Path)
    p.add_argument("--upstream-root",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--z-start-mm",default=14.4,type=float)
    p.add_argument("--n-depth",default=32,type=int)
    a=p.parse_args()
    result=rebeam(a.acquisition,a.upstream_root,a.output,a.z_start_mm,a.n_depth)
    print(json.dumps({k:result[k] for k in ("seed","motion","broad_shape","core_replay_relative_errors")}))
