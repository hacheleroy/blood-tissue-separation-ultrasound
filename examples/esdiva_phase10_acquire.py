"""Acquire phase-10 source-resolved eSDIva RF with a guarded tissue core.

The eSDIva Example 22 cloud is kept at a declared particle fraction. Blood
travels along its original tilted vessel at a controlled velocity multiplier;
tissue follows an explicitly prescribed axial/lateral trajectory in millimetres.
The calibration frames have stationary tissue. The core training annotation
uses the fixed 5.5 mm vessel guard, giving 40 tissue pixels. Isolated
sources remain evaluator-only. A source failure is reported and archived
without erasing the acquired IQ; the scorer refuses to filter it.
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

from esdiva_angle_calibration import angular_diagnostics, calibration_gate
from esdiva_example22_pilot import (GRID_MM, UPSTREAM_COMMIT, digest,
                                    numerical_source_hash, relative_error, view_rf)

VIEW_INDICES = (4, 3, 5)


def tissue_shift_mm(positions, eval_time_fraction, motion):
    """Return (x,z) scatterer shifts in mm at one slow-time frame."""
    if motion == "none" or eval_time_fraction < 0:
        return np.zeros((len(positions), 2))
    phase = 2*np.pi*eval_time_fraction
    dx = .30*np.sin(2*phase)*np.ones(len(positions))
    dz = .06*np.sin(phase)*np.ones(len(positions))
    if motion == "smooth":
        dx *= 1 + .20*(positions[:, 2]-16.)/3.
        dz *= 1 + .25*positions[:, 0]/8.
    elif motion != "rigid":
        raise ValueError("Motion must be none, rigid, or smooth")
    return np.column_stack((dx, dz))


def make_trajectories(example, n_events, emission_prf, seed, fraction,
                      calibration_frames, evaluation_frames, motion, flow_scale):
    if not 0 <= flow_scale <= 1 or not 0 < fraction <= 1:
        raise ValueError("Invalid flow scale or particle fraction")
    old_seed = example.SEED
    try:
        example.SEED = seed
        positions, amps, truth = example.build_phantom(n_events, emission_prf)
    finally:
        example.SEED = old_seed
    source_tissue = truth["n_tissue"]
    if fraction < 1:
        rng = np.random.default_rng(seed+1)
        keep_t = np.sort(rng.choice(source_tissue, max(1, round(source_tissue*fraction)), replace=False))
        keep_b = np.sort(rng.choice(truth["n_blood"], max(1, round(truth["n_blood"]*fraction)), replace=False))
    else:
        keep_t = np.arange(source_tissue)
        keep_b = np.arange(truth["n_blood"])
    selected = np.r_[keep_t, source_tissue+keep_b]
    positions = positions[:, selected].copy()
    amps = amps[selected]
    n_tissue = len(keep_t)
    if flow_scale != 1:
        p0 = positions[0, n_tissue:].astype(float)
        relative = p0-example.VESSEL_CENTER_MM
        s0 = relative@example.FLOW_DIR
        transverse = relative-s0[:, None]*example.FLOW_DIR
        t = np.arange(n_events)/emission_prf
        s = s0[None]+1e3*flow_scale*t[:, None]*truth["v_blood"][keep_b][None]
        s = (s+example.VESSEL_LENGTH_MM/2) % example.VESSEL_LENGTH_MM-example.VESSEL_LENGTH_MM/2
        positions[:, n_tissue:] = (example.VESSEL_CENTER_MM[None, None]+
                                   s[:, :, None]*example.FLOW_DIR[None, None]+transverse[None]).astype(np.float32)
    if motion != "none":
        base = positions[0, :n_tissue].copy()
        for event in range(n_events):
            frame = event//len(VIEW_INDICES)-calibration_frames
            fraction_of_test = frame/evaluation_frames if frame >= 0 else -1
            shift = tissue_shift_mm(base, fraction_of_test, motion)
            positions[event, :n_tissue, 0] = base[:, 0]+shift[:, 0]
            positions[event, :n_tissue, 2] = base[:, 2]+shift[:, 1]
    return positions, amps, n_tissue, truth["v_blood"][keep_b]


def acquire(upstream_root, output, seed=2026, scatterer_fraction=.1,
            calibration_frames=16, evaluation_frames=96, motion="none", flow_scale=1.):
    upstream_root, output = Path(upstream_root).resolve(), Path(output).resolve()
    script = upstream_root/"examples/example22_flow_doppler/step1_define_flow_phantom.py"
    head = subprocess.check_output(["git", "-C", str(upstream_root), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_COMMIT or not script.is_file():
        raise ValueError("Expected pinned upstream Example 22")
    if (calibration_frames < 8 or evaluation_frames < 80 or
            evaluation_frames % 16 or motion not in ("none", "rigid", "smooth")):
        raise ValueError("Need disjoint calibration and at least 80 whole evaluation frames")
    if (output/"report.json").exists():
        raise FileExistsError(f"Refusing to overwrite completed acquisition {output}")
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(script.parent))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq
    from esdiva.utilities import add_noise, rf_rms

    n_views = len(VIEW_INDICES)
    n_frames = calibration_frames+evaluation_frames
    n_events = n_frames*n_views
    emission_prf = example.PRF_DOPPLER*n_views
    sim, tx, _ = example.build_simulator()
    _, nine_events = example.compound_events(tx, 1)
    frame_events = [nine_events[i] for i in VIEW_INDICES]
    events = frame_events*n_frames
    positions, amps, n_tissue, v_blood = make_trajectories(
        example, n_events, emission_prf, seed, scatterer_fraction,
        calibration_frames, evaluation_frames, motion, flow_scale)
    print(f"Phase 10: {n_events} events, {n_tissue} tissue + {len(v_blood)} blood scatterers, "
          f"motion={motion}, flow scale={flow_scale}", flush=True)
    amps_t = amps.copy();amps_t[n_tissue:] = 0.
    amps_b = amps.copy();amps_b[:n_tissue] = 0.
    rf_t, ct = sim.sequence_rf(positions, amps_t, events, out_path=output/"rf_tissue")
    rf_b, cb = sim.sequence_rf(positions, amps_b, events, out_path=output/"rf_blood")
    if (rf_t.shape != rf_b.shape or ct["dt"] != cb["dt"] or
            not np.array_equal(ct["t0_per_event"], cb["t0_per_event"])):
        raise AssertionError("Component RF event grids differ")
    check_event = calibration_frames*n_views+n_views
    rf_check, cc = sim.sequence_rf(positions[check_event:check_event+1], amps,
                                   events[check_event:check_event+1])
    if (rf_check.shape[-1] > rf_t.shape[-1] or
            not np.array_equal(cc["t0_per_event"], ct["t0_per_event"][check_event:check_event+1])):
        raise AssertionError("Combined RF event grid differs")
    rf_error = relative_error(rf_t[check_event:check_event+1, :, :rf_check.shape[-1]]+
                              rf_b[check_event:check_event+1, :, :rf_check.shape[-1]], rf_check)
    if rf_error > 2e-4:
        raise AssertionError(f"RF source additivity failed: {rf_error}")
    rf_clean = rf_t+rf_b
    rf_mix = add_noise(rf_clean, 6., reference=rf_rms(rf_clean),
                       rng=np.random.default_rng(seed+101))
    rf_noise = rf_mix-rf_clean
    stop_cal = calibration_frames*n_views
    beam_cal, z_mm, x_mm = view_rf(rf_t[:stop_cal]+rf_noise[:stop_cal],
                                  {"dt": ct["dt"], "t0_per_event": ct["t0_per_event"][:stop_cal]},
                                  frame_events, tx, GRID_MM, example.C)
    coord_eval = {"dt": ct["dt"], "t0_per_event": ct["t0_per_event"][stop_cal:]}
    beam_mix, _, _ = view_rf(rf_mix[stop_cal:], coord_eval, frame_events, tx, GRID_MM, example.C)
    beam_t, _, _ = view_rf(rf_t[stop_cal:], coord_eval, frame_events, tx, GRID_MM, example.C)
    beam_b, _, _ = view_rf(rf_b[stop_cal:], coord_eval, frame_events, tx, GRID_MM, example.C)
    beam_clean_first, _, _ = view_rf(rf_clean[stop_cal:stop_cal+3],
                                    {"dt": ct["dt"],
                                     "t0_per_event": ct["t0_per_event"][stop_cal:stop_cal+3]},
                                    frame_events, tx, GRID_MM, example.C)
    beam_error = relative_error(beam_t[:1]+beam_b[:1], beam_clean_first)
    if beam_error > 2e-4:
        raise AssertionError(f"Beamformer additivity failed: {beam_error}")
    fs_depth = example.C/(2*GRID_MM["dz"]*1e-3)
    f_echo = echo_center_frequency(np.moveaxis(beam_cal, -2, -1),1/fs_depth,
                                   band=(.3*example.FC,1.7*example.FC))
    cal = rf2iq(beam_cal, fs=fs_depth, fc=f_echo, axis=-2)
    mixture = rf2iq(beam_mix, fs=fs_depth, fc=f_echo, axis=-2)
    tissue = rf2iq(beam_t, fs=fs_depth, fc=f_echo, axis=-2)
    blood = rf2iq(beam_b, fs=fs_depth, fc=f_echo, axis=-2)
    iq_error = relative_error(tissue[:1]+blood[:1],
                              rf2iq(beam_clean_first, fs=fs_depth, fc=f_echo, axis=-2))
    if iq_error > 2e-4:
        raise AssertionError(f"IQ additivity failed: {iq_error}")
    first_z = int(np.argmin(abs(z_mm-15.6)))
    rows = slice(first_z,first_z+8)
    z_mm = z_mm[rows]
    cal,mixture,tissue,blood = (a[:, :, rows, :] for a in (cal,mixture,tissue,blood))
    if mixture.shape != (evaluation_frames,n_views,8,16):
        raise AssertionError(f"Unexpected cropped IQ shape {mixture.shape}")
    noise = mixture-tissue-blood
    training_roi = (abs(x_mm[None]-np.tan(np.deg2rad(60.))*
                        (z_mm[:,None]-16.)) >= 5.5)
    if training_roi.shape != (8,16) or training_roi.sum() != 40:
        raise AssertionError("Changed phase-10 anatomical tissue core")
    flow_roi = np.broadcast_to(abs(x_mm[None]) <= 2.,training_roi.shape).copy()
    blood_to_tissue = float(np.mean(abs(blood[:, :, training_roi])**2)/
                            max(np.mean(abs(tissue[:, :, training_roi])**2),1e-30))
    legacy_roi = np.broadcast_to(abs(x_mm[None]) >= 6.,training_roi.shape)
    legacy_blood_to_tissue = float(np.mean(abs(blood[:, :, legacy_roi])**2)/
                                   max(np.mean(abs(tissue[:, :, legacy_roi])**2),1e-30))
    source_pass = blood_to_tissue <= .1
    angle_rows = angular_diagnostics(cal[:calibration_frames//2],
                                     cal[calibration_frames//2:],training_roi,VIEW_INDICES)
    angle_spec, angle_pass = calibration_gate(angle_rows,VIEW_INDICES)
    geometry = {"prf_hz": example.PRF_DOPPLER,"f0_hz":float(f_echo),
                "sound_speed_mps":example.C,"phase_sign":-1,
                "tx_rad":np.deg2rad([example.ANGLES_DEG[i] for i in VIEW_INDICES]).tolist(),
                "rx_rad":[0.]*n_views,"view_offsets_s":(np.arange(n_views)/emission_prf).tolist()}
    # This geometric displacement refers to scatterer transport, not an exact
    # beamformed-IQ warp. It is never included in PublicInput or a gate.
    known_motion = np.zeros((evaluation_frames,2,8,16))
    for frame in range(evaluation_frames):
        zz,xx=np.meshgrid(z_mm,x_mm,indexing="ij")
        shift=tissue_shift_mm(np.column_stack((xx.ravel(),np.zeros(xx.size),zz.ravel())),
                              frame/evaluation_frames,motion)
        known_motion[frame,0]=shift[:,1].reshape(8,16)/GRID_MM["dz"]
        known_motion[frame,1]=shift[:,0].reshape(8,16)/GRID_MM["dx"]
    np.savez_compressed(output/"iq_sources.npz", views=mixture, calibration_views=cal,
                        blood=blood,tissue=tissue,noise=noise,z_mm=z_mm,x_mm=x_mm,
                        training_roi=training_roi,flow_roi=flow_roi,
                        known_scatterer_displacement_pixels=known_motion)
    report = {"status":("phase-10 source-resolved development acquisition" if source_pass
                        else "phase-10 source eligibility failed; no filter scores"),
              "upstream_commit":head,"seed":seed,"view_indices":VIEW_INDICES,
              "scatterer_fraction":scatterer_fraction,"calibration_frames":calibration_frames,
              "evaluation_frames":evaluation_frames,"motion":motion,"flow_scale":flow_scale,
              "n_tissue":n_tissue,"n_blood":len(v_blood),
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "pilot_adapter_sha256":hashlib.sha256(Path(__file__).with_name("esdiva_example22_pilot.py").read_bytes()).hexdigest(),
              "benchmark_source_sha256":numerical_source_hash(Path(__file__).resolve().parents[1]),
              "environment":{"python":platform.python_version(),"esdiva":importlib.metadata.version("esdiva")},
              "grid_mm":GRID_MM,"geometry":geometry,"carrier_hz":float(f_echo),
              "rf_shape":rf_t.shape,"rf_tissue_sha256":digest(rf_t),"rf_blood_sha256":digest(rf_b),
              "rf_t0_sha256":digest(ct["t0_per_event"]),"rf_noise_sha256":digest(rf_noise),
              "rf_source_additivity_relative_error":rf_error,"beam_source_additivity_relative_error":beam_error,
              "iq_source_additivity_relative_error":iq_error,
              "iq_shape":mixture.shape,"mixture_sha256":digest(mixture),"calibration_sha256":digest(cal),
              "blood_sha256":digest(blood),"tissue_sha256":digest(tissue),"noise_sha256":digest(noise),
              "known_motion_sha256":digest(known_motion),"training_roi_sha256":digest(training_roi),
              "training_blood_to_tissue_power_ratio":blood_to_tissue,
              "legacy_core_blood_to_tissue_power_evaluator_only":legacy_blood_to_tissue,
              "source_eligibility_pass":source_pass,"source_eligibility_limit":.1,
              "calibration_angle_diagnostics":angle_rows,"calibration_gate":angle_spec,
              "calibration_gate_pass":angle_pass,
              "calibration_gate_pass_is_legacy_diagnostic_only":True,
              "source_and_known_motion_are_evaluator_only":True,
              "limitation":"Particle subset and one scatterer seed; analytic scatterer motion is not exact beamformed-IQ warp"}
    (output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False,
                                 default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+"\n")
    return report


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--upstream-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--seed",type=int,default=2026)
    p.add_argument("--scatterer-fraction",type=float,default=.1)
    p.add_argument("--calibration-frames",type=int,default=16)
    p.add_argument("--evaluation-frames",type=int,default=80)
    p.add_argument("--motion",choices=("none","rigid","smooth"),default="none")
    p.add_argument("--flow-scale",type=float,default=1.)
    a=p.parse_args()
    r=acquire(a.upstream_root,a.output,a.seed,a.scatterer_fraction,
              a.calibration_frames,a.evaluation_frames,a.motion,a.flow_scale)
    print(json.dumps({k:r[k] for k in ("status","rf_shape","iq_shape",
                                      "calibration_gate_pass","rf_source_additivity_relative_error")},indent=2))
