"""Angle-resolved eSDIva Example 22 RF pilot; optional upstream dependency.

Usage:
  python examples/esdiva_example22_pilot.py --upstream-root /path/to/eSDIva \
      --output /path/to/new/pilot

The default 16 calibration + 32 evaluation frames and a 10% scatterer subset
exercise the three-view processing contract cheaply. They are NOT a
long-ensemble benchmark, realistic speckle assessment, or held-out study.
RF checkpoints and IQ arrays stay in the output directory, outside Git.
"""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np


UPSTREAM_COMMIT = "ab7eaef6d23428762ab5b9079ab03a5aa99af085"
VIEW_INDICES = (1, 4, 7)  # -15, 0, +15 degrees in the nine-angle source sequence
GRID_MM = {
    # The oblique 4-mm vessel occupies much of a narrow lateral tile. Include
    # outer tissue pixels so the geometric training ROI has an actual guard.
    "x_extent": [-8.0, 8.0], "y_extent": [-0.05, 0.05],
    # Keep enough fast-time samples for the received-frequency estimate and
    # Hilbert transform; crop the common eight-row score tile afterward.
    "z_extent": [13.0, 19.0], "dx": 1.0, "dy": 0.1, "dz": 0.1,
}


def digest(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def numerical_source_hash(root):
    h = hashlib.sha256()
    for p in sorted((root / "src/btbench").glob("*.py")):
        h.update(p.name.encode()); h.update(p.read_bytes())
    return h.hexdigest()


def relative_error(actual, expected):
    return float(np.linalg.norm((actual - expected).ravel()) /
                 max(np.linalg.norm(expected.ravel()), 1e-30))


def view_rf(rf, coords, frame_events, probe, grid, c):
    """Beamform each emission separately, preserving (time,view,z,x)."""
    from esdiva.beamforming import das_volume

    n_view = len(frame_events)
    if rf.ndim != 3 or len(rf) % n_view or np.shape(coords["t0_per_event"]) != (len(rf),):
        raise ValueError("RF, event origins and complete view cycles must match")
    frames = []
    for start in range(0, len(rf), n_view):
        stop = start + n_view
        # The upstream DAS sums events; one call PER event retains the view axis.
        block = []
        for j in range(start, stop):
            v, axes = das_volume(
                rf[j:j+1],
                {"dt": coords["dt"], "t0_per_event": coords["t0_per_event"][j:j+1]},
                [frame_events[j-start]], probe, grid, c=c, fnum=1.0,
            )
            block.append(v[:, 0, :].T)  # DAS returns (x,y,z)
        frames.append(np.stack(block))
    return np.stack(frames), np.asarray(axes["z_mm"]), np.asarray(axes["x_mm"])


def make_pilot(upstream_root, output, calibration_frames=16, evaluation_frames=32,
               scatterer_fraction=.1):
    upstream_root = Path(upstream_root).resolve()
    output = Path(output).resolve()
    script = upstream_root / "examples/example22_flow_doppler/step1_define_flow_phantom.py"
    if not script.is_file():
        raise FileNotFoundError(f"Missing upstream Example 22 at {script}")
    head = subprocess.check_output(["git", "-C", str(upstream_root), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_COMMIT:
        raise ValueError(f"Upstream commit {head} differs from pinned {UPSTREAM_COMMIT}")
    if calibration_frames < 8 or evaluation_frames < 32 or evaluation_frames % 16:
        raise ValueError("Pilot requires >=8 calibration and >=32 whole evaluation frames")
    if not 0 < scatterer_fraction <= 1:
        raise ValueError("scatterer_fraction must be in (0,1]")
    if (output / "report.json").exists():
        raise FileExistsError(f"Refusing to replace completed pilot in {output}")
    output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(script.parent))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq
    from esdiva.utilities import add_noise, rf_rms

    n_views = len(VIEW_INDICES)
    n_frames = calibration_frames + evaluation_frames
    n_events = n_frames * n_views
    emission_prf = example.PRF_DOPPLER * n_views
    sim, tx, _ = example.build_simulator()
    _, nine_events = example.compound_events(tx, 1)
    frame_events = [nine_events[i] for i in VIEW_INDICES]
    events = frame_events * n_frames
    positions, amps, truth = example.build_phantom(n_events, emission_prf)
    n_tissue = truth["n_tissue"]
    if scatterer_fraction < 1:
        rng = np.random.default_rng(2027)
        keep_t = rng.choice(n_tissue, max(1, round(n_tissue*scatterer_fraction)), replace=False)
        keep_b = n_tissue + rng.choice(truth["n_blood"],
                                        max(1, round(truth["n_blood"]*scatterer_fraction)), replace=False)
        selected = np.sort(np.r_[keep_t, keep_b])
        positions = positions[:, selected]
        amps = amps[selected]
        n_tissue = len(keep_t)
    print(f"Pilot: {n_events} events, {n_tissue} tissue + {len(amps)-n_tissue} blood scatterers "
          f"({scatterer_fraction:.3g} of upstream; not a speckle-fidelity study)", flush=True)
    amps_t = amps.copy(); amps_t[n_tissue:] = 0.
    amps_b = amps.copy(); amps_b[:n_tissue] = 0.
    rf_t, ct = sim.sequence_rf(positions, amps_t, events, out_path=output / "rf_tissue")
    rf_b, cb = sim.sequence_rf(positions, amps_b, events, out_path=output / "rf_blood")
    if (rf_t.shape != rf_b.shape or ct["dt"] != cb["dt"] or
            not np.array_equal(ct["t0_per_event"], cb["t0_per_event"])):
        raise AssertionError("Component RF grids differ")
    # Independently check the simulator's source superposition on an event.
    rf_check, cc = sim.sequence_rf(positions[:1], amps, events[:1])
    if not np.array_equal(cc["t0_per_event"], ct["t0_per_event"][:1]):
        raise AssertionError("Combined RF time origin differs")
    if rf_t.shape[-1] < rf_check.shape[-1]:
        raise AssertionError("Full RF trace is shorter than the one-event control")
    rf_error = relative_error(rf_t[:1, :, :rf_check.shape[-1]] +
                              rf_b[:1, :, :rf_check.shape[-1]], rf_check)
    if rf_error > 2e-4:
        raise AssertionError(f"RF source additivity failed: {rf_error:.3g}")

    rf_clean = rf_t + rf_b
    # One absolute receiver-noise scale and one draw over calibration and test.
    rf_mix = add_noise(rf_clean, 6.0, reference=rf_rms(rf_clean), rng=np.random.default_rng(11))
    rf_noise = rf_mix - rf_clean
    stop_cal = calibration_frames * n_views
    beam_cal, z_mm, x_mm = view_rf(rf_t[:stop_cal] + rf_noise[:stop_cal],
                                  {"dt": ct["dt"], "t0_per_event": ct["t0_per_event"][:stop_cal]},
                                  frame_events, tx, GRID_MM, example.C)
    coords_eval = {"dt": ct["dt"], "t0_per_event": ct["t0_per_event"][stop_cal:]}
    beam_mix, _, _ = view_rf(rf_mix[stop_cal:], coords_eval, frame_events, tx, GRID_MM, example.C)
    beam_t, _, _ = view_rf(rf_t[stop_cal:], coords_eval, frame_events, tx, GRID_MM, example.C)
    beam_b, _, _ = view_rf(rf_b[stop_cal:], coords_eval, frame_events, tx, GRID_MM, example.C)
    beam_clean, _, _ = view_rf(rf_clean[stop_cal:], coords_eval, frame_events, tx, GRID_MM, example.C)
    beam_error = relative_error(beam_t + beam_b, beam_clean)
    if beam_error > 2e-4:
        raise AssertionError(f"Beamformer source additivity failed: {beam_error:.3g}")
    if beam_mix.shape[:2] != (evaluation_frames, n_views) or beam_mix.shape[-1] != 16:
        raise AssertionError(f"Unexpected view/grid layout: {beam_mix.shape}")
    fs_depth = example.C / (2 * GRID_MM["dz"] * 1e-3)
    f_echo = echo_center_frequency(np.moveaxis(beam_cal, -2, -1), 1/fs_depth,
                                   band=(.3*example.FC, 1.7*example.FC))
    # A single calibration-derived carrier and identical operation on all RF.
    cal = rf2iq(beam_cal, fs=fs_depth, fc=f_echo, axis=-2)
    mixture = rf2iq(beam_mix, fs=fs_depth, fc=f_echo, axis=-2)
    tissue = rf2iq(beam_t, fs=fs_depth, fc=f_echo, axis=-2)
    blood = rf2iq(beam_b, fs=fs_depth, fc=f_echo, axis=-2)
    iq_clean = rf2iq(beam_clean, fs=fs_depth, fc=f_echo, axis=-2)
    iq_error = relative_error(tissue + blood, iq_clean)
    if iq_error > 2e-4:
        raise AssertionError(f"IQ source additivity failed: {iq_error:.3g}")
    first_z = int(np.argmin(abs(z_mm - 15.6)))
    rows = slice(first_z, first_z + 8)
    if len(z_mm[rows]) != 8:
        raise AssertionError("Score tile extends beyond the beamformed grid")
    z_mm = z_mm[rows]
    cal = cal[:, :, rows, :]
    mixture = mixture[:, :, rows, :]
    tissue = tissue[:, :, rows, :]
    blood = blood[:, :, rows, :]
    if mixture.shape != (evaluation_frames, n_views, 8, 16):
        raise AssertionError(f"Unexpected cropped IQ layout: {mixture.shape}")
    noise = mixture - tissue - blood
    x = x_mm[None, :]
    training_roi = np.broadcast_to(abs(x) >= 6., (len(z_mm), len(x_mm))).copy()
    flow_roi = np.broadcast_to(abs(x) <= 2., training_roi.shape).copy()
    if not training_roi.any() or not flow_roi.any() or np.any(training_roi & flow_roi):
        raise AssertionError("Geometric ROIs overlap or are empty")
    # Truth is used only by the evaluator to reject a contaminated annotation;
    # no source component is passed into PublicInput or the filtering API.
    training_blood_to_tissue = float(np.mean(abs(blood[:, :, training_roi])**2) /
                                     max(np.mean(abs(tissue[:, :, training_roi])**2), 1e-30))
    if training_blood_to_tissue > .1:
        raise AssertionError(f"Training ROI receives too much blood: {training_blood_to_tissue:.3g}")
    flow_blood_to_tissue = float(np.mean(abs(blood[:, :, flow_roi])**2) /
                                 max(np.mean(abs(tissue[:, :, flow_roi])**2), 1e-30))
    tissue_lag_phase = np.angle(np.sum(tissue[1:] * tissue[:-1].conj(), axis=(0, 2, 3)))
    blood_lag_phase = np.angle(np.sum(blood[1:] * blood[:-1].conj(), axis=(0, 2, 3)))
    angles = np.deg2rad([example.ANGLES_DEG[i] for i in VIEW_INDICES])
    offsets = np.arange(n_views) / emission_prf
    geometry = {"prf_hz": example.PRF_DOPPLER, "f0_hz": float(f_echo),
                "sound_speed_mps": example.C, "phase_sign": -1,
                "tx_rad": angles.tolist(), "rx_rad": np.zeros(n_views).tolist(),
                "view_offsets_s": offsets.tolist()}
    from btbench.shared import PublicInput, SharedContext, fit_shared

    public = PublicInput(mixture, cal, training_roi, geometry)
    context = SharedContext(public, packet=16, score_start=16, motion="none")
    specs = {
        "identity": {"method": "identity", "horizon": 16},
        "svd_cut1": {"method": "svd", "horizon": 16, "discard": 1},
        "axis1_cut1": {"method": "axis1", "horizon": 16, "discard": 1, "rho": .75},
        "axis2": {"method": "axis2", "horizon": 16, "rank": 4, "gain_cap": 1.},
        "axis3_orthogonal": {"method": "stap", "horizon": 16, "adaptive": False},
        "axis3_loaded": {"method": "stap", "horizon": 16, "adaptive": True, "loading_rel": .01},
        "lee_aer": {"method": "lee_aer", "horizon": 16, "domain": "temporal", "epsilon_abs": 26.},
    }
    records = {}
    for name, spec in specs.items():
        try:
            result = fit_shared(context, spec)
            b = result.apply_frozen(blood)
            t = result.apply_frozen(tissue)
            score = slice(context.score_start, None)
            records[name] = {
                "blood_iq_nrmse": relative_error(result.views, blood[score]),
                "blood_frozen_power_ratio": float(np.mean(abs(b)**2) / max(np.mean(abs(blood[score])**2), 1e-30)),
                "tissue_frozen_power_ratio": float(np.mean(abs(t)**2) / max(np.mean(abs(tissue[score])**2), 1e-30)),
                "fitted_result_finite": bool(np.isfinite(result.views).all()),
                "diagnostics": result.diagnostics,
            }
        except (ValueError, np.linalg.LinAlgError) as exc:
            records[name] = {"status": "failed", "error": str(exc)}

    np.savez_compressed(output / "pilot_iq.npz", views=mixture, calibration_views=cal,
                        blood=blood, tissue=tissue, noise=noise, z_mm=z_mm, x_mm=x_mm,
                        training_roi=training_roi, flow_roi=flow_roi)
    report = {
        "status": "development plumbing pilot; no method ranking or held-out validation",
        "upstream_commit": head, "upstream_seed": example.SEED,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "benchmark_source_sha256": numerical_source_hash(Path(__file__).resolve().parents[1]),
        "environment": {"python": platform.python_version(),
                        "packages": {name: importlib.metadata.version(name) for name in
                                     ("esdiva", "numpy", "scipy", "numba", "h5py", "pyvista")}},
        "selected_angles_deg": np.rad2deg(angles).tolist(), "grid_mm": GRID_MM,
        "scatterer_fraction": scatterer_fraction,
        "n_tissue": n_tissue, "n_blood": len(amps)-n_tissue,
        "calibration_frames": calibration_frames, "evaluation_frames": evaluation_frames,
        "scored_frames": evaluation_frames-context.score_start,
        "geometry": geometry, "rf_shape": list(rf_t.shape), "view_shape": list(mixture.shape),
        "rf_tissue_sha256": digest(rf_t), "rf_blood_sha256": digest(rf_b),
        "rf_t0_sha256": digest(ct["t0_per_event"]), "rf_sample_period_s": float(ct["dt"]),
        "rf_noise_sha256": digest(rf_noise), "rf_source_additivity_relative_error": rf_error,
        "beam_source_additivity_relative_error": beam_error,
        "iq_source_additivity_relative_error": iq_error,
        "iq_carrier_hz_from_calibration": float(f_echo),
        "training_pixels": int(training_roi.sum()), "flow_pixels": int(flow_roi.sum()),
        "training_blood_to_tissue_power_ratio": training_blood_to_tissue,
        "flow_blood_to_tissue_power_ratio": flow_blood_to_tissue,
        "tissue_lag1_phase_rad_by_view": tissue_lag_phase.tolist(),
        "blood_lag1_phase_rad_by_view": blood_lag_phase.tolist(),
        "calibration_sha256": digest(cal), "mixture_sha256": digest(mixture),
        "calibration": context.calibration, "methods": records,
        "limitation": "reduced scatterer density, stationary tissue, fast flow, short history, one scene/seed; no filter ranking, motion or slow-flow claim",
    }
    (output / "report.json").write_text(json.dumps(
        report, indent=2, allow_nan=False,
        default=lambda value: value.tolist() if isinstance(value, np.ndarray) else value.item()
    ) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--calibration-frames", type=int, default=16)
    parser.add_argument("--evaluation-frames", type=int, default=32)
    parser.add_argument("--scatterer-fraction", type=float, default=.1)
    args = parser.parse_args()
    result = make_pilot(args.upstream_root, args.output, args.calibration_frames,
                        args.evaluation_frames, args.scatterer_fraction)
    print(json.dumps({k: result[k] for k in ("status", "rf_shape", "view_shape",
                                               "rf_source_additivity_relative_error",
                                               "iq_carrier_hz_from_calibration")}, indent=2))


if __name__ == "__main__":
    main()
