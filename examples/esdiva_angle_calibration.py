"""Measure cross-angle tissue calibration on disjoint eSDIva RF frames.

This is a development diagnostic. It uses tissue-only channel RF and no blood
truth to compare the regression gain with a phase plus RMS gain. The latter
preserves each view's measured amplitude when angular speckle decorrelates.

Usage:
  python examples/esdiva_angle_calibration.py --upstream-root /path/to/eSDIva \
      --output /path/to/new/control.json --view-indices 3 4 5
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, view_rf


def _hash(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def angular_diagnostics(calibration, validation, roi, indices):
    """Audit tissue decorrelation and scaling without fitting to validation."""
    reference = int(np.flatnonzero(np.asarray(indices) == 4)[0])
    if not (roi.any() and calibration.shape[1] == len(indices) and
            validation.shape[1] == len(indices)):
        raise ValueError("Require a 0-degree reference and nonempty tissue ROI")
    a = calibration[:, :, roi].transpose(1, 0, 2).reshape(len(indices), -1)
    b = validation[:, :, roi].transpose(1, 0, 2).reshape(len(indices), -1)
    ref = a[reference]
    ref_power = np.vdot(ref, ref).real
    rows = []
    for j, index in enumerate(indices):
        cross = np.vdot(ref, a[j])
        regression = cross / ref_power
        coherence = abs(cross) / np.sqrt(ref_power * np.vdot(a[j], a[j]).real)
        rms = np.sqrt(np.vdot(a[j], a[j]).real / ref_power)
        bounded = rms * np.exp(1j * np.angle(cross))
        val_cross = np.vdot(b[reference], b[j])
        val_coherence = abs(val_cross) / np.sqrt(np.vdot(b[reference], b[reference]).real *
                                                 np.vdot(b[j], b[j]).real)
        rows.append({"angle_deg": float(-20 + 5 * index),
                     "gain_regression_magnitude": float(abs(regression)),
                     "gain_rms_magnitude": float(rms),
                     "coherence_train": float(coherence),
                     "coherence_validation": float(val_coherence),
                     "relative_phase_train_rad": float(np.angle(cross)),
                     "relative_phase_validation_rad": float(np.angle(val_cross)),
                     "regression_input_noise_power_gain": float(1 / abs(regression)**2),
                     "bounded_input_noise_power_gain": float(1 / abs(bounded)**2),
                     "regression_residual_validation": float(np.linalg.norm(b[j] - regression*b[reference])**2 /
                                                           np.vdot(b[j], b[j]).real),
                     "bounded_residual_validation": float(np.linalg.norm(b[j] - bounded*b[reference])**2 /
                                                        np.vdot(b[j], b[j]).real)})
    return rows


def calibration_gate(rows, indices, coherence_min=.7, max_noise_power_gain=4.):
    """Reject unsupported reference geometry before comparing blood filters.

    The bound applies to scalar regression normalization, upstream of the
    filter's own noise gain. It is a development acceptance rule, not an
    assertion that the output will preserve blood or reject tissue.
    """
    specification = {"coherence_min": coherence_min,
                     "max_input_noise_power_gain": max_noise_power_gain,
                     "reference_view_must_be_zero_deg": True}
    passes = bool(tuple(indices)[0] == 4 and all(
        row["coherence_train"] >= coherence_min and
        row["coherence_validation"] >= coherence_min and
        row["regression_input_noise_power_gain"] <= max_noise_power_gain
        for row in rows))
    return specification, passes


def run(upstream_root, output, indices=(3, 4, 5), scatterer_fraction=1.,
        calibration_frames=8, validation_frames=8):
    upstream_root, output = Path(upstream_root).resolve(), Path(output).resolve()
    script = upstream_root / "examples/example22_flow_doppler/step1_define_flow_phantom.py"
    head = subprocess.check_output(["git", "-C", str(upstream_root), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_COMMIT or not script.is_file():
        raise ValueError("Expected the pinned eSDIva Example 22 source")
    indices = tuple(indices)
    if (len(indices) < 3 or len(set(indices)) != len(indices) or
            4 not in indices or any(i < 0 or i > 8 for i in indices) or
            min(calibration_frames, validation_frames) < 4 or
            not 0 < scatterer_fraction <= 1):
        raise ValueError("Invalid view, frame or density control")
    if output.exists():
        raise FileExistsError(f"Refusing to replace {output}")
    sys.path.insert(0, str(script.parent))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq

    start = time.monotonic()
    sim, tx, _ = example.build_simulator()
    _, nine_events = example.compound_events(tx, 1)
    frame_events = [nine_events[i] for i in indices]
    n_frames = calibration_frames + validation_frames
    positions, amps, truth = example.build_phantom(n_frames*len(indices),
                                                   example.PRF_DOPPLER*len(indices))
    n_tissue = truth["n_tissue"]
    if scatterer_fraction < 1:
        rng = np.random.default_rng(2027)
        keep = np.sort(rng.choice(n_tissue, max(1, round(n_tissue*scatterer_fraction)), replace=False))
    else:
        keep = np.arange(n_tissue)
    positions = positions[:, keep]
    amps = amps[keep]
    print(f"Calibration: {n_frames*len(indices)} tissue RF events, {len(keep)} scatterers", flush=True)
    # Zero blood rows are not needed because this diagnostic acquires tissue only.
    rf, coords = sim.sequence_rf(positions, amps, frame_events*n_frames)
    beam, z_mm, x_mm = view_rf(rf, coords, frame_events, tx, GRID_MM, example.C)
    fs_depth = example.C / (2*GRID_MM["dz"]*1e-3)
    f_echo = echo_center_frequency(np.moveaxis(beam[:calibration_frames], -2, -1),
                                   1/fs_depth, band=(.3*example.FC, 1.7*example.FC))
    iq = rf2iq(beam, fs=fs_depth, fc=f_echo, axis=-2)
    first_z = int(np.argmin(abs(z_mm-15.6)))
    iq = iq[:, :, first_z:first_z+8, :]
    roi = np.broadcast_to(abs(x_mm[None, :]) >= 6., iq.shape[2:])
    rows = angular_diagnostics(iq[:calibration_frames], iq[calibration_frames:], roi, indices)
    result = {"status": "development tissue-only calibration control; no filter ranking",
              "upstream_commit": head, "upstream_seed": example.SEED,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "esdiva_version": importlib.metadata.version("esdiva"),
              "view_indices": indices, "angles_deg": [float(example.ANGLES_DEG[i]) for i in indices],
              "scatterer_fraction": scatterer_fraction, "n_tissue": len(keep),
              "calibration_frames": calibration_frames, "validation_frames": validation_frames,
              "carrier_hz": float(f_echo), "rf_shape": rf.shape, "rf_sha256": _hash(rf),
              "iq_shape": iq.shape, "iq_sha256": _hash(iq),
              "t0_sha256": _hash(coords["t0_per_event"]),
              "tissue_roi_pixels": int(roi.sum()), "diagnostics": rows,
              "validation_shares_stationary_tissue_scatterers": True,
              "calibration_gate": calibration_gate(rows, indices)[0],
              "calibration_gate_pass": calibration_gate(rows, indices)[1],
              "elapsed_s": time.monotonic()-start}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False,
                                 default=lambda x: list(x) if isinstance(x, tuple) else x.tolist())+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--view-indices", nargs="+", type=int, default=(3, 4, 5))
    parser.add_argument("--scatterer-fraction", type=float, default=1.)
    parser.add_argument("--calibration-frames", type=int, default=8)
    parser.add_argument("--validation-frames", type=int, default=8)
    args = parser.parse_args()
    report = run(args.upstream_root, args.output, args.view_indices,
                 args.scatterer_fraction, args.calibration_frames, args.validation_frames)
    print(json.dumps({"elapsed_s": report["elapsed_s"], "diagnostics": report["diagnostics"]}, indent=2))
