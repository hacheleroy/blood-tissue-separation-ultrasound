"""Survey independent eSDIva scatterer clouds before selecting new views.

Each angle is physically acquired as channel RF. Stationary tissue and flowing
blood use separate RF simulations; one common 6 dB receiver-noise draw is
scaled by their summed RF, as in phase 2. The tissue-plus-noise calibration is
split into disjoint eight-frame train and validation blocks. This development
survey reports input geometry only and never computes filter scores.
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

from esdiva_angle_calibration import angular_diagnostics
from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, digest, view_rf


def survey(upstream_root, output, seed, fraction=.1, indices=(4, 0, 1, 2, 3, 5, 6, 7, 8),
           frames=16):
    """Measure spatially guarded angular coherence at declared particle density."""
    upstream_root, output = Path(upstream_root).resolve(), Path(output).resolve()
    if (seed in (2026, 2031, 2053) or frames < 8 or frames % 2 or
            not 0 < fraction <= 1 or len(set(indices)) != len(indices) or
            indices[0] != 4 or any(i not in range(9) for i in indices)):
        raise ValueError("Require new development seed, 0-degree reference and even frames")
    if output.exists():
        raise FileExistsError(output)
    head = subprocess.check_output(["git", "-C", str(upstream_root), "rev-parse", "HEAD"], text=True).strip()
    source = upstream_root / "examples/example22_flow_doppler/step1_define_flow_phantom.py"
    if head != UPSTREAM_COMMIT or not source.is_file():
        raise ValueError("Expected pinned eSDIva Example 22")
    sys.path.insert(0, str(source.parent))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq
    from esdiva.utilities import add_noise, rf_rms

    t0 = time.monotonic()
    sim, tx, _ = example.build_simulator()
    _, nine_events = example.compound_events(tx, 1)
    events = [nine_events[i] for i in indices]
    old_seed = example.SEED
    try:
        example.SEED = seed
        positions, amps, truth = example.build_phantom(frames*len(indices),
                                                       example.PRF_DOPPLER*len(indices))
    finally:
        example.SEED = old_seed
    source_tissue = truth["n_tissue"]
    rng = np.random.default_rng(seed+1)
    if fraction < 1:
        keep_t = np.sort(rng.choice(source_tissue, round(source_tissue*fraction), replace=False))
        keep_b = np.sort(rng.choice(truth["n_blood"], round(truth["n_blood"]*fraction), replace=False))
    else:
        keep_t, keep_b = np.arange(source_tissue), np.arange(truth["n_blood"])
    selected = np.r_[keep_t, source_tissue+keep_b]
    positions, amps = positions[:, selected], amps[selected]
    n_tissue = len(keep_t)
    print(f"Development seed {seed}: {frames*len(indices)} RF events, "
          f"{n_tissue} tissue + {len(keep_b)} blood scatterers", flush=True)
    amps_t = amps.copy(); amps_t[n_tissue:] = 0
    amps_b = amps.copy(); amps_b[:n_tissue] = 0
    rf_t, ct = sim.sequence_rf(positions, amps_t, events*frames)
    rf_b, cb = sim.sequence_rf(positions, amps_b, events*frames)
    if rf_t.shape != rf_b.shape or ct["dt"] != cb["dt"] or not np.array_equal(
            ct["t0_per_event"], cb["t0_per_event"]):
        raise AssertionError("Source RF time grids differ")
    rf_clean = rf_t+rf_b
    rf_cal = rf_t+add_noise(rf_clean, 6., reference=rf_rms(rf_clean),
                            rng=np.random.default_rng(seed+101))-rf_clean
    coords = {"dt": ct["dt"], "t0_per_event": ct["t0_per_event"]}
    beam_cal, z_mm, x_mm = view_rf(rf_cal, coords, events, tx, GRID_MM, example.C)
    beam_b, _, _ = view_rf(rf_b, coords, events, tx, GRID_MM, example.C)
    beam_t, _, _ = view_rf(rf_t, coords, events, tx, GRID_MM, example.C)
    fs_depth = example.C/(2*GRID_MM["dz"]*1e-3)
    carrier = echo_center_frequency(np.moveaxis(beam_cal, -2, -1), 1/fs_depth,
                                    band=(.3*example.FC, 1.7*example.FC))
    first_z = int(np.argmin(abs(z_mm-15.6)))
    rows = slice(first_z, first_z+8)
    cal = rf2iq(beam_cal, fs=fs_depth, fc=carrier, axis=-2)[:, :, rows, :]
    blood = rf2iq(beam_b, fs=fs_depth, fc=carrier, axis=-2)[:, :, rows, :]
    tissue = rf2iq(beam_t, fs=fs_depth, fc=carrier, axis=-2)[:, :, rows, :]
    rois = {"outer_both": abs(x_mm) >= 6.,
            "outer_left": x_mm <= -6., "outer_right": x_mm >= 6.,
            "wide_both": abs(x_mm) >= 5.}
    diagnostics = {}
    for name, lateral in rois.items():
        roi = np.broadcast_to(lateral[None, :], cal.shape[2:])
        diagnostics[name] = {
            "pixels": int(roi.sum()),
            "blood_to_tissue_power": float(np.mean(abs(blood[:, :, roi])**2)/
                                           max(np.mean(abs(tissue[:, :, roi])**2), 1e-30)),
            "angles": angular_diagnostics(cal[:frames//2], cal[frames//2:], roi, indices),
        }
    result = {"status": "new-seed development angle survey; no filter ranking",
              "seed": seed, "upstream_commit": head, "script_sha256": hashlib.sha256(
                  Path(__file__).read_bytes()).hexdigest(),
              "esdiva_version": importlib.metadata.version("esdiva"),
              "view_indices": list(indices), "angles_deg": [float(example.ANGLES_DEG[i]) for i in indices],
              "scatterer_fraction": fraction, "n_tissue": n_tissue, "n_blood": len(keep_b),
              "frames": frames, "calibration_train_frames": frames//2,
              "calibration_validation_frames": frames//2,
              "carrier_hz": float(carrier), "rf_shape": rf_t.shape,
              "rf_tissue_sha256": digest(rf_t), "rf_blood_sha256": digest(rf_b),
              "rf_t0_sha256": digest(ct["t0_per_event"]), "calibration_iq_sha256": digest(cal),
              "tissue_iq_sha256": digest(tissue), "blood_iq_sha256": digest(blood),
              "calibration_iq_shape": cal.shape, "roi_diagnostics": diagnostics,
              "validation_shares_stationary_tissue_scatterers": True,
              "note": "Nine-angle survey uses nine emissions per frame; three-angle timing and noise change slightly and require a separate confirmation.",
              "elapsed_s": time.monotonic()-t0}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--upstream-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--scatterer-fraction", type=float, default=.1)
    p.add_argument("--view-indices", type=int, nargs="+", default=(4, 0, 1, 2, 3, 5, 6, 7, 8))
    p.add_argument("--frames", type=int, default=16)
    a = p.parse_args()
    result = survey(a.upstream_root, a.output, a.seed, a.scatterer_fraction,
                    tuple(a.view_indices), a.frames)
    print(json.dumps({"seed": result["seed"], "elapsed_s": result["elapsed_s"]}))
