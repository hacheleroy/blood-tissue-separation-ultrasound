"""Diagnose an ineligible phase-8 cloud from existing RF; never score it.

This script was written after the phase-8 freeze and cannot alter its input
gate, selection rule, or all-scene decision. RF checkpoints are evaluator-only.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, view_rf
from esdiva_phase5_rebeam import load_rf
from esdiva_phase8_selector import choose_wide_guard


def rf_samples(folder):
    length = 0
    for path in sorted(folder.glob("rf_event_*.npz")):
        with np.load(path) as packet:
            length = max(length, packet["rf"].shape[1])
    if len(list(folder.glob("rf_event_*.npz"))) != 288:
        raise ValueError("Incomplete RF source")
    return length


def diagnose(acquisitions, upstream, output, seed):
    acquisitions, upstream, output = map(lambda p: Path(p).resolve(),
                                          (acquisitions, upstream, output))
    if output.exists():
        raise FileExistsError(output)
    head = subprocess.check_output(
        ["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_COMMIT:
        raise ValueError("Changed Example 22 source")
    sys.path.insert(0, str(upstream/"examples/example22_flow_doppler"))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq
    from esdiva.utilities import add_noise, rf_rms
    _, tx, _ = example.build_simulator()
    _, nine = example.compound_events(tx, 1)
    events = [nine[i] for i in (4, 3, 5)]
    rows = []
    for scene in ("stationary_fast", "rigid_fast", "smooth_slow", "smooth_zero"):
        folder = acquisitions/str(seed)/scene
        length = max(rf_samples(folder/"rf_tissue"),
                     rf_samples(folder/"rf_blood"))
        tissue, coords = load_rf(folder/"rf_tissue", 288, 64, length)
        blood, bcoords = load_rf(folder/"rf_blood", 288, 64, length)
        if (not np.array_equal(coords["t0_per_event"], bcoords["t0_per_event"]) or
                coords["dt"] != bcoords["dt"]):
            raise ValueError("RF source coordinates differ")
        clean = tissue+blood
        noise = add_noise(clean, 6., reference=rf_rms(clean),
                          rng=np.random.default_rng(seed+101))-clean
        cal, z, x = view_rf(
            tissue[:48]+noise[:48],
            {"dt":coords["dt"], "t0_per_event":coords["t0_per_event"][:48]},
            events, tx, GRID_MM, example.C)
        eval_coords = {"dt":coords["dt"],
                       "t0_per_event":coords["t0_per_event"][48:]}
        b, _, _ = view_rf(blood[48:], eval_coords, events, tx, GRID_MM, example.C)
        t, _, _ = view_rf(tissue[48:], eval_coords, events, tx, GRID_MM, example.C)
        fs = example.C/(2*GRID_MM["dz"]*1e-3)
        carrier = echo_center_frequency(
            np.moveaxis(cal, -2, -1), 1/fs,
            band=(.3*example.FC, 1.7*example.FC))
        cal, b, t = (rf2iq(a, fs=fs, fc=carrier, axis=-2)
                     for a in (cal, b, t))
        start = int(np.argmin(abs(z-15.6)))
        core = slice(start, start+8)
        mask = np.broadcast_to(abs(x[None]) >= 6., (8, len(x)))
        power = float(np.mean(abs(b[:, :, core, :][:, :, mask])**2)/
                      np.mean(abs(t[:, :, core, :][:, :, mask])**2))
        wide_start = int(np.argmin(abs(z-14.4)))
        wide = slice(wide_start, wide_start+32)
        _, gate = choose_wide_guard(cal[:, :, wide, :], z[wide], x)
        report_path = folder/"report.json"
        if report_path.exists():
            report = json.loads(report_path.read_text())
            if (abs(power-report["training_blood_to_tissue_power_ratio"]) > 1e-12 or
                    abs(carrier-report["carrier_hz"]) > 1e-6):
                raise AssertionError("Post-freeze replay differs from acquisition")
        rows.append({
            "scene":scene, "report_exists":report_path.exists(),
            "core_training_blood_to_tissue_power":power,
            "frozen_acquisition_limit":.1,
            "postfreeze_wide_calibration":gate,
            "carrier_hz":float(carrier),
        })
    result = {"status":"post-freeze RF diagnosis only; no scoring or retuning",
              "seed":seed,"upstream_commit":head,
              "acquisition_script_sha256":hashlib.sha256(
                  Path(__file__).with_name("esdiva_phase2_acquire.py").read_bytes()).hexdigest(),
              "diagnostic_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "rows":rows}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisitions", type=Path, required=True)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=2729)
    args = parser.parse_args()
    report = diagnose(args.acquisitions, args.upstream, args.output, args.seed)
    print(json.dumps([{key: row[key] for key in
                       ("scene", "report_exists", "core_training_blood_to_tissue_power")}
                      for row in report["rows"]]))
