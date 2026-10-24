"""Verify and summarize an unmodified eSDIva Example 22 reproduction.

Run this after the upstream acquisition, Doppler and figure steps:
  python examples/audit_esdiva_example22.py --upstream-root /path/to/eSDIva \
      --output /path/to/reproduction.json

The output is a small manifest, not a copy of the channel RF or IQ arrays.
"""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path

import numpy as np


UPSTREAM_COMMIT = "ab7eaef6d23428762ab5b9079ab03a5aa99af085"
EXPECTED_EVENTS = {"A": 576, "B": 72, "C": 288}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit(root):
    root = Path(root).resolve()
    head = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_COMMIT:
        raise ValueError(f"Expected {UPSTREAM_COMMIT}, got {head}")
    changed = subprocess.check_output(
        ["git", "-C", str(root), "diff", "--name-only"], text=True
    ).splitlines()
    figure_names = ("ex22_bercoff_arms.png", "ex22_bercoff_tradeoff.png")
    allowed = {f"docs/examples/assets/{name}" for name in figure_names}
    if set(changed) - allowed:
        raise ValueError(f"Tracked upstream source changed: {set(changed)-allowed}")
    base = root / "examples/example22_flow_doppler"
    out = base / "out"
    arms = {}
    for arm, expected in EXPECTED_EVENTS.items():
        folder = out / f"RF_{arm}"
        manifest = folder / "contents.json"
        record = json.loads(manifest.read_text())
        events = record["events"]
        if set(events) != {str(i) for i in range(expected)}:
            raise AssertionError(f"{arm} has {len(events)}/{expected} RF events")
        hashes = []
        for i in range(expected):
            event = events[str(i)]
            path = folder / event["file"]
            digest = sha256(path)
            if digest != event["sha256"]:
                raise AssertionError(f"{arm} event {i} checksum mismatch")
            with np.load(path) as rf_file:
                if (list(rf_file["rf"].shape) != event["shape"] or
                        float(rf_file["t0"]) != event["t0"] or
                        float(rf_file["dt"]) != event["dt"]):
                    raise AssertionError(f"{arm} event {i} RF coordinates mismatch")
            hashes.append(digest)
        dts = {events[str(i)]["dt"] for i in range(expected)}
        if len(dts) != 1:
            raise AssertionError(f"{arm} RF sample periods differ across events")
        origins = [events[str(i)]["t0"] for i in range(expected)]
        arms[arm] = {
            "events": expected, "rf_sample_period_s": dts.pop(),
            "t0_min_s": min(origins), "t0_max_s": max(origins),
            "rf_shape_range": [min(events[str(i)]["shape"][1] for i in range(expected)),
                               max(events[str(i)]["shape"][1] for i in range(expected))],
            "manifest_sha256": sha256(manifest),
            "event_sha256_root": hashlib.sha256("\n".join(hashes).encode()).hexdigest(),
        }
    metrics_file = out / "metrics.json"
    iq_file = out / "doppler_results.npz"
    figures = {name: sha256(root / "docs/examples/assets" / name)
               for name in figure_names}
    metrics = json.loads(metrics_file.read_text())
    if set(metrics["arms"]) != set(EXPECTED_EVENTS):
        raise AssertionError("Missing Doppler metrics")
    for arm, expected in EXPECTED_EVENTS.items():
        if metrics["arms"][arm]["n_emissions"] != expected:
            raise AssertionError(f"{arm} Doppler emission count mismatch")
    return {
        "status": "upstream Example 22 steps 2-4 completed; all RF event checksums verified",
        "upstream_commit": head, "upstream_seed": 2026,
        "source_file_sha256": {p.name: sha256(p) for p in
                               sorted(base.glob("step[1234]_*.py"))},
        "environment": {"python": platform.python_version(),
                        "packages": {name: importlib.metadata.version(name) for name in
                                     ("esdiva", "numpy", "scipy", "numba", "h5py", "pyvista")}},
        "arms": arms, "metrics": metrics,
        "metrics_sha256": sha256(metrics_file), "iq_npz_sha256": sha256(iq_file),
        "figure_sha256": figures,
        "boundary": "Upstream nine-angle compound Doppler validates its own flow example, not our three-view separation filters",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.upstream_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": report["status"], "arms": report["arms"]}, indent=2))


if __name__ == "__main__":
    main()
