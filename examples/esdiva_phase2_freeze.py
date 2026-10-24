"""Verify the committed eSDIva phase-2 settings before opening held-out seeds."""

import hashlib
import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_freeze(root, manifest_path, seed, motion, flow_scale, report=None):
    root, manifest_path = Path(root).resolve(), Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    for relative, expected in manifest["files_sha256"].items():
        actual = sha256(root/relative)
        if actual != expected:
            raise ValueError(f"Frozen file changed: {relative}: {actual} != {expected}")
    actual = numerical_source_hash(root)
    if actual != manifest["numerical_source_sha256"]:
        raise ValueError("Numerical source changed after freeze")
    scenarios = json.loads((root/"configs/esdiva_phase2_scenarios.json").read_text())
    if seed not in scenarios["reserved_heldout_seeds"]:
        raise ValueError("Seed is not in the frozen held-out partition")
    if not any(s["motion"] == motion and s["flow_scale"] == flow_scale for s in scenarios["scenes"]):
        raise ValueError("Motion/flow scene is not frozen")
    if scenarios["upstream_commit"] != manifest["upstream_commit"]:
        raise ValueError("Upstream pin changed after freeze")
    if report is not None:
        if (report["seed"] != seed or report["motion"] != motion or
                report["flow_scale"] != flow_scale or
                report["upstream_commit"] != scenarios["upstream_commit"] or
                report["script_sha256"] != manifest["files_sha256"]["examples/esdiva_phase2_acquire.py"] or
                report["benchmark_source_sha256"] != actual or
                report["view_indices"] != scenarios["view_indices"] or
                report["scatterer_fraction"] != scenarios["particle_fraction"] or
                report["calibration_frames"] != scenarios["calibration_frames"] or
                report["evaluation_frames"] != scenarios["evaluation_frames"]):
            raise ValueError("Acquisition disagrees with frozen held-out protocol")
    return manifest
