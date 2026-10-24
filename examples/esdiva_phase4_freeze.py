"""Verify the phase-4 pre-held-out source, rules, scenes and acquisition."""

import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash
from esdiva_phase4_score import sha


def verify_freeze(root, manifest_path, seed, scene_name, report=None):
    root, manifest_path = Path(root).resolve(), Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    config = json.loads((root/"configs/esdiva_phase4_candidate.json").read_text())
    if manifest["status"] != "phase-4 pre-held-out comparison freeze":
        raise ValueError("Invalid phase-4 freeze status")
    for relative, expected in manifest["files_sha256"].items():
        if sha(root/relative) != expected:
            raise ValueError(f"Frozen file changed: {relative}")
    if numerical_source_hash(root) != manifest["numerical_source_sha256"]:
        raise ValueError("Numerical source changed after freeze")
    if (seed not in config["reserved_heldout_seeds"] or
        scene_name not in config["scenes"] or
        config["upstream_commit"] != manifest["upstream_commit"]):
        raise ValueError("Seed or scene is outside frozen reserve")
    if report is not None:
        scene = config["scenes"][scene_name]
        if (report["seed"] != seed or report["motion"] != scene["motion"] or
            report["flow_scale"] != scene["flow_scale"] or
            report["upstream_commit"] != config["upstream_commit"] or
            report["script_sha256"] != manifest["files_sha256"]["examples/esdiva_phase2_acquire.py"] or
            report["benchmark_source_sha256"] != manifest["numerical_source_sha256"] or
            report["view_indices"] != [4,3,5] or
            report["scatterer_fraction"] != config["scatterer_fraction"] or
            report["calibration_frames"] != config["calibration_frames"] or
            report["evaluation_frames"] != config["evaluation_frames"]):
            raise ValueError("Acquisition differs from frozen contract")
    return manifest
