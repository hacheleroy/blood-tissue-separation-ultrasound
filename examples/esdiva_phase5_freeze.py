"""Verify the published phase-5 development pass before any reserved run."""

import hashlib
import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_freeze(root, manifest_path, seed, scene_name, report=None):
    root,manifest_path=Path(root).resolve(),Path(manifest_path).resolve()
    manifest=json.loads(manifest_path.read_text())
    config=json.loads((root/"configs/esdiva_phase5_candidate.json").read_text())
    if manifest["status"]!="phase-5 pre-held-out source and results freeze":
        raise ValueError("The complete phase-5 development gate has not passed")
    for name,expected in manifest["files_sha256"].items():
        if sha(root/name)!=expected:raise ValueError(f"Frozen file changed: {name}")
    if (numerical_source_hash(root)!=manifest["numerical_source_sha256"] or
            seed not in config["reserved_heldout_seeds"] or
            scene_name not in config["scenes"] or
            config["upstream_commit"]!=manifest["upstream_commit"]):
        raise ValueError("Source, scene or reserved seed differs from the freeze")
    if report is not None:
        scene=config["scenes"][scene_name]
        if (report["seed"]!=seed or report["motion"]!=scene["motion"] or
                report["flow_scale"]!=scene["flow_scale"] or
                report["upstream_commit"]!=config["upstream_commit"] or
                report["script_sha256"]!=manifest["files_sha256"]["examples/esdiva_phase2_acquire.py"] or
                report["benchmark_source_sha256"]!=manifest["numerical_source_sha256"] or
                report["view_indices"]!=[4,3,5] or
                report["scatterer_fraction"]!=config["scatterer_fraction"] or
                report["calibration_frames"]!=config["calibration_frames"] or
                report["evaluation_frames"]!=config["evaluation_frames"]):
            raise ValueError("Reserved acquisition violates frozen protocol")
    return manifest
