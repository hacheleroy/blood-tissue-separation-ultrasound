"""Check the sealed Axis 1 source and prospective simulation archives."""
from hashlib import sha256
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/axis1_motion_order_v1"
HELDOUT_SHA256 = "3d03b3ec35b0afd9b62620f93dc9b1a4dbb3a0a3642793625ff1a0da1666a6d1"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify():
    freeze = json.loads((RUN / "freeze.json").read_text())
    config = json.loads((ROOT / "configs/axis1_motion_order_v1.json").read_text())
    dev = json.loads((RUN / "development.json").read_text())
    held = json.loads((RUN / "heldout.json").read_text())
    if (digest(RUN / "development.json") != freeze["development_sha256"]
            or digest(RUN / "heldout.json") != HELDOUT_SHA256
            or not freeze["locked_before_heldout"]
            or freeze["heldout_seeds"] != config["heldout_seeds"]
            or set(config["development_seeds"]) & set(config["heldout_seeds"])):
        raise ValueError("Split, freeze, or archive digest changed")
    if not freeze["sources"] or freeze["sources"] != dev["source_sha256"] != held["source_sha256"]:
        raise ValueError("Source lease changed")
    for name, expected in freeze["sources"].items():
        path = Path(name)
        if path.is_absolute() or ".." in path.parts or digest(ROOT / path) != expected:
            raise ValueError(f"Frozen source changed: {name}")
    names = None
    for split, archive in (("development", dev), ("heldout", held)):
        if archive["split"] != split or len(archive["cases"]) != 8:
            raise ValueError("Incomplete split")
        pairs = [(c["seed"], c["condition"]) for c in archive["cases"]]
        if pairs != [(seed, condition) for seed in config[f"{split}_seeds"]
                    for condition in config["conditions"]]:
            raise ValueError("Changed seed/condition matrix")
        for case in archive["cases"]:
            methods = [r["method"] for r in case["results"]]
            if len(methods) != 16 or len(set(methods)) != 16:
                raise ValueError("Incomplete or duplicated method matrix")
            if names is None:names = methods
            elif methods != names:raise ValueError("Methods changed across scenes")
            for row in case["results"]:
                metrics = (row["added_blood_iq_error"], row["frozen_blood_iq_error"],
                    row["frozen_tissue_power_retention"],
                    row["refitted_background_power_remaining"])
                if not all(math.isfinite(q) and q >= 0 for q in metrics):
                    raise ValueError("Nonfinite or negative metrics")
                if row["frozen_closure_absolute"] > 1e-8:
                    raise ValueError("Frozen operator replay mismatch")
    return {"status":"sealed Axis 1 motion-order archive verified",
            "cases":len(held["cases"]),"methods":len(names),
            "heldout_sha256":HELDOUT_SHA256}


if __name__ == "__main__":print(json.dumps(verify(), indent=2))
