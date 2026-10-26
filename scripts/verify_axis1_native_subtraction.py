"""Verify the source lease and complete native-grid Axis 1 result archive."""
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/axis1_native_subtraction_v1"
HELDOUT_SHA256 = "0879ab23172f31db7587ccf99680976327186ebba6f42609930b8612a7bb21fd"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    config=json.loads((ROOT/"configs/axis1_native_subtraction_v1.json").read_text())
    freeze=json.loads((RUN/"freeze.json").read_text())
    archives={name:json.loads((RUN/f"{name}.json").read_text())
              for name in ("seen_diagnostic","development","heldout")}
    if (not freeze["locked_before_heldout"] or
        freeze["heldout_seeds"]!=config["heldout_seeds"] or
        digest(RUN/"heldout.json")!=HELDOUT_SHA256 or
        any(digest(RUN/f"{name}.json")!=freeze[f"{name}_sha256"]
            for name in ("seen_diagnostic","development"))):
        raise ValueError("Freeze or archive digest changed")
    sources=freeze["sources"]
    if not sources or any(x["source_sha256"]!=sources for x in archives.values()):
        raise ValueError("Source lease differs between splits")
    for name,expected in sources.items():
        path=Path(name)
        if path.is_absolute() or ".." in path.parts or digest(ROOT/path)!=expected:
            raise ValueError(f"Numerical source changed: {name}")
    seeds=[set(config[f"{split}_seeds"]) for split in archives]
    if any(seeds[i]&seeds[j] for i in range(3) for j in range(i+1,3)):
        raise ValueError("Split seeds overlap")
    for split,archive in archives.items():
        declared=(config["seen_conditions"] if split=="seen_diagnostic"
                  else list(config["conditions"]))
        expected=[(seed,condition) for seed in config[f"{split}_seeds"]
                  for condition in declared]
        if archive["split"]!=split or [(c["seed"],c["condition"]) for c in archive["cases"]]!=expected:
            raise ValueError("Incomplete split matrix")
        for case in archive["cases"]:
            if set(case["results"])!={"raw_svd","registered_return","native_subtraction"}:
                raise ValueError("Missing method")
            mechanism=case["mechanism"]
            if max(mechanism["output_identity_nrmse"],mechanism["blood_identity_nrmse"])>1e-10:
                raise ValueError("Native-grid identity failed")
            for result in case["results"].values():
                for key in ("added_blood_iq_error","frozen_blood_iq_error",
                    "frozen_tissue_power_retention","refitted_background_power_remaining"):
                    value=result[key]
                    if not math.isfinite(value) or value<0:
                        raise ValueError("Invalid primary measure")
                if max(result["frozen_source_closure_nrmse"],
                       result["frozen_mixture_replay_nrmse"])>1e-10:
                    raise ValueError("Frozen source replay failed")
        if split=="seen_diagnostic":
            if archive["decision"] is not None:raise ValueError("Seen diagnostic given a decision")
        elif archive["decision"]["joint_criterion_passed"]:
            raise ValueError("Archived negative joint decision changed")
    return {"status":"frozen native-grid Axis 1 study verified",
        "seen_cases":len(archives["seen_diagnostic"]["cases"]),
        "development_cases":len(archives["development"]["cases"]),
        "heldout_cases":len(archives["heldout"]["cases"]),
        "heldout_sha256":HELDOUT_SHA256}


if __name__=="__main__":print(json.dumps(verify(),indent=2))
