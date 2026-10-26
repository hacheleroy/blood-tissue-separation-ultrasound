"""Verify the sealed Axis 1 motion-stability simulation and its result matrix."""
from hashlib import sha256
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/axis1_stable_tissue_v1"
HELDOUT_SHA256 = "27af33a904e6f2bd7757a5e70136670624ecd5d13266f7d592a917cffe30974f"
METHODS = {"raw_svd", "registered_return", "native_subtraction",
           "smooth_full_field", "smooth_bright_quartile"}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify():
    config=json.loads((ROOT/"configs/axis1_stable_tissue_v1.json").read_text())
    freeze=json.loads((RUN/"freeze.json").read_text())
    archive={name:json.loads((RUN/f"{name}.json").read_text())
             for name in ("seen_diagnostic","development","heldout")}
    if (not freeze["locked_before_heldout"] or
        config["primary_method"]!="smooth_full_field" or
        config["heldout_seeds"]!=freeze["heldout_seeds"] or
        digest(RUN/"heldout.json")!=HELDOUT_SHA256 or
        any(digest(RUN/f"{name}.json")!=freeze[f"{name}_sha256"]
            for name in ("seen_diagnostic","development"))):
        raise ValueError("Study freeze or result archive changed")
    sources=freeze["sources"]
    if not sources or any(item["source_sha256"]!=sources for item in archive.values()):
        raise ValueError("Source lease changed")
    for name,expected in sources.items():
        path=Path(name)
        if path.is_absolute() or ".." in path.parts or digest(ROOT/path)!=expected:
            raise ValueError(f"Numerical source changed: {name}")
    seeds=[set(config[f"{name}_seeds"]) for name in archive]
    if any(seeds[i]&seeds[j] for i in range(3) for j in range(i+1,3)):
        raise ValueError("Split seed overlap")
    for split,item in archive.items():
        expected=[(seed,condition) for seed in config[f"{split}_seeds"]
                  for condition in config["conditions"]]
        if (item["split"]!=split or
            [(c["seed"],c["condition"]) for c in item["cases"]]!=expected):
            raise ValueError("Incomplete split matrix")
        for case in item["cases"]:
            if set(case["results"])!=METHODS:
                raise ValueError("Incomplete method matrix")
            for result in case["results"].values():
                for key in ("added_blood_iq_error","frozen_blood_iq_error",
                    "frozen_tissue_power_retention","refitted_background_power_remaining",
                    "same_motion_added_blood_iq_error","motion_refit_background_change_iq_norm"):
                    if not math.isfinite(result[key]) or result[key]<0:
                        raise ValueError("Invalid result measure")
                if max(result["frozen_source_closure_nrmse"],
                       result["frozen_mixture_replay_nrmse"])>1e-10:
                    raise ValueError("Frozen source replay failed")
            if not math.isfinite(case["mechanism"]["stable_motion_refit_rms_px"]):
                raise ValueError("Invalid motion difference")
        if split=="seen_diagnostic":
            if item["decision"] is not None:
                raise ValueError("Seen bank assigned a decision")
            continue
        decision=item["decision"]
        if decision["primary_method"]!=config["primary_method"] or decision["joint_criterion_passed"]:
            raise ValueError("Negative joint decision changed")
        rows={ (row["seed"],row["condition"]):row for row in decision["case_gates"] }
        if set(rows)!=set(expected):raise ValueError("Incomplete decision rows")
        thresholds=config["decision"]
        for case in item["cases"]:
            is_moving=case["scene"]["motion"]!="none"
            r=case["results"]
            actual=rows[(case["seed"],case["condition"])]
            limit=(thresholds["moving_added_blood_error_margin_vs_raw"] if is_moving
                   else thresholds["stationary_added_blood_error_margin_vs_raw"])
            computed={
                "blood_noninferior_to_raw":r[config["primary_method"]]["added_blood_iq_error"] <= r["raw_svd"]["added_blood_iq_error"]+limit,
                "tissue_below_limit":not is_moving or r[config["primary_method"]]["refitted_background_power_remaining"] <= thresholds["moving_background_power_max"],
                "registration_failure_limit":case["mechanism"]["registration_failure_frames"]<=thresholds["maximum_registration_failures_per_case"]}
            if any(actual[key]!=value for key,value in computed.items()):
                raise ValueError("Decision gate inconsistent with stored measures")
    return {"status":"frozen stable-tissue Axis 1 study verified",
            "seen_cases":len(archive["seen_diagnostic"]["cases"]),
            "development_cases":len(archive["development"]["cases"]),
            "heldout_cases":len(archive["heldout"]["cases"]),
            "heldout_sha256":HELDOUT_SHA256}


if __name__=="__main__":print(json.dumps(verify(),indent=2))
