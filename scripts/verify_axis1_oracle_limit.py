"""Verify source lease, nested targets and full oracle-limit case matrix."""
from hashlib import sha256
import json
import math
from pathlib import Path
from statistics import median


ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/"runs/axis1_oracle_limit_v1"
HELDOUT_SHA256="3d3c282c5aafa04cdbccd33bd80fca617e818886faa811f957bf401c4950f794"


def digest(path):return sha256(path.read_bytes()).hexdigest()


def verify():
    config=json.loads((ROOT/"configs/axis1_oracle_limit_v1.json").read_text())
    freeze=json.loads((RUN/"freeze.json").read_text())
    archives={name:json.loads((RUN/f"{name}.json").read_text()) for name in
              ("seen_diagnostic","development","heldout")}
    if (not freeze["locked_before_heldout"] or
        freeze["heldout_seeds"]!=config["heldout_seeds"] or
        digest(RUN/"heldout.json")!=HELDOUT_SHA256 or
        any(digest(RUN/f"{name}.json")!=freeze[f"{name}_sha256"] for name in
            ("seen_diagnostic","development"))):
        raise ValueError("Archive or freeze changed")
    sources=freeze["sources"]
    if not sources or any(a["source_sha256"]!=sources for a in archives.values()):
        raise ValueError("Source lease differs across splits")
    for name,expected in sources.items():
        p=Path(name)
        if p.is_absolute() or ".." in p.parts or digest(ROOT/p)!=expected:
            raise ValueError(f"Numerical source changed: {name}")
    sets=[set(config[f"{name}_seeds"]) for name in archives]
    if any(sets[i]&sets[j] for i in range(3) for j in range(i+1,3)):
        raise ValueError("Seed split overlap")
    for split,a in archives.items():
        expected=[(seed,condition,flow) for seed in config[f"{split}_seeds"]
                  for condition in config["conditions"] for flow in config["flow_mps"]]
        if a["split"]!=split or [(c["seed"],c["condition"],c["flow_mps"])
                                  for c in a["cases"]]!=expected:
            raise ValueError("Incomplete case matrix")
        target_groups={}
        for c in a["cases"]:
            rows=c["lengths"]
            if [r["ensemble_frames"] for r in rows]!=config["ensemble_frames"]:
                raise ValueError("Missing ensemble length")
            target=[r["target_sha256"] for r in rows]
            if not all(t==target[0] for t in target) or set(target[0])!={"tissue","blood","noise"}:
                raise ValueError("Different scored frames across histories")
            target_groups.setdefault((c["seed"],c["condition"]),[]).append(target[0])
            for r in rows:
                if set(r["controls"])!={"raw_svd","known_motion_mixture_svd"} or set(r["oracles"])!={"1","2"}:
                    raise ValueError("Missing comparator")
                if abs(r["duration_ms"]-1000*r["ensemble_frames"]/config["prf_hz"])>1e-10:
                    raise ValueError("Duration inconsistency")
                for m in r["controls"].values():
                    if any(not math.isfinite(v) or v<0 for v in m.values()):
                        raise ValueError("Invalid control score")
                for k,v in r["oracles"].items():
                    if v["frozen_source_closure_nrmse"]>1e-10:
                        raise ValueError("Frozen-source closure failed")
                    for field in ("blood_iq_error_at_full_subtraction",
                        "background_power_at_full_subtraction",
                        "aligned_blood_projection_norm_fraction",
                        "tissue_training_energy_fraction"):
                        if not math.isfinite(v[field]) or v[field]<0:
                            raise ValueError("Invalid oracle score")
                    gain=v["minimum_gain_for_tissue_limit"]
                    blood=v["minimum_blood_error_at_tissue_limit"]
                    if (gain is None)!=(blood is None):
                        raise ValueError("Inconsistent gain feasibility")
                    if gain is not None:
                        if not 0<=gain<=1 or abs(blood-gain*v["blood_iq_error_at_full_subtraction"])>1e-10:
                            raise ValueError("Invalid minimum-gain blood error")
                    if v["joint_ambitious_goal_feasible"]!=(blood is not None and blood<=config["ambitious_blood_iq_error_limit"]):
                        raise ValueError("Ambitious-gate mismatch")
                    if v["joint_vs_raw_feasible"]!=(blood is not None and blood<=r["controls"]["raw_svd"]["added_blood_iq_error"]+config["raw_svd_blood_error_margin"]):
                        raise ValueError("Raw-gate mismatch")
        for targets in target_groups.values():
            if (len(targets)!=len(config["flow_mps"]) or
                any(v["tissue"]!=targets[0]["tissue"] or v["noise"]!=targets[0]["noise"]
                    for v in targets)):
                raise ValueError("Speed comparison changed tissue or noise target")
        for length in config["ensemble_frames"]:
            rows=[next(r for r in c["lengths"] if r["ensemble_frames"]==length)["oracles"]["1"]
                  for c in a["cases"] if c["motion"]!="none"]
            summary=a["summary"][str(length)]
            check={"moving_cases":len(rows),
                "median_oracle_blood_error_full":median(r["blood_iq_error_at_full_subtraction"] for r in rows),
                "median_oracle_background_power_full":median(r["background_power_at_full_subtraction"] for r in rows),
                "rank1_tissue_limit_feasible_cases":sum(r["minimum_gain_for_tissue_limit"] is not None for r in rows),
                "rank1_joint_ambitious_cases":sum(r["joint_ambitious_goal_feasible"] for r in rows),
                "rank1_joint_vs_raw_cases":sum(r["joint_vs_raw_feasible"] for r in rows)}
            if any(abs(summary[k]-v)>1e-10 for k,v in check.items()):
                raise ValueError("Summary differs from cases")
    return {"status":"Axis 1 oracle-limit archives verified",
            "seen_cases":len(archives["seen_diagnostic"]["cases"]),
            "development_cases":len(archives["development"]["cases"]),
            "heldout_cases":len(archives["heldout"]["cases"]),
            "heldout_sha256":HELDOUT_SHA256}


if __name__=="__main__":print(json.dumps(verify(),indent=2))
