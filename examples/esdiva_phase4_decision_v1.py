"""Apply the declared phase-4 development decision without changing settings."""

import argparse
import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash
from esdiva_phase4_score import sha


def decide(config_path, records):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text())
    root = Path(__file__).resolve().parents[1]
    expected = {"config_sha256": sha(config_path),
                "script_sha256": sha(root/"examples/esdiva_phase4_score.py"),
                "phase_context_sha256": sha(root/"examples/esdiva_phase4_differential_context.py"),
                "bounded_context_sha256": sha(root/"examples/esdiva_phase3_bounded_context.py"),
                "input_gate_sha256": sha(root/"examples/esdiva_phase3_input_gate.py"),
                "numerical_source_sha256": numerical_source_hash(root)}
    values = {}
    for path in records:
        record = json.loads(Path(path).read_text())
        key = (record["seed"], record["scene"])
        if key in values:
            raise ValueError(f"Duplicate seed/scene: {key}")
        if (record["seed"] not in (config["design_seed"], config["fresh_development_seed"]) or
            record["scene"] not in config["scenes"] or
            any(record[k] != value for k, value in expected.items())):
            raise ValueError(f"Unexpected seed, scene, or source/config hashes: {path}")
        values[key] = (record, Path(path))
    required = {(seed, scene) for seed in (config["design_seed"],
                                          config["fresh_development_seed"])
                for scene in config["scenes"]}
    if set(values) != required:
        raise ValueError(f"Incomplete phase-4 development: missing {sorted(required-set(values))}")
    limits = config["comparison"]
    findings = {}
    all_pass = True
    for (seed, scene), (record, path) in sorted(values.items()):
        checks = {"input_gate": record["status"] == "phase-4 paired source-resolved comparison"
                  and record["input_gate"]["input_gate_pass"]}
        if checks["input_gate"]:
            variants = record["variants"]
            for name in ("none", "differential_phase"):
                method = variants[name]["methods"]["roundtrip"]
                checks[f"roundtrip_{name}"] = method["paired_added_blood_iq_nrmse"] <= limits["roundtrip_nrmse_max"]
            for method in ("svd_cut1", "axis1_cut1"):
                base = variants["none"]["methods"][method]
                phase = variants["differential_phase"]["methods"][method]
                checks[f"blood_{method}"] = (phase["paired_added_blood_iq_nrmse"]-
                    base["paired_added_blood_iq_nrmse"] <= limits["blood_nrmse_increase_max_motion"])
                if scene in ("rigid_fast", "smooth_slow"):
                    checks[f"leakage_{method}"] = (phase["refitted_background_output_to_input_tissue_power"] <=
                        limits["moving_tissue_leakage_ratio_max"]*
                        base["refitted_background_output_to_input_tissue_power"])
            if scene == "smooth_slow":
                slow = variants["none"]["methods"]["axis2_slow_r6_cap045"]
                checks["slow_blood"] = slow["paired_added_blood_iq_nrmse"] <= limits["slow_axis2_blood_nrmse_max"]
                checks["slow_leakage"] = (slow["refitted_background_output_to_input_tissue_power"] <=
                                           limits["slow_axis2_tissue_leakage_max"])
        findings[f"{seed}/{scene}"] = {"checks":checks, "passed":all(checks.values()),
                                       "score_sha256":sha(path),
                                       "acquisition_sha256":record["acquisition_sha256"]}
        all_pass &= findings[f"{seed}/{scene}"]["passed"]
    return {"status":"development criteria pass; eligible for a pre-held-out freeze"
            if all_pass else "development criteria fail; reserved seeds remain unopened",
            "all_pass":bool(all_pass),"config_sha256":sha(config_path),
            "decision_script_sha256":sha(__file__),"checks":findings,
            "reserved_seeds_unopened_at_decision":[2081,2093]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("records",nargs="+",type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = decide(args.config,args.records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(result["status"])
