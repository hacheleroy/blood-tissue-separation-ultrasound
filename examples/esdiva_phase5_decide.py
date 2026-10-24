"""Apply the fixed phase-5 paired gates without selecting a candidate."""

import argparse
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decide(config_path, input_folder, output):
    config_path,input_folder,output=map(lambda p:Path(p).resolve(),
                                         (config_path,input_folder,output))
    if output.exists():raise FileExistsError(output)
    config=json.loads(config_path.read_text());rule=config["comparison"]
    seed_list=[*config["design_seeds"],config["new_development_seed"]]
    failures=[];rows=[];inputs=[]
    for seed in seed_list:
        for scene in config["scenes"]:
            filename=input_folder/f"phase5_candidate_{scene}_s{seed}.json"
            if not filename.exists():
                failures.append(f"missing: {scene} seed {seed}")
                continue
            item=json.loads(filename.read_text())
            inputs.append({"scene":scene,"seed":seed,"sha256":sha(filename)})
            scripts=Path(__file__).parent
            if (item["config_sha256"]!=sha(config_path) or
                    item["script_sha256"]!=sha(scripts/"esdiva_phase5_score.py") or
                    item["local_phase_context_sha256"]!=sha(scripts/"esdiva_phase5_local_phase.py") or
                    item["rebeam_sha256"]!=sha(scripts/"esdiva_phase5_rebeam.py") or
                    item["bounded_context_sha256"]!=sha(scripts/"esdiva_phase3_bounded_context.py") or
                    item["input_gate_sha256"]!=sha(scripts/"esdiva_phase3_input_gate.py") or
                    item["numerical_source_sha256"]!=config["numerical_source_sha256"] or
                    item["seed"]!=seed or
                    item["scene"]!=scene or item["status"]!="phase-5 paired source-resolved comparison" or
                    not item["input_gate_core"]["input_gate_pass"] or
                    not item["input_gate_wide"]["pass"] or
                    item["blood_to_tissue_wide_roi_power_evaluator_only"]>
                    config["evaluator_only_blood_to_tissue_roi_power_max"]):
                failures.append(f"input/source/score: {scene} seed {seed}")
                continue
            none=item["variants"]["none"]["methods"]
            phase=item["variants"]["local_phase"]["methods"]
            for variant,methods in (("none",none),("local_phase",phase)):
                if methods["roundtrip"]["paired_added_blood_iq_nrmse"]>rule["roundtrip_nrmse_max"]:
                    failures.append(f"roundtrip {variant}: {scene} seed {seed}")
            for name in ("svd_cut1","axis1_cut1"):
                blood_delta=(phase[name]["paired_added_blood_iq_nrmse"]-
                             none[name]["paired_added_blood_iq_nrmse"])
                leakage_ratio=(phase[name]["refitted_background_output_to_input_tissue_power"]/
                               none[name]["refitted_background_output_to_input_tissue_power"])
                rows.append({"seed":seed,"scene":scene,"method":name,
                             "blood_nrmse_increase":blood_delta,"tissue_leakage_ratio":leakage_ratio})
                if blood_delta>rule["blood_nrmse_increase_max_motion"]:
                    failures.append(f"blood {name}: {scene} seed {seed}")
                if (scene in ("rigid_fast","smooth_slow") and
                        leakage_ratio>rule["moving_tissue_leakage_ratio_max"]):
                    failures.append(f"moving tissue {name}: {scene} seed {seed}")
            if scene=="smooth_slow":
                slow=none["axis2_slow_r4_cap040"]
                if (slow["paired_added_blood_iq_nrmse"]>rule["slow_axis2_blood_nrmse_max"] or
                        slow["refitted_background_output_to_input_tissue_power"]>
                        rule["slow_axis2_tissue_leakage_max"]):
                    failures.append(f"slow Axis 2: {scene} seed {seed}")
    result={"status":"pass" if not failures else "fail",
            "decision":"Eligible for a pre-held-out source/config/results freeze" if not failures
            else "Reserved held-out seeds remain unopened",
            "config_sha256":sha(config_path),"decision_script_sha256":sha(__file__),
            "input_scores":inputs,"rows":rows,"failures":failures,
            "reserved_heldout_seeds":config["reserved_heldout_seeds"]}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--input-folder",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();r=decide(a.config,a.input_folder,a.output)
    print(json.dumps({"status":r["status"],"failures":r["failures"]}))
