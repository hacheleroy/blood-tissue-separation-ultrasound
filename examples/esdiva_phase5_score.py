"""Predeclared paired score for broad-support local phase registration."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_angle_calibration import angular_diagnostics
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase3_input_gate import choose_input
from esdiva_phase5_local_phase import LocalPhaseContext, tissue_mask


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(a,b):return float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))


def evaluate(acquisition, broad_path, config_path, output):
    acquisition,broad_path,config_path,output=map(lambda p:Path(p).resolve(),
                                                 (acquisition,broad_path,config_path,output))
    if output.exists():raise FileExistsError(output)
    config=json.loads(config_path.read_text());report=json.loads((acquisition/"report.json").read_text())
    audit=json.loads(broad_path.with_suffix(".json").read_text())
    if (report["upstream_commit"]!=config["upstream_commit"] or
            report["seed"] not in [*config["design_seeds"],config["new_development_seed"]] or
            report["scatterer_fraction"]!=config["scatterer_fraction"] or
            report["calibration_frames"]!=config["calibration_frames"] or
            report["evaluation_frames"]!=config["evaluation_frames"] or
            report["view_indices"]!=[4,3,5] or
            audit["source_report_sha256"]!=sha(acquisition/"report.json") or
            audit["script_sha256"]!=sha(Path(__file__).with_name("esdiva_phase5_rebeam.py")) or
            report["benchmark_source_sha256"]!=config["numerical_source_sha256"] or
            not audit["rf_sources_and_noise_verified"]):
        raise ValueError("Acquisition or RF replay violates the fixed source contract")
    scenes=[key for key,entry in config["scenes"].items() if
            entry["motion"]==report["motion"] and entry["flow_scale"]==report["flow_scale"]]
    if len(scenes)!=1:raise ValueError("Unknown scene")
    with np.load(acquisition/"iq_sources.npz") as old,np.load(broad_path) as wide:
        core={k:old[k] for k in ("views","calibration_views","tissue","blood","noise","z_mm","x_mm")}
        large={k:wide[k] for k in ("views","calibration_views","tissue","blood","noise","z_mm","x_mm","core_rows")}
    for key,array in (("mixture_sha256",core["views"]),("calibration_sha256",core["calibration_views"]),
                      ("tissue_sha256",core["tissue"]),("blood_sha256",core["blood"]),
                      ("noise_sha256",core["noise"])):
        if digest(array)!=report[key]:raise ValueError(f"IQ source digest mismatch: {key}")
    for key in ("views","calibration_views","tissue","blood","noise"):
        if (digest(large[key])!=audit["broad_iq_sha256"][key] or
                not np.array_equal(large[key][:,:,large["core_rows"]],core[key])):
            raise ValueError(f"Wider source IQ replay mismatch: {key}")
    if (not np.array_equal(large["z_mm"][large["core_rows"]],core["z_mm"]) or
            not np.array_equal(large["x_mm"],core["x_mm"]) or
            relative_error(core["tissue"]+core["blood"]+core["noise"],core["views"])>1e-8):
        raise ValueError("Core shape or source additivity mismatch")
    core_roi,core_gate=choose_input(core["calibration_views"],core["x_mm"],
                                    coherence_min=config["input_coherence_min"],
                                    max_noise_power_gain=config["input_noise_power_gain_max"],
                                    min_roi_pixels=config["minimum_core_roi_pixels"])
    wide_roi=tissue_mask(large["z_mm"],large["x_mm"],config["local_phase"]["vessel_guard_mm"])
    rows=angular_diagnostics(large["calibration_views"][:8],large["calibration_views"][8:],wide_roi,[4,3,5])
    wide_gate={"pixels":int(wide_roi.sum()),
               "minimum_coherence_train":min(r["coherence_train"] for r in rows[1:]),
               "minimum_coherence_validation":min(r["coherence_validation"] for r in rows[1:]),
               "maximum_bounded_noise_power_gain_train":max(r["bounded_input_noise_power_gain"] for r in rows[1:])}
    wide_gate["pass"]=(wide_gate["pixels"]>=config["minimum_wide_roi_pixels"] and
                       min(wide_gate["minimum_coherence_train"],wide_gate["minimum_coherence_validation"])
                       >=config["input_coherence_min"] and
                       wide_gate["maximum_bounded_noise_power_gain_train"]<=
                       config["input_noise_power_gain_max"])
    contamination=ratio(large["blood"][:,:,wide_roi],large["tissue"][:,:,wide_roi])
    result={"status":"input gate failed; no filter scores","seed":report["seed"],"scene":scenes[0],
            "acquisition_sha256":sha(acquisition/"report.json"),
            "broad_audit_sha256":sha(broad_path.with_suffix(".json")),
            "config_sha256":sha(config_path),"script_sha256":sha(__file__),
            "local_phase_context_sha256":sha(Path(__file__).with_name("esdiva_phase5_local_phase.py")),
            "rebeam_sha256":sha(Path(__file__).with_name("esdiva_phase5_rebeam.py")),
            "bounded_context_sha256":sha(Path(__file__).with_name("esdiva_phase3_bounded_context.py")),
            "input_gate_sha256":sha(Path(__file__).with_name("esdiva_phase3_input_gate.py")),
            "numerical_source_sha256":numerical_source_hash(Path(__file__).resolve().parents[1]),
            "input_gate_core":core_gate,"input_gate_wide":wide_gate,
            "blood_to_tissue_wide_roi_power_evaluator_only":contamination,
            "source_truth_evaluator_only":True,"score_frames":len(core["views"])-config["score_start"]}
    if (not core_gate["input_gate_pass"] or not wide_gate["pass"] or
            contamination>config["evaluator_only_blood_to_tissue_roi_power_max"]):
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
        return result
    geo=report["geometry"]
    observations={"mixture":(core["views"],large["views"]),
                  "background":(core["tissue"]+core["noise"],
                                large["tissue"]+large["noise"])}
    inputs={key:(PublicInput(c,core["calibration_views"],core_roi,geo),
                 PublicInput(w,large["calibration_views"],wide_roi,geo))
            for key,(c,w) in observations.items()}
    variants={}
    for variant in ("none","local_phase"):
        ctx={}
        for label,(ci,wi) in inputs.items():
            ctx[label]=(Phase3Context(ci,packet=config["packet"],score_start=config["score_start"],motion="none")
                        if variant=="none" else LocalPhaseContext(ci,wi,large["core_rows"],
                            large["z_mm"],large["x_mm"],packet=config["packet"],
                            score_start=config["score_start"],phase_options=config["local_phase"]))
        if ctx["mixture"].calibration["calibration_sha256"]!=ctx["background"].calibration["calibration_sha256"]:
            raise AssertionError("Paired calibration differs")
        methods={}
        names=config["methods"] if variant=="none" else {
            key:config["methods"][key] for key in ("roundtrip","svd_cut1","axis1_cut1")}
        for method,spec in names.items():
            mixed,bg=(fit_shared(ctx[label],spec) for label in ("mixture","background"))
            methods[method]={"paired_added_blood_iq_nrmse":relative_error(
                                 mixed.views-bg.views,core["blood"][config["score_start"]:]),
                             "refitted_background_output_to_input_tissue_power":ratio(
                                 bg.views,core["tissue"][config["score_start"]:]),
                             "frozen_blood_power_ratio":ratio(
                                 mixed.apply_frozen(core["blood"]),core["blood"][config["score_start"]:])}
        variants[variant]={"motion_fit":ctx["mixture"].motion_summary,
                           "background_motion_fit":ctx["background"].motion_summary,
                           "methods":methods}
    result["status"]="phase-5 paired source-resolved comparison"
    result["variants"]=variants
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",type=Path,required=True)
    p.add_argument("--broad",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();r=evaluate(a.acquisition,a.broad,a.config,a.output)
    print(json.dumps({k:r[k] for k in ("seed","scene","status")}))
