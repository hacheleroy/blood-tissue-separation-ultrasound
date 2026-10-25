"""Reserved-bank phase-10 scorer with the unchanged numerical comparison."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_phase10_support import (
    build_phase10_context, choose_phase10_mode, phase10_support,
)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(a,b):return float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))


def evaluate(acquisition, broad_path, config_path, output):
    acquisition,broad_path,config_path,output=map(lambda p:Path(p).resolve(),
                                                 (acquisition,broad_path,config_path,output))
    if output.exists():raise FileExistsError(output)
    config=json.loads(config_path.read_text());report=json.loads((acquisition/"report.json").read_text())
    audit=json.loads(broad_path.with_suffix(".json").read_text())
    if (config["development_manifest_sha256"]!=sha(Path(__file__).resolve().parents[1]/
                                           "runs/esdiva_phase10_design/manifest.json") or
            config["development_decision_sha256"]!=sha(Path(__file__).resolve().parents[1]/
                                           "runs/esdiva_phase10_design/decision.json") or
            config["phase10_reserve_split_sha256"]!=sha(Path(__file__).resolve().parents[1]/
                                           "configs/esdiva_phase10_reserve_split.json") or
            config["prospective_decision_sha256"]!=sha(Path(__file__).resolve().parents[1]/
                                           "runs/esdiva_phase10_prospective/decision.json") or
            config["prospective_manifest_sha256"]!=sha(Path(__file__).resolve().parents[1]/
                                           "runs/esdiva_phase10_prospective/manifest.json") or
            config["support_source_sha256"]!=sha(Path(__file__).with_name(
                                           "esdiva_phase10_support.py")) or
            config["acquisition_source_sha256"]!=sha(Path(__file__).with_name(
                                           "esdiva_phase10_acquire.py")) or
            config["selector_source_sha256"]!=sha(Path(__file__).with_name(
                                           "esdiva_phase8_selector.py")) or
            config["transport_source_sha256"]!=sha(Path(__file__).with_name(
                                           "esdiva_phase7_reference_transport.py")) or
            config["magnitude_shift_fit_sha256"]!=sha(Path(__file__).with_name(
                                           "esdiva_phase6_unit_shift.py")) or
            config["score_script_sha256"]!=sha(__file__) or
            config["admitted_reserved_seeds"]!=[2401,2411] or
            report["upstream_commit"]!=config["upstream_commit"] or
            report["seed"] not in config["admitted_reserved_seeds"] or
            report["script_sha256"]!=config["acquisition_source_sha256"] or
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
        core={k:old[k] for k in ("views","calibration_views","tissue","blood","noise","z_mm","x_mm","training_roi")}
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
    wide_roi,core_roi,calibration_gate=phase10_support(
        large["calibration_views"],large["z_mm"],large["x_mm"],large["core_rows"])
    if (not np.array_equal(core_roi,core["training_roi"]) or
            digest(core["training_roi"])!=report["training_roi_sha256"] or
            not np.array_equal(large["calibration_views"][:,:,large["core_rows"]],
                               core["calibration_views"])):
        raise ValueError("Changed source-defined core tissue support")
    contamination=ratio(large["blood"][:,:,wide_roi],large["tissue"][:,:,wide_roi])
    core_contamination=ratio(core["blood"][:,:,core_roi],core["tissue"][:,:,core_roi])
    if abs(core_contamination-report["training_blood_to_tissue_power_ratio"])>1e-12:
        raise ValueError("Acquisition source-eligibility statistic changed")
    result={"status":"input gate failed; no filter scores","seed":report["seed"],"scene":scenes[0],
            "acquisition_sha256":sha(acquisition/"report.json"),
            "broad_audit_sha256":sha(broad_path.with_suffix(".json")),
            "config_sha256":sha(config_path),"script_sha256":sha(__file__),
            "local_phase_context_sha256":sha(Path(__file__).with_name("esdiva_phase5_local_phase.py")),
            "selector_context_sha256":sha(Path(__file__).with_name("esdiva_phase8_selector.py")),
            "support_context_sha256":sha(Path(__file__).with_name("esdiva_phase10_support.py")),
            "acquisition_script_sha256":sha(Path(__file__).with_name("esdiva_phase10_acquire.py")),
            "transport_context_sha256":sha(Path(__file__).with_name("esdiva_phase7_reference_transport.py")),
            "magnitude_shift_fit_sha256":sha(Path(__file__).with_name("esdiva_phase6_unit_shift.py")),
            "rebeam_sha256":sha(Path(__file__).with_name("esdiva_phase5_rebeam.py")),
            "bounded_context_sha256":sha(Path(__file__).with_name("esdiva_phase3_bounded_context.py")),
            "numerical_source_sha256":numerical_source_hash(Path(__file__).resolve().parents[1]),
            "calibration_support":calibration_gate,
            "source_eligibility_pass":report["source_eligibility_pass"],
            "blood_to_tissue_core_roi_power_evaluator_only":core_contamination,
            "blood_to_tissue_wide_roi_power_evaluator_only":contamination,
            "source_truth_evaluator_only":True,"score_frames":len(core["views"])-config["score_start"]}
    if (not calibration_gate["pass"] or not report["source_eligibility_pass"] or
            core_contamination>config["evaluator_only_blood_to_tissue_roi_power_max"] or
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
    selected,selection=choose_phase10_mode(inputs["mixture"][0],inputs["mixture"][1],
                          large["core_rows"],large["z_mm"],large["x_mm"])
    result["selection"]=selection
    result["selected_mode"]=selected
    for variant in ("none","selected"):
        ctx={}
        for label,(ci,wi) in inputs.items():
            ctx[label]=build_phase10_context("none" if variant=="none" else selected,
                       ci,wi,large["core_rows"],large["z_mm"],large["x_mm"])
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
    result["status"]="phase-10 paired source-resolved reserved diagnostic"
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
