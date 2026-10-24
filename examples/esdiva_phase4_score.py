"""Paired source score for fixed differential-phase and slow-flow candidates.

Only the evaluator reads isolated IQ or known motion. The filter receives
mixture or separately fitted tissue-plus-identical-noise background, tissue
calibration, a geometrically selected tissue ROI and acquisition geometry.
This runner writes input failure records without invoking a filter.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase3_input_gate import choose_input
from esdiva_phase4_differential_context import DifferentialPhaseContext


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ratio(a, b):
    return float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2), 1e-30))


def evaluate(acquisition, config_path, output):
    acquisition, config_path, output = map(lambda p: Path(p).resolve(),
                                           (acquisition, config_path, output))
    if output.exists():
        raise FileExistsError(output)
    config = json.loads(config_path.read_text())
    report = json.loads((acquisition/"report.json").read_text())
    with np.load(acquisition/"iq_sources.npz") as f:
        mixture, cal, tissue, blood, noise, x = (f[k] for k in
            ("views", "calibration_views", "tissue", "blood", "noise", "x_mm"))
    for key, value in (("mixture_sha256", mixture), ("calibration_sha256", cal),
                       ("tissue_sha256", tissue), ("blood_sha256", blood),
                       ("noise_sha256", noise)):
        if digest(value) != report[key]:
            raise ValueError(f"Acquisition digest mismatch: {key}")
    if (relative_error(tissue+blood+noise, mixture) > 1e-8 or
        report["upstream_commit"] != config["upstream_commit"] or
        report["view_indices"] != [4, 3, 5] or
        report["scatterer_fraction"] != config["scatterer_fraction"] or
        report["calibration_frames"] != config["calibration_frames"] or
        report["evaluation_frames"] != config["evaluation_frames"] or
        report["seed"] not in [*config["design_seeds"], config["fresh_development_seed"],
                             *config["reserved_heldout_seeds"]]):
        raise ValueError("Acquisition is outside declared design or source contract")
    scenes = [name for name, entry in config["scenes"].items()
              if entry["motion"] == report["motion"] and
                 entry["flow_scale"] == report["flow_scale"]]
    if len(scenes) != 1:
        raise ValueError("Acquisition motion/flow does not match a frozen scene")
    roi, gate = choose_input(cal, x,
                             coherence_min=config["input_coherence_min"],
                             max_noise_power_gain=config["input_noise_power_gain_max"],
                             min_roi_pixels=config["minimum_roi_pixels"])
    contamination = ratio(blood[:,:,roi], tissue[:,:,roi]) if roi is not None else None
    result = {"status": "input gate failed; no filter scores",
              "seed": report["seed"], "scene": scenes[0], "motion": report["motion"],
              "flow_scale": report["flow_scale"],
              "scatterer_fraction": report["scatterer_fraction"],
              "acquisition_sha256": sha(acquisition/"report.json"),
              "config_sha256": sha(config_path), "script_sha256": sha(__file__),
              "phase_context_sha256": sha(Path(__file__).with_name("esdiva_phase4_differential_context.py")),
              "bounded_context_sha256": sha(Path(__file__).with_name("esdiva_phase3_bounded_context.py")),
              "input_gate_sha256": sha(Path(__file__).with_name("esdiva_phase3_input_gate.py")),
              "numerical_source_sha256": numerical_source_hash(Path(__file__).resolve().parents[1]),
              "selected_roi_sha256": digest(roi) if roi is not None else None,
              "input_gate": gate,
              "blood_to_tissue_roi_power_evaluator_only": contamination,
              "source_truth_evaluator_only": True,
              "score_frames": len(mixture)-config["score_start"]}
    if gate["input_gate_pass"] and contamination <= config["evaluator_only_blood_to_tissue_roi_power_max"]:
        public = PublicInput(mixture, cal, roi, report["geometry"])
        background = PublicInput(tissue+noise, cal, roi, report["geometry"])
        records = {}
        for name, cls in (("none", Phase3Context),
                          ("differential_phase", DifferentialPhaseContext)):
            kwargs = {"packet": config["packet"], "score_start": config["score_start"]}
            if name == "none":
                kwargs["motion"] = "none"
            else:
                kwargs["phase_options"] = config["phase_gate"]
            mixed, bg = (cls(p, **kwargs) for p in (public, background))
            if mixed.calibration["calibration_sha256"] != bg.calibration["calibration_sha256"]:
                raise AssertionError("Paired calibration differs")
            methods = {}
            names = config["methods"] if name == "none" else {
                k:config["methods"][k] for k in ("roundtrip", "svd_cut1", "axis1_cut1")}
            for method, spec in names.items():
                fit, control = fit_shared(mixed, spec), fit_shared(bg, spec)
                methods[method] = {
                    "paired_added_blood_iq_nrmse": relative_error(
                        fit.views-control.views, blood[config["score_start"]:]),
                    "refitted_background_output_to_input_tissue_power": ratio(
                        control.views, tissue[config["score_start"]:]),
                    "frozen_blood_power_ratio": ratio(
                        fit.apply_frozen(blood), blood[config["score_start"]:])}
            records[name] = {"motion_fit": mixed.motion_summary,
                             "background_motion_fit": bg.motion_summary,
                             "methods": methods}
        result["status"] = "phase-4 paired source-resolved comparison"
        result["variants"] = records
    elif gate["input_gate_pass"]:
        result["status"] = "evaluator-only blood contamination gate failed; no filter scores"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    r = evaluate(a.acquisition, a.config, a.output)
    print(json.dumps({k:r[k] for k in ("seed", "motion", "flow_scale", "status")}))
