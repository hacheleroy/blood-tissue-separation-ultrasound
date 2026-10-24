"""Seal the declared comparison and eight passing development records."""

import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash
from esdiva_phase4_decision import decide
from esdiva_phase4_score import sha


SOURCE_FILES = (
    "examples/esdiva_angle_calibration.py",
    "examples/esdiva_example22_pilot.py",
    "examples/esdiva_phase2_acquire.py",
    "examples/esdiva_phase3_bounded_context.py",
    "examples/esdiva_phase3_input_gate.py",
    "examples/esdiva_phase4_differential_context.py",
    "examples/esdiva_phase4_differential_context_v1.py",
    "examples/esdiva_phase4_motion_context.py",
    "examples/esdiva_phase4_motion_study.py",
    "examples/esdiva_phase4_phase_context.py",
    "examples/esdiva_phase4_phase_study.py",
    "examples/esdiva_phase4_score.py",
    "examples/esdiva_phase4_score_v1.py",
    "examples/esdiva_phase4_decision.py",
    "examples/esdiva_phase4_decision_v1.py",
    "examples/esdiva_phase4_plot.py",
    "examples/esdiva_phase4_freeze.py",
    "examples/esdiva_phase4_run_heldout.py",
    "examples/esdiva_phase4_make_freeze.py",
    "configs/esdiva_phase4_candidate.json",
    "configs/esdiva_phase4_candidate_v1.json",
    "tests/test_esdiva_phase4_differential.py",
    "tests/test_esdiva_phase4_phase.py",
    "docs/RESULTS_ESDIVA_PHASE4_DEVELOPMENT.md",
    "docs/figures/esdiva_phase4_motion_paired.png",
    "docs/figures/esdiva_phase4_slow_flow.png",
)


def make_freeze(root):
    root = Path(root).resolve()
    path = root/"configs/esdiva_phase4_freeze.json"
    if path.exists():
        raise FileExistsError("Refusing to replace the phase-4 freeze")
    config = json.loads((root/"configs/esdiva_phase4_candidate.json").read_text())
    if (set((*config["design_seeds"],config["fresh_development_seed"])) &
        set(config["reserved_heldout_seeds"])):
        raise ValueError("Development and reserved seeds overlap")
    files = {relative:sha(root/relative) for relative in SOURCE_FILES}
    paths = []
    for seed in (*config["design_seeds"],config["fresh_development_seed"]):
        for scene,spec in config["scenes"].items():
            folder = f"runs/esdiva_phase4_development/{seed}/{scene}"
            report_path,score_path = root/folder/"report.json",root/folder/"score.json"
            report,score = json.loads(report_path.read_text()),json.loads(score_path.read_text())
            if (report["seed"] != seed or report["motion"] != spec["motion"] or
                report["flow_scale"] != spec["flow_scale"] or
                report["upstream_commit"] != config["upstream_commit"] or
                report["script_sha256"] != files["examples/esdiva_phase2_acquire.py"] or
                report["benchmark_source_sha256"] != numerical_source_hash(root) or
                report["view_indices"] != [4,3,5] or
                score["acquisition_sha256"] != sha(report_path) or
                score["numerical_source_sha256"] != numerical_source_hash(root)):
                raise ValueError(f"Development source or acquisition mismatch: {folder}")
            files[f"{folder}/report.json"] = sha(report_path)
            files[f"{folder}/score.json"] = sha(score_path)
            paths.append(score_path)
    decision = decide(root/"configs/esdiva_phase4_candidate.json", paths)
    if not decision["all_pass"]:
        raise ValueError("Development criteria failed; no held-out freeze")
    summary_path = root/"runs/esdiva_phase4_development/decision.json"
    if json.loads(summary_path.read_text()) != decision:
        raise ValueError("Development decision differs from committed summary")
    files["runs/esdiva_phase4_development/decision.json"] = sha(summary_path)
    for name in ("motion_rigid", "phase_stationary", "phase_rigid",
                 "phase_smooth_slow", "phase_smooth_zero"):
        relative = f"runs/esdiva_phase4_development/2083_exploratory/{name}.json"
        record = json.loads((root/relative).read_text())
        if record["seed"] != config["design_seeds"][0]:
            raise ValueError(f"Exploratory diagnostic is not a design seed: {relative}")
        files[relative] = sha(root/relative)
    v1 = root/"runs/esdiva_phase4_development/v1_decision.json"
    if (json.loads(v1.read_text())["all_pass"] or
        json.loads(v1.read_text())["config_sha256"] != files["configs/esdiva_phase4_candidate_v1.json"] or
        json.loads(v1.read_text())["decision_script_sha256"] != files["examples/esdiva_phase4_decision_v1.py"]):
        raise ValueError("Prior failed decision must be preserved")
    files["runs/esdiva_phase4_development/v1_decision.json"] = sha(v1)
    for seed in config["design_seeds"]:
        for scene in config["scenes"]:
            folder = f"runs/esdiva_phase4_development/v1_{seed}/{scene}"
            report_path, score_path = root/folder/"report.json", root/folder/"score.json"
            report, score = json.loads(report_path.read_text()), json.loads(score_path.read_text())
            if (report["seed"] != seed or
                score["acquisition_sha256"] != sha(report_path) or
                score["config_sha256"] != files["configs/esdiva_phase4_candidate_v1.json"] or
                score["script_sha256"] != files["examples/esdiva_phase4_score_v1.py"] or
                score["phase_context_sha256"] != files["examples/esdiva_phase4_differential_context_v1.py"]):
                raise ValueError(f"Failed first candidate record mismatch: {folder}")
            files[f"{folder}/report.json"] = sha(report_path)
            files[f"{folder}/score.json"] = sha(score_path)
    manifest = {"status":"phase-4 pre-held-out comparison freeze",
                "upstream_commit":config["upstream_commit"],
                "numerical_source_sha256":numerical_source_hash(root),
                "files_sha256":files,
                "decision":"No reselection or parameter change after opening reserved seeds; failed input gates produce no filter scores. Preserve negative controls and the method tradeoff.",
                "limitation":"Only two independent scatterer clouds at 10% density with four correlated scenes and 16 scored frames each; this is not phantom or in vivo validation."}
    path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return manifest


if __name__ == "__main__":
    result = make_freeze(Path(__file__).resolve().parents[1])
    print(json.dumps({"status":result["status"],"file_count":len(result["files_sha256"])}))
