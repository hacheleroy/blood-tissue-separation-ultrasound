"""Seal source, settings and all five development records before held-out RF."""

import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash
from esdiva_phase2_freeze import sha256


SOURCE_FILES = (
    "examples/esdiva_angle_calibration.py",
    "examples/esdiva_example22_pilot.py",
    "examples/esdiva_gated_stap.py",
    "examples/esdiva_phase2_acquire.py",
    "examples/esdiva_phase2_score.py",
    "examples/esdiva_phase2_freeze.py",
    "examples/esdiva_phase2_run_heldout.py",
    "examples/esdiva_phase2_plot.py",
    "examples/esdiva_phase2_make_freeze.py",
    "configs/esdiva_phase2_development.json",
    "configs/esdiva_phase2_motion_relaxation_dev.json",
    "configs/esdiva_phase2_scenarios.json",
)


def make_freeze(root):
    root=Path(root).resolve()
    path=root/"configs/esdiva_phase2_freeze.json"
    if path.exists():
        raise FileExistsError("Refusing to replace the freeze manifest")
    scenarios=json.loads((root/"configs/esdiva_phase2_scenarios.json").read_text())
    if (len(scenarios["development_seeds"])!=1 or len(scenarios["reserved_heldout_seeds"])<2 or
            set(scenarios["development_seeds"]) & set(scenarios["reserved_heldout_seeds"])):
        raise ValueError("Need disjoint development and reserved seeds")
    files={relative:sha256(root/relative) for relative in SOURCE_FILES}
    config_sha=files["configs/esdiva_phase2_development.json"]
    for scene in scenarios["scenes"]:
        folder=f"runs/esdiva_phase2_development/{scene['id']}"
        r_path,s_path=root/folder/"report.json",root/folder/"score.json"
        report,score=json.loads(r_path.read_text()),json.loads(s_path.read_text())
        if (not report["calibration_gate_pass"] or
                report["seed"]!=scenarios["development_seeds"][0] or
                report["motion"]!=scene["motion"] or report["flow_scale"]!=scene["flow_scale"] or
                report["upstream_commit"]!=scenarios["upstream_commit"] or
                report["script_sha256"]!=files["examples/esdiva_phase2_acquire.py"] or
                report["benchmark_source_sha256"]!=numerical_source_hash(root) or
                score["seed"]!=report["seed"] or score["status"]!="development" or
                score["config_sha256"]!=config_sha or
                score["script_sha256"]!=files["examples/esdiva_phase2_score.py"] or
                score["acquisition_report_sha256"]!=sha256(r_path) or
                any(method["status"]!="ok" for mode in score["results"].values()
                    for method in mode["methods"].values())):
            raise ValueError(f"Development evidence invalid: {scene['id']}")
        files[f"{folder}/report.json"]=sha256(r_path)
        files[f"{folder}/score.json"]=sha256(s_path)
    for name in ("stationary_fast","rigid_fast"):
        relative=f"runs/esdiva_phase2_development/{name}/motion_relaxation_score.json"
        record=json.loads((root/relative).read_text())
        if (record["status"]!="development" or
                record["config_sha256"]!=files["configs/esdiva_phase2_motion_relaxation_dev.json"] or
                record["script_sha256"]!=files["examples/esdiva_phase2_score.py"] or
                record["acquisition_report_sha256"]!=files[f"runs/esdiva_phase2_development/{name}/report.json"]):
            raise ValueError(f"Invalid motion-relaxation diagnostic: {name}")
        files[relative]=sha256(root/relative)
    manifest={"status":"pre-held-out freeze; no reserved seed opened",
              "upstream_commit":scenarios["upstream_commit"],
              "numerical_source_sha256":numerical_source_hash(root),
              "files_sha256":files,
              "decision":"The loaded STAP and motion-correction arms remain diagnostics; all baselines and negative controls are retained. No tuning or promotion after reserved seeds are opened.",
              "limitation":"Nominal guarded samples are correlated; five reduced-density development scenes share one scatterer seed and score 16 frames each."}
    path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return manifest


if __name__=="__main__":
    root=Path(__file__).resolve().parents[1]
    result=make_freeze(root)
    print(json.dumps({"status":result["status"],"file_count":len(result["files_sha256"])},indent=2))
