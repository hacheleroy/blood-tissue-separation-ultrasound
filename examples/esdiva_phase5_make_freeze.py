"""Seal passing phase-5 development results; refuse a failed/missing scene."""

import json
from pathlib import Path

from esdiva_example22_pilot import numerical_source_hash
from esdiva_phase5_freeze import sha


def make_freeze(root):
    root=Path(root).resolve();path=root/"configs/esdiva_phase5_freeze.json"
    if path.exists():raise FileExistsError(path)
    pre=json.loads((root/"configs/esdiva_phase5_design_freeze.json").read_text())
    config=json.loads((root/"configs/esdiva_phase5_candidate.json").read_text())
    if set((*config["design_seeds"],config["new_development_seed"]))&set(config["reserved_heldout_seeds"]):
        raise ValueError("Development and held-out banks overlap")
    for name,expected in pre["files_sha256"].items():
        if sha(root/name)!=expected:raise ValueError(f"Pre-new-cloud file changed: {name}")
    if numerical_source_hash(root)!=pre["numerical_source_sha256"]:
        raise ValueError("Numerical source changed")
    decision_path=root/"runs/esdiva_phase5_development/decision.json"
    decision=json.loads(decision_path.read_text())
    if (decision["status"]!="pass" or decision["failures"] or
            decision["config_sha256"]!=sha(root/"configs/esdiva_phase5_candidate.json") or
            decision["decision_script_sha256"]!=sha(root/"examples/esdiva_phase5_decide.py") or
            len(decision["input_scores"])!=16):
        raise ValueError("Full development decision did not pass")
    files=dict(pre["files_sha256"])
    fresh=config["new_development_seed"]
    for scene,spec in config["scenes"].items():
        prefix=f"runs/esdiva_phase5_development/{fresh}/{scene}"
        for leaf in ("report.json","wide_audit.json","score.json"):
            name=f"{prefix}/{leaf}"; files[name]=sha(root/name)
        report=json.loads((root/f"{prefix}/report.json").read_text())
        score=json.loads((root/f"{prefix}/score.json").read_text())
        audit=json.loads((root/f"{prefix}/wide_audit.json").read_text())
        if (report["seed"]!=fresh or report["motion"]!=spec["motion"] or
                report["flow_scale"]!=spec["flow_scale"] or
                report["script_sha256"]!=files["examples/esdiva_phase2_acquire.py"] or
                report["benchmark_source_sha256"]!=pre["numerical_source_sha256"] or
                score["acquisition_sha256"]!=files[f"{prefix}/report.json"] or
                score["broad_audit_sha256"]!=files[f"{prefix}/wide_audit.json"] or
                audit["source_report_sha256"]!=files[f"{prefix}/report.json"]):
            raise ValueError(f"New-cloud source mismatch: {scene}")
        recorded=next((i for i in decision["input_scores"] if i["seed"]==fresh and i["scene"]==scene),None)
        if recorded is None or recorded["sha256"]!=files[f"{prefix}/score.json"]:
            raise ValueError(f"New-cloud score changed: {scene}")
    files["runs/esdiva_phase5_development/decision.json"]=sha(decision_path)
    for name in ("examples/esdiva_phase5_plot.py","examples/esdiva_phase5_freeze.py",
                 "examples/esdiva_phase5_make_freeze.py","README.md",
                 "docs/RESULTS_ESDIVA_PHASE5_DEVELOPMENT.md",
                 "docs/figures/esdiva_phase5_motion_paired.png"):
        files[name]=sha(root/name)
    manifest={"status":"phase-5 pre-held-out source and results freeze",
              "upstream_commit":config["upstream_commit"],
              "numerical_source_sha256":pre["numerical_source_sha256"],
              "files_sha256":files,
              "reserved_heldout_seeds":config["reserved_heldout_seeds"],
              "decision":"The full new-cloud paired rule passed; reserved seeds remain unopened. No reselection or retuning is allowed after this freeze.",
              "limitation":"Four independent 10%-density scatterer clouds, four correlated scenes per cloud, 16 scored frames per scene; RF software simulation is not phantom validation."}
    path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return manifest


if __name__=="__main__":
    r=make_freeze(Path(__file__).resolve().parents[1])
    print(json.dumps({"status":r["status"],"file_count":len(r["files_sha256"])}))
