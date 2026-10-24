"""Audit the fixed phase-3 data-only ROI candidate on RF survey reports."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from esdiva_phase3_input_gate import select_diagnostics


def summarize(config_path, reports, output, figure):
    config = json.loads(Path(config_path).read_text())
    rows = []
    for path in reports:
        source = json.loads(Path(path).read_text())
        if (source["upstream_commit"] != config["upstream_commit"] or
                source["seed"] in config["spent_heldout_seeds_excluded"] or
                source["view_indices"] != config["view_indices"] or
                source["frames"] != config["train_frames"]+config["validation_frames"]):
            raise ValueError(f"Survey source mismatch: {path}")
        decision = select_diagnostics(source["roi_diagnostics"],
                                      min_roi_pixels=config["minimum_roi_pixels"])
        chosen = decision["selected_roi"]
        blood_ratio = (source["roi_diagnostics"][chosen]["blood_to_tissue_power"]
                       if chosen else None)
        row = {"seed": source["seed"], "scatterer_fraction": source["scatterer_fraction"],
               "n_tissue": source["n_tissue"], "n_blood": source["n_blood"],
               "source_report": Path(path).name,
               "source_report_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
               "input_gate": decision, "blood_to_tissue_power_evaluator_only": blood_ratio,
               "passed": bool(decision["input_gate_pass"] and blood_ratio is not None and
                              blood_ratio <= config["max_blood_to_tissue_power_selected_roi_evaluator_only"])}
        rows.append(row)
    required = [(seed,.1) for seed in config["prospective_development_seeds"]]+[(2063,1.)]
    for key in required:
        if sum((r["seed"],r["scatterer_fraction"]) == key for r in rows) != 1:
            raise ValueError(f"Missing or duplicate prospective condition {key}")
    prospective = [r for r in rows if (r["seed"],r["scatterer_fraction"]) in required]
    summary = {"status":"input-only development; no filter or held-out scores",
               "config_sha256":hashlib.sha256(Path(config_path).read_bytes()).hexdigest(),
               "gate_script_sha256":hashlib.sha256(Path(__file__).with_name(
                   "esdiva_phase3_input_gate.py").read_bytes()).hexdigest(),
               "required_conditions_pass":bool(all(r["passed"] for r in prospective)),
               "required_conditions":required,"records":rows}
    output = Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    ordered = sorted(rows,key=lambda r:(r["seed"],r["scatterer_fraction"]))
    labels=[f"{r['seed']}\n{r['scatterer_fraction']:g}" for r in ordered]
    val=[r["input_gate"]["candidates"][r["input_gate"]["selected_roi"]]
         if r["input_gate"]["selected_roi"] else None for r in ordered]
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,layout="constrained")
    x=np.arange(len(ordered))
    palette={"outer_left":"#207f8e", "wide_both":"#d38326",
             "outer_both":"#7657a3", "outer_right":"#ab5a62", None:"#777777"}
    colors=[palette[r["input_gate"]["selected_roi"]] for r in ordered]
    axes[0].bar(x,[v["min_coherence_validation"] if v else 0 for v in val],color=colors)
    axes[1].bar(x,[v["max_bounded_noise_power_gain_train"] if v else 0 for v in val],color=colors)
    axes[0].axhline(.7,color="black",ls="--",lw=1)
    axes[1].axhline(4.,color="black",ls="--",lw=1)
    axes[0].set_ylabel("Selected ROI: validation coherence")
    axes[1].set_ylabel("Selected ROI: bounded noise gain")
    axes[1].set_xticks(x,labels)
    axes[1].set_xlabel("Development scatterer seed and particle fraction")
    axes[0].set_title("Three-view input candidate; ROI selected only from first eight frames")
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(color=palette["outer_left"],label="left outer"),
                             Patch(color=palette["wide_both"],label="symmetric wide"),
                             Patch(color=palette["outer_both"],label="original symmetric")],
                   loc="lower left")
    figure=Path(figure);figure.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure,dpi=180);plt.close(fig)
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--figure",type=Path,required=True)
    p.add_argument("reports",type=Path,nargs="+")
    a=p.parse_args()
    result=summarize(a.config,a.reports,a.output,a.figure)
    print(json.dumps({"required_conditions_pass":result["required_conditions_pass"]}))
