"""Summarize the phase-3 development input gate without scoring blood filters."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def summarize(config_path, reports, summary_path, figure_path):
    """Plot angle quality and check prospective development conditions."""
    config = json.loads(Path(config_path).read_text())
    records = []
    for path in reports:
        report = json.loads(Path(path).read_text())
        if (report["upstream_commit"] != config["upstream_commit"] or
                report["seed"] in config["spent_heldout_seeds_excluded"] or
                report["frames"] != config["calibration_frames"] or
                report["view_indices"][0] != 4):
            raise ValueError(f"Survey provenance mismatch: {path}")
        for name in ("outer_both", "outer_left", "wide_both"):
            roi = report["roi_diagnostics"][name]
            angles = {row["angle_deg"]: row for row in roi["angles"]}
            if not all(angle in angles for angle in (-5, 5)):
                raise ValueError("Expected both nearby views")
            rows = [angles[-5], angles[5]]
            coherence = min(min(r["coherence_train"], r["coherence_validation"])
                            for r in rows)
            gain = max(r["regression_input_noise_power_gain"] for r in rows)
            contamination = roi["blood_to_tissue_power"]
            records.append({"seed": report["seed"], "fraction": report["scatterer_fraction"],
                            "n_views": len(report["view_indices"]), "roi": name,
                            "roi_pixels": roi["pixels"], "min_coherence": coherence,
                            "max_noise_power_gain": gain,
                            "blood_to_tissue_power": contamination,
                            "gate_pass": bool(coherence >= config["min_coherence_per_angle_and_half"]
                                              and gain <= config["max_regression_input_noise_power_gain_per_angle"]
                                              and contamination <= config["max_blood_to_tissue_power_in_roi"]),
                            "report": Path(path).name})
    required = [(seed, .1, 3) for seed in config["prospective_development_seeds"]]
    required.append((2063, 1., 3))
    selected = []
    for seed, fraction, views in required:
        matches = [r for r in records if (r["seed"], r["fraction"], r["n_views"], r["roi"]) ==
                   (seed, fraction, views, "outer_left")]
        if len(matches) != 1:
            raise ValueError(f"Missing or duplicate required development condition {seed}/{fraction}")
        selected.extend(matches)
    decision = bool(all(r["gate_pass"] for r in selected))
    summary = {"status": "development input diagnostics only; no filter ranking",
               "selected_roi": config["training_roi"], "required_conditions_pass": decision,
               "required_records": selected, "all_records": records,
               "next_step": ("develop a full 80-frame scene and paired blood-response controls before any new held-out freeze"
                             if decision else "repair input reconstruction on new development data; do not open new held-out seeds")}
    summary_path = Path(summary_path)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")

    # The plot keeps discovery and prospective checks distinct. Bars show the
    # weaker of the two ±5-degree views, not an average that hides a bad view.
    groups = [r for r in records if r["n_views"] == 3 and r["roi"] in ("outer_both", "outer_left")]
    keys = sorted(set((r["seed"], r["fraction"]) for r in groups))
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True, layout="constrained")
    xpos = np.arange(len(keys)); width = .34
    for j, roi in enumerate(("outer_both", "outer_left")):
        subset = {(r["seed"],r["fraction"]): r for r in groups if r["roi"] == roi}
        x = xpos + (j-.5)*width
        label = "Original |x| ≥ 6 mm" if j == 0 else "Candidate x ≤ −6 mm"
        ax[0].bar(x, [subset[k]["min_coherence"] for k in keys], width, label=label)
        ax[1].bar(x, [subset[k]["max_noise_power_gain"] for k in keys], width, label=label)
    ax[0].axhline(config["min_coherence_per_angle_and_half"], color="black", ls="--", lw=1)
    ax[1].axhline(config["max_regression_input_noise_power_gain_per_angle"], color="black", ls="--", lw=1)
    ax[0].set_ylabel("Minimum ±5° coherence")
    ax[1].set_ylabel("Maximum scalar noise-power gain")
    ax[1].set_xlabel("Development scatterer seed and particle fraction")
    ax[1].set_xticks(xpos, [f"{s}\n{f:g}" for s,f in keys])
    ax[0].set_ylim(0, 1)
    ax[0].legend(loc="lower left")
    ax[0].set_title("Independent RF development controls; two disjoint eight-frame calibration halves")
    figure_path = Path(figure_path)
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--figure", type=Path, required=True)
    p.add_argument("reports", type=Path, nargs="+")
    a = p.parse_args()
    result = summarize(a.config, a.reports, a.summary, a.figure)
    print(json.dumps({"required_conditions_pass": result["required_conditions_pass"],
                      "next_step": result["next_step"]}))
