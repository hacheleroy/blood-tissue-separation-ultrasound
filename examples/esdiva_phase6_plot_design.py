"""Plot source-resolved phase-6 development gates from archived compact scores."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


def render(archive, output):
    archive, output = Path(archive), Path(output)
    decision = json.loads((archive / "decision.json").read_text())
    if decision["status"] != "no design candidate passes" or len(decision["trials"]) != 12:
        raise ValueError("Expected the complete negative phase-6 design decision")
    colors = {0.5: "#0072B2", 0.65: "#009E73", 0.8: "#D55E00", 1.: "#CC79A7"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), sharex=True, sharey=True)
    for row, seed in enumerate((2309, 2311)):
        for col, scene in enumerate(("rigid_fast", "smooth_slow")):
            ax = axes[row, col]
            ax.add_patch(Rectangle((-.09, .57), .11, .33,
                                   facecolor="#e9f5ee", edgecolor="none", zorder=0))
            ax.axhline(.9, color="#315c42", linestyle="--", linewidth=1)
            ax.axvline(.02, color="#315c42", linestyle="--", linewidth=1)
            n = 0
            for record in decision["trials"]:
                name = record["trial"]
                if not name.startswith("wide_supported_core_diagnostic_local_phase_"):
                    continue
                config = json.loads((archive / "configs" / name).read_text())
                score = json.loads((archive / str(seed) / scene / "scores" / name).read_text())
                if score["status"] != "phase-6 paired source-resolved development comparison":
                    continue
                none = score["variants"]["none"]["methods"]
                moved = score["variants"]["local_phase"]["methods"]
                fraction = config["local_phase"]["phase_fraction"]
                narrow = config["local_phase"]["spatial_sigma_pixels"] == 2
                for method, marker in (("svd_cut1", "o"), ("axis1_cut1", "s")):
                    blood = (moved[method]["paired_added_blood_iq_nrmse"] -
                             none[method]["paired_added_blood_iq_nrmse"])
                    leakage = (moved[method]["refitted_background_output_to_input_tissue_power"] /
                               none[method]["refitted_background_output_to_input_tissue_power"])
                    ax.scatter(blood, leakage, s=56, marker=marker,
                               facecolor=colors[fraction] if narrow else "white",
                               edgecolor=colors[fraction], linewidth=1.5, zorder=3)
                    n += 1
            ax.set_title(f"Seed {seed}  ·  {scene.replace('_', ' ')}", loc="left", fontsize=11)
            ax.set_xlim(-.09, .19)
            ax.set_ylim(.57, 1.08)
            ax.grid(alpha=.2)
            if n == 0:
                ax.text(.5, .5, "Input gate blocked all phase trials", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9)
    for ax in axes[-1]: ax.set_xlabel("Added blood IQ NRMSE versus no motion")
    fig.supylabel("Refitted tissue leakage ratio versus no motion", x=.045)
    handles = [Line2D([], [], marker="o", linestyle="", color=c, label=f"phase {f:g}")
               for f, c in colors.items()]
    handles += [Line2D([], [], marker=m, linestyle="", color="#333", label=label)
                for m, label in (("o", "SVD cut 1"), ("s", "Axis 1 cut 1"))]
    handles += [Line2D([], [], marker="o", linestyle="", markerfacecolor="white",
                       markeredgecolor="#333", label="width 4; filled = width 2")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(.5, .01), fontsize=8)
    fig.suptitle("Phase-6 design clouds: paired blood response and moving-tissue leakage",
                 y=.97, fontsize=13)
    fig.subplots_adjust(left=.11, right=.98, top=.87, bottom=.18,
                        hspace=.22, wspace=.15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    render(args.archive, args.output)
