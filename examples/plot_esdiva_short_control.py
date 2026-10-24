"""Plot the paired blood response and background output of one RF control."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt


LABELS = {
    "identity": "Identity", "svd_cut1": "SVD", "baranger_ssm": "SSM",
    "axis1_cut1": "Axis 1", "axis2": "Axis 2", "axis3_orthogonal": "Axis 3 orth",
    "axis3_loaded": "Axis 3 loaded", "lee_aer": "Lee AER",
    "axis3_loaded_gated": "Axis 3 gated",
}


def plot(replay, pilot, output):
    replay = json.loads(Path(replay).read_text())
    pilot = json.loads(Path(pilot).read_text())
    fig, ax = plt.subplots(figsize=(7.3, 5.2), constrained_layout=True)
    points = []
    for name, row in replay["methods"].items():
        if "blood_marginal_iq_nrmse" not in row:
            continue
        x = row["refitted_background_output_power_to_input_tissue"]
        y = row["blood_marginal_iq_nrmse"]
        gated = pilot["methods"].get(name, {}).get("adaptive_stap_support_gate")
        color = "tab:red" if gated is not None and not gated["passed"] else "tab:blue"
        duplicate = next((p for p in points if abs(p[0]-x) < 1e-10 and
                          abs(p[1]-y) < 1e-10 and p[2] == color), None)
        if duplicate is not None:
            duplicate[3].append(LABELS.get(name, name))
            continue
        points.append([x, y, color, [LABELS.get(name, name)]])
    for x, y, color, labels in points:
        ax.scatter(x, y, s=70, color=color, marker="x" if color == "tab:red" else "o")
        offset = (5, -14) if labels == ["Axis 2"] else (5, 4)
        ax.annotate(" / ".join(labels), (x, y), xytext=offset,
                    textcoords="offset points", fontsize=9)
    ax.set(xscale="log", xlabel="Refitted tissue + noise output / input tissue power",
           ylabel="Added-blood response IQ error (NRMSE)",
           title="Example 22: one full-density, stationary fast-flow development case")
    ax.grid(alpha=.25, which="both")
    ax.text(.02, .98, "Lower on both axes is preferred. Red: adaptive support gate failed.",
            transform=ax.transAxes, va="top", fontsize=8)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.replay, args.pilot, args.output)
