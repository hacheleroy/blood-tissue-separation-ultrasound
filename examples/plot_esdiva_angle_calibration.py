"""Plot tissue angle coherence and scalar normalization gain from RF controls."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt


def plot(controls, output):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), constrained_layout=True)
    for path in controls:
        report = json.loads(Path(path).read_text())
        rows = sorted(report["diagnostics"], key=lambda row: row["angle_deg"])
        angle = [row["angle_deg"] for row in rows]
        label = f"{report['scatterer_fraction']:.0%} scatterers"
        axes[0].plot(angle, [row["coherence_train"] for row in rows], "o-", label=label)
        axes[1].plot(angle, [row["regression_input_noise_power_gain"] for row in rows],
                     "o-", label=label)
    axes[0].axhline(.7, color="0.35", linestyle="--", linewidth=1, label="gate")
    axes[1].axhline(4, color="0.35", linestyle="--", linewidth=1, label="gate")
    axes[0].set(ylabel="Tissue coherence to 0°", ylim=(0, 1.05))
    axes[1].set(ylabel="Normalization noise power gain", yscale="log", ylim=(.5, 100))
    for ax in axes:
        ax.set(xlabel="Transmit angle (degrees)")
        ax.grid(alpha=.25)
        ax.legend(frameon=False)
    fig.suptitle("Example 22 tissue-only RF: angular calibration development control")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.controls, args.output)
