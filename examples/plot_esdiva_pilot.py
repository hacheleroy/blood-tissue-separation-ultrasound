"""Plot evaluator-only source and phase diagnostics from the RF pilot.

Usage: python examples/plot_esdiva_pilot.py --pilot /path/to/output \
    --output docs/figures/esdiva_pilot_v08.png
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads((args.pilot / "report.json").read_text())
    with np.load(args.pilot / "pilot_iq.npz") as d:
        x = d["x_mm"]
        t = np.mean(abs(d["tissue"])**2, axis=(0, 1, 2))
        b = np.mean(abs(d["blood"])**2, axis=(0, 1, 2))
        train = d["training_roi"].any(axis=0)
        flow = d["flow_roi"].any(axis=0)

    ref = max(np.max(t), np.max(b))
    angles = np.asarray(report["selected_angles_deg"])
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.8), constrained_layout=True)
    ax = axes[0]
    ax.plot(x, 10*np.log10(t/ref + 1e-15), "o-", label="Tissue only", color="#276e85")
    ax.plot(x, 10*np.log10(b/ref + 1e-15), "o-", label="Blood only", color="#be6135")
    ax.scatter(x[train], np.full(train.sum(), -47), marker="s", color="#276e85",
               label="Tissue training annotation")
    ax.scatter(x[flow], np.full(flow.sum(), -44), marker="s", color="#be6135",
               label="Flow annotation")
    ax.set(xlabel="Lateral position x (mm)", ylabel="Mean isolated-source power (dB)",
           ylim=(-51, 5), title="Source separation and guarded training pixels")
    ax.grid(alpha=.2)
    ax.legend(ncol=2, fontsize=8, loc="upper center")
    ax = axes[1]
    ax.plot(angles, report["tissue_lag1_phase_rad_by_view"], "o-",
            label="Stationary tissue", color="#276e85")
    ax.plot(angles, report["blood_lag1_phase_rad_by_view"], "o-",
            label="Receding blood", color="#be6135")
    ax.axhline(0, color="gray", lw=.8)
    ax.set(xlabel="Transmit angle (degrees)", ylabel="Same-view lag-one IQ phase (rad)",
           xticks=angles, title="Phase survives independent per-view beamforming")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    fig.suptitle("eSDIva Example 22 three-view wiring pilot: 10% scatterer density, one seed",
                 fontsize=11)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=170)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
