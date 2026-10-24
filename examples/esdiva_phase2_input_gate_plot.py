"""Visualize frozen held-out angular eligibility without tuning its threshold."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCENES=("stationary_fast","rigid_fast","smooth_fast","smooth_slow","smooth_zero")
SEEDS=(2031,2053)


def plot(root,output):
    root,output=Path(root),Path(output)
    fig,axes=plt.subplots(1,2,figsize=(10,3.9),sharey=True,layout="constrained")
    for ax,j in zip(axes,(1,2)):
        for seed,color,marker,offset in ((2031,"#bd495e","o",-.08),(2053,"#177a99","s",.08)):
            values=[]
            for scene in SCENES:
                report=json.loads((root/str(seed)/scene/"report.json").read_text())
                row=report["calibration_angle_diagnostics"][j]
                values.append(min(row["coherence_train"],row["coherence_validation"]))
            ax.plot(np.arange(len(SCENES))+offset,values,marker=marker,color=color,
                    linewidth=1.2,markersize=6,label=f"seed {seed}")
        ax.axhline(.7,ls="--",lw=1.2,color="black",label="frozen 0.70 gate")
        ax.set_title(f"{[-5,+5][j-1]:+d}° view")
        ax.set_xticks(range(len(SCENES)),["Static\nfast","Rigid\nfast","Smooth\nfast",
                                          "Smooth\nslow","Smooth\nzero"])
        ax.grid(axis="y",alpha=.2)
        ax.set_ylim(0,1.02)
    axes[0].set_ylabel("Minimum calibration-half tissue coherence")
    axes[1].legend(loc="upper right",fontsize=8)
    fig.suptitle("eSDIva held-out angular input gate (before filter scoring)")
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=190,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    plot(args.reports_root,args.output)
