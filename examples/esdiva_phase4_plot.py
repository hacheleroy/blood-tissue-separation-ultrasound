"""Plot exact paired development scores for motion and slow-flow candidates."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCENES = ("stationary_fast", "rigid_fast", "smooth_slow", "smooth_zero")
LABELS = {"stationary_fast":"Stationary tissue · fast flow",
          "rigid_fast":"Rigid tissue · fast flow",
          "smooth_slow":"Smooth tissue · slow flow",
          "smooth_zero":"Smooth tissue · zero flow"}
METRICS = ("paired_added_blood_iq_nrmse",
           "refitted_background_output_to_input_tissue_power")


def plot(records, motion_figure, flow_figure):
    rows = [json.loads(Path(p).read_text()) for p in records]
    samples = {(r["seed"], r["scene"]):r for r in rows}
    if len(rows) != len(samples) or len(samples) != 12 or any(
        r["status"] != "phase-4 paired source-resolved comparison" for r in rows):
        raise ValueError("Need twelve distinct eligible development scenes")
    seeds = sorted({seed for seed,_ in samples})
    if {s for _,s in samples} != set(SCENES):
        raise ValueError("Incomplete scene set")
    colors = {seeds[0]:"#246a9b", seeds[1]:"#bc7b24", seeds[2]:"#a44b65"}
    fig,axes = plt.subplots(2,2,figsize=(11,6.5),sharex=True,layout="constrained")
    for row, method in enumerate(("svd_cut1","axis1_cut1")):
        for col in range(2):
            ax = axes[row,col]
            for offset,seed in zip((-.14,0.,.14),seeds):
                base = [samples[seed,s]["variants"]["none"]["methods"][method]
                        for s in SCENES]
                corrected = [samples[seed,s]["variants"]["differential_phase"]["methods"][method]
                             for s in SCENES]
                y = ([b[METRICS[0]]-a[METRICS[0]] for a,b in zip(base,corrected)]
                     if col==0 else [b[METRICS[1]]/a[METRICS[1]] for a,b in zip(base,corrected)])
                ax.scatter(np.arange(4)+offset,y,c=colors[seed],s=42,label=str(seed))
            ax.axhline(.02 if col==0 else .9,color="#555555",lw=1,ls="--")
            ax.set_xticks(range(4),["Stationary\nfast","Rigid\nfast",
                                     "Smooth\nslow","Smooth\nzero"],fontsize=8)
            ax.grid(axis="y",alpha=.2)
            if row==0:
                ax.set_title("Added-blood NRMSE increase" if col==0 else
                             "Tissue leakage ratio: corrected / no motion",fontsize=10)
            if col==0:
                ax.set_ylabel(("SVD discard 1" if row==0 else
                              "Axis 1: Motion-compensated virtual long ensemble")+"\nLower is better")
            if col==1:
                ax.axvspan(.5,1.5,color="#aac9d0",alpha=.12)
                ax.axvspan(1.5,2.5,color="#aac9d0",alpha=.12)
    axes[0,0].legend(fontsize=8,ncol=3,loc="upper left",title="Scatterer seed")
    fig.suptitle("Differential phase: predeclared blood (+0.020) and moving-leakage (0.90) gates")
    fig.savefig(motion_figure,dpi=180)
    plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(13,4.8),layout="constrained")
    methods=("svd_cut1","baranger_ssm","axis2_default","axis2_slow_r4_cap040","lee_aer")
    names=("SVD","Baranger SSM","Axis 2 default","Axis 2 slow cap","Stephen Lee AER")
    palette=("#266a9f","#557e78","#7b788b","#a74c48","#765a9d")
    markers=("o","s","^","D","X")
    for ax,seed in zip(axes,seeds):
        record=samples[seed,"smooth_slow"]["variants"]["none"]["methods"]
        x=[record[m][METRICS[0]] for m in methods]
        y=[record[m][METRICS[1]] for m in methods]
        for label,xx,yy,color,marker in zip(names,x,y,palette,markers):
            ax.scatter(xx,yy,s=85,marker=marker,color=color,label=label)
        ax.axvline(.7,color="gray",lw=.8,linestyle="--")
        ax.axhline(.18,color="gray",lw=.8,linestyle="--")
        ax.set_xlim(-.03,max(1.18,max(x)+.08))
        ax.set_ylim(-.02,max(.42,max(y)+.03))
        ax.set_title(f"Development seed {seed}")
        ax.set_xlabel("Added-blood IQ NRMSE ↓")
        ax.grid(alpha=.18)
    axes[0].set_ylabel("Refitted background / input tissue power ↓")
    fig.suptitle("Smooth tissue, 2% blood velocity: preservation and clutter tradeoff")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",ncol=5,bbox_to_anchor=(.5,-.05),fontsize=8)
    fig.savefig(flow_figure,dpi=180,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--motion-figure",required=True,type=Path)
    p.add_argument("--flow-figure",required=True,type=Path)
    p.add_argument("records",nargs="+",type=Path)
    a=p.parse_args()
    plot(a.records,a.motion_figure,a.flow_figure)
