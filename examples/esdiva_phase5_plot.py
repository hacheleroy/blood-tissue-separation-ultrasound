"""Plot the locked paired motion gates from phase-5 score records."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(folder, output):
    folder,output=Path(folder),Path(output)
    if output.exists():raise FileExistsError(output)
    seeds=(2083,2097,2101,2111)
    scenes=("rigid_fast","smooth_slow")
    methods={"svd_cut1":("Global SVD", "#175d9d","o"),
             "axis1_cut1":("Axis 1: Motion-compensated virtual long ensemble","#c05b24","s")}
    fig,axes=plt.subplots(2,2,figsize=(11.0,7.0),sharex="col",layout="constrained")
    for col,scene in enumerate(scenes):
        for method,(name,color,marker) in methods.items():
            delta=[];leak=[]
            for seed in seeds:
                r=json.loads((folder/f"phase5_candidate_{scene}_s{seed}.json").read_text())
                before=r["variants"]["none"]["methods"][method]
                after=r["variants"]["local_phase"]["methods"][method]
                delta.append(after["paired_added_blood_iq_nrmse"]-
                             before["paired_added_blood_iq_nrmse"])
                leak.append(after["refitted_background_output_to_input_tissue_power"]/
                            before["refitted_background_output_to_input_tissue_power"])
            for row,data in enumerate((delta,leak)):
                ax=axes[row,col]
                ax.plot(range(4),data,color=color,marker=marker,linewidth=1.8,markersize=6,
                        label=name if row==0 and col==0 else None)
        axes[0,col].axhline(.02,color="#a52929",linestyle="--",linewidth=1,label="Gate" if col==0 else None)
        axes[1,col].axhline(.9,color="#a52929",linestyle="--",linewidth=1)
        axes[0,col].set_title(scene.replace("_"," ").title())
        axes[1,col].set_xticks(range(4),[str(s) for s in seeds])
        axes[1,col].set_xlabel("Independent tissue scatterer cloud")
        axes[0,col].grid(alpha=.2);axes[1,col].grid(alpha=.2)
    axes[0,0].set_ylabel("Paired added-blood IQ NRMSE change\n(local phase − no motion)")
    axes[1,0].set_ylabel("Refitted tissue leakage ratio\n(local phase / no motion)")
    axes[0,0].legend(loc="upper left",fontsize=8)
    fig.suptitle("Wider-support, disjoint-validated local phase registration",fontsize=13)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=170)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scores",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();plot(a.scores,a.output)
