"""Plot the paired blood/background tradeoff across frozen eSDIva scenes."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = {
    "svd_cut1": ("SVD discard 1", "#177a99"),
    "baranger_ssm": ("Baranger SSM", "#4a9d6b"),
    "axis1_cut1": ("Axis 1", "#4354a5"),
    "axis2": ("Axis 2", "#d8852b"),
    "axis3_orthogonal": ("Axis 3 orthogonal", "#8b5eb2"),
    "axis3_loaded_gated": ("Axis 3 guarded loaded", "#bd495e"),
    "lee_aer": ("Stephen Lee AER", "#686868"),
}
SCENES = [("none",1.,"Static\nfast"),("rigid",1.,"Rigid\nfast"),
          ("smooth",1.,"Smooth\nfast"),("smooth",.02,"Smooth\nslow"),
          ("smooth",0.,"Smooth\nzero")]


def plot(scores, output):
    records=[json.loads(Path(p).read_text()) for p in scores]
    if (len({r["status"] for r in records}) != 1 or
            len({r["config_sha256"] for r in records}) != 1 or
            len({r["script_sha256"] for r in records}) != 1 or
            len({(r["seed"],r["motion_condition"],r["flow_scale"]) for r in records}) != len(records)):
        raise ValueError("Require one frozen configuration, script and distinct scenes/seeds")
    seeds={r["seed"] for r in records}
    if any({r["seed"] for r in records if
            r["motion_condition"]==motion and r["flow_scale"]==flow} != seeds
           for motion,flow,_ in SCENES):
        raise ValueError("Every plotted scene needs the same set of seeds")
    fig,ax=plt.subplots(2,1,figsize=(11,7.2),sharex=True,layout="constrained")
    for method,(label,color) in METHODS.items():
        for row,metric in enumerate(("blood_marginal_iq_nrmse",
                                      "refitted_background_output_power_to_input_tissue")):
            values=[];low=[];high=[]
            for motion,flow,_ in SCENES:
                selected=[r for r in records if r["motion_condition"]==motion and r["flow_scale"]==flow]
                mode="none" if motion=="none" else "smooth"
                v=np.array([r["results"][mode]["methods"][method][metric] for r in selected])
                values.append(np.median(v));low.append(np.min(v));high.append(np.max(v))
            x=np.arange(len(SCENES))
            ax[row].plot(x,values,marker="o",lw=1.7,ms=5,color=color,label=label)
            if len(seeds)>1:
                ax[row].fill_between(x,low,high,color=color,alpha=.09)
    ax[0].set_ylabel("Added-blood IQ NRMSE  ↓")
    ax[1].set_ylabel("Refitted background power / input tissue  ↓")
    ax[1].set_yscale("log")
    ax[1].set_xticks(range(len(SCENES)),[s[2] for s in SCENES])
    for axis in ax:
        axis.grid(axis="y",alpha=.25)
        axis.set_xlim(-.15,len(SCENES)-.85)
    ax[0].legend(ncol=4,fontsize=8,loc="upper center",bbox_to_anchor=(.5,1.28))
    fig.suptitle(f"eSDIva source-resolved {records[0]['status']}: "
                 f"{len(seeds)} seed{'s' if len(seeds)!=1 else ''}, "
                 "16 scored frames per scene/seed",fontsize=12)
    output=Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=180,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores",type=Path,nargs="+",required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    plot(args.scores,args.output)
