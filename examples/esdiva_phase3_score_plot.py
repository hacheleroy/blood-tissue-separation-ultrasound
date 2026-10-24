"""Draw paired source-response development outcomes for the bounded candidate."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


METHODS = (("svd_cut1", "SVD discard 1"),
           ("baranger_ssm", "Baranger SSM"),
           ("axis1_cut1", "Axis 1: Motion-compensated virtual long ensemble"),
           ("axis2", "Axis 2: Motion-demodulated multi-lag separation"),
           ("axis3_orthogonal", "Axis 3: Local ultrasound STAP, orthogonal"),
           ("axis3_loaded_gated", "Axis 3: Local ultrasound STAP, guarded loaded"),
           ("lee_aer", "Stephen Lee AER"))


def plot(paths, output):
    records=[json.loads(Path(path).read_text()) for path in paths]
    names={"none":"Stationary fast", "rigid":"Rigid fast", "smooth_slow":"Smooth slow",
           "smooth_zero":"Smooth zero"}
    keyed={("none" if d["motion"]=="none" else
            "rigid" if d["motion"]=="rigid" else
            "smooth_zero" if d["flow_scale"]==0 else "smooth_slow"):d for d in records}
    if set(keyed)!=set(names) or any(d["seed"]!=2083 or not d["input_gate"]["input_gate_pass"]
                                or d["input_gate"]["min_roi_pixels"]!=40 for d in records):
        raise ValueError("Require four comparable supported-ROI development scenes")
    matrices=[]
    for metric in ("added_blood_iq_nrmse","background_output_to_input_tissue_power"):
        values=np.array([[keyed[k]["results"]["phase_times_rms"]["methods"][method][metric]
                          for k in names] for method,_ in METHODS])
        if not np.isfinite(values).all():
            raise ValueError("Missing finite paired score")
        matrices.append(values)
    fig,axes=plt.subplots(1,2,figsize=(16,6),layout="constrained")
    for ax,metric,title,cmap,norm in zip(axes,matrices,
            ("Added-blood IQ NRMSE ↓","Refitted background / input tissue power ↓"),
            ("YlOrRd","YlOrRd"),(None,LogNorm(vmin=.0007,vmax=1.))):
        im=ax.imshow(metric,aspect="auto",cmap=cmap,norm=norm,vmin=None if norm else 0,
                     vmax=None if norm else 1.3)
        for i in range(metric.shape[0]):
            for j in range(metric.shape[1]):
                ax.text(j,i,f"{metric[i,j]:.3f}",ha="center",va="center",fontsize=8,
                        color="white" if (metric[i,j]>.72 if norm is None else metric[i,j]>.15)
                        else "black")
        ax.set_xticks(range(4),list(names.values()),rotation=40,ha="right")
        ax.set_yticks(range(len(METHODS)),[label for _,label in METHODS],fontsize=8)
        ax.set_title(title)
        fig.colorbar(im,ax=ax,fraction=.046,pad=.03)
    fig.suptitle("eSDIva development seed 2083; bounded calibration, 16 scored frames per scene")
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=180)
    plt.close(fig)
    return matrices


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("scores",type=Path,nargs="+")
    a=p.parse_args()
    plot(a.scores,a.output)
