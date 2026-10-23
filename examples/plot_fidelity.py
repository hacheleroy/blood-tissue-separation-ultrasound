"""Plot all three seeds for four illustrative, explicitly labeled conditions."""
import argparse
import json
import gzip
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser();p.add_argument("results");p.add_argument("output")
    args=p.parse_args()
    raw=gzip.decompress(Path(args.results).read_bytes()) if args.results.endswith(".gz") else Path(args.results).read_bytes()
    report=json.loads(raw)
    choices=[("svd_fixed_uncentered","SVD, short"),("svd_fixed_centered","SVD, centered"),
        ("svd_long_uncentered","SVD, long"),("svd_baranger_uncentered","Baranger SSM"),
        ("butterworth_hp_0.03_prf","Butterworth, 30 Hz"),
        ("aer_temporal_alpha0.5_uncentered","Temporal AER, raw gain"),
        ("axis1_uncentered","Axis 1"),("axis2_with_memory_uncentered","Axis 2, memory"),
        ("ica_circular_complex","ICA, guarded")]
    panels=[("slow_overlap","Statistical IQ: slow spectral overlap"),
        ("weak_decorrelating","Statistical IQ: weak, decorrelating blood"),
        ("psf_deforming_tissue","PSF transport: deforming surrounding tissue"),
        ("psf_equal_static","PSF transport: both sources static")]
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout="constrained")
    colors=plt.get_cmap("tab10").colors
    for ax,(condition,title) in zip(axes.flat,panels):
        for j,(name,label) in enumerate(choices):
            rows=[r for c in report["records"] if c["condition"]["label"]==condition
                for r in c["results"] if r["id"]==name and r["status"]=="ok"]
            values=[(r["metrics"].get("frozen_blood_iq_nrmse"),r["metrics"].get("frozen_tissue_power_retention")) for r in rows]
            values=[(max(x,1e-6),max(y,1e-12)) for x,y in values if x is not None and y is not None]
            if values:
                xy=np.asarray(values);ax.scatter(xy[:,0],xy[:,1],s=37,color=colors[j],label=label,alpha=.85,edgecolor="white",linewidth=.4)
        ax.set(xscale="log",yscale="log",title=title,xlabel="Isolated blood IQ distortion εB",ylabel="Tissue power retention rT")
        ax.grid(True,which="major",alpha=.2)
        ax.axvline(1,color=".7",lw=.8,ls="--")
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=3,frameon=False)
    fig.suptitle("Tissue rejection and blood fidelity · each point is one seed\nLower is better on both axes; no gain or phase fitted to truth",fontsize=13)
    fig.get_layout_engine().set(rect=(0,.12,1,.88))
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output.with_suffix(".png"),dpi=180,bbox_inches="tight")
    svg=output.with_suffix(".svg")
    fig.savefig(svg,bbox_inches="tight")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    plt.close(fig)


if __name__=="__main__":main()
