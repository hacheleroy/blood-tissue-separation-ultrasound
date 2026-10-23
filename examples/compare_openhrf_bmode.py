"""Provider B-mode structural check; no quantitative phase/velocity claim."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from btbench.io import _open_h5,load_acquisition


def main():
    p=argparse.ArgumentParser();p.add_argument("dataset");p.add_argument("output")
    p.add_argument("--cached-reference",help="Optional locally cached image_db/coordinates NPZ")
    args=p.parse_args();a=load_acquisition(args.dataset)
    if args.cached_reference:
        with np.load(args.cached_reference) as r:reference,coords=r["image_db"],r["coordinates"]
    else:
        with _open_h5(a.metadata["source"]) as f:
            root="tracks/track_0/data/image/"
            reference=f[root+"values"][a.metadata["frames_start_stop"][0]]
            coords=f[root+"coordinates"][()]
    zz,xx=np.meshgrid(a.z_m,a.x_m,indexing="ij")
    ref=RegularGridInterpolator((coords[:,0,2],coords[0,:,0]),reference)(np.stack((zz,xx),axis=-1))
    own=20*np.log10(np.maximum(abs(a.iq[0]),1e-30));own-=own.max()
    corr=float(np.corrcoef(own.ravel(),ref.ravel())[0,1])
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps({"source":a.metadata["source"],"frame":a.metadata["frames_start_stop"][0],
        "log_amplitude_spatial_correlation":corr,"interpolation":"provider dB image to exported Cartesian grid",
        "interpretation":"structure check; different receive steering/apodization; no phase/velocity validation"},indent=2)+"\n")
    fig,axes=plt.subplots(1,2,figsize=(9,4.7),layout="constrained")
    extent=[a.x_m[0]*1e3,a.x_m[-1]*1e3,a.z_m[-1]*1e3,a.z_m[0]*1e3]
    for ax,img,title in zip(axes,[ref,own],["Provider B-mode","CPU DAS, coherent view mean"]):
        im=ax.imshow(img,extent=extent,aspect="equal",cmap="gray",vmin=-60,vmax=0)
        ax.set(title=title,xlabel="Lateral position (mm)",ylabel="Depth (mm)")
    fig.colorbar(im,ax=axes,label="Display amplitude (dB)",shrink=.8)
    fig.suptitle(f"Waterloo CCA · first frame · structural correlation {corr:.3f}\nImage structure alone does not validate Doppler phase",fontsize=12)
    fig.savefig(out/"bmode_comparison.png",dpi=180)
    plt.close(fig)
    print(json.dumps({"correlation":corr,"output":str(out)}))


if __name__=="__main__":main()
