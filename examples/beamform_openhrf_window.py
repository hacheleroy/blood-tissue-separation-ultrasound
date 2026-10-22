"""Bounded local/HTTP RF window -> unfiltered IQ. Read docs/DATA_CONTRACT.md first.

Usage: python examples/beamform_openhrf_window.py source.hdf5 geometry.json out.npz
The geometry file is an acquisition-specific, manually verified contract.
This example deliberately refuses guessed timing/multiplexed reconstructions.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from btbench.beamform import plane_wave_das
from btbench.io import _open_h5,save_acquisition
from btbench.core import Acquisition


def main():
    import hdf5plugin  # register Blosc filters used in Waterloo HDF5 files
    p=argparse.ArgumentParser()
    p.add_argument("source");p.add_argument("geometry");p.add_argument("output")
    args=p.parse_args()
    g=json.loads(Path(args.geometry).read_text())
    if not g.get("time_reference_verified") or not g.get("full_receive_aperture_verified"):
        raise ValueError("Verify transmit time origin and receive multiplexing before RF reconstruction")
    start,stop=g["frames_start_stop"]
    if stop-start>2048 or stop<=start:
        raise ValueError("This CPU example is limited to 1..2048 contiguous frames")
    tx=g["transmit_indices"]
    root=g.get("root","").strip("/")
    def key(s):return root+"/"+s if root else s
    with _open_h5(args.source) as f:
        ds=f[key("data/raw_data")]
        estimate=(stop-start)*len(tx)*np.prod(ds.shape[2:])*ds.dtype.itemsize
        if estimate>g.get("max_raw_bytes",512*1024**2):
            raise MemoryError("RF window exceeds configured memory budget")
        raw=np.stack([ds[start:stop,i,...] for i in tx],axis=1)
        if raw.shape[-1]!=1:
            raise ValueError("This example accepts real RF (n_ch=1) only; DDC-IQ requires its own phase model")
        raw=raw[...,0]
        elements=f["probe/probe_geometry"][()][:,[2,0]]
    z=np.asarray(g["z_m"]);x=np.asarray(g["x_m"])
    views=plane_wave_das(raw,elements,z,x,g["tx_angles_rad"],g["fs_hz"],g["f0_hz"],
        g["initial_times_s"],g["tx_reference_times_s"],g.get("tx_origins_zx_m"),
        g.get("sound_speed_mps",1540.),g.get("f_number",1.5))
    # Preserve view data for STAP/vector Doppler; no angle averaging is forced.
    view_path=Path(args.output).with_suffix(".views.npy")
    np.save(view_path,views)
    a=Acquisition(views[:,0],g["revisit_prf_hz"],g["f0_hz"],
        sound_speed_mps=g.get("sound_speed_mps",1540.),phase_sign=g["phase_sign"],
        timestamps_s=np.asarray(g["frame_times_s"]),z_m=z,x_m=x,
        metadata={"source":args.source,"unfiltered_iq":True,"contrast_status":"unconfirmed",
        "single_view_index":tx[0],"views_file":view_path.name,"independent_truth_available":False,
        "beamformer":"bounded CPU far-field plane-wave DAS; provider comparison required"})
    save_acquisition(args.output,a)


if __name__=="__main__":main()
