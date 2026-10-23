"""Bounded RF -> unfiltered multi-view IQ. See docs/DATA_CONTRACT.md.

Usage: python examples/beamform_openhrf_window.py source.hdf5 geometry.json out.npz
"""
import argparse
import json
from pathlib import Path
from btbench.io import export_openhrf_rf_window


def main():
    p=argparse.ArgumentParser()
    p.add_argument("source");p.add_argument("geometry");p.add_argument("output")
    args=p.parse_args()
    a=export_openhrf_rf_window(args.source,json.loads(Path(args.geometry).read_text()),args.output)
    print(json.dumps({"iq_shape":list(a.iq.shape),"view_shape":list(a.iq_views.shape),
        "revisit_prf_hz":a.prf_hz,"view_offsets_s":a.view_offsets_s.tolist(),"output":args.output}))


if __name__=="__main__":main()
