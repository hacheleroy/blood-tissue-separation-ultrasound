import argparse
import json
from btbench.io import load_acquisition
from btbench.experimental import run_openhrf_pilot


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("dataset");p.add_argument("output")
    p.add_argument("--calibration-frames",type=int,default=32)
    p.add_argument("--spatial-stride",type=int,default=2)
    args=p.parse_args()
    r=run_openhrf_pilot(load_acquisition(args.dataset),args.output,args.calibration_frames,args.spatial_stride)
    print(json.dumps({"outputs":len(r["results"]),"ok":sum(x["status"]=="ok" for x in r["results"]),
        "code_sha256":r["code_sha256"],"output":args.output}))
