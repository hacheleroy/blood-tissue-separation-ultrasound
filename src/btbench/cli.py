import argparse
import json
from pathlib import Path
from .io import load_acquisition,save_acquisition,inspect_openhrf
from .synthetic import poiseuille
from .benchmark import run_benchmark,standard_methods


def main(argv=None):
    parser=argparse.ArgumentParser(description="Blood/tissue complex-IQ research benchmark")
    sub=parser.add_subparsers(dest="cmd",required=True)
    p=sub.add_parser("demo");p.add_argument("--output",required=True);p.add_argument("--seed",type=int,default=0)
    p=sub.add_parser("simulate");p.add_argument("output");p.add_argument("--seed",type=int,default=0)
    p=sub.add_parser("run");p.add_argument("dataset");p.add_argument("--config",required=True);p.add_argument("--output",required=True)
    p=sub.add_parser("inspect-openhrf");p.add_argument("source");p.add_argument("--output")
    p=sub.add_parser("write-config");p.add_argument("output")
    args=parser.parse_args(argv)
    if args.cmd=="simulate":
        save_acquisition(args.output,poiseuille(seed=args.seed))
    elif args.cmd=="write-config":
        Path(args.output).parent.mkdir(parents=True,exist_ok=True)
        Path(args.output).write_text(json.dumps({"blas_threads":1,"methods":standard_methods()},indent=2))
    elif args.cmd=="run":
        r=run_benchmark(load_acquisition(args.dataset),json.loads(Path(args.config).read_text()),args.output)
        print(json.dumps({"methods":len(r["results"]),"ok":sum(x["status"]=="ok" for x in r["results"]),"output":args.output}))
    elif args.cmd=="demo":
        a=poiseuille(seed=args.seed)
        cfg={"blas_threads":1,"save_iq":False,"methods":standard_methods()}
        for centered in (False,True):
            cfg["methods"].append({"id":f"oracle_svd_velocity_centered{centered}","method":"svd",
                "params":{"centered":centered,"window":64},"oracle_grid":{"discard":list(range(9))},"objective":"velocity_nrmse"})
        cfg["methods"].append({"id":"oracle_aer_velocity","method":"aer","params":{"centered":False},
            "oracle_grid":{"epsilon_rel":[1e-5,1e-3,0.1],"alpha":[0.25,0.5]},"objective":"velocity_nrmse"})
        r=run_benchmark(a,cfg,args.output)
        save_acquisition(Path(args.output)/"synthetic_iq.npz",a)
        print(json.dumps({"methods":len(r["results"]),"ok":sum(x["status"]=="ok" for x in r["results"]),"elapsed_s":r["runtime"]["elapsed_s"],"output":args.output}))
    elif args.cmd=="inspect-openhrf":
        data=inspect_openhrf(args.source)
        result=json.dumps(data,indent=2)
        if args.output: Path(args.output).write_text(result)
        else: print(result)


if __name__=="__main__":
    main()
