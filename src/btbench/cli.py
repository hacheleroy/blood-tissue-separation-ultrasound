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
    p.add_argument("--prf-hz",type=float,default=1000.)
    p.add_argument("--packet",type=int,default=64);p.add_argument("--long-window",type=int,default=256)
    p=sub.add_parser("export-openhrf-rf");p.add_argument("source");p.add_argument("geometry");p.add_argument("output")
    p=sub.add_parser("separation-study");p.add_argument("--config");p.add_argument("--output",required=True)
    p=sub.add_parser("stap-study");p.add_argument("--config");p.add_argument("--output",required=True)
    p=sub.add_parser("literature-study");p.add_argument("--config");p.add_argument("--output",required=True)
    p=sub.add_parser("fidelity-followup");p.add_argument("--config");p.add_argument("--output",required=True)
    p=sub.add_parser("shared-study");p.add_argument("--config");p.add_argument("--output",required=True)
    p=sub.add_parser("openhrf-followup");p.add_argument("datasets",nargs='+');p.add_argument("--config",required=True);p.add_argument("--output",required=True)
    args=parser.parse_args(argv)
    if args.cmd=="simulate":
        save_acquisition(args.output,poiseuille(seed=args.seed))
    elif args.cmd=="write-config":
        Path(args.output).parent.mkdir(parents=True,exist_ok=True)
        Path(args.output).write_text(json.dumps({"blas_threads":1,
            "methods":standard_methods(args.packet,args.long_window,args.prf_hz)},indent=2))
    elif args.cmd=="run":
        r=run_benchmark(load_acquisition(args.dataset),json.loads(Path(args.config).read_text()),args.output)
        print(json.dumps({"methods":len(r["results"]),"ok":sum(x["status"]=="ok" for x in r["results"]),"output":args.output}))
    elif args.cmd=="demo":
        a=poiseuille(seed=args.seed)
        cfg={"blas_threads":1,"save_iq":False,"methods":standard_methods(prf_hz=a.prf_hz)}
        for centered in (False,True):
            cfg["methods"].append({"id":f"oracle_svd_preservation_centered{centered}","method":"svd",
                "params":{"centered":centered,"window":64},"oracle_grid":{"discard":list(range(9))},
                "objective":"frozen_tissue_power_retention",
                "constraints":{"frozen_blood_iq_nrmse":{"max":.5}}})
        cfg["methods"].append({"id":"oracle_aer_iq","method":"aer","params":{"centered":False},
            "oracle_grid":{"epsilon_rel":[1e-5,1e-3,0.1],"alpha":[0.25,0.5]},"objective":"blood_iq_nrmse"})
        r=run_benchmark(a,cfg,args.output)
        save_acquisition(Path(args.output)/"synthetic_iq.npz",a)
        print(json.dumps({"methods":len(r["results"]),"ok":sum(x["status"]=="ok" for x in r["results"]),"elapsed_s":r["runtime"]["elapsed_s"],"output":args.output}))
    elif args.cmd=="inspect-openhrf":
        data=inspect_openhrf(args.source)
        result=json.dumps(data,indent=2)
        if args.output: Path(args.output).write_text(result)
        else: print(result)
    elif args.cmd=="export-openhrf-rf":
        from .io import export_openhrf_rf_window
        a=export_openhrf_rf_window(args.source,json.loads(Path(args.geometry).read_text()),args.output)
        print(json.dumps({"iq_shape":list(a.iq.shape),"views":list(a.iq_views.shape),"output":args.output}))
    elif args.cmd=="separation-study":
        from .separation_experiment import run_separation_study,default_separation_config
        config=json.loads(Path(args.config).read_text()) if args.config else default_separation_config()
        result=run_separation_study(config,args.output)
        print(json.dumps({"cases":len(result["records"]),"output":args.output}))
    elif args.cmd=="stap-study":
        from .stap_experiment import run_stap_study,default_stap_config
        config=json.loads(Path(args.config).read_text()) if args.config else default_stap_config()
        result=run_stap_study(config,args.output)
        print(json.dumps({"cases":len(result["records"]),"output":args.output,"elapsed_s":result["runtime"]["elapsed_s"]}))
    elif args.cmd=="literature-study":
        from .literature_experiment import run_literature_study,default_literature_config
        config=json.loads(Path(args.config).read_text()) if args.config else default_literature_config()
        result=run_literature_study(config,args.output)
        print(json.dumps({"cases":len(result["records"]),"output":args.output}))
    elif args.cmd=="fidelity-followup":
        from .followup_experiment import default_followup_config,run_followup_study
        config=json.loads(Path(args.config).read_text()) if args.config else default_followup_config()
        result=run_followup_study(config,args.output)
        print(json.dumps({"cases":len(result['records']),"output":args.output}))
    elif args.cmd=="shared-study":
        from .shared_experiment import default_shared_config,run_shared_study
        config=json.loads(Path(args.config).read_text()) if args.config else default_shared_config()
        result=run_shared_study(config,args.output)
        print(json.dumps({"cases":len(result['records']),"output":args.output}))
    elif args.cmd=="openhrf-followup":
        from .followup_experiment import run_openhrf_followup
        result=run_openhrf_followup([load_acquisition(p) for p in args.datasets],json.loads(Path(args.config).read_text()),args.output)
        print(json.dumps({"acquisitions":len(result['records']),"output":args.output}))


if __name__=="__main__":
    main()
