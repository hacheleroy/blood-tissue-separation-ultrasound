"""Small prospectively configured Cartesian experiment bank (statistical IQ).

Usage: python examples/run_phase_diagram.py configs/phase_diagram_small.json runs/grid
All simulator axes and method settings are fixed before reading outcomes.
This is an implementation test, not a realistic acoustic identifiability map.
"""
import argparse
import itertools
import json
from pathlib import Path
from btbench.synthetic import poiseuille
from btbench.benchmark import run_benchmark


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("config");parser.add_argument("output")
    args=parser.parse_args()
    config=json.loads(Path(args.config).read_text())
    root=Path(args.output)
    if root.exists():raise FileExistsError("Choose a new sweep output directory")
    root.mkdir(parents=True)
    (root/"frozen_config.json").write_text(json.dumps(config,indent=2))
    keys=list(config["simulation_grid"])
    summary=[]
    for index,values in enumerate(itertools.product(*(config["simulation_grid"][k] for k in keys))):
        parameters={**config.get("simulation_base",{}),**dict(zip(keys,values))}
        a=poiseuille(**parameters)
        a.metadata["simulation_parameters"]=parameters
        case=f"case_{index:04d}"
        report=run_benchmark(a,config["benchmark"],root/case)
        strict_report=json.loads((root/case/"results.json").read_text())
        summary.append({"case":case,"parameters":parameters,"input_sha256":report["input_sha256"],
            "results":[{k:r[k] for k in ("id","status","metrics") if k in r} for r in strict_report["results"]]})
    (root/"phase_diagram.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps({"cases":len(summary),"output":str(root)}))


if __name__=="__main__":main()
