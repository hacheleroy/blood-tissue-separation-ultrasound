"""Acquire reserved eSDIva RF only after verifying the published phase-4 freeze."""

import argparse
import json
from pathlib import Path

from esdiva_phase2_acquire import acquire
from esdiva_phase4_freeze import verify_freeze
from esdiva_phase4_score import evaluate


def run(upstream_root, output, seed, scene_name, freeze_path):
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root/"configs/esdiva_phase4_candidate.json").read_text())
    verify_freeze(root, freeze_path, seed, scene_name)
    output = Path(output).resolve()
    score_path = output/"phase4_score.json"
    if score_path.exists():
        raise FileExistsError("Refusing to overwrite reserved score")
    report_path = output/"report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text())
    else:
        s = config["scenes"][scene_name]
        report = acquire(upstream_root, output, seed,
                         config["scatterer_fraction"], config["calibration_frames"],
                         config["evaluation_frames"], s["motion"], s["flow_scale"])
    verify_freeze(root, freeze_path, seed, scene_name, report)
    score = evaluate(output, root/"configs/esdiva_phase4_candidate.json", score_path)
    if score["scene"] != scene_name:
        raise AssertionError("Scored a different frozen scene")
    return score


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--seed",type=int,required=True)
    parser.add_argument("--scene",required=True)
    parser.add_argument("--freeze",type=Path,required=True)
    args = parser.parse_args()
    result = run(args.upstream_root,args.output,args.seed,args.scene,args.freeze)
    print(json.dumps({k:result[k] for k in ("status","seed","scene","score_frames")}))
