"""Open a reserved eSDIva seed only after verifying the frozen protocol.

The source-resolved RF/IQ stays outside Git; its report and compact score can
be reviewed and then committed without exposing arrays to the filters.
"""

import argparse
import json
from pathlib import Path

from esdiva_phase2_acquire import acquire
from esdiva_phase2_freeze import verify_freeze
from esdiva_phase2_score import score


def run(upstream_root, output, seed, scene_id, freeze_path):
    root = Path(__file__).resolve().parents[1]
    scenarios = json.loads((root/"configs/esdiva_phase2_scenarios.json").read_text())
    scene = next((s for s in scenarios["scenes"] if s["id"] == scene_id), None)
    if scene is None:
        raise ValueError("Unknown frozen scene")
    verify_freeze(root, freeze_path, seed, scene["motion"], scene["flow_scale"])
    output = Path(output).resolve()
    if (output/"score.json").exists():
        raise FileExistsError("Refusing to overwrite held-out scores")
    report_path = output/"report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text())
        verify_freeze(root, freeze_path, seed, scene["motion"], scene["flow_scale"],report)
    else:
        report = acquire(upstream_root,output,seed,
                         scenarios["particle_fraction"],scenarios["calibration_frames"],
                         scenarios["evaluation_frames"],scene["motion"],scene["flow_scale"])
    if not report["calibration_gate_pass"]:
        raise ValueError("Held-out input calibration gate failed; no filtering score")
    return score(output,root/"configs/esdiva_phase2_development.json",
                 output/"score.json",phase="heldout",freeze_path=freeze_path)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--seed",type=int,required=True)
    parser.add_argument("--scene",required=True)
    parser.add_argument("--freeze",type=Path,required=True)
    args=parser.parse_args()
    result=run(args.upstream_root,args.output,args.seed,args.scene,args.freeze)
    print(json.dumps({k:result[k] for k in ("status","seed","motion_condition",
                                             "flow_scale","scored_frames")},indent=2))
