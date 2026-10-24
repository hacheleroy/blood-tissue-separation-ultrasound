"""Pair blood-addition and background fits after the phase-3 input gate.

Only the evaluator opens isolated tissue/blood/noise and known transport. All
filters receive mixture or independently refitted background IQ, the selected
geometric ROI, and the same tissue calibration. This is development scoring;
there is no held-out claim or parameter selection from a filter score.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_gated_stap import fit_gated_loaded_stap
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase3_input_gate import choose_input


def score(acquisition, config_path, output):
    acquisition,config_path,output=map(lambda p:Path(p).resolve(),
                                       (acquisition,config_path,output))
    if output.exists():
        raise FileExistsError(output)
    report=json.loads((acquisition/"report.json").read_text())
    config=json.loads(config_path.read_text())
    with np.load(acquisition/"iq_sources.npz") as data:
        mixture,cal,tissue,blood,noise=(data[k] for k in
                                          ("views","calibration_views","tissue","blood","noise"))
        x_mm=data["x_mm"]
    for key,value in (("mixture_sha256",mixture),("calibration_sha256",cal),
                      ("blood_sha256",blood),("tissue_sha256",tissue),
                      ("noise_sha256",noise)):
        if digest(value)!=report[key]:
            raise ValueError(f"IQ acquisition hash mismatch: {key}")
    if (relative_error(blood+tissue+noise,mixture)>1e-8 or
            report["view_indices"]!=[4,3,5] or len(cal)!=16 or len(mixture)!=80):
        raise ValueError("Unsupported source or scene geometry")
    roi,input_audit=choose_input(cal,x_mm)
    if not input_audit["input_gate_pass"]:
        raise ValueError("Phase-3 calibration gate failed; no filter scores")
    ratio=lambda a,b:float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))
    blood_roi=ratio(blood[:,:,roi],tissue[:,:,roi])
    if blood_roi>.1:
        raise ValueError("Selected tissue ROI contaminated by blood; no filter scores")
    public=PublicInput(mixture,cal,roi,report["geometry"])
    background=PublicInput(tissue+noise,cal,roi,report["geometry"])
    score_slice=slice(config["score_start"],None)
    records={}
    modes=["none"] if report["motion"]=="none" else ["smooth"]
    for mode in modes:
        options=config["motion_options"] if mode!="none" else None
        for gain_model in ("regression8","phase_times_rms"):
            start=time.monotonic()
            mixed=Phase3Context(public,packet=config["packet"],
                                score_start=config["score_start"],motion=mode,
                                motion_options=options,gain_model=gain_model)
            bg=Phase3Context(background,packet=config["packet"],
                             score_start=config["score_start"],motion=mode,
                             motion_options=options,gain_model=gain_model)
            if mixed.calibration!=bg.calibration:
                raise AssertionError("Paired calibration differs")
            methods={}
            for item in config["methods"]:
                name=item["id"]
                spec={k:v for k,v in item.items() if k!="id"}
                t=time.monotonic()
                try:
                    if name=="axis3_loaded_gated":
                        fit,gate=fit_gated_loaded_stap(mixed,spec,
                            config["stap_min_snapshots_per_dimension"],
                            config["stap_max_operator_gain"])
                        control,bg_gate=fit_gated_loaded_stap(bg,spec,
                            config["stap_min_snapshots_per_dimension"],
                            config["stap_max_operator_gain"])
                    else:
                        fit,control=fit_shared(mixed,spec),fit_shared(bg,spec)
                        gate=bg_gate=None
                    marginal=fit.views-control.views
                    methods[name]={"status":"ok", "spec":spec,
                        "added_blood_iq_nrmse":relative_error(marginal,blood[score_slice]),
                        "background_output_to_input_tissue_power":ratio(
                            control.views,tissue[score_slice]),
                        "frozen_blood_power_ratio":ratio(fit.apply_frozen(blood),blood[score_slice]),
                        "mixture_gate":gate,"background_gate":bg_gate,
                        "elapsed_s":time.monotonic()-t}
                except (ValueError,np.linalg.LinAlgError) as exc:
                    methods[name]={"status":"failed","error":str(exc),"spec":spec}
            records[gain_model]={"calibration":mixed.calibration,
                                 "motion_fit":mixed.motion_summary,
                                 "background_motion_fit":bg.motion_summary,
                                 "context_elapsed_s":time.monotonic()-start,
                                 "methods":methods}
    outcome={"status":"phase-3 development source-resolved paired score",
             "seed":report["seed"],"motion":report["motion"],"flow_scale":report["flow_scale"],
             "scatterer_fraction":report["scatterer_fraction"],
             "acquisition_sha256":hashlib.sha256((acquisition/"report.json").read_bytes()).hexdigest(),
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "input_gate_sha256":hashlib.sha256(Path(__file__).with_name(
                 "esdiva_phase3_input_gate.py").read_bytes()).hexdigest(),
             "bounded_context_sha256":hashlib.sha256(Path(__file__).with_name(
                 "esdiva_phase3_bounded_context.py").read_bytes()).hexdigest(),
             "config_sha256":hashlib.sha256(config_path.read_bytes()).hexdigest(),
             "numerical_source_sha256":numerical_source_hash(Path(__file__).resolve().parents[1]),
             "selected_training_roi_sha256":digest(roi),"input_gate":input_audit,
             "blood_to_tissue_power_selected_roi_evaluator_only":blood_roi,
             "scored_frames":len(mixture)-config["score_start"],
             "source_and_known_motion_are_evaluator_only":True,
             "results":records}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(outcome,indent=2,allow_nan=False)+"\n")
    return outcome


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=score(a.acquisition,a.config,a.output)
    print(json.dumps({"seed":result["seed"],"roi":result["input_gate"]["selected_roi"],
                      "scored_frames":result["scored_frames"]}))
