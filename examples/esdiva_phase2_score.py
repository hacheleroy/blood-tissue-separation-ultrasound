"""Score longer eSDIva development/held-out IQ with frozen, data-only fits.

The mixture and separately acquired background arms share disjoint calibration
and a geometric tissue ROI. The evaluator holds isolated sources and prescribed
scatterer displacement; neither source labels nor known motion enter a filter
or gate. Scores on a new seed are valid only if code/settings were frozen first.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from btbench.shared import PublicInput, SharedContext, fit_shared
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_phase2_freeze import sha256, verify_freeze
from esdiva_gated_stap import fit_gated_loaded_stap


def score(acquisition, config_path, output, phase="development", freeze_path=None):
    acquisition, config_path, output = map(lambda p:Path(p).resolve(),
                                            (acquisition,config_path,output))
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    report = json.loads((acquisition/"report.json").read_text())
    config = json.loads(config_path.read_text())
    if phase not in ("development", "heldout"):
        raise ValueError("Unknown evaluation phase")
    if phase == "heldout":
        if freeze_path is None:
            raise ValueError("Held-out scoring requires the frozen manifest")
        manifest = verify_freeze(Path(__file__).resolve().parents[1], freeze_path,
                                 report["seed"], report["motion"], report["flow_scale"],report)
        if sha256(config_path) != manifest["files_sha256"]["configs/esdiva_phase2_development.json"]:
            raise ValueError("Scorer config differs from the frozen configuration")
    if not report["calibration_gate_pass"]:
        raise ValueError("Acquisition failed its input calibration gate")
    with np.load(acquisition/"iq_sources.npz") as arrays:
        mixture, cal = arrays["views"], arrays["calibration_views"]
        tissue, blood, noise = arrays["tissue"], arrays["blood"], arrays["noise"]
        roi = arrays["training_roi"]
        known = arrays["known_scatterer_displacement_pixels"]
    for key, value in (("mixture_sha256",mixture),("calibration_sha256",cal),
                       ("blood_sha256",blood),("tissue_sha256",tissue),
                       ("noise_sha256",noise),("known_motion_sha256",known),
                       ("training_roi_sha256",roi)):
        if digest(value) != report[key]:
            raise ValueError(f"Acquisition array hash mismatch: {key}")
    if relative_error(blood+tissue+noise,mixture)>1e-8:
        raise ValueError("IQ source additivity failed")
    packet,score_start,horizon=(config[k] for k in ("packet","score_start","horizon"))
    if (mixture.shape[0] < score_start+packet or horizon < packet or
            horizon%packet or score_start%packet):
        raise ValueError("Invalid score/history support")
    public = PublicInput(mixture,cal,roi,report["geometry"])
    background = PublicInput(tissue+noise,cal,roi,report["geometry"])
    results = {}
    contexts = {}
    for mode in config["motion_modes"]:
        t0 = time.monotonic()
        options=config.get("motion_options") if mode!="none" else None
        mixed = SharedContext(public,packet=packet,score_start=score_start,motion=mode,
                              motion_options=options)
        bg = SharedContext(background,packet=packet,score_start=score_start,motion=mode,
                           motion_options=options)
        if mixed.calibration != bg.calibration:
            raise AssertionError("Counterfactual calibration changed")
        contexts[mode] = (mixed,bg)
        estimated = mixed.displacement
        motion_rmse = float(np.sqrt(np.mean((estimated-known)**2)))
        zero_rmse = float(np.sqrt(np.mean(known**2)))
        motion_record = {"registration_mode":mode,"fit":mixed.motion_summary,
                         "background_fit":bg.motion_summary,
                         "estimated_nonzero_frame_fraction":float(np.mean(np.any(estimated!=0,axis=(1,2,3)))),
                         "scatterer_displacement_rmse_pixels_evaluation_only":motion_rmse,
                         "zero_displacement_rmse_pixels_evaluation_only":zero_rmse,
                         "context_elapsed_s":time.monotonic()-t0}
        results[mode] = {"motion":motion_record,"methods":{}}
    score_slice=slice(score_start,None)
    def ratio(value, reference):
        return float(np.mean(abs(value)**2)/max(np.mean(abs(reference)**2),1e-30))
    for mode,(mixed,bg) in contexts.items():
        for method in config["methods"]:
            name=method["id"]
            if mode!="none" and name=="identity":
                continue
            spec={k:v for k,v in method.items() if k!="id"}
            start=time.monotonic()
            try:
                if name=="axis3_loaded_gated":
                    fitted,gate=fit_gated_loaded_stap(mixed,spec,
                        config["stap_min_snapshots_per_dimension"],
                        config["stap_max_operator_gain"])
                    control,background_gate=fit_gated_loaded_stap(bg,spec,
                        config["stap_min_snapshots_per_dimension"],
                        config["stap_max_operator_gain"])
                else:
                    fitted=fit_shared(mixed,spec)
                    control=fit_shared(bg,spec)
                    gate=background_gate=None
                b=fitted.apply_frozen(blood)
                t=fitted.apply_frozen(tissue)
                n=fitted.apply_frozen(noise)
                marginal=fitted.views-control.views
                blocks=fitted.diagnostics["blocks"]
                results[mode]["methods"][name]={
                    "status":"ok","spec":spec,"mixture_gate":gate,
                    "background_gate":background_gate,
                    "blood_marginal_iq_nrmse":relative_error(marginal,blood[score_slice]),
                    "blood_marginal_power_ratio":ratio(marginal,blood[score_slice]),
                    "refitted_background_output_power_to_input_tissue":ratio(control.views,tissue[score_slice]),
                    "mixture_blood_iq_nrmse":relative_error(fitted.views,blood[score_slice]),
                    "frozen_blood_power_ratio":ratio(b,blood[score_slice]),
                    "frozen_tissue_power_ratio":ratio(t,tissue[score_slice]),
                    "frozen_noise_power_ratio":ratio(n,noise[score_slice]),
                    "recomposition_relative_error":relative_error(b+t+n,fitted.views),
                    "max_operator_gain":float(max((x.get("operator_gain_max",1.) for x in blocks),default=1.)),
                    "training_snapshots_min":min((x.get("training_snapshots_min",0) for x in blocks),default=0)
                      if spec["method"]=="stap" else None,
                    "elapsed_s":time.monotonic()-start}
            except (ValueError,np.linalg.LinAlgError) as exc:
                results[mode]["methods"][name]={"status":"failed","spec":spec,"error":str(exc)}
    output.parent.mkdir(parents=True,exist_ok=True)
    outcome={"status":"development" if phase=="development" else "held-out evaluation",
             "acquisition_report_sha256":hashlib.sha256((acquisition/"report.json").read_bytes()).hexdigest(),
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "gate_wrapper_sha256":hashlib.sha256(Path(__file__).with_name("esdiva_gated_stap.py").read_bytes()).hexdigest(),
             "numerical_source_sha256":numerical_source_hash(Path(__file__).resolve().parents[1]),
             "config_sha256":hashlib.sha256(config_path.read_bytes()).hexdigest(),
             "seed":report["seed"],"motion_condition":report["motion"],
             "flow_scale":report["flow_scale"],"particle_fraction":report["scatterer_fraction"],
             "scored_frames":mixture.shape[0]-score_start,
             "source_and_known_motion_are_evaluator_only":True,"results":results,
             "limitation":"known scatterer displacement is a geometric reference, not an exact IQ warp; particle density and motion/flow strata must be stated"}
    output.write_text(json.dumps(outcome,indent=2,allow_nan=False)+"\n")
    return outcome


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--phase",choices=("development","heldout"),default="development")
    p.add_argument("--freeze",type=Path,help="Committed freeze manifest, required for held-out scoring")
    a=p.parse_args()
    result=score(a.acquisition,a.config,a.output,a.phase,a.freeze)
    print(json.dumps({"status":result["status"],"seed":result["seed"],
                      "scored_frames":result["scored_frames"]},indent=2))
