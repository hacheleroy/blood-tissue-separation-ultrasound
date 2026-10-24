"""Score data-only RF motion candidates against paired source controls.

Development input is the same verified eSDIva three-view IQ archive used by
phase 3. Isolated blood/tissue and prescribed scatterer shifts enter only the
evaluator. The latter are a geometric reference, not an exact IQ image warp.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, relative_error
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase3_input_gate import choose_input
from esdiva_phase4_motion_context import EnvelopeMotionContext


METHODS = {
    "roundtrip": {"method":"roundtrip", "horizon":64},
    "svd_cut1": {"method":"svd", "horizon":64, "discard":1},
    "axis1_cut1": {"method":"axis1", "horizon":64, "discard":1, "rho":.75},
}


def evaluate(acquisition, output, variants=("none","complex_smooth","envelope_rigid",
                                              "envelope_smooth")):
    acquisition,output=Path(acquisition).resolve(),Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    report=json.loads((acquisition/"report.json").read_text())
    with np.load(acquisition/"iq_sources.npz") as f:
        mixture,cal,tissue,blood,noise,known,x=(f[k] for k in (
            "views","calibration_views","tissue","blood","noise",
            "known_scatterer_displacement_pixels","x_mm"))
    for key,value in (("mixture_sha256",mixture),("calibration_sha256",cal),
                      ("tissue_sha256",tissue),("blood_sha256",blood),
                      ("noise_sha256",noise),("known_motion_sha256",known)):
        if digest(value)!=report[key]:
            raise ValueError(f"RF acquisition hash mismatch: {key}")
    roi,gate=choose_input(cal,x,min_roi_pixels=40)
    if not gate["input_gate_pass"] or relative_error(tissue+blood+noise,mixture)>1e-8:
        raise ValueError("Source or input gate failed; no motion score")
    ratio=lambda a,b:float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))
    contamination=ratio(blood[:,:,roi],tissue[:,:,roi])
    if contamination>.1:
        raise ValueError("Training region is blood-contaminated")
    public=PublicInput(mixture,cal,roi,report["geometry"])
    bg_public=PublicInput(tissue+noise,cal,roi,report["geometry"])
    result={}
    for variant in variants:
        if variant=="none":
            cls,mode,opts=Phase3Context,"none",None
        elif variant=="complex_smooth":
            cls,mode,opts=Phase3Context,"smooth",{"maximum_validation_error":.12}
        elif variant in ("envelope_rigid","envelope_smooth"):
            cls,mode,opts=EnvelopeMotionContext,variant.split("_")[1],{
                "maximum_validation_error":.12}
        else:
            raise ValueError(f"Unknown candidate: {variant}")
        mixed=cls(public,motion=mode,motion_options=opts)
        bg=cls(bg_public,motion=mode,motion_options=opts)
        if mixed.calibration!=bg.calibration:
            raise AssertionError("Counterfactual calibration changed")
        methods={}
        for name,spec in METHODS.items():
            a=fit_shared(mixed,spec);b=fit_shared(bg,spec)
            methods[name]={"added_blood_iq_nrmse":relative_error(a.views-b.views,blood[64:]),
                           "background_to_input_tissue_power":ratio(b.views,tissue[64:]),
                           "mixture_to_blood_iq_nrmse":relative_error(a.views,blood[64:])}
        result[variant]={"mixture_motion":mixed.motion_summary,
                         "background_motion":bg.motion_summary,
                         "motion_rmse_pixels_evaluation_only":float(np.sqrt(np.mean(
                             (mixed.displacement-known)**2))),
                         "zero_map_rmse_pixels_evaluation_only":float(np.sqrt(np.mean(known**2))),
                         "methods":methods}
    outcome={"status":"development motion candidate; no held-out scores",
             "seed":report["seed"],"motion":report["motion"],
             "flow_scale":report["flow_scale"],"particle_fraction":report["scatterer_fraction"],
             "acquisition_sha256":hashlib.sha256((acquisition/"report.json").read_bytes()).hexdigest(),
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "motion_context_sha256":hashlib.sha256(Path(__file__).with_name(
                 "esdiva_phase4_motion_context.py").read_bytes()).hexdigest(),
             "training_roi_sha256":digest(roi),"input_gate":gate,
             "blood_to_tissue_roi_power_evaluator_only":contamination,
             "source_and_motion_truth_are_evaluator_only":True,
             "variants":result}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(outcome,indent=2,allow_nan=False)+"\n")
    return outcome


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--variants",nargs="+",default=("none","complex_smooth",
                                                 "envelope_rigid","envelope_smooth"))
    a=p.parse_args()
    r=evaluate(a.acquisition,a.output,a.variants)
    print(json.dumps({"seed":r["seed"],"motion":r["motion"],
                      "variants":list(r["variants"])}))
