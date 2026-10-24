"""Development paired source score for phase-only tissue tracking and slow flow."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, relative_error
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase3_input_gate import choose_input
from esdiva_phase4_phase_context import PhaseOnlyContext


METHODS={
    "roundtrip":{"method":"roundtrip","horizon":64},
    "svd_cut1":{"method":"svd","horizon":64,"discard":1},
    "axis1_cut1":{"method":"axis1","horizon":64,"discard":1,"rho":.75},
    "axis2_default":{"method":"axis2","horizon":64,"rank":4,"gain_cap":1.},
    "axis2_r4_cap045":{"method":"axis2","horizon":64,"rank":4,"gain_cap":.45},
    "axis2_r6_cap045":{"method":"axis2","horizon":64,"rank":6,"gain_cap":.45},
}


def evaluate(acquisition, output, variants=("none","phase_gate","phase_no_gate")):
    acquisition,output=Path(acquisition).resolve(),Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    report=json.loads((acquisition/"report.json").read_text())
    with np.load(acquisition/"iq_sources.npz") as f:
        mixture,cal,tissue,blood,noise,x=(f[k] for k in
            ("views","calibration_views","tissue","blood","noise","x_mm"))
    for key,value in (("mixture_sha256",mixture),("calibration_sha256",cal),
                      ("tissue_sha256",tissue),("blood_sha256",blood),("noise_sha256",noise)):
        if digest(value)!=report[key]:
            raise ValueError(f"Acquisition hash mismatch: {key}")
    roi,gate=choose_input(cal,x,min_roi_pixels=40)
    if not gate["input_gate_pass"] or relative_error(tissue+blood+noise,mixture)>1e-8:
        raise ValueError("Source or input gate failed")
    ratio=lambda a,b:float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))
    if ratio(blood[:,:,roi],tissue[:,:,roi])>.1:
        raise ValueError("Blood contamination in training ROI")
    public=PublicInput(mixture,cal,roi,report["geometry"])
    background=PublicInput(tissue+noise,cal,roi,report["geometry"])
    rows={}
    for name in variants:
        if name=="none":
            a,b=(Phase3Context(p,motion="none") for p in (public,background))
        elif name=="phase_gate":
            a,b=(PhaseOnlyContext(p) for p in (public,background))
        elif name=="phase_no_gate":
            settings={"minimum_phase_rad":0.,"maximum_phase_disagreement_rad":np.pi,
                      "minimum_validation_coherence":0.,
                      "minimum_validation_nrmse_improvement":-10.}
            a,b=(PhaseOnlyContext(p,phase_options=settings) for p in (public,background))
        else:
            raise ValueError(f"Unknown candidate {name}")
        if a.calibration["calibration_sha256"]!=b.calibration["calibration_sha256"]:
            raise AssertionError("Counterfactual calibration changed")
        methods={}
        for method,spec in METHODS.items():
            fit,control=fit_shared(a,spec),fit_shared(b,spec)
            methods[method]={"spec":spec,
                             "added_blood_iq_nrmse":relative_error(fit.views-control.views,blood[64:]),
                             "background_to_input_tissue_power":ratio(control.views,tissue[64:]),
                             "frozen_blood_power_ratio":ratio(fit.apply_frozen(blood),blood[64:])}
        rows[name]={"motion_fit":a.motion_summary,
                    "background_motion_fit":b.motion_summary,"methods":methods}
    result={"status":"development paired RF tissue phase and slow-flow control",
            "seed":report["seed"],"motion":report["motion"],"flow_scale":report["flow_scale"],
            "scatterer_fraction":report["scatterer_fraction"],
            "acquisition_sha256":hashlib.sha256((acquisition/"report.json").read_bytes()).hexdigest(),
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "phase_context_sha256":hashlib.sha256(Path(__file__).with_name(
                "esdiva_phase4_phase_context.py").read_bytes()).hexdigest(),
            "input_gate":gate,"selected_roi_sha256":digest(roi),
            "source_truth_evaluator_only":True,"scored_frames":16,"variants":rows}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquisition",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--variants",nargs="+",default=("none","phase_gate","phase_no_gate"))
    a=p.parse_args()
    r=evaluate(a.acquisition,a.output,a.variants)
    print(json.dumps({"seed":r["seed"],"motion":r["motion"],"variants":list(r["variants"])}))
