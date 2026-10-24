"""Replay one eSDIva development IQ acquisition with paired source controls.

The mixture and separately acquired tissue-plus-noise condition use identical
calibration, timestamps, ROI, packet, history and receiver-noise realization.
Their refitted output difference measures the response to adding blood while
allowing each method to re-estimate its operator. Source IQ remains evaluator
data; the PublicInput of every fit receives observations only.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from btbench.shared import PublicInput, SharedContext, fit_shared
from esdiva_example22_pilot import digest, numerical_source_hash, relative_error
from esdiva_gated_stap import fit_gated_loaded_stap


SPECS = {
    "identity": {"method": "identity", "horizon": 16},
    "svd_cut1": {"method": "svd", "horizon": 16, "discard": 1},
    "baranger_ssm": {"method": "baranger", "horizon": 16},
    "axis1_cut1": {"method": "axis1", "horizon": 16, "discard": 1, "rho": .75},
    "axis2": {"method": "axis2", "horizon": 16, "rank": 4, "gain_cap": 1.},
    "axis3_orthogonal": {"method": "stap", "horizon": 16, "adaptive": False},
    "axis3_loaded": {"method": "stap", "horizon": 16, "adaptive": True, "loading_rel": .01},
    "axis3_loaded_gated": {"method": "stap", "horizon": 16, "adaptive": True, "loading_rel": .01},
    "lee_aer": {"method": "lee_aer", "horizon": 16, "domain": "temporal", "epsilon_abs": 26.},
}


def run(pilot, output):
    pilot, output = Path(pilot).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to replace {output}")
    report = json.loads((pilot / "report.json").read_text())
    if not report.get("calibration_gate_pass"):
        raise ValueError("Calibration gate failed in source acquisition")
    arrays = np.load(pilot / "pilot_iq.npz")
    mixture, cal = arrays["views"], arrays["calibration_views"]
    blood, tissue, noise = arrays["blood"], arrays["tissue"], arrays["noise"]
    roi = arrays["training_roi"]
    if (digest(mixture) != report["mixture_sha256"] or
            digest(cal) != report["calibration_sha256"] or
            not np.allclose(mixture, blood+tissue+noise, atol=1e-7)):
        raise ValueError("Pilot IQ identity or additive sources failed")
    geometry = report["geometry"]
    public = PublicInput(mixture, cal, roi, geometry)
    background = PublicInput(tissue+noise, cal, roi, geometry)
    mixed = SharedContext(public, packet=16, score_start=16, motion="none")
    bg = SharedContext(background, packet=16, score_start=16, motion="none")
    assert mixed.calibration == bg.calibration
    score = slice(16, None)
    results = {}
    for name, spec in SPECS.items():
        start = time.monotonic()
        try:
            if name == "axis3_loaded_gated":
                fitted, gate = fit_gated_loaded_stap(mixed, spec)
                control, bg_gate = fit_gated_loaded_stap(bg, spec)
            else:
                fitted = fit_shared(mixed, spec)
                control = fit_shared(bg, spec)
                gate = bg_gate = None
            marginal = fitted.views - control.views
            b = fitted.apply_frozen(blood)
            t = fitted.apply_frozen(tissue)
            n = fitted.apply_frozen(noise)
            control_t = control.apply_frozen(tissue)
            control_n = control.apply_frozen(noise)
            results[name] = {
                "spec": spec,
                "mixture_gate": gate, "background_gate": bg_gate,
                "blood_marginal_iq_nrmse": relative_error(marginal, blood[score]),
                "blood_marginal_power_ratio": float(np.mean(abs(marginal)**2) /
                                                    np.mean(abs(blood[score])**2)),
                "refitted_background_output_power_to_input_tissue": float(
                    np.mean(abs(control.views)**2)/np.mean(abs(tissue[score])**2)),
                "refitted_tissue_frozen_power_ratio": float(np.mean(abs(control_t)**2) /
                                                              np.mean(abs(tissue[score])**2)),
                "refitted_noise_frozen_power_ratio": float(np.mean(abs(control_n)**2) /
                                                             np.mean(abs(noise[score])**2)),
                "mixture_blood_iq_nrmse": relative_error(fitted.views, blood[score]),
                "mixture_frozen_blood_power_ratio": float(np.mean(abs(b)**2) /
                                                           np.mean(abs(blood[score])**2)),
                "mixture_frozen_tissue_power_ratio": float(np.mean(abs(t)**2) /
                                                            np.mean(abs(tissue[score])**2)),
                "mixture_frozen_noise_power_ratio": float(np.mean(abs(n)**2) /
                                                           np.mean(abs(noise[score])**2)),
                "mixture_recomposition_relative_error": relative_error(b+t+n, fitted.views),
                "background_recomposition_relative_error": relative_error(
                    control_t+control_n, control.views),
                "max_operator_gain": float(max((r.get("operator_gain_max", 1.)
                                                 for r in fitted.diagnostics["blocks"]), default=1.)),
                "elapsed_s": time.monotonic()-start,
            }
        except (ValueError, np.linalg.LinAlgError) as exc:
            results[name] = {"status": "failed", "error": str(exc), "spec": spec}
    outcome = {
        "status": "one development acquisition; no settings selection or held-out ranking",
        "pilot_report_sha256": hashlib.sha256((pilot/"report.json").read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "gate_wrapper_sha256": hashlib.sha256(
            Path(__file__).with_name("esdiva_gated_stap.py").read_bytes()).hexdigest(),
        "numerical_source_sha256": numerical_source_hash(Path(__file__).resolve().parents[1]),
        "mixture_sha256": digest(mixture), "calibration_sha256": digest(cal),
        "angle_gate_pass": True, "calibration_identical_in_counterfactual": True,
        "source_arrays_not_in_filter_api": True, "methods": results,
        "limitation": "stationary fast-flow, one scatterer seed, 16-frame history; AER regularization imported without retuning; no motion or slow-flow validation",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(outcome, indent=2, allow_nan=False)+"\n")
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.pilot, args.output)
    for name, row in result["methods"].items():
        print(name, row.get("blood_marginal_iq_nrmse", row.get("error")))
