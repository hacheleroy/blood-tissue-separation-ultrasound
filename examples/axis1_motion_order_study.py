"""Prospective Axis 1 simulation study; freeze source/config before heldout.

Run development first, inspect it, then --freeze. Heldout requires an exact
source/config lease and disjoint seeds. No source truth reaches any fitting API.
"""
import argparse
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from btbench.filters import svd_filter
from btbench.motion_first import estimate_bmode_motion, registration_order_filter, motion_windows
from btbench.shared import SharedContext, fit_shared, shared_scene


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/axis1_motion_order_v1.json"
SOURCE = ROOT / "src/btbench/motion_first.py"
DEFAULT_OUT = ROOT / "runs/axis1_motion_order_v1"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run_case(seed, condition, config):
    a = shared_scene(seed=seed, n_frames=config["n_frames"],
        calibration_frames=config["calibration_frames"], shape=tuple(config["shape_zx"]),
        flow_mps=config["flow_mps"], decorrelation=config["decorrelation"],
        motion=condition, amplitude_pixels=config["amplitude_pixels"])
    mixture = a.public.views
    background = a.tissue+a.noise
    s = config["score_start"]
    crop = config["common_crop_pixels"]
    region = (slice(crop, -crop), slice(crop, -crop))
    b_ref = a.blood[s:].mean(axis=1)[:, region[0], region[1]]
    t_ref = a.tissue[s:].mean(axis=1)[:, region[0], region[1]]
    t_den = max(np.sum(abs(t_ref)**2), 1e-30)
    b_den = max(np.sum(abs(b_ref)**2), 1e-30)
    inputs = {"mixture":mixture,"background":background}
    fit_time = {}
    fits = {}
    for model in ("translation", "auto_affine", "auto_smooth"):
        tic = perf_counter()
        fits[model] = {name:estimate_bmode_motion(x, model=model,
            max_shift_px=config["max_shift_px"]) for name,x in inputs.items()}
        fit_time[model] = perf_counter()-tic
    legacy = {}
    tic = perf_counter()
    for name, x in inputs.items():
        public = replace(a.public, views=x)
        legacy[name] = SharedContext(public, packet=16, score_start=s, motion="smooth")
    fit_time["existing_axis1"] = perf_counter()-tic
    dz = .003/(config["shape_zx"][0]-1)
    dx = .003/(config["shape_zx"][1]-1)
    def get_result(method, name):
        x = inputs[name]
        discard = 2 if method.endswith("discard2") else config["fixed_svd_discard"]
        if method.startswith("raw_svd"):
            return svd_filter(x, discard=discard)
        if method == "existing_axis1":
            return fit_shared(legacy[name], config["legacy_axis1"])
        if method in ("fixed_window_svd", "adaptive_window_svd",
                      "fixed_window_maresca", "adaptive_window_maresca"):
            selector = "maresca" if "maresca" in method else "fixed"
            return motion_windows(x, fits["translation"][name], a.public.geometry["prf_hz"],
                a.public.geometry["f0_hz"], dz, dx,
                window=config["window_max"], min_window=config["window_min"],
                hop=config["window_hop"], adaptive=method.startswith("adaptive"),
                selector=selector, discard=config["fixed_svd_discard"],
                offset=config["maresca_offset"])
        parts = method.split("_")
        order = "before" if method.startswith("filter_before_") else "after"
        selector = "baranger" if "baranger" in method else "fixed"
        model = "auto_smooth" if "smooth" in method else ("auto_affine" if "affine" in method else "translation")
        phase = "phase" in parts
        return registration_order_filter(x, fits[model][name], order=order,
            selector=selector, discard=discard, phase=phase)
    methods = ["raw_svd", "filter_before_svd", "registration_before_svd",
        "registration_before_svd_phase", "registration_before_affine_svd_phase",
        "registration_before_smooth_svd_phase", "filter_before_baranger",
        "registration_before_baranger_phase", "existing_axis1",
        "fixed_window_svd", "adaptive_window_svd",
        "fixed_window_maresca", "adaptive_window_maresca",
        "raw_svd_discard2", "filter_before_svd_discard2",
        "registration_before_svd_phase_discard2"]
    scored = []
    for method in methods:
        tic = perf_counter()
        mix, bg = get_result(method,"mixture"), get_result(method,"background")
        def scored_iq(value):
            return value[s:].mean(axis=1)[:, region[0], region[1]]
        # The legacy fit returns only the scored block and includes 3 views.
        if method == "existing_axis1":
            def legacy_iq(value):return value.mean(axis=1)[:, region[0], region[1]]
            mixed = legacy_iq(mix.views);base = legacy_iq(bg.views)
            frozen_b = legacy_iq(mix.apply_frozen(a.blood))
            frozen_t = legacy_iq(mix.apply_frozen(a.tissue))
        else:
            mixed = scored_iq(mix.iq);base = scored_iq(bg.iq)
            frozen_b = scored_iq(mix.apply_frozen(a.blood))
            frozen_t = scored_iq(mix.apply_frozen(a.tissue))
        added = mixed-base
        closure = float(np.linalg.norm(mixed-scored_iq(mix.apply_frozen(mixture)))) if method != "existing_axis1" else float(np.linalg.norm(mixed-legacy_iq(mix.apply_frozen(mixture))))
        if method == "existing_axis1":
            fit_label = "existing_axis1"
            registration = legacy["mixture"].motion_summary
            failures = registration["frames"]-registration["accepted_frames"]
            phase_coherence = None
            roundtrip = None
        else:
            fit_label = ("auto_smooth" if "smooth" in method else
                         "auto_affine" if "affine" in method else "translation")
            fit = fits[fit_label]["mixture"]
            registration = fit.diagnostics
            failures = registration["failed_frames"]
            phase_coherence = registration["phase_coherence_min"]
            roundtrip = float(np.linalg.norm(scored_iq(fit.restore(fit.align(mixture)))-
                                             scored_iq(mixture))/max(np.linalg.norm(scored_iq(mixture)),1e-30))
        windows = mix.diagnostics.get("windows", []) if method != "existing_axis1" else []
        scored.append({"method":method,
            "added_blood_iq_error":float(np.linalg.norm(added-b_ref)/np.sqrt(b_den)),
            "added_blood_power_retention":float(np.sum(abs(added)**2)/b_den),
            "frozen_blood_iq_error":float(np.linalg.norm(frozen_b-b_ref)/np.sqrt(b_den)),
            "frozen_tissue_power_retention":float(np.sum(abs(frozen_t)**2)/t_den),
            "refitted_background_power_remaining":float(np.sum(abs(base)**2)/t_den),
            "frozen_closure_absolute":closure,"roundtrip_iq_error":roundtrip,
            "registration_failure_frames":failures,
            "registration_failure_definition":("not accepted by complex-reference gate, including identity"
                if method == "existing_axis1" else "envelope held-out NCC below 0.5"),
            "minimum_phase_coherence":phase_coherence,
            "registered_model_counts":({k:sum(r["model"]==k for r in registration["frames"])
                for k in ("identity","translation","affine","smooth","identity_fallback")}
                if method != "existing_axis1" else None),
            "window_lengths":([w["stop"]-w["start"] for w in windows]
                 if method.startswith(("fixed_window", "adaptive_window")) else None),
            "window_discards":([w["svd"]["discard"] for w in windows]
                 if method.startswith(("fixed_window", "adaptive_window")) else None),
            "method_runtime_s":perf_counter()-tic,
            "fit_runtime_s":fit_time[fit_label],
            "parameters":({"selector":"maresca","fixed_offset":config["maresca_offset"]}
                if "maresca" in method else {"discard":2 if method.endswith("discard2")
                    else config["fixed_svd_discard"]} if "baranger" not in method else None)})
    return {"seed":seed,"condition":condition,"blood_input_power":float(b_den),
        "tissue_input_power":float(t_den),"results":scored,
        "registration_check":{k:{"failures":v["mixture"].diagnostics["failed_frames"],
            "models":[r["model"] for r in v["mixture"].diagnostics["frames"]]}
            for k,v in fits.items()},
        "legacy_registration":legacy["mixture"].motion_summary}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--split", choices=["development","heldout"], required=True)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--freeze", action="store_true", help="seal development and source lease")
    args = p.parse_args()
    if args.freeze and args.split != "development":p.error("Freeze only after development")
    c = json.loads(CONFIG.read_text())
    assert not set(c["development_seeds"]) & set(c["heldout_seeds"])
    args.out.mkdir(parents=True, exist_ok=True)
    freeze = args.out / "freeze.json"
    dependencies = ["src/btbench/core.py", "src/btbench/filters.py",
        "src/btbench/shared.py", "src/btbench/registration.py",
        "src/btbench/svd_literature.py", "src/btbench/paper_filters.py",
        "src/btbench/stap.py", "src/btbench/axes.py"]
    sources = {str(path.relative_to(ROOT)):digest(path) for path in
               (CONFIG,SOURCE,Path(__file__),*(ROOT / name for name in dependencies))}
    if args.split == "heldout":
        if not freeze.exists():p.error("Development freeze missing")
        manifest = json.loads(freeze.read_text())
        if sources != manifest["sources"] or digest(args.out / "development.json") != manifest["development_sha256"]:
            p.error("Source/config/development mismatch; heldout lease refused")
    seeds = c[f"{args.split}_seeds"]
    tic = perf_counter()
    rows = [run_case(seed, condition, c) for seed in seeds for condition in c["conditions"]]
    result = {"study":c["study"],"split":args.split,"source_sha256":sources,
        "score_region":"last 32 frames, common 4-pixel spatial crop; coherent mean of all 3 views",
        "known_truth_use":"scoring only, never motion or filter fit",
        "motion_window_status":"Zhang-inspired automatic ROI; fixed Maresca offset, not published operator calibration",
        "elapsed_s":perf_counter()-tic,"cases":rows}
    target = args.out / f"{args.split}.json"
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    if args.freeze:
        freeze.write_text(json.dumps({"sources":sources,"development_sha256":digest(target),
            "heldout_seeds":c["heldout_seeds"],"locked_before_heldout":True},indent=2)+"\n")
    print(json.dumps({"split":args.split,"cases":len(rows),"elapsed_s":result["elapsed_s"],
        "archive":str(target),"sha256":digest(target),"freeze":freeze.exists()}))


if __name__ == "__main__":main()
