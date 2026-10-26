"""Seen-bank mechanism diagnosis followed by frozen development/heldout study."""
import argparse
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import ndimage

from btbench.filters import svd_filter
from btbench.motion_first import estimate_bmode_motion, registration_order_filter
from btbench.native_subtraction import native_tissue_subtraction
from btbench.shared import shared_scene


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/axis1_native_subtraction_v1.json"
OUT = ROOT / "runs/axis1_native_subtraction_v1"
SOURCE_PATHS = [CONFIG, Path(__file__), ROOT / "src/btbench/native_subtraction.py",
    ROOT / "src/btbench/motion_first.py", ROOT / "src/btbench/filters.py",
    ROOT / "src/btbench/core.py", ROOT / "src/btbench/shared.py",
    ROOT / "src/btbench/registration.py", ROOT / "src/btbench/svd_literature.py",
    ROOT / "src/btbench/paper_filters.py", ROOT / "src/btbench/stap.py",
    ROOT / "src/btbench/axes.py"]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def turnover_tissue(scene, seed, fraction):
    """Evaluation-only source construction, not an elevational acoustic model."""
    if fraction == 0:
        return scene
    rng = np.random.default_rng(seed + 876432)
    n, views, nz, nx = scene.tissue.shape
    original = scene.tissue
    amplitude = np.sqrt(np.mean(abs(original)**2, axis=(0, 2, 3)))
    innovations = rng.normal(size=(n, nz, nx))+1j*rng.normal(size=(n, nz, nx))
    innovations = ndimage.gaussian_filter(innovations.real, (0, .65, .65)) + 1j * ndimage.gaussian_filter(innovations.imag, (0, .65, .65))
    innovations /= max(np.sqrt(np.mean(abs(innovations)**2)), 1e-30)
    state = np.zeros((nz, nx), complex)
    fields = np.empty_like(innovations)
    rho = .96
    for t in range(n):
        state = rho*state + np.sqrt(1-rho*rho)*innovations[t]
        fields[t] = state
    change = fraction*fields[:, None]*amplitude[None, :, None, None]
    tissue = original+change
    scene.tissue = tissue
    scene.public = replace(scene.public, views=scene.blood+tissue+scene.noise)
    scene.metadata["tissue_speckle_turnover_surrogate"] = {"fraction":fraction,"rho":rho,
        "description":"correlated complex texture innovation added to tissue; not 3-D acoustics"}
    return scene


def scene_for(seed, condition, config, split):
    if split == "seen_diagnostic":
        params = {"motion":condition, "flow_mps":.0012, "amplitude_pixels":1.1,
                  "decorrelation":.99, "tissue_turnover":0.}
    else:
        params = config["conditions"][condition]
    scene = shared_scene(seed=seed, n_frames=config["n_frames"],
        calibration_frames=config["calibration_frames"],
        shape=tuple(config["shape_zx"]), motion=params["motion"],
        flow_mps=params["flow_mps"], amplitude_pixels=params["amplitude_pixels"],
        decorrelation=params["decorrelation"])
    return turnover_tissue(scene, seed, params["tissue_turnover"]), params


def run_case(seed, condition, config, split):
    scene, params = scene_for(seed, condition, config, split)
    mixture = scene.public.views
    inputs = {"mixture":mixture, "background":scene.tissue+scene.noise}
    start = config["score_start"]
    crop = config["common_crop_pixels"]
    region = (slice(crop, -crop), slice(crop, -crop))
    def score_domain(x):
        return x[start:].mean(axis=1)[:, region[0], region[1]]
    truth_b = score_domain(scene.blood)
    truth_t = score_domain(scene.tissue)
    power = lambda x:float(np.sum(abs(x)**2))
    norm = lambda x:float(np.linalg.norm(x))
    bden = max(power(truth_b), 1e-30)
    tden = max(power(truth_t), 1e-30)
    tic = perf_counter()
    motion = {name:estimate_bmode_motion(x,model="translation",
        max_shift_px=config["max_shift_px"]) for name,x in inputs.items()}
    fit_time = perf_counter()-tic
    def fit(name, method):
        x = inputs[name]
        if method == "raw_svd":
            return svd_filter(x,discard=config["svd_discard"])
        if method == "registered_return":
            return registration_order_filter(x,motion[name],order="after",selector="fixed",
                discard=config["svd_discard"],phase=config["phase_alignment"])
        if method == "native_subtraction":
            return native_tissue_subtraction(x,motion[name],
                discard=config["svd_discard"],phase=config["phase_alignment"])
        raise ValueError(method)
    results = {}
    fits = {}
    for method in ("raw_svd","registered_return","native_subtraction"):
        tic = perf_counter()
        mix,bg=fit("mixture",method),fit("background",method)
        fits[method]=mix
        output=score_domain(mix.iq);base=score_domain(bg.iq)
        added=output-base
        fb=score_domain(mix.apply_frozen(scene.blood))
        ft=score_domain(mix.apply_frozen(scene.tissue))
        fn=score_domain(mix.apply_frozen(scene.noise))
        results[method]={
            "added_blood_iq_error":norm(added-truth_b)/np.sqrt(bden),
            "added_blood_power_retention":power(added)/bden,
            "frozen_blood_iq_error":norm(fb-truth_b)/np.sqrt(bden),
            "frozen_tissue_power_retention":power(ft)/tden,
            "refitted_background_power_remaining":power(base)/tden,
            "frozen_source_closure_nrmse":norm(output-fb-ft-fn)/max(norm(score_domain(mixture)),1e-30),
            "frozen_mixture_replay_nrmse":norm(output-score_domain(mix.apply_frozen(mixture)))/max(norm(score_domain(mixture)),1e-30),
            "elapsed_s":perf_counter()-tic}
    aligned = motion["mixture"].align(mixture,phase=config["phase_alignment"])
    roundtrip = motion["mixture"].restore(aligned,phase=config["phase_alignment"])
    bround = motion["mixture"].restore(motion["mixture"].align(scene.blood,
        phase=config["phase_alignment"]),phase=config["phase_alignment"])
    native = fits["native_subtraction"]
    registered = fits["registered_return"]
    identity_error = norm(score_domain(native.iq-registered.iq)-score_domain(mixture-roundtrip))/max(norm(score_domain(mixture)),1e-30)
    if identity_error > 1e-10:
        raise ValueError("Native-grid identity failed")
    source_identity_error = norm(score_domain(native.apply_frozen(scene.blood)-registered.apply_frozen(scene.blood))-score_domain(scene.blood-bround))/np.sqrt(bden)
    if source_identity_error > 1e-10:
        raise ValueError("Frozen blood identity failed")
    return {"seed":seed,"condition":condition,"scene":params,
        "blood_input_power":bden,"tissue_input_power":tden,
        "results":results,
        "mechanism":{
            "mixture_roundtrip_iq_error":norm(score_domain(mixture-roundtrip))/max(norm(score_domain(mixture)),1e-30),
            "blood_roundtrip_iq_error":norm(score_domain(scene.blood-bround))/np.sqrt(bden),
            "blood_removed_as_tissue_iq_error":norm(score_domain(scene.blood-native.apply_frozen(scene.blood)))/np.sqrt(bden),
            "output_identity_nrmse":identity_error,"blood_identity_nrmse":source_identity_error,
            "motion_fit_seconds":fit_time,
            "registration_failure_frames":motion["mixture"].diagnostics["failed_frames"],
            "background_registration_failure_frames":motion["background"].diagnostics["failed_frames"],
            "model_counts":{model:sum(r["model"]==model for r in motion["mixture"].diagnostics["frames"])
                for model in ("identity","translation","identity_fallback")},
            "minimum_phase_coherence":motion["mixture"].diagnostics["phase_coherence_min"]}}


def decide(cases, config):
    d=config["decision"]
    moving=[c for c in cases if c["scene"]["motion"]!="none"]
    stationary=[c for c in cases if c["scene"]["motion"]=="none"]
    def gate(c, stationary_case=False):
        r=c["results"]
        blood=r["native_subtraction"]["added_blood_iq_error"]
        raw=r["raw_svd"]["added_blood_iq_error"]
        tissue=r["native_subtraction"]["refitted_background_power_remaining"]
        return {"seed":c["seed"],"condition":c["condition"],
            "blood_noninferior_to_raw":bool(blood <= raw + d["stationary_added_blood_error_margin_vs_raw" if stationary_case else "moving_added_blood_error_margin_vs_raw"]),
            "tissue_below_limit":bool(stationary_case or tissue <= d["moving_background_power_max"]),
            "registration_failure_limit":bool(c["mechanism"]["registration_failure_frames"]<=d["maximum_registration_failures_per_case"]),
            "native_better_than_registered_return":bool(blood<r["registered_return"]["added_blood_iq_error"])}
    rows=[gate(c,True) for c in stationary]+[gate(c) for c in moving]
    eligible=(len(moving)==d["required_moving_cases"] and len(stationary)==d["required_stationary_cases"]
        and all(row["blood_noninferior_to_raw"] and row["tissue_below_limit"] and row["registration_failure_limit"] for row in rows))
    return {"joint_criterion_passed":bool(eligible),"case_gates":rows,
        "criterion":"every moving case: added-blood error <= raw SVD +0.05; background <=0.10; zero registration failures; stationary blood <= raw+0.05"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--split",choices=["seen_diagnostic","development","heldout"],required=True)
    p.add_argument("--freeze",action="store_true")
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    if args.freeze and args.split!="development":p.error("Only development can be frozen")
    config=json.loads(CONFIG.read_text())
    assert not set(config["heldout_seeds"]) & (set(config["development_seeds"]) | set(config["seen_diagnostic_seeds"]))
    sources={str(path.relative_to(ROOT)):digest(path) for path in SOURCE_PATHS}
    args.out.mkdir(parents=True,exist_ok=True)
    target=args.out/f"{args.split}.json"
    if args.split=="heldout":
        freeze=json.loads((args.out/"freeze.json").read_text())
        if target.exists() or sources!=freeze["sources"] or any(
            digest(args.out/f"{name}.json")!=freeze[f"{name}_sha256"]
            for name in ("development","seen_diagnostic")):
            p.error("Held-out lease invalid or target already scored")
    if args.split=="development" and not (args.out/"seen_diagnostic.json").exists():
        p.error("Run seen-bank diagnosis first")
    if args.split=="seen_diagnostic":
        seeds=config["seen_diagnostic_seeds"];conditions=config["seen_conditions"]
    else:
        seeds=config[f"{args.split}_seeds"];conditions=list(config["conditions"])
    tic=perf_counter()
    cases=[run_case(seed,cond,config,args.split) for seed in seeds for cond in conditions]
    result={"study":config["study"],"split":args.split,"source_sha256":sources,
        "cases":cases,"decision":decide(cases,config) if args.split!="seen_diagnostic" else None,
        "elapsed_s":perf_counter()-tic,
        "truth_use":"source truth scores outputs after fitting; fit APIs see IQ only",
        "physical_status":"simulation only; turnover is specified speckle surrogate"}
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if args.freeze:
        (args.out/"freeze.json").write_text(json.dumps({"sources":sources,
            "seen_diagnostic_sha256":digest(args.out/"seen_diagnostic.json"),
            "development_sha256":digest(target),"heldout_seeds":config["heldout_seeds"],
            "locked_before_heldout":True},indent=2)+"\n")
    print(json.dumps({"split":args.split,"cases":len(cases),"elapsed_s":result["elapsed_s"],
        "decision":result["decision"]["joint_criterion_passed"] if result["decision"] else None,
        "sha256":digest(target)}))


if __name__=="__main__":main()
