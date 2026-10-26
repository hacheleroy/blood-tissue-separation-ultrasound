"""Seen diagnosis, fresh development, then source-frozen prospective simulation."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from axis1_native_subtraction_study import scene_for
from btbench.filters import svd_filter
from btbench.motion_first import estimate_bmode_motion, registration_order_filter
from btbench.native_subtraction import native_tissue_subtraction
from btbench.stable_tissue import bright_reference_subtraction, smooth_motion_path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/axis1_stable_tissue_v1.json"
OUT = ROOT / "runs/axis1_stable_tissue_v1"
SOURCE_PATHS = [CONFIG, Path(__file__),
    ROOT / "examples/axis1_native_subtraction_study.py",
    ROOT / "src/btbench/stable_tissue.py", ROOT / "src/btbench/native_subtraction.py",
    ROOT / "src/btbench/motion_first.py", ROOT / "src/btbench/filters.py",
    ROOT / "src/btbench/core.py", ROOT / "src/btbench/shared.py",
    ROOT / "src/btbench/registration.py", ROOT / "src/btbench/svd_literature.py",
    ROOT / "src/btbench/paper_filters.py", ROOT / "src/btbench/stap.py",
    ROOT / "src/btbench/axes.py"]
METHODS = ("raw_svd", "registered_return", "native_subtraction",
           "smooth_full_field", "smooth_bright_quartile")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run_case(seed, condition, config, split):
    scene, params = scene_for(seed, condition, config, "development")
    mixture = scene.public.views
    background = scene.tissue + scene.noise
    inputs = {"mixture":mixture,"background":background}
    crop, start = config["common_crop_pixels"], config["score_start"]
    def score_domain(y):
        return y[start:].mean(axis=1)[:,crop:-crop,crop:-crop]
    truth_b,truth_t = score_domain(scene.blood),score_domain(scene.tissue)
    norm = np.linalg.norm
    power = lambda y:float(np.sum(abs(y)**2))
    bden,tden = max(power(truth_b),1e-30),max(power(truth_t),1e-30)
    tic = perf_counter()
    base_motion = {key:estimate_bmode_motion(x,model="translation",
        max_shift_px=config["max_shift_px"]) for key,x in inputs.items()}
    stable_motion = {key:smooth_motion_path(m,
        sigma_frames=config["motion_smoothing_sigma_frames"])
        for key,m in base_motion.items()}
    motion_seconds = perf_counter()-tic

    def fit(name, method, motion_key=None):
        x = inputs[name]
        mk = motion_key or name
        if method == "raw_svd":
            return svd_filter(x,discard=config["svd_discard"])
        if method == "registered_return":
            return registration_order_filter(x,base_motion[mk],order="after",
                selector="fixed",discard=config["svd_discard"],
                phase=config["phase_alignment"])
        if method == "native_subtraction":
            return native_tissue_subtraction(x,base_motion[mk],
                discard=config["svd_discard"],phase=config["phase_alignment"])
        fraction = 1.0 if method == "smooth_full_field" else config["bright_training_fraction"]
        if method not in ("smooth_full_field","smooth_bright_quartile"):
            raise ValueError(method)
        return bright_reference_subtraction(x,stable_motion[mk],
            fraction=fraction,phase=config["phase_alignment"])

    results={}
    for method in METHODS:
        tic = perf_counter()
        mix,bg = fit("mixture",method),fit("background",method)
        output,base = score_domain(mix.iq),score_domain(bg.iq)
        fb=score_domain(mix.apply_frozen(scene.blood))
        ft=score_domain(mix.apply_frozen(scene.tissue))
        fn=score_domain(mix.apply_frozen(scene.noise))
        shared_motion_bg = score_domain(fit("background",method,"mixture").iq)
        results[method]={
            "added_blood_iq_error":float(norm(output-base-truth_b)/np.sqrt(bden)),
            "added_blood_power_retention":power(output-base)/bden,
            "frozen_blood_iq_error":float(norm(fb-truth_b)/np.sqrt(bden)),
            "frozen_tissue_power_retention":power(ft)/tden,
            "refitted_background_power_remaining":power(base)/tden,
            "same_motion_background_power_remaining":power(shared_motion_bg)/tden,
            "same_motion_added_blood_iq_error":float(norm(output-shared_motion_bg-truth_b)/np.sqrt(bden)),
            "motion_refit_background_change_iq_norm":float(norm(shared_motion_bg-base)/np.sqrt(bden)),
            "frozen_source_closure_nrmse":float(norm(output-fb-ft-fn)/max(norm(score_domain(mixture)),1e-30)),
            "frozen_mixture_replay_nrmse":float(norm(output-score_domain(mix.apply_frozen(mixture)))/max(norm(score_domain(mixture)),1e-30)),
            "elapsed_s":perf_counter()-tic}
    diff=base_motion["mixture"].forward-base_motion["background"].forward
    stable_diff=stable_motion["mixture"].forward-stable_motion["background"].forward
    return {"seed":seed,"condition":condition,"scene":params,
            "blood_input_power":bden,"tissue_input_power":tden,
            "results":results,"mechanism":{
                "motion_fit_seconds":motion_seconds,
                "base_motion_refit_rms_px":float(np.sqrt(np.mean(diff*diff))),
                "stable_motion_refit_rms_px":float(np.sqrt(np.mean(stable_diff*stable_diff))),
                "registration_failure_frames":base_motion["mixture"].diagnostics["failed_frames"],
                "background_registration_failure_frames":base_motion["background"].diagnostics["failed_frames"]}}


def decide(cases, config):
    d=config["decision"]
    moving=[c for c in cases if c["scene"]["motion"]!="none"]
    stationary=[c for c in cases if c["scene"]["motion"]=="none"]
    rows=[]
    for c in stationary+moving:
        active=c in moving
        blood=c["results"][config["primary_method"]]["added_blood_iq_error"]
        raw=c["results"]["raw_svd"]["added_blood_iq_error"]
        tissue=c["results"][config["primary_method"]]["refitted_background_power_remaining"]
        rows.append({"seed":c["seed"],"condition":c["condition"],
            "blood_noninferior_to_raw":bool(blood<=raw+d["moving_added_blood_error_margin_vs_raw" if active else "stationary_added_blood_error_margin_vs_raw"]),
            "tissue_below_limit":bool(not active or tissue<=d["moving_background_power_max"]),
            "registration_failure_limit":bool(c["mechanism"]["registration_failure_frames"]<=d["maximum_registration_failures_per_case"])})
    eligible=(len(moving)==d["required_moving_cases"] and len(stationary)==d["required_stationary_cases"]
              and all(all(row[key] for key in ("blood_noninferior_to_raw","tissue_below_limit","registration_failure_limit")) for row in rows))
    return {"primary_method":config["primary_method"],
       "joint_criterion_passed":bool(eligible),"case_gates":rows,
       "criterion":"all moving: blood error <= raw SVD+0.05, background <=0.10, zero registration failures; stationary blood <= raw+0.05"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--split",choices=("seen_diagnostic","development","heldout"),required=True)
    p.add_argument("--freeze",action="store_true")
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    if args.freeze and args.split!="development":p.error("Only development can be frozen")
    config=json.loads(CONFIG.read_text())
    all_seeds=[set(config[f"{name}_seeds"]) for name in ("seen_diagnostic","development","heldout")]
    if any(all_seeds[i]&all_seeds[j] for i in range(3) for j in range(i+1,3)):
        p.error("Seed overlap")
    sources={str(path.relative_to(ROOT)):digest(path) for path in SOURCE_PATHS}
    args.out.mkdir(parents=True,exist_ok=True)
    target=args.out/f"{args.split}.json"
    if target.exists():p.error("Split target already scored; choose a new output directory")
    if args.split=="heldout":
        freeze=json.loads((args.out/"freeze.json").read_text())
        if sources!=freeze["sources"] or freeze["heldout_seeds"]!=config["heldout_seeds"] or any(
            digest(args.out/f"{name}.json")!=freeze[f"{name}_sha256"]
            for name in ("seen_diagnostic","development")):
            p.error("Held-out lease invalid")
    elif args.split=="development" and not (args.out/"seen_diagnostic.json").exists():
        p.error("Run seen-bank diagnosis first")
    tic=perf_counter()
    cases=[run_case(seed,condition,config,args.split)
           for seed in config[f"{args.split}_seeds"] for condition in config["conditions"]]
    result={"study":config["study"],"split":args.split,"source_sha256":sources,
        "cases":cases,"decision":None if args.split=="seen_diagnostic" else decide(cases,config),
        "elapsed_s":perf_counter()-tic,
        "truth_use":"source truth is evaluator-only; filter sees observed mixture IQ",
        "physical_status":"source-labeled simulation; no physical measurements"}
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
