"""Evaluation-only tissue-mode oracle; never a deployable clutter filter.

Each scene is generated once at a fixed 256-frame duration. Shorter ensembles
are nested histories ending on the same final 64 scored frames. Known motion
and isolated tissue are used only here, for an explanatory ceiling test.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from axis1_native_subtraction_study import turnover_tissue
from btbench.core import casorati, frames
from btbench.filters import svd_filter
from btbench.motion_first import MotionFit, _inverse_coordinates
from btbench.native_subtraction import native_tissue_subtraction
from btbench.shared import array_hash, shared_scene


ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/"configs/axis1_oracle_limit_v1.json"
OUT=ROOT/"runs/axis1_oracle_limit_v1"
SOURCE_PATHS=[CONFIG,Path(__file__),ROOT/"examples/axis1_native_subtraction_study.py",
    ROOT/"src/btbench/shared.py",ROOT/"src/btbench/core.py",
    ROOT/"src/btbench/filters.py",ROOT/"src/btbench/motion_first.py",
    ROOT/"src/btbench/native_subtraction.py",ROOT/"src/btbench/registration.py",
    ROOT/"src/btbench/svd_literature.py",ROOT/"src/btbench/paper_filters.py",
    ROOT/"src/btbench/stap.py",ROOT/"src/btbench/axes.py"]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def oracle_motion(displacement, views):
    """Known geometric maps; explicitly restricted to this evaluator."""
    grid=np.indices(displacement.shape[2:],dtype=float)
    forward=grid[None]+displacement
    inverse=np.stack([_inverse_coordinates(q) for q in forward])
    return MotionFit(forward,inverse,np.zeros((len(displacement),views)),
        {"evaluation_only":"imposed simulator displacement; unavailable to filter"},False)


def minimum_gain_for_tissue_limit(background, removed, tissue_power, limit):
    """Exact smallest gain in [0,1] meeting tissue power for a fixed projector.

    For alpha in [0,1], C(alpha) = ||background-alpha*removed||^2/tissue_power.
    Because blood error is alpha*||removed_blood||/||blood||, the smallest
    feasible alpha also gives the best blood preservation within this ray.
    This says nothing about projectors in other directions or nonlinear rules.
    """
    if not 0 < limit < 1 or tissue_power <= 0:
        raise ValueError("Need positive tissue power and limit in (0,1)")
    a=float(np.vdot(removed,removed).real)
    b=float(np.vdot(removed,background).real)
    c=float(np.vdot(background,background).real-limit*tissue_power)
    if c<=0:return 0.0
    if a<=0:return None
    disc=b*b-a*c
    if disc<0:return None
    lower=(b-np.sqrt(disc))/a
    upper=(b+np.sqrt(disc))/a
    gain=max(0.0,float(lower))
    if gain>min(1.0,float(upper))+1e-12:return None
    if np.linalg.norm(background-gain*removed)**2 > limit*tissue_power*(1+1e-10):
        raise ValueError("Quadratic tissue-gain solution failed")
    return gain


def fit_tissue_modes(aligned_tissue, max_rank):
    """One source-truth SVD, shared by both oracle ranks."""
    _,singular,vh=np.linalg.svd(casorati(aligned_tissue),full_matrices=False)
    modes=vh[:max_rank].conj().T
    return modes, singular


def run_scene(seed,condition,flow_mps,config):
    params=config["conditions"][condition]
    scene=shared_scene(seed=seed,n_frames=config["total_frames"],
        calibration_frames=config["calibration_frames_unused"],
        shape=tuple(config["shape_zx"]),prf_hz=config["prf_hz"],
        motion=params["motion"],amplitude_pixels=params["amplitude_pixels"],
        flow_mps=flow_mps,decorrelation=config["flow_decorrelation"])
    scene=turnover_tissue(scene,seed,params["tissue_turnover"])
    x=scene.public.views
    t,b,n=scene.tissue,scene.blood,scene.noise
    background=t+n
    full_motion=oracle_motion(scene.known_displacement,x.shape[1])
    score_n,crop=config["score_last_frames"],config["common_crop_pixels"]
    def score(z):return z[-score_n:].mean(axis=1)[:,crop:-crop,crop:-crop]
    truth_b,truth_t=score(b),score(t)
    bden=max(float(np.linalg.norm(truth_b)**2),1e-30)
    tden=max(float(np.linalg.norm(truth_t)**2),1e-30)
    target_hashes={source:array_hash(score(value)) for source,value in
        (("blood",b),("tissue",t),("noise",n))}
    rows=[]
    for length in config["ensemble_frames"]:
        start=config["total_frames"]-length
        xs,ts,bs,ns=x[start:],t[start:],b[start:],n[start:]
        bg=background[start:]
        motion=MotionFit(full_motion.forward[start:],full_motion.inverse[start:],
            full_motion.phase[start:],full_motion.diagnostics,False)
        aligned_t=motion.align(ts,phase=False)
        aligned_b=motion.align(bs,phase=False)
        modes,sigma=fit_tissue_modes(aligned_t,max(config["oracle_tissue_ranks"]))
        def removed(y,rank):
            aligned=motion.align(y,phase=False)
            v=modes[:,:rank]
            projected=(casorati(aligned)@v)@v.conj().T
            return motion.restore(frames(projected,y.shape),phase=False)
        raw_mix,raw_bg=svd_filter(xs,discard=1),svd_filter(bg,discard=1)
        oracle_motion_mix=native_tissue_subtraction(xs,motion,phase=False)
        oracle_motion_bg=native_tissue_subtraction(bg,motion,phase=False)
        def metrics(output,base):
            increment=score(output-base)
            return {"added_blood_iq_error":float(np.linalg.norm(increment-truth_b)/np.sqrt(bden)),
                "background_power_remaining":float(np.linalg.norm(score(base))**2/tden)}
        controls={"raw_svd":metrics(raw_mix.iq,raw_bg.iq),
            "known_motion_mixture_svd":metrics(oracle_motion_mix.iq,oracle_motion_bg.iq)}
        oracles={}
        for rank in config["oracle_tissue_ranks"]:
            rem_b,rem_bg=removed(bs,rank),removed(bg,rank)
            output_b=bs-rem_b
            output_bg=bg-rem_bg
            target_rem_b,target_rem_bg=score(rem_b),score(rem_bg)
            gain=minimum_gain_for_tissue_limit(score(bg),target_rem_bg,tden,
                config["tissue_power_limit"])
            blood_error=float(np.linalg.norm(target_rem_b)/np.sqrt(bden))
            aligned_v=modes[:,:rank]
            aligned_removed_b=(casorati(aligned_b)@aligned_v)@aligned_v.conj().T
            aligned_overlap=float(np.linalg.norm(aligned_removed_b)/max(np.linalg.norm(aligned_b),1e-30))
            oracles[str(rank)]={
                "blood_iq_error_at_full_subtraction":blood_error,
                "background_power_at_full_subtraction":float(np.linalg.norm(score(output_bg))**2/tden),
                "aligned_blood_projection_norm_fraction":aligned_overlap,
                "tissue_training_energy_fraction":float(np.sum(sigma[:rank]**2)/max(np.sum(sigma*sigma),1e-30)),
                "minimum_gain_for_tissue_limit":gain,
                "minimum_blood_error_at_tissue_limit":None if gain is None else float(gain*blood_error),
                "joint_ambitious_goal_feasible":bool(gain is not None and gain*blood_error<=config["ambitious_blood_iq_error_limit"]),
                "joint_vs_raw_feasible":bool(gain is not None and gain*blood_error<=controls["raw_svd"]["added_blood_iq_error"]+config["raw_svd_blood_error_margin"]),
                "frozen_source_closure_nrmse":float(np.linalg.norm(score(xs-removed(xs,rank))-score(output_bg+output_b))/
                    max(np.linalg.norm(score(xs)),1e-30))}
        roundtrip_b=motion.restore(aligned_b,phase=False)
        rows.append({"ensemble_frames":length,"duration_ms":1000*length/config["prf_hz"],
            "target_sha256":target_hashes,"controls":controls,"oracles":oracles,
            "blood_roundtrip_iq_error":float(np.linalg.norm(score(bs-roundtrip_b))/np.sqrt(bden))})
    return {"seed":seed,"condition":condition,"flow_mps":flow_mps,
        "motion":params["motion"],"tissue_turnover":params["tissue_turnover"],
        "blood_input_power":bden,"tissue_input_power":tden,"lengths":rows}


def summarize(cases,config):
    by_length={}
    for length in config["ensemble_frames"]:
        subset=[next(x for x in c["lengths"] if x["ensemble_frames"]==length)
                for c in cases if c["motion"]!="none"]
        rank1=[x["oracles"]["1"] for x in subset]
        by_length[str(length)]={"moving_cases":len(subset),
            "median_oracle_blood_error_full":float(np.median([x["blood_iq_error_at_full_subtraction"] for x in rank1])),
            "median_oracle_background_power_full":float(np.median([x["background_power_at_full_subtraction"] for x in rank1])),
            "rank1_tissue_limit_feasible_cases":sum(x["minimum_gain_for_tissue_limit"] is not None for x in rank1),
            "rank1_joint_ambitious_cases":sum(x["joint_ambitious_goal_feasible"] for x in rank1),
            "rank1_joint_vs_raw_cases":sum(x["joint_vs_raw_feasible"] for x in rank1)}
    return by_length


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--split",choices=("seen_diagnostic","development","heldout"),required=True)
    p.add_argument("--freeze",action="store_true")
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    if args.freeze and args.split!="development":p.error("Only development can be frozen")
    config=json.loads(CONFIG.read_text())
    if config["total_frames"]!=max(config["ensemble_frames"]) or config["score_last_frames"]>min(config["ensemble_frames"]):
        p.error("Invalid nested-history design")
    sets=[set(config[f"{name}_seeds"]) for name in ("seen_diagnostic","development","heldout")]
    if any(sets[i]&sets[j] for i in range(3) for j in range(i+1,3)):
        p.error("Seed split overlap")
    sources={str(path.relative_to(ROOT)):digest(path) for path in SOURCE_PATHS}
    args.out.mkdir(parents=True,exist_ok=True)
    target=args.out/f"{args.split}.json"
    if target.exists():p.error("Split target already exists")
    if args.split=="development" and not (args.out/"seen_diagnostic.json").exists():
        p.error("Run seen diagnosis first")
    if args.split=="heldout":
        freeze=json.loads((args.out/"freeze.json").read_text())
        if (sources!=freeze["sources"] or config["heldout_seeds"]!=freeze["heldout_seeds"] or
            any(digest(args.out/f"{name}.json")!=freeze[f"{name}_sha256"] for name in
                ("seen_diagnostic","development"))):
            p.error("Held-out lease invalid")
    tic=perf_counter()
    cases=[run_scene(seed,condition,flow,config)
           for seed in config[f"{args.split}_seeds"]
           for condition in config["conditions"] for flow in config["flow_mps"]]
    result={"study":config["study"],"split":args.split,"source_sha256":sources,
        "cases":cases,"summary":summarize(cases,config),"elapsed_s":perf_counter()-tic,
        "truth_use":"known motion and true tissue modes are evaluator-only; this is not a proposed filter",
        "scope":"one fixed matched-tissue temporal projection ray per rank, not a universal separation bound",
        "physical_status":"specified IQ software simulation only"}
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    if args.freeze:
        (args.out/"freeze.json").write_text(json.dumps({"sources":sources,
            "seen_diagnostic_sha256":digest(args.out/"seen_diagnostic.json"),
            "development_sha256":digest(target),"heldout_seeds":config["heldout_seeds"],
            "locked_before_heldout":True},indent=2)+"\n")
    print(json.dumps({"split":args.split,"cases":len(cases),"elapsed_s":result["elapsed_s"],
        "summary":result["summary"],"sha256":digest(target)}))


if __name__=="__main__":main()
