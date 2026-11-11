"""Frozen development/evaluation of three source-blind complex-IQ RPCA arms.

Run from the repository root with PYTHONPATH=src and BLAS threads limited to one:
  python scripts/benchmark_rpca.py development
  python scripts/benchmark_rpca.py evaluation
The evaluation command requires the development selection and verifies that
neither the grid nor the solver source changed after selection.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from time import perf_counter

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

from btbench.filters import rpca_filter, svd_filter
from btbench.measurements import doppler_velocity
from btbench.rpca_comparators import stable_rpca, cauchy_rpca_core, batch_irls_rpca
from btbench.synthetic import poiseuille, advected_scatterers

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/rpca_comparison_v1.json"
SOURCE = ROOT / "configs/separation_fidelity.json"
SOLVER = ROOT / "src/btbench/rpca_comparators.py"
DEFAULT_OUTPUT = ROOT / "results/rpca_comparison_v1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scenes(config, names, seeds):
    src = json.loads(SOURCE.read_text())
    args = {"n_frames": config["frames"], "shape": tuple(config["shape"]),
            "prf_hz": src["generator"]["prf_hz"], "f0_hz": src["generator"]["f0_hz"]}
    bank = [("statistical", x) for x in src["statistical_conditions"]]
    bank += [("psf_scatterer", x) for x in src["scatterer_conditions"]]
    wanted = set(names)
    if wanted != {x["label"] for _, x in bank if x["label"] in wanted}:
        raise ValueError("Unknown or duplicate condition label")
    for family, row in bank:
        if row["label"] not in wanted:
            continue
        kwargs = dict(row)
        name = kwargs.pop("label")
        for seed in seeds:
            generator = poiseuille if family == "statistical" else advected_scatterers
            yield name, family, seed, generator(seed=seed, **args, **kwargs)


def run(method, iq, params, max_iter, tol):
    start = perf_counter()
    if method == "stable_rpca":
        result, diag = stable_rpca(iq, max_iter=max_iter, tol=tol, **params)
    elif method == "cauchy_rpca_core":
        result, diag = cauchy_rpca_core(iq, max_iter=max_iter, tol=tol, **params)
    elif method == "batch_irls_reconstruction":
        result, diag = batch_irls_rpca(iq, max_iter=max_iter, tol=tol, **params)
    elif method == "svd_discard_1":
        fit = svd_filter(iq, discard=1, window=None)
        result, diag = fit.iq, {"iterations": 1, "converged": True}
    elif method == "rpca_pcp":
        fit = rpca_filter(iq, max_iter=300, tol=1e-7)
        result, diag = fit.iq, fit.diagnostics
    elif method == "unfiltered":
        result, diag = iq.copy(), {"iterations": 0, "converged": True}
    else:
        raise ValueError(method)
    return result, diag, perf_counter() - start


def nrmse(output, reference, mask):
    b = reference[:, mask]
    den = np.linalg.norm(b)
    return float(np.linalg.norm(output[:, mask] - b) / den) if den else None


def power(a, mask):
    return float(np.mean(np.abs(a[:, mask])**2))


def ratio(num, den):
    return float(num / den) if den > 0 else None


def evaluate(acq, output, paired, background, last):
    blood, tissue = acq.blood_true, acq.tissue_true
    flow, wall, tissue_roi = acq.flow_mask, acq.near_wall_mask, acq.tissue_roi
    observed = output[-last:]
    response = paired[-last:]
    bg = background[-last:]
    ref = blood[-last:]
    t_ref = tissue[-last:]
    truth_power = np.mean(np.abs(ref)**2, axis=0)
    observed_power = np.mean(np.abs(observed)**2, axis=0)
    # Report velocity only where the output has measurable power. The gate is
    # fixed from the blood-only truth; it is used solely during evaluation.
    valid_power = observed_power >= 0.01 * truth_power
    v_out, coh = doppler_velocity(observed, acq.prf_hz, acq.f0_hz,
                                  phase_sign=acq.phase_sign)
    v_ref, _ = doppler_velocity(ref, acq.prf_hz, acq.f0_hz,
                                phase_sign=acq.phase_sign)
    valid = wall & valid_power & np.isfinite(v_out) & np.isfinite(v_ref)
    result = {
        "observed_flow_iq_nrmse": nrmse(observed, ref, flow),
        "observed_near_wall_iq_nrmse": nrmse(observed, ref, wall),
        "paired_flow_iq_nrmse": nrmse(response, ref, flow),
        "paired_near_wall_iq_nrmse": nrmse(response, ref, wall),
        "tissue_only_background_power_ratio": ratio(power(bg, tissue_roi), power(t_ref, tissue_roi)),
        "tissue_only_flow_power_to_blood_ratio": ratio(power(bg, flow), power(ref, flow)),
        "near_wall_output_power_to_blood_ratio": ratio(power(observed, wall), power(ref, wall)),
        "near_wall_velocity_coverage": float(np.mean(valid[wall])),
        "near_wall_velocity_mae_vs_blood_iq_mps": float(np.mean(np.abs(v_out[valid] - v_ref[valid]))) if valid.any() else None,
        "near_wall_velocity_mae_vs_prescribed_mps": float(np.mean(np.abs(v_out[valid] - acq.velocity_true_mps[valid]))) if valid.any() else None,
        "blood_only_near_wall_velocity_mae_vs_prescribed_mps": float(np.mean(np.abs(v_ref[wall] - acq.velocity_true_mps[wall]))),
        "near_wall_pixel_count": int(wall.sum()),
        "velocity_gate_power_fraction": 0.01,
        "mean_valid_near_wall_coherence": float(np.mean(coh[valid])) if valid.any() else None
    }
    if not all(v is None or not isinstance(v, float) or np.isfinite(v) for v in result.values()):
        raise FloatingPointError("Nonfinite benchmark metric")
    return result


def one_case(acq, method, params, config):
    noise = acq.iq - acq.blood_true - acq.tissue_true
    bg_input = acq.tissue_true + noise
    output, diag, seconds = run(method, acq.iq, params, config["max_iterations"], config["tolerance"])
    background, _, _ = run(method, bg_input, params, config["max_iterations"], config["tolerance"])
    metrics = evaluate(acq, output, output - background, background, config["score_last_frames"])
    return {"metrics": metrics, "diagnostics": diag, "fit_seconds": seconds}


def selection_score(metrics):
    a = metrics["observed_near_wall_iq_nrmse"]
    b = metrics["paired_near_wall_iq_nrmse"]
    c = metrics["tissue_only_background_power_ratio"]
    return math.log1p(a) + 0.5 * math.log1p(b) + 0.2 * math.log1p(100*c)


def development(config, out):
    cases = list(scenes(config, config["development_conditions"], [config["development_seed"]]))
    selected = {}
    grid = []
    for entry in config["methods_in_order"]:
        method = entry["id"]
        for index, params in enumerate(entry["candidates"]):
            rows = []
            for name, family, seed, acq in cases:
                fit = one_case(acq, method, params, config)
                rows.append({"condition": name, "family": family, "seed": seed, **fit,
                             "selection_score": selection_score(fit["metrics"])})
            score = float(np.mean([r["selection_score"] for r in rows]))
            grid.append({"method": method, "candidate_index": index, "params": params,
                         "mean_selection_score": score, "cases": rows})
            print(f"DEV {method} {index}: {score:.4f}", flush=True)
        best = min((g for g in grid if g["method"] == method),
                   key=lambda g: (g["mean_selection_score"], g["candidate_index"]))
        selected[method] = {"candidate_index": best["candidate_index"], "params": best["params"],
                            "mean_selection_score": best["mean_selection_score"]}
        print(f"SELECTED {method}: candidate {best['candidate_index']}", flush=True)
    payload = {"config_sha256": digest(CONFIG), "source_config_sha256": digest(SOURCE),
               "solver_sha256": digest(SOLVER), "selected": selected, "grid": grid}
    (out / "development.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    (out / "selected.json").write_text(json.dumps({k: payload[k] for k in
        ("config_sha256", "source_config_sha256", "solver_sha256", "selected")},
        indent=2, allow_nan=False) + "\n")


def evaluation(config, out):
    selected = json.loads((out / "selected.json").read_text())
    for key, path in [("config_sha256", CONFIG), ("source_config_sha256", SOURCE), ("solver_sha256", SOLVER)]:
        if selected[key] != digest(path):
            raise ValueError(f"Evaluation refused: {path} changed since development")
    src = json.loads(SOURCE.read_text())
    names = [x["label"] for x in src["statistical_conditions"] + src["scatterer_conditions"]]
    rows = []
    methods = [x["id"] for x in config["methods_in_order"]] + config["reference_methods"]
    for name, family, seed, acq in scenes(config, names, config["evaluation_seeds"]):
        for method in methods:
            params = selected["selected"][method]["params"] if method in selected["selected"] else {}
            fit = one_case(acq, method, params, config)
            rows.append({"method": method, "condition": name, "family": family, "seed": seed,
                         "params": params, **fit})
        print(f"EVAL {family} {name} seed={seed}: {len(methods)} methods", flush=True)
    metrics = ["observed_flow_iq_nrmse", "observed_near_wall_iq_nrmse",
               "paired_near_wall_iq_nrmse", "tissue_only_background_power_ratio",
               "tissue_only_flow_power_to_blood_ratio", "near_wall_velocity_coverage",
               "near_wall_velocity_mae_vs_blood_iq_mps", "near_wall_velocity_mae_vs_prescribed_mps"]
    summary = []
    for method in methods:
        for family in ["statistical", "psf_scatterer", "all"]:
            subset = [r for r in rows if r["method"] == method and
                      (family == "all" or r["family"] == family)]
            medians = {m: float(np.median([r["metrics"][m] for r in subset
                                if r["metrics"][m] is not None])) for m in metrics}
            summary.append({"method": method, "family": family, "cases": len(subset),
                            "median": medians,
                            "converged_fraction": float(np.mean([r["diagnostics"].get("converged", True) for r in subset])),
                            "median_fit_seconds": float(np.median([r["fit_seconds"] for r in subset]))})
    payload = {"config_sha256": selected["config_sha256"], "source_config_sha256": selected["source_config_sha256"],
               "solver_sha256": selected["solver_sha256"], "selected": selected["selected"],
               "rows": rows, "summary": summary,
               "note": "Each nonlinear counterfactual was independently refit to tissue plus noise; last 16 frames scored. Seeds 61,73,89 were excluded from tuning."}
    (out / "evaluation.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["development", "evaluation"])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    cli = parser.parse_args()
    cli.output.mkdir(parents=True, exist_ok=True)
    specification = json.loads(CONFIG.read_text())
    (development if cli.stage == "development" else evaluation)(specification, cli.output)
