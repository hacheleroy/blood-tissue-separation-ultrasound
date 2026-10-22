"""Prospectively configured comparison of the literature SVD/EVD baselines."""
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
from .benchmark import standard_methods,run_benchmark
from .synthetic import poiseuille
from .stap_experiment import _clean


def default_literature_config():
    ids=["svd_fixed_uncentered","svd_fixed_centered","svd_energy_uncentered",
        "svd_baranger_uncentered","svd_baranger_centered","svd_local_uncentered",
        "svd_long_uncentered","song_local_reconstruction_uncentered",
        "svd_temporal_frequency_uncentered","evd_fixed_uncentered","song_rsvd_uncentered",
        "song_rsvd_downsampled_uncentered","huang_evd_median_uncentered",
        "huang_evd_ncc_uncentered","huang_evd_sigmoid_uncentered","huang_gaussian_iq_extension_uncentered",
        "butterworth_hp_0.01_prf","butterworth_hp_0.03_prf","butterworth_hp_0.08_prf",
        "butterworth_hp_causal_0.03_prf"]
    registry={s["id"]:s for s in standard_methods(packet=32,long_window=96,prf_hz=1000.)}
    return {"version":1,"blas_threads":1,"seeds":[0,1,2],
        "generator":{"n_frames":96,"shape":[8,13],"prf_hz":1000.,"f0_hz":5e6},
        "factorial":{"umax_mps":[.003,.025],"tissue_acceleration_mps2":[0.,.03],"clutter_db":[20.,40.]},
        "measurement_window":32,"measurement_hop":32,"methods":[registry[i] for i in ids],
        "scope":"Statistical IQ, explicit uncentered default, independently frozen parameters; no test-based tuning"}


def run_literature_study(config,output):
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError("Choose a new literature output directory")
    out.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(_clean(config),sort_keys=True,indent=2,allow_nan=False)
    (out/"frozen_config.json").write_text(encoded)
    keys=list(config["factorial"]);records=[]
    for values,seed in itertools.product(itertools.product(*(config["factorial"][k] for k in keys)),config["seeds"]):
        condition=dict(zip(keys,values))
        a=poiseuille(**{**config["generator"],**condition,"seed":seed})
        cfg={"blas_threads":config.get("blas_threads",1),"methods":config["methods"],
            "measurement_window":config["measurement_window"],"measurement_hop":config["measurement_hop"]}
        report=run_benchmark(a,cfg,out/f"case_{len(records):04d}")
        # Full per-window readouts remain reproducible in the case directories.
        # The portable aggregate retains every method's metrics, fitted cutoffs,
        # convergence/fallback diagnostics, status, input hash and code hash.
        rows=[{k:v for k,v in row.items() if k not in ("readout","spec")} for row in report["results"]]
        record={"condition":condition,"seed":seed,"input_sha256":report["input_sha256"],
            "code_sha256":report["code_sha256"],"results":rows,"runtime":report["runtime"]}
        records.append(record)
        print(json.dumps({"literature_cases_completed":len(records),"seed":seed}),flush=True)
    result={"config_sha256":hashlib.sha256(encoded.encode()).hexdigest(),"records":records,
        "parameters_frozen_before_execution":True,"experimental_data":False,
        "provisional_arms":["song_local_reconstruction_uncentered","huang_gaussian_iq_extension_uncentered",
                            "song_rsvd_downsampled_uncentered"]}
    (out/"results.json").write_text(json.dumps(_clean(result),separators=(',',':'),allow_nan=False))
    lines=["# Literature SVD/EVD comparison","",
        "Frozen statistical-IQ cases. See SVD_LITERATURE_AUDIT.md for source fidelity and provisional arms.","",
        "| Method | Completed | Median velocity NRMSE | Median WSS NRMSE | Median blood retention | Median tissue retention |",
        "|---|---:|---:|---:|---:|---:|"]
    for method in config["methods"]:
        rows=[r for c in records for r in c["results"] if r["id"]==method["id"] and r["status"]=="ok"]
        def median(key):
            vals=[r.get("metrics",{}).get(key) for r in rows]
            vals=[v for v in vals if v is not None]
            return 'n/a' if not vals else f'{np.median(vals):.4g}'
        lines.append(f"| {method['id']} | {len(rows)}/{len(records)} | {median('velocity_nrmse')} | {median('wss_nrmse')} | {median('frozen_blood_power_retention')} | {median('frozen_tissue_power_retention')} |")
    lines.extend(["","Every failed or skipped method is retained in results.json. Pooled medians do not establish a universal ranking.",
        "Gaussian-localized IQ overlap-add is an extension; paper power aggregation is reported separately.",
        "Song thresholds are a documented reconstruction, pending original-paper/software parity."])
    (out/"SUMMARY.md").write_text('\n'.join(lines)+'\n')
    return result
