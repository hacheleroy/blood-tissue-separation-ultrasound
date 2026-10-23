"""Prospectively configured source-fidelity studies for Axes 1 and 2."""
from dataclasses import replace
from pathlib import Path
import copy
import hashlib
import json
import numpy as np
from .benchmark import run_benchmark,standard_methods
from .synthetic import poiseuille,advected_scatterers
from .stap_experiment import _clean


def default_separation_config():
    ids=["svd_fixed_uncentered","svd_fixed_centered","svd_baranger_uncentered","svd_baranger_centered",
        "svd_local_uncentered","svd_local_centered","svd_long_uncentered","svd_long_centered",
        "butterworth_hp_0.01_prf","butterworth_hp_0.03_prf","butterworth_hp_0.08_prf",
        "aer_temporal_alpha0.5_uncentered","aer_temporal_alpha0.25_uncentered",
        "aer_temporal_alpha0.5_centered","aer_temporal_alpha0.25_centered",
        "aer_spatial_alpha0.5_uncentered","aer_spatial_alpha0.5_centered",
        "axis1_uncentered","axis1_centered","axis1_rigid_uncentered","axis1_rigid_centered",
        "axis2_uncentered","axis2_centered","axis2_with_memory_uncentered","axis2_with_memory_centered",
        "ica_circular_complex","rpca_pcp"]
    registry={m["id"]:m for m in standard_methods(packet=32,long_window=96,prf_hz=1000.)}
    methods=[{"id":"unfiltered","method":"identity"}]+[registry[i] for i in ids]
    for name,extra in [("ica_uncentered_retained_300",{"centered":False,"ica_nonconvergence":"retain"}),
                       ("ica_centered_retained_300",{"centered":True,"ica_nonconvergence":"retain"}),
                       ("ica_centered_retained_1500",{"centered":True,"ica_nonconvergence":"retain","ica_options":{"max_iter":1500}})]:
        spec=copy.deepcopy(registry["ica_circular_complex"])
        spec["id"]=name;spec["params"].update(extra);methods.append(spec)
    for original,name,extra in [
        ("axis1_uncentered","axis1_no_memory",{"rho":0.}),
        ("axis1_uncentered","axis1_no_phase",{"phase_compensation":False}),
        ("axis2_with_memory_uncentered","axis2_no_phase",{"phase_compensation":False}),
        ("axis1_uncentered","axis1_phase_error",{"motion_phase_error_rad_per_frame":.02}),
        ("axis2_with_memory_uncentered","axis2_phase_error",{"motion_phase_error_rad_per_frame":.02})]:
        spec=copy.deepcopy(registry[original]);spec["id"]=name;spec["params"].update(extra);methods.append(spec)
    for objective,constraints,name in [
        ("blood_iq_nrmse",{},"oracle_svd_iq"),
        ("frozen_tissue_power_retention",{"frozen_blood_iq_nrmse":{"max":.25}},"oracle_svd_distortion_025"),
        ("frozen_tissue_power_retention",{"frozen_blood_iq_nrmse":{"max":.5}},"oracle_svd_distortion_050")]:
        methods.append({"id":name,"method":"svd","params":{"window":32},
            "oracle_grid":{"discard":[0,1,2,3,4,8,16,32]},"objective":objective,"constraints":constraints})
    statistical=[
        {"label":"slow_overlap","umax_mps":.003,"tissue_velocity_mps":.002,"tissue_acceleration_mps2":0.,"clutter_db":40.,"decorrelation":.995},
        {"label":"fast_accelerating","umax_mps":.025,"tissue_acceleration_mps2":.03,"clutter_db":30.},
        {"label":"texture_jump","umax_mps":.01,"motion_event":True,"clutter_db":30.},
        {"label":"weak_decorrelating","umax_mps":.003,"snr_db":0.,"decorrelation":.8,"clutter_db":40.},
        {"label":"opposed_flow","umax_mps":-.006,"tissue_velocity_mps":.002,"tissue_acceleration_mps2":0.}]
    scatterers=[
        {"label":"psf_stationary_tissue","umax_mps":.02,"tissue_velocity_mps":0.},
        {"label":"psf_bulk_motion","umax_mps":.006,"tissue_velocity_mps":.002},
        {"label":"psf_deforming_tissue","umax_mps":.003,"tissue_velocity_mps":.001,"strain_motion_m":.00008,"lateral_motion_m":.00006},
        {"label":"psf_opposed_flow","umax_mps":-.006,"tissue_velocity_mps":.002},
        {"label":"psf_equal_static","umax_mps":0.,"tissue_velocity_mps":0.,"uniform_flow":True}]
    return {"version":3,"blas_threads":1,"seeds":[17,29,43],
        "generator":{"n_frames":96,"shape":[8,13],"prf_hz":1000.,"f0_hz":5e6},
        "statistical_conditions":statistical,"scatterer_conditions":scatterers,"methods":methods,
        "measurement_window":32,"measurement_hop":32,"injection_amplitudes":[.1,1.],
        "filter_roi":"estimated from mixture; evaluation flow mask is never passed to filters",
        "scope":"statistical IQ and linear 2-D PSF transport; no RF, RBC or fluid-structure validation"}


def run_separation_study(config,output):
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError("Choose a new separation output directory")
    out.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(_clean(config),sort_keys=True,indent=2,allow_nan=False)
    (out/"frozen_config.json").write_text(encoded)
    records=[]
    for family,conditions,generator in [("statistical",config["statistical_conditions"],poiseuille),
                                        ("psf_transport",config["scatterer_conditions"],advected_scatterers)]:
        for condition in conditions:
            for seed in config["seeds"]:
                params={**config["generator"],**{k:v for k,v in condition.items() if k!="label"},"seed":seed}
                a=generator(**params)
                # ROI used by motion fitting is inferred from IQ; no anatomical
                # mask from the forward model is handed to a data-only filter.
                a=replace(a,tissue_roi=None)
                cfg={k:config[k] for k in ("blas_threads","methods","measurement_window","measurement_hop","injection_amplitudes")}
                report=run_benchmark(a,cfg,out/f"case_{len(records):04d}")
                rows=[{k:v for k,v in r.items() if k not in ("readout","spec")} for r in report["results"]]
                records.append({"family":family,"condition":condition,"seed":seed,"scene_parameters":params,
                    "input_sha256":report["input_sha256"],"code_sha256":report["code_sha256"],
                    "runtime":report["runtime"],"data_only_pareto_ids":report["data_only_pareto_ids"],"results":rows})
                print(json.dumps({"separation_cases_completed":len(records),"family":family,"condition":condition["label"],"seed":seed}),flush=True)
    report={"config_sha256":hashlib.sha256(encoded.encode()).hexdigest(),"records":records,
        "parameters_frozen_before_execution":True,"experimental_data":False,
        "primary_endpoints":["frozen_blood_iq_nrmse","frozen_tissue_power_retention"],
        "scope":config["scope"]}
    (out/"results.json").write_text(json.dumps(_clean(report),separators=(',',':'),allow_nan=False))
    lines=["# Separation fidelity study","",
        "Primary comparison: isolated blood complex-IQ distortion versus tissue power retention.",
        "No fitted gain/phase correction. Filter parameters and case bank are frozen before execution.",
        "The equal-static case is an identifiability control; complete separation is not expected.","",
        "| Family | Method | Completed | Median blood distortion | Median tissue retention | Median spectral TV | ICA convergence |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for family in ("statistical","psf_transport"):
        cases=[c for c in records if c["family"]==family]
        for method in config["methods"]:
            rows=[r for c in cases for r in c["results"] if r["id"]==method["id"] and r["status"]=="ok"]
            def median(key):
                vals=[r["metrics"].get(key) for r in rows];vals=[v for v in vals if v is not None]
                return "n/a" if not vals else f"{np.median(vals):.4g}"
            fits=[p for r in rows for p in r.get("diagnostics",{}).get("packets",[]) if "ica_converged" in p]
            convergence="n/a" if not fits else f"{sum(p['ica_converged'] for p in fits)}/{len(fits)}"
            lines.append(f"| {family} | {method['id']} | {len(rows)}/{len(cases)} | {median('frozen_blood_iq_nrmse')} | {median('frozen_tissue_power_retention')} | {median('frozen_blood_spectral_tv')} | {convergence} |")
    lines.extend(["","Use per-condition Pareto sets; pooled medians are descriptive, not a universal ranking.",
        "Infeasible oracle searches retain their full candidate curves. RPCA uses labeled refitted injection responses.",
        "The 2-D PSF model includes transport and boundary mixing but omits channel propagation and RBC physics."])
    (out/"SUMMARY.md").write_text("\n".join(lines)+"\n")
    return report
