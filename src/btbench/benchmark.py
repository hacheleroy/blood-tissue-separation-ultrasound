"""Explicit method registry; filters receive IQ/settings, never reference fields."""
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

from .filters import svd_filter,aer_filter,local_filter,rpca_filter,temporal_highpass,sliding_filter
from .axes import motion_packet_filter,multilag_filter,ica_filter
from .measurements import evaluate,frozen_component_metrics,evaluate_windows


def run_filter(iq,spec,tissue_roi=None,views=None):
    p=dict(spec.get("params",{}))
    method=spec["method"]
    if method=="stap":
        from .stap import local_stap
        if views is None:
            raise ValueError("STAP requires uncompounded IQ views")
        if not np.allclose(iq,np.mean(views,axis=1)):
            raise ValueError("STAP comparison needs iq to be the same coherent view mean")
        return local_stap(views,**p)
    registry={"svd":svd_filter,"aer":aer_filter,"local":local_filter,"rpca":rpca_filter,
              "axis1":motion_packet_filter,"axis2":multilag_filter,"ica":ica_filter,
              "highpass":temporal_highpass,"sliding":sliding_filter}
    if method not in registry:
        raise ValueError("Unknown IQ method "+method)
    if method in ("axis1","axis2","ica"):
        p["tissue_roi"]=tissue_roi
    return registry[method](iq,**p)


def grid_specs(base,grid):
    keys=list(grid)
    for values in itertools.product(*(grid[k] for k in keys)):
        spec={**base,"params":{**base.get("params",{}),**dict(zip(keys,values))}}
        yield spec


def oracle(acquisition,base,grid,objective="velocity_nrmse",geometry_factor=2.):
    """Same-case truth-assisted bound over the DECLARED candidate grid only.

    Never used as a deployable setting. Every candidate is logged. Separate
    objectives are needed for velocity, IQ, and WSS; no implicit equivalence.
    """
    if acquisition.velocity_true_mps is None and objective!="blood_iq_nrmse":
        raise ValueError("Oracle requires independent velocity/WSS truth")
    if objective=="blood_iq_nrmse" and acquisition.blood_true is None:
        raise ValueError("IQ oracle requires independent blood-IQ truth")
    grid=dict(grid)
    if grid.get("discard")=="all":
        p=int(np.prod(acquisition.iq.shape[1:]))
        t=min(base.get("params",{}).get("window") or len(acquisition.iq),len(acquisition.iq))
        grid["discard"]=list(range(min(p,t)))
    best=None;curve=[]
    for spec in grid_specs(base,grid):
        try:
            result=run_filter(acquisition.iq,spec,acquisition.tissue_roi,acquisition.iq_views)
        except (ValueError,MemoryError,np.linalg.LinAlgError) as exc:
            curve.append({"params":spec["params"],"loss":None,"status":"invalid","reason":str(exc)})
            continue
        metrics=evaluate(result.iq,acquisition,geometry_factor)
        loss=metrics.get(objective)
        if loss is None:
            raise ValueError("Oracle objective unavailable: "+objective)
        curve.append({"params":spec["params"],"loss":float(loss)})
        if best is None or loss<best[0]:
            best=(loss,result,spec,metrics)
    if best is None:
        raise ValueError("Oracle grid has no valid candidates")
    return best[1],{"oracle":True,"selection_uses_test_truth":True,
                    "objective":objective,"selected":best[2],"curve":curve}


def standard_methods(packet=64,long_window=256):
    """Factorial controls; short window and long window report distinct latency."""
    methods=[]
    for centered in (False,True):
        tag="centered" if centered else "uncentered"
        for sel in ("fixed","energy","ssm_surrogate"):
            methods.append({"id":f"svd_{sel}_{tag}","method":"svd","params":{
                "centered":centered,"window":packet,"selector":sel,"discard":1}})
        methods.append({"id":f"svd_long_{tag}","method":"sliding","params":{
            "method":"svd","centered":centered,"window":long_window,"hop":packet,"discard":1}})
        methods.append({"id":f"svd_local_{tag}","method":"local","params":{
            "method":"svd","centered":centered,"patch_shape":[6,6],"stride":[3,3],"window":packet,"discard":1}})
        for domain in ("temporal","spatial"):
            for alpha in (0.5,0.25):
                methods.append({"id":f"aer_{domain}_alpha{alpha}_{tag}","method":"aer",
                    "params":{"domain":domain,"centered":centered,"epsilon_rel":1e-3,"alpha":alpha,"window":long_window}})
        methods.append({"id":f"aer_noise_floor_{tag}","method":"aer","params":{
            "centered":centered,"domain":"temporal","epsilon_rule":"noise_floor","window":long_window}})
        methods.append({"id":f"axis1_{tag}","method":"axis1","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3]}})
        methods.append({"id":f"axis1_rigid_{tag}","method":"axis1","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],"rigid":True}})
        methods.append({"id":f"axis2_{tag}","method":"axis2","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],
            "rank":6,"lags":[1,2,4,8,16],"use_memory":False}})
        methods.append({"id":f"axis2_with_memory_{tag}","method":"axis2","params":{
            "packet":packet,"centered":centered,"patch_shape":[6,6],"stride":[3,3],
            "rank":6,"lags":[1,2,4,8,16],"use_memory":True}})
    methods += [
        {"id":"ica_circular_complex","method":"ica","params":{"packet":packet,"patch_shape":[6,6],"stride":[3,3],"rank":6}},
        {"id":"rpca_pcp","method":"rpca","params":{}},
    ]
    return methods


def run_benchmark(acquisition,config,output):
    out=Path(output)
    if (out/"results.json").exists():
        raise FileExistsError("Output already contains results; choose a new directory")
    out.mkdir(parents=True,exist_ok=True)
    a=acquisition
    results=[]
    start_all=time.perf_counter()
    threads=int(config.get("blas_threads",1))
    with threadpool_limits(limits=threads):
        for spec in config["methods"]:
            t0=time.perf_counter()
            row={"id":spec["id"],"method":spec["method"],"oracle":False,"spec":spec}
            try:
                if "oracle_grid" in spec:
                    result,diag=oracle(a,spec,spec["oracle_grid"],spec.get("objective","velocity_nrmse"),config.get("geometry_factor",2.))
                    row.update(diag)
                else:
                    result=run_filter(a.iq,spec,a.tissue_roi,a.iq_views)
                metrics=evaluate(result.iq,a,config.get("geometry_factor",2.))
                metrics.update(frozen_component_metrics(result,a))
                row.update(status="ok",metrics=metrics,diagnostics=result.diagnostics)
                row["readout"]=evaluate_windows(result.iq,a,config.get("measurement_window",64),
                    config.get("measurement_hop",64),config.get("geometry_factor",2.))
                if config.get("save_iq",False):
                    # IDs are labels, never paths supplied by external datasets.
                    name=hashlib.sha256(spec["id"].encode()).hexdigest()[:16]
                    np.save(out/f"filtered_{name}.npy",result.iq)
            except (ValueError,MemoryError,np.linalg.LinAlgError) as exc:
                row.update(status="skipped",reason=str(exc))
            row["wall_time_s"]=time.perf_counter()-t0
            try:
                import resource,sys
                peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                row["process_high_water_rss_bytes"]=int(peak if sys.platform=="darwin" else peak*1024)
            except ImportError:
                row["process_high_water_rss_bytes"]=None
            results.append(row)
    code=hashlib.sha256()
    for module in sorted(Path(__file__).parent.glob("*.py")):
        code.update(module.name.encode());code.update(module.read_bytes())
    report={"config":config,"results":results,"acquisition":a.metadata,"code_sha256":code.hexdigest(),
            "iq_shape":list(a.iq.shape),"input_sha256":hashlib.sha256(a.iq.tobytes()).hexdigest(),
            "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
            "platform":platform.platform(),"blas_threads":threads,"elapsed_s":time.perf_counter()-start_all},
            "interpretation":"synthetic smoke test unless independent experimental truth is provided"}
    # Strict JSON: invalid numerical estimates become null with explicit fractions
    # and candidate selection records, not nonstandard Infinity tokens.
    def clean(v):
        if isinstance(v,dict): return {k:clean(w) for k,w in v.items()}
        if isinstance(v,(tuple,list)): return [clean(w) for w in v]
        if isinstance(v,(float,np.floating)) and not np.isfinite(v): return None
        if isinstance(v,np.generic): return v.item()
        return v
    (out/"results.json").write_text(json.dumps(clean(report),indent=2,allow_nan=False))
    lines=["# Benchmark run", "", "This run is a software check; it does not establish clinical or acoustic performance.","",
           "| Method | Status | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE | Time (s) |",
           "|---|---|---:|---:|---:|---:|"]
    for row in results:
        m=row.get("metrics",{})
        def fmt(key):
            v=m.get(key)
            return "—" if v is None else f"{v:.3g}"
        lines.append(f"| {row['id']} | {row['status']} | {fmt('velocity_nrmse')} | {fmt('near_wall_velocity_nrmse')} | {fmt('wss_nrmse')} | {row['wall_time_s']:.3g} |")
    (out/"SUMMARY.md").write_text("\n".join(lines)+"\n")
    return report
