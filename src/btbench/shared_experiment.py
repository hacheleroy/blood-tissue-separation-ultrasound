"""Frozen common-input study and global training-only operating-point selection."""
from dataclasses import replace, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from pathlib import Path
import gzip
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from threadpoolctl import threadpool_limits
from .shared import AXIS_NAMES, SharedContext, fit_shared, shared_scene, array_hash
from .registration import MotionGate
from .measurements import relative_norm, blood_fidelity
from .followup_experiment import numerical_source_hash
from .stap_experiment import _clean


def default_shared_config():
    methods=[]
    def add(method, motion, horizon, suffix='', family=None, **params):
        methods.append(dict(id=f'{method}_{motion}_h{horizon}'+('_'+suffix if suffix else ''),
                            family=family or method, method=method, motion=motion,
                            horizon=horizon, evaluation_only=motion=='known', **params))
    for motion in ('none','smooth'):
        for k in (1,4,'all'):
            add('svd',motion,64,f'k{k}',discard=k)
        for method in ('weighted_svd','axis1'):
            for k in (1,4):add(method,motion,64,f'k{k}',discard=k,rho=.75)
        for method in ('svd','weighted_svd','axis1'):
            add(method,motion,64,'centered_k1',discard=1,centered=True,rho=.75)
        add('baranger',motion,64)
        for eps in (1.,26.,1000.):
            add('lee_aer',motion,64,f'temporal_abs{eps:g}',family='lee_temporal_raw',domain='temporal',epsilon_abs=eps)
        add('lee_aer',motion,64,'temporal_relative',family='lee_temporal_raw',domain='temporal',epsilon_rel=.001)
        add('lee_aer',motion,64,'spatial_abs26',family='lee_spatial_raw',domain='spatial',epsilon_abs=26.)
        add('lee_aer',motion,64,'spatial_relative',family='lee_spatial_raw',domain='spatial',epsilon_rel=.001)
        add('lee_aer',motion,64,'temporal_unit',family='lee_unit_extension',domain='temporal',epsilon_rel=.001,unit_operator=True)
        for rank in (4,8):
            for width in (.1,.3):
                add('axis2',motion,64,f'r{rank}_p{width:g}',rank=rank,phase_width=width,rho=.75,
                    ridge_rel=.001,epsilon_rel=1e-6,lags=[1,2,4,8],gain_cap=1.)
        add('stap',motion,64,'orthogonal',family='stap_orthogonal',adaptive=False)
        for loading in (.01,.1):
            add('stap',motion,64,f'loaded{loading:g}',family='stap_loaded',adaptive=True,loading_rel=loading)
    for k in (1,4,'all'):add('svd','none',16,f'k{k}',discard=k)
    add('weighted_svd','none',16,'k1',discard=1,rho=.75)
    add('axis1','none',16,'k1',discard=1,rho=.75)
    add('axis2','none',16,'r4_p0.1',rank=4,phase_width=.1,rho=.75,ridge_rel=.001,
        epsilon_rel=1e-6,lags=[1,2,4,8],gain_cap=1.)
    add('stap','none',16,'orthogonal',family='stap_orthogonal',adaptive=False)
    add('stap','none',16,'loaded0.01',family='stap_loaded',adaptive=True,loading_rel=.01)
    add('lee_aer','none',16,'temporal_abs26',family='lee_temporal_raw',domain='temporal',epsilon_abs=26.)
    for motion in ('rigid','known'):
        add('svd',motion,64,'k1',discard=1)
        add('axis1',motion,64,'k1',discard=1,rho=.75)
        add('axis2',motion,64,'r4_p0.1',rank=4,phase_width=.1,rho=.75,ridge_rel=.001,
            epsilon_rel=1e-6,lags=[1,2,4,8],gain_cap=1.)
        add('stap',motion,64,'loaded0.01',family='stap_loaded',adaptive=True,loading_rel=.01)
        add('lee_aer',motion,64,'temporal_abs26',family='lee_temporal_raw',domain='temporal',epsilon_abs=26.)
    add('identity','none',16)
    for motion in ('rigid','smooth','known'):add('roundtrip',motion,16)
    return dict(study='shared_axes_v05',version='0.4.0',blas_threads=1,case_workers=4,packet_frames=16,score_start=64,
        generator=dict(n_frames=128,calibration_frames=64,shape=[8,16],prf_hz=1000.,f0_hz=5e6,
                       clutter_db=30.,snr_db=20.,decorrelation=.985),
        training_seeds=[101,211],evaluation_seeds=[307,419,523],
        background_targets=[.01,.0001],selection='global pooled training means, grouped by family/motion/support',
        conditions=[
            dict(id='slow_overlap',flow_mps=.003,clutter_db=40.),
            dict(id='separated_low_rank',flow_mps=.02,low_rank=True),
            dict(id='subpixel_translation',flow_mps=.006,motion='translation',amplitude_pixels=.8),
            dict(id='smooth_deformation',flow_mps=.006,motion='smooth',amplitude_pixels=.8),
            dict(id='unmodeled_deformation',flow_mps=.006,motion='unmodeled',amplitude_pixels=.8),
            dict(id='calibration_drift',flow_mps=.006,acceleration_mps2=.02,view_phase_drift_rad=.3),
            dict(id='opposed_flow',flow_mps=-.008,phase_sign=-1),
            dict(id='identical_statistics',flow_mps=0.,identical_statistics=True)],methods=methods)


def default_motion_repair_config():
    """Predeclared matched-arm repair study; no evaluation-seed tuning."""
    base=default_shared_config()
    base.update(study='motion_repair_v06',version='0.5.0',
                training_seeds=[601,787,919],evaluation_seeds=[1601,1907,2203,2503,2801],
                motion_options=asdict(MotionGate()))
    templates=[m for m in base['methods'] if m['motion']=='none' and m['horizon']==64 and m['id'] in {
        'svd_none_h64_k1','weighted_svd_none_h64_k1','axis1_none_h64_k1',
        'axis2_none_h64_r4_p0.1','stap_none_h64_loaded0.01',
        'lee_aer_none_h64_temporal_abs26','lee_aer_none_h64_spatial_abs26','baranger_none_h64'}]
    methods=[]
    for mode in ('none','rigid','smooth','known'):
        for template in templates:
            m=dict(template);m['motion']=mode;m['id']=template['id'].replace('_none_',f'_{mode}_')
            m['evaluation_only']=mode=='known';m['diagnostic_only']=False
            methods.append(m)
    for mode in ('none','rigid','smooth','known'):
        methods.append(dict(id=f'roundtrip_{mode}_h16',family='roundtrip',method='roundtrip',motion=mode,horizon=16,
                            evaluation_only=mode=='known',diagnostic_only=False))
    base['methods']=methods
    return base


def validate_shared_config(config):
    if set(config['training_seeds'])&set(config['evaluation_seeds']):raise ValueError('Seed splits overlap')
    if not config['training_seeds'] or not config['evaluation_seeds']:raise ValueError('Both seed splits are required')
    ids=[m['id'] for m in config['methods']]
    if len(set(ids))!=len(ids) or not ids:raise ValueError('Method IDs must be unique and nonempty')
    if len({c['id'] for c in config['conditions']})!=len(config['conditions']):raise ValueError('Duplicate condition')
    packet=config['packet_frames'];n=config['generator']['n_frames'];start=config['score_start']
    if packet<2 or n%packet or start%packet or not 0<=start<n:raise ValueError('Invalid common scoring interval')
    for spec in config['methods']:
        if spec.get('evaluation_only',False)!=(spec['motion']=='known'):
            raise ValueError('Only the explicit known-motion path is evaluation-only')
        if spec.get('diagnostic_only',False) and not spec['motion'].endswith('_ungated'):
            raise ValueError('Only ungated ablations may be diagnostic-only')
        if spec['horizon']>start or spec['horizon']<packet or spec['horizon']%packet:
            raise ValueError('The warm-up must cover every finite history')
        if spec['method']=='axis2' and max(spec.get('lags',[1,2,4,8]))>=packet:
            raise ValueError('Lag pairs must fit inside each packet')
    return config


def select_global_operating_points(records, methods, targets):
    """One setting per family, motion and support; no condition-specific tuning."""
    groups={};by_id={m['id']:m for m in methods}
    for m in methods:
        if not m.get('evaluation_only') and not m.get('diagnostic_only') and m['method'] not in ('identity','roundtrip'):
            key=(m['family'],m['motion'],m['horizon'])
            groups.setdefault(key,[]).append(m['id'])
    training=[r for r in records if r['split']=='training'];evaluation=[r for r in records if r['split']=='evaluation']
    if not training or not evaluation:raise ValueError('Training and held-out records are required')
    lookup=lambda case,mid:next(r for r in case['results'] if r['id']==mid)
    selections=[]
    for (family,motion,horizon),mids in sorted(groups.items()):
        candidates=[]
        for mid in mids:
            rows=[lookup(c,mid) for c in training]
            valid=all(r['status']=='ok' and all(np.isfinite(r['metrics'][k]) for k in
                      ('blood_iq_nrmse','refitted_background_to_tissue_power')) for r in rows)
            candidates.append(dict(id=mid,valid=valid,training_mean_iq_nrmse=float(np.mean([r['metrics']['blood_iq_nrmse'] for r in rows])) if valid else None,
                training_mean_background=float(np.mean([r['metrics']['refitted_background_to_tissue_power'] for r in rows])) if valid else None))
        for target in [None]+list(targets):
            eligible=[c for c in candidates if c['valid'] and (target is None or c['training_mean_background']<=target)]
            item=dict(family=family,motion=motion,horizon_frames=horizon,background_target=target,
                      target_scope='arithmetic mean over all training conditions and seeds',candidates=candidates,
                      evaluation_truth_used_for_selection=False)
            if not eligible:item.update(status='infeasible_on_training',selected_id=None,evaluation=[])
            else:
                winner=min(eligible,key=lambda c:(c['training_mean_iq_nrmse'],c['id']))
                mid=winner['id'];held=[]
                for case in evaluation:
                    r=lookup(case,mid)
                    held.append(dict(condition=case['condition'],seed=case['seed'],status=r['status'],metrics=r['metrics'],
                        meets_target=None if target is None else bool(r['status']=='ok' and r['metrics']['refitted_background_to_tissue_power']<=target)))
                item.update(status='selected',selected_id=mid,selected_spec=by_id[mid],evaluation=held)
            selections.append(item)
    return selections


def _diagnostics(context, scene):
    error=context.displacement-scene.known_displacement
    return dict(motion_fit=context.motion_diagnostics,calibration=context.calibration,
        input_sha256=array_hash(context.public.views),
        evaluation_only_displacement_rmse_pixels=float(np.sqrt(np.mean(error**2))),
        evaluation_only_interior_displacement_rmse_pixels=float(np.sqrt(np.mean(error[:,:,1:-1,1:-1]**2))))


def _evaluate_result(actual, background, weak, scene, start):
    truth=scene.blood[start:];tissue=scene.tissue[start:]
    power=lambda x:float(np.sum(abs(x)**2))
    error=lambda x:relative_norm(x-truth,truth)
    response1=actual.views-background.views;response01=(weak.views-background.views)/.1
    b=actual.apply_frozen(scene.blood);t=actual.apply_frozen(scene.tissue);n=actual.apply_frozen(scene.noise)
    metrics=dict(blood_iq_nrmse=error(actual.views),
        output_to_blood_power_ratio=power(actual.views)/max(power(truth),1e-30),
        coherent_mean_iq_nrmse=relative_norm(actual.views.mean(1)-truth.mean(1),truth.mean(1)),
        refitted_background_to_tissue_power=power(background.views)/max(power(tissue),1e-30),
        refitted_blood_response_nrmse_alpha1=error(response1),
        refitted_blood_response_nrmse_alpha0_1=error(response01),
        incremental_power_ratio_alpha1=power(response1)/max(power(truth),1e-30),
        incremental_power_ratio_alpha0_1=power(response01)/max(power(truth),1e-30),
        frozen_component_closure_nrmse=float(np.linalg.norm(b+t+n-actual.views)/max(np.linalg.norm(actual.views),1e-30)),
        frozen_blood_iq_nrmse=error(b),frozen_blood_power_retention=power(b)/max(power(truth),1e-30),
        frozen_tissue_power_retention=power(t)/max(power(tissue),1e-30),
        frozen_noise_power_retention=power(n)/max(power(scene.noise[start:]),1e-30))
    metrics.update({k.replace('frozen_','actual_',1):v for k,v in blood_fidelity(actual.views,truth).items()})
    if actual.source_diagnostic is not None:
        diagnostic=actual.source_diagnostic(scene.blood)
        if diagnostic:metrics['axis2_basis_diagnostic']=diagnostic
    return metrics


def run_shared_case(config, condition, seed, split):
    params={**config['generator'],**{k:v for k,v in condition.items() if k!='id'},'seed':seed}
    scene=shared_scene(**params);start=config['score_start'];packet=config['packet_frames'];results=[];contexts={}
    for motion in dict.fromkeys(m['motion'] for m in config['methods']):
        known={'known_motion':scene.known_displacement} if motion=='known' else {}
        background=scene.tissue+scene.noise
        publics=[scene.public,replace(scene.public,views=background),replace(scene.public,views=background+.1*scene.blood)]
        fits=[SharedContext(p,packet,start,motion,motion_options=config.get('motion_options'),**known) for p in publics]
        contexts[motion]={name:_diagnostics(ctx,scene) for name,ctx in zip(('actual','background','weak_injection'),fits)}
        for spec in [m for m in config['methods'] if m['motion']==motion]:
            t0=time.perf_counter()
            actual,bg,weak=[fit_shared(c,spec) for c in fits]
            metrics=_evaluate_result(actual,bg,weak,scene,start)
            results.append(dict(id=spec['id'],family=spec['family'],status='ok',evaluation_only=spec['evaluation_only'],
                metrics=metrics,diagnostics=actual.diagnostics,
                counterfactual_blocks=dict(background=bg.diagnostics['blocks'],weak_injection=weak.diagnostics['blocks']),
                elapsed_s=time.perf_counter()-t0))
            del actual,bg,weak
        del fits
    order={m['id']:i for i,m in enumerate(config['methods'])};results.sort(key=lambda r:order[r['id']])
    return dict(condition=condition['id'],seed=seed,split=split,source_metadata=scene.metadata,
        input_sha256=array_hash(scene.public.views),calibration_sha256=array_hash(scene.public.calibration_views),
        blood_truth_sha256=array_hash(scene.blood),tissue_truth_sha256=array_hash(scene.tissue),
        contexts=contexts,results=results)


def _case_worker(config, condition, seed, split, code_hash, config_hash):
    if numerical_source_hash()!=code_hash:raise RuntimeError('Worker source differs from frozen parent')
    with threadpool_limits(limits=config.get('blas_threads',1)):
        record=run_shared_case(config,condition,seed,split)
    record.update(code_sha256=code_hash,config_sha256=config_hash)
    return _clean(record)


def run_shared_study(config, output):
    validate_shared_config(config);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    frozen=(json.dumps(config,indent=2)+'\n').encode();code_hash=numerical_source_hash()
    config_path=output/'frozen_config.json'
    if config_path.exists() and config_path.read_bytes()!=frozen:raise FileExistsError('Frozen study configuration changed')
    config_path.write_bytes(frozen);checkpoint=output/'checkpoints';checkpoint.mkdir(exist_ok=True)
    config_hash=hashlib.sha256(frozen).hexdigest();jobs=[];t0=time.perf_counter()
    for condition in config['conditions']:
        for split,key in [('training','training_seeds'),('evaluation','evaluation_seeds')]:
            for seed in config[key]:
                jobs.append((condition,seed,split,checkpoint/f'{condition["id"]}_{split}_{seed}.json.gz'))
    records=[None]*len(jobs);pending=[]
    for i,(condition,seed,split,path) in enumerate(jobs):
        if path.exists():
            with gzip.open(path,'rt') as f:record=json.load(f)
            if (record['code_sha256']!=code_hash or record['config_sha256']!=config_hash
                    or (record['condition'],record['seed'],record['split'])!=(condition['id'],seed,split)):
                raise ValueError('Checkpoint identity/source/config mismatch')
            records[i]=record
        else:pending.append(i)
    def accept(i, record):
        path=jobs[i][3];temporary=path.with_suffix('.tmp')
        temporary.write_bytes(gzip.compress((json.dumps(record,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        temporary.replace(path);records[i]=record
        print(json.dumps(dict(completed=sum(r is not None for r in records),condition=record['condition'],
            split=record['split'],seed=record['seed'],outputs=len(record['results']),
            elapsed_s=round(time.perf_counter()-t0,1))),flush=True)
    workers=config.get('case_workers',1)
    if workers<1:raise ValueError('Need at least one case worker')
    if workers==1:
        for i in pending:accept(i,_case_worker(config,*jobs[i][:3],code_hash,config_hash))
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures={pool.submit(_case_worker,config,*jobs[i][:3],code_hash,config_hash):i for i in pending}
            for future in as_completed(futures):accept(futures[future],future.result())
    if numerical_source_hash()!=code_hash:raise RuntimeError('Numerical source changed during study')
    report=dict(study=config['study'],version=config['version'],axis_names=AXIS_NAMES,
        code_sha256=code_hash,config_sha256=hashlib.sha256(frozen).hexdigest(),records=records,
        operating_points=select_global_operating_points(records,config['methods'],config['background_targets']),
        evaluation_seed_truth_used_for_selection=False,
        protocol=dict(primary_domain='uncompounded complex views, lab-frame amplitude and phase',
            background='independently refitted background-only output; not actual-mixture tissue attribution',
            injection='independently refitted finite incremental response; not additive attribution',
            known_motion_and_best_diagonal_weights='evaluation-only diagnostics, excluded from selection',
            frozen_sources='representation-dependent conditional diagnostics',
            calibration='same disjoint-time packet for all fits, with a shared tissue texture'),
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,elapsed_s=time.perf_counter()-t0,
                     blas_threads=config.get('blas_threads',1),case_workers=workers))
    report=_clean(report);raw=(json.dumps(report,separators=(',',':'),allow_nan=False)+'\n').encode()
    target=output/'results.json.gz'
    if target.exists():raise FileExistsError('Completed result already exists; choose a new output directory')
    target.write_bytes(gzip.compress(raw,mtime=0))
    return report
