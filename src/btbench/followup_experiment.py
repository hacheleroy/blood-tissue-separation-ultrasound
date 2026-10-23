"""Frozen operating-point calibration, held-out source tests and real-IQ transfer.

Source truth selects operating points only on the declared calibration seeds.
This is a calibration-assisted research comparison, not a data-only selector.
Every adaptive transform is fitted to its evaluation mixture without truth.
"""
from dataclasses import replace
from pathlib import Path
import copy
import gzip
import hashlib
import json
import platform
import time
import numpy as np
from threadpoolctl import threadpool_limits
from .benchmark import run_filter
from .core import Acquisition
from .measurements import frozen_component_metrics,relative_norm,blood_fidelity
from .synthetic import poiseuille,advected_scatterers
from .stap_experiment import _clean


def numerical_source_hash():
    h=hashlib.sha256()
    for p in sorted(Path(__file__).parent.glob('*.py')):
        h.update(p.name.encode());h.update(p.read_bytes())
    return h.hexdigest()


def default_followup_config():
    methods=[]
    def add(name,family,method,params=None):
        methods.append({'id':name,'family':family,'method':method,'params':params or {}})
    add('unfiltered','identity','identity')
    for window,tag in [(32,'short'),(256,'long')]:
        for centered in (False,True):
            for discard in [0,1,4,16,min(window,96)]:
                add(f'svd_{tag}_c{int(centered)}_k{discard}',f'svd_{tag}','svd',
                    {'window':window,'centered':centered,'discard':discard})
            for discard in [1,4,16]:
                add(f'local_svd_{tag}_c{int(centered)}_k{discard}',f'local_svd_{tag}','local',
                    {'method':'svd','window':window,'centered':centered,'discard':discard,
                     'patch_shape':[4,6],'stride':[4,6]})
    for memory in (False,True):
        for rigid in (False,True):
            for centered in (False,True):
                family='axis1_'+('rigid_' if rigid else '')+('memory' if memory else 'no_memory')
                for discard in [1,4,16]:
                    add(f'{family}_c{int(centered)}_k{discard}',family,'axis1',
                        {'packet':32,'rho':.9 if memory else 0.,'centered':centered,'rigid':rigid,
                         'discard':discard,'q_min':.3,'patch_shape':[4,6],'stride':[4,6]})
    for centered in (False,True):
        add(f'axis1_no_phase_c{int(centered)}','axis1_no_phase','axis1',
            {'packet':32,'rho':.9,'centered':centered,'discard':1,'q_min':.3,
             'patch_shape':[4,6],'stride':[4,6],'phase_compensation':False})
        base={'packet':32,'centered':centered,'patch_shape':[4,6],'stride':[4,6],
              'rank':6,'lags':[1,2,4,8,16],'q_min':.3}
        for memory in (False,True):
            add(f'axis2_original_c{int(centered)}_m{int(memory)}','axis2_original','axis2',
                {**base,'use_memory':memory})
        add(f'axis2_cap_c{int(centered)}','axis2_gain_capped','axis2',
            {**base,'use_memory':True,'operator_gain_cap':1.})
        add(f'axis2_ridge_c{int(centered)}','axis2_ridge','axis2',
            {**base,'use_memory':True,'whitening_regularization_rel':1e-3})
        for cutoff in [1e-4,1e-2]:
            add(f'axis2_ridge_cap_c{int(centered)}_e{cutoff:g}','axis2_ridge_capped','axis2',
                {**base,'use_memory':True,'epsilon_rel':cutoff,
                 'whitening_regularization_rel':1e-3,'operator_gain_cap':1.})
    for domain in ['temporal','spatial']:
        for window in [32,256]:
            for eps in [1e-6,1e-3,.1]:
                add(f'lee_{domain}_w{window}_relative{eps:g}',f'lee_{domain}_raw_w{window}',
                    'lee_aer',{'domain':domain,'window':window,'epsilon_rel':eps})
                if window==256:
                    add(f'lee_{domain}_unit_relative{eps:g}',f'lee_{domain}_unit_extension',
                        'lee_aer',{'domain':domain,'window':window,'epsilon_rel':eps,
                                   'gain_normalization':'unit_operator'})
        for eps in [1.,26.,1000.]:
            add(f'lee_{domain}_absolute{eps:g}',f'lee_{domain}_raw_w256','lee_aer',
                {'domain':domain,'epsilon_abs':eps})
    for window in [32,128,256]:
        add(f'baranger_w{window}','baranger','svd',{'selector':'baranger','window':window})
    for offset in [0,8,26]:
        add(f'maresca_offset{offset}','maresca','sliding',{'method':'svd','window':60,'hop':45,
            'selector':'maresca','selector_options':{'offset':offset}})
    for normalization in ['paper','none']:
        add('song_pdf_'+normalization,'song_pdf_'+normalization,'song_pdf_local',
            {'window':60,'patch_shape':[8,8],'stride':[4,4],'normalization':normalization,
             'selector_options':{'prf_hz':1000.,'tissue_cutoff_hz':10.,'slope_threshold':.1,
                                 'noise_tail_fraction':.25,'tail_tolerance_decades':.05}})
    for cutoff in [10.,30.,80.]:
        add(f'butterworth_{cutoff:g}Hz','butterworth','highpass',{'cutoff_hz':cutoff,'prf_hz':1000.})
    for window in [32,256]:
        add(f'svd_spatial_equivalence_w{window}','svd_spatial_equivalence','svd',
            {'window':window,'discard':1,'projection_side':'spatial'})
    for domain in ['temporal','spatial']:
        add(f'aer_common_equivalence_{domain}','aer_common_equivalence','aer',
            {'centered':True,'domain':domain,'normalization':'common','epsilon_rel':1e-3})
    return {'version':4,'date':'2026-10-23','blas_threads':1,
        'calibration_seeds':[11,23],'evaluation_seeds':[61,73,89],
        'generator':{'n_frames':256,'shape':[8,12],'prf_hz':1000.,'f0_hz':5e6},
        'conditions':[
            {'label':'slow_overlap','generator':'statistical','umax_mps':.003,'tissue_velocity_mps':.002,
             'tissue_acceleration_mps2':0.,'clutter_db':40.,'decorrelation':.995},
            {'label':'accelerating','generator':'statistical','umax_mps':.025,'tissue_acceleration_mps2':.03,'clutter_db':30.},
            {'label':'texture_jump','generator':'statistical','umax_mps':.006,'motion_event':True,'clutter_db':30.},
            {'label':'opposed_flow','generator':'statistical','umax_mps':-.006,'tissue_velocity_mps':.002,
             'tissue_acceleration_mps2':0.,'clutter_db':30.},
            {'label':'translated_texture','generator':'translated','umax_mps':.006,'clutter_db':30.,'translation_amplitude_pixels':2},
            {'label':'psf_deforming','generator':'psf','umax_mps':.003,'tissue_velocity_mps':.001,
             'strain_motion_m':.00015,'lateral_motion_m':.00035,'motion_frequency_hz':8.,'n_blood':192,'n_tissue':192},
            {'label':'equal_static_failure_control','generator':'psf','umax_mps':0.,'tissue_velocity_mps':0.,'uniform_flow':True,
             'n_blood':192,'n_tissue':192}],
        'methods':methods,'matched_tissue_retention_targets':[1e-2,1e-4,1e-6],
        'selection':'minimum mean calibration output-to-blood IQ error subject to every calibration refitted background/tissue power <= target',
        'source_truth_use':'calibration-assisted operating points; evaluation truth never selects settings',
        'primary_rejection_metric':'refitted_background_to_tissue_power',
        'primary_loss_metric':'blood_iq_nrmse',
        'refitted_injection_amplitudes':[.1,1.],
        'frozen_attribution':'secondary, representation-dependent diagnostic; left/right equivalent outputs can have different component attributions',
        'comparison_scope':'declared discrete grids, two calibration seeds and three held-out seeds per condition; no universal ranking',
        'phase_control':'uncentered Axis 1 is exactly invariant to global phase demodulation/restoration; spatial registration is separate'}


def make_followup_scene(config,condition,seed):
    params={**config['generator'],**{k:v for k,v in condition.items() if k not in ('label','generator')},'seed':seed}
    kind=condition['generator']
    if kind=='translated':
        amplitude=int(params.pop('translation_amplitude_pixels'))
        shape=tuple(params.pop('shape'));pad=amplitude+1
        source=poiseuille(**params,shape=(shape[0]+2*pad,shape[1]+2*pad))
        n=len(source.iq);phase=2*np.pi*np.arange(n)/n
        shifts=np.rint(amplitude*np.column_stack((np.sin(phase),np.sin(2*phase)))).astype(int)
        def crop(value):
            return np.stack([value[i,pad+dz:pad+dz+shape[0],pad+dx:pad+dx+shape[1]] for i,(dz,dx) in enumerate(shifts)])
        blood,tissue,iq=crop(source.blood_true),crop(source.tissue_true),crop(source.iq)
        flow=np.any(crop(np.broadcast_to(source.flow_mask,source.iq.shape)),axis=0)
        a=Acquisition(iq,source.prf_hz,source.f0_hz,blood_true=blood,tissue_true=tissue,flow_mask=flow,
            z_m=source.z_m[pad:pad+shape[0]],x_m=source.x_m[pad:pad+shape[1]],
            metadata={**source.metadata,'generator':'padded integer-translated statistical IQ',
                      'evaluation_only_crop_shifts_pixels':shifts.tolist(),'boundary_model':'crop from padded scene; no wrapping'})
    elif kind=='psf':a=advected_scatterers(**params)
    else:a=poiseuille(**params)
    return replace(a,tissue_roi=None)


def compact_diagnostics(diag):
    d=copy.deepcopy(diag)
    for key in ['motion_phase_rad','motion_shifts_pixels']:
        if key in d:
            a=np.asarray(d.pop(key));d[key+'_sha256']=hashlib.sha256(a.tobytes()).hexdigest()
            d[key+'_range']=[float(a.min()),float(a.max())]
    return d


def calibrate_operating_points(records,methods,targets,
        rejection_key='refitted_background_to_tissue_power',loss_key='blood_iq_nrmse'):
    """Select using calibration rows only; never interpolate unmatched settings."""
    families=sorted({s['family'] for s in methods});selections=[]
    conditions=sorted({r['condition'] for r in records})
    for condition in conditions:
        calibration=[r for r in records if r['condition']==condition and r['split']=='calibration']
        evaluation=[r for r in records if r['condition']==condition and r['split']=='evaluation']
        for family in families:
            ids=[s['id'] for s in methods if s['family']==family]
            for target in targets:
                candidates=[]
                for mid in ids:
                    rows=[next(x for x in r['results'] if x['id']==mid) for r in calibration]
                    if not rows or any(r['status']!='ok' for r in rows):continue
                    m=[r['metrics'] for r in rows]
                    if any(v.get(rejection_key) is None or v.get(loss_key) is None for v in m):continue
                    worst=max(v[rejection_key] for v in m)
                    if worst<=target:
                        candidates.append((float(np.mean([v[loss_key] for v in m])),mid,worst))
                row={'condition':condition,'family':family,'target':target,'evaluation_only_truth_matching':False,
                     'calibration_assisted':True,'calibration_candidates':len(candidates),
                     'rejection_metric':rejection_key,'loss_metric':loss_key}
                if not candidates:
                    selections.append({**row,'status':'infeasible_on_calibration','selected_id':None,'evaluation':[]});continue
                loss,mid,worst=min(candidates)
                held=[]
                for case in evaluation:
                    measured=next(x for x in case['results'] if x['id']==mid)
                    metrics=measured.get('metrics',{})
                    held.append({'seed':case['seed'],'status':measured['status'],'metrics':metrics,
                                 'meets_target':measured['status']=='ok' and metrics.get(rejection_key,np.inf)<=target})
                selections.append({**row,'status':'selected_on_calibration','selected_id':mid,
                    'calibration_mean_loss':loss,'calibration_max_rejection_metric':worst,'evaluation':held})
    return selections


def refitted_response_metrics(a,spec,mixture_result,amplitudes=(.1,1.)):
    """Representation-independent counterfactuals, not additive attribution.

    F(T+N) is a separately fitted background-only run. [F(T+N+alpha B)-
    F(T+N)]/alpha is an incremental response. None is labeled actual tissue
    leakage inside F(X). The output-to-blood error directly scores F(X).
    """
    background=a.iq-a.blood_true
    base=run_filter(background,spec,tissue_roi=None,prf_hz=a.prf_hz).iq
    power=lambda v:float(np.sum(abs(v)**2))
    metrics={'refitted_background_to_tissue_power':power(base)/max(power(a.tissue_true),1e-30),
             'refitted_background_to_blood_power':power(base)/max(power(a.blood_true),1e-30)}
    for alpha in amplitudes:
        injected=mixture_result.iq if alpha==1 else run_filter(background+alpha*a.blood_true,spec,tissue_roi=None,prf_hz=a.prf_hz).iq
        response=(injected-base)/alpha;tag=f'_alpha{alpha:g}'
        metrics['refitted_blood_response_nrmse'+tag]=relative_norm(response-a.blood_true,a.blood_true)
        metrics['refitted_blood_response_power_retention'+tag]=power(response)/max(power(a.blood_true),1e-30)
        if alpha==1:
            companion=blood_fidelity(response,a.blood_true,a.flow_mask)
            metrics.update({k.replace('frozen_','refitted_response_'):v for k,v in companion.items()})
    return metrics


def run_followup_study(config,output):
    if set(config['calibration_seeds'])&set(config['evaluation_seeds']):
        raise ValueError('Calibration and evaluation seeds must be disjoint')
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose a new follow-up output directory')
    out.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(config,sort_keys=True,indent=2,allow_nan=False)+'\n'
    (out/'frozen_config.json').write_text(encoded)
    code=numerical_source_hash();records=[];began=time.perf_counter()
    with threadpool_limits(limits=config.get('blas_threads',1)):
        for condition in config['conditions']:
            for split,seeds in [('calibration',config['calibration_seeds']),('evaluation',config['evaluation_seeds'])]:
                for seed in seeds:
                    a=make_followup_scene(config,condition,seed);rows=[]
                    for spec in config['methods']:
                        start=time.perf_counter()
                        row={'id':spec['id'],'family':spec['family']}
                        try:
                            result=run_filter(a.iq,spec,tissue_roi=None,prf_hz=a.prf_hz)
                            metrics=frozen_component_metrics(result,a)
                            metrics['blood_iq_nrmse']=relative_norm(result.iq-a.blood_true,a.blood_true)
                            metrics.update(refitted_response_metrics(a,spec,result,config['refitted_injection_amplitudes']))
                            row.update(status='ok',metrics=metrics,diagnostics=compact_diagnostics(result.diagnostics))
                        except Exception as e:
                            row.update(status='error',error=f'{type(e).__name__}: {e}')
                        row['wall_time_s']=time.perf_counter()-start;rows.append(row)
                    records.append({'condition':condition['label'],'split':split,'seed':seed,
                        'input_sha256':hashlib.sha256(a.iq.tobytes()).hexdigest(),'code_sha256':code,'results':rows})
                    (out/'checkpoint.json.gz').write_bytes(gzip.compress(json.dumps(_clean(records),separators=(',',':')).encode(),mtime=0))
                    print(json.dumps({'completed_cases':len(records),'condition':condition['label'],'split':split,'seed':seed,
                                      'errors':sum(r['status']!='ok' for r in rows)}),flush=True)
    selections=calibrate_operating_points(records,config['methods'],config['matched_tissue_retention_targets'])
    frozen_selections=calibrate_operating_points(records,config['methods'],config['matched_tissue_retention_targets'],
        rejection_key='frozen_tissue_power_retention',loss_key='frozen_blood_iq_nrmse')
    report={'code_sha256':code,'config_sha256':hashlib.sha256(encoded.encode()).hexdigest(),'records':records,
        'matched_operating_points':selections,'parameters_frozen_before_execution':True,
        'secondary_representation_dependent_frozen_matching':frozen_selections,
        'calibration_seed_truth_used_for_operating_points':True,'evaluation_seed_truth_used_for_selection':False,
        'runtime':{'python':platform.python_version(),'numpy':np.__version__,'blas_threads':config.get('blas_threads',1),
                   'elapsed_s':time.perf_counter()-began},'scope':config['comparison_scope']}
    (out/'results.json.gz').write_bytes(gzip.compress(json.dumps(_clean(report),separators=(',',':'),allow_nan=False).encode(),mtime=0))
    (out/'checkpoint.json.gz').unlink()
    return report


def default_openhrf_followup_config():
    """Bounded 1024-frame, two-subject protocol fixed before new RF processing."""
    methods=[{'id':'unfiltered','method':'identity'}]
    for window in [32,256,768]:
        for k in [1,4]:
            methods.append({'id':f'svd_w{window}_k{k}','method':'svd','params':{'window':window,'discard':k}})
    methods.append({'id':'baranger_w256','method':'svd','params':{'selector':'baranger','window':256}})
    methods.append({'id':'maresca_k26_w60','method':'sliding','params':{'method':'svd','window':60,'hop':45,
        'selector':'maresca','selector_options':{'offset':26}}})
    for normalization in ['paper','none']:
        methods.append({'id':'song_pdf_'+normalization,'method':'song_pdf_local','params':{'window':60,
            'patch_shape':[8,8],'stride':[8,8],'normalization':normalization,
            'selector_options':{'tissue_cutoff_hz':10.,'slope_threshold':.1,'noise_tail_fraction':.25,'tail_tolerance_decades':.05}}})
    for domain in ['temporal','spatial']:
        for window in [256,768]:
            methods.append({'id':f'lee_{domain}_w{window}','method':'lee_aer',
                'params':{'domain':domain,'window':window,'epsilon_rel':1e-3}})
        methods.append({'id':f'lee_{domain}_absolute26','method':'lee_aer',
            'params':{'domain':domain,'epsilon_abs':26.}})
        methods.append({'id':f'lee_{domain}_unit_extension','method':'lee_aer',
            'params':{'domain':domain,'epsilon_rel':1e-3,'gain_normalization':'unit_operator'}})
    for memory,rigid in [(False,False),(True,False),(True,True)]:
        methods.append({'id':f'axis1_memory{int(memory)}_rigid{int(rigid)}','method':'axis1',
            'params':{'packet':32,'rho':.9 if memory else 0.,'rigid':rigid,'patch_shape':[6,6],'stride':[6,6],'q_min':.3}})
    for centered in [False,True]:
        for guarded in [False,True]:
            params={'packet':32,'centered':centered,'rank':6,'lags':[1,2,4,8,16],
                    'use_memory':True,'patch_shape':[6,6],'stride':[6,6],'q_min':.3}
            if guarded:params.update(epsilon_rel=1e-3,whitening_regularization_rel=1e-3,operator_gain_cap=1.)
            methods.append({'id':f'axis2_c{int(centered)}_guard{int(guarded)}','method':'axis2','params':params})
    methods.append({'id':'view_svd_matched','method':'view_svd','params':{'window':4,'hop':2}})
    for adaptive in [False,True]:
        methods.append({'id':'stap_covariance' if adaptive else 'stap_orthogonal','method':'stap',
                        'params':{'taps':4,'hop':2,'adaptive':adaptive,'training_radius':3,'guard':1}})
    return {'version':4,'date':'2026-10-23','blas_threads':1,
        'dataset_revision':'05234aa084b27c8e4d0096c3a9adc33fbbc123d8',
        'acquisitions':[{'name':'waterloo_acq0','file':'Acq0.hdf5','subject_id':'1'},
                        {'name':'waterloo_acq12','file':'Acq12.hdf5','subject_id':'3'}],
        'frames_start_stop':[0,1024],'frame_batch_size':16,'reconstruction_grid_stride':2,
        'calibration_frames':256,'evaluation_frames':768,'evaluation_spatial_stride':2,
        'methods':methods,'measurement_window':64,'measurement_hop':64,
        'operating_point_selection':'fixed a priori; no experimental output used to choose a method or parameter',
        'calibration':'first 256 frames fit a nuisance ROI and bulk velocity; subsequent 768 frames evaluate',
        'scope':'two short CCA baseline records from different subjects; no independent blood/tissue or phase truth; no cardiac-cycle generalization',
        'contrast_status':'unconfirmed; native blood inferred from provider description'}


def run_openhrf_followup(acquisitions,config,output):
    from .benchmark import run_benchmark
    from .stap import estimate_tissue_velocity
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose a new experimental output directory')
    out.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(config,sort_keys=True,indent=2,allow_nan=False)+'\n'
    (out/'frozen_config.json').write_text(encoded)
    records=[];seen=set()
    for source,a in zip(config['acquisitions'],acquisitions,strict=True):
        subject=a.metadata.get('subject_id')
        if subject!=source['subject_id'] or subject in seen:
            raise ValueError('Experimental subjects must match the declared distinct-subject split')
        seen.add(subject)
        ncal=config['calibration_frames'];stride=config['evaluation_spatial_stride']
        views=a.iq_views[:,:,::stride,::stride]
        if len(views)!=ncal+config['evaluation_frames']:raise ValueError('Unexpected experimental window length')
        cal=views[:ncal];test=views[ncal:]
        power=np.mean(abs(cal)**2,axis=(0,1));roi=power>=np.quantile(power,.75)
        v=estimate_tissue_velocity(cal,roi,a.tx_angles_rad,a.rx_angles_rad,a.prf_hz,a.f0_hz,a.sound_speed_mps,a.phase_sign)
        clock=np.arange(len(a.iq))/a.prf_hz if a.timestamps_s is None else a.timestamps_s
        meta={**a.metadata,'independent_truth_available':False,'calibration_frames':ncal,
              'evaluation_frames':len(test),'estimated_tissue_velocity_zx_mps':v.tolist(),
              'calibration_sha256':hashlib.sha256(cal.tobytes()).hexdigest(),
              'evaluation_time_origin_s':float(clock[ncal]),
              'calibration_roi':'brightest calibration-power quartile; not verified tissue-only'}
        evaluation=replace(a,iq=test.mean(axis=1),iq_views=test,timestamps_s=clock[ncal:],
            z_m=a.z_m[::stride],x_m=a.x_m[::stride],tissue_roi=roi,flow_mask=None,near_wall_mask=None,
            blood_true=None,tissue_true=None,blood_views_true=None,tissue_views_true=None,velocity_true_mps=None,wall=None,metadata=meta)
        cfg=copy.deepcopy(config)
        for method in cfg['methods']:
            if method['method']=='stap':
                method['params']['tissue_velocities_zx_mps']=[v.tolist()]
                if method['params']['adaptive']:method['params']['training_roi']=roi.tolist()
        case=out/source['name'];case.mkdir()
        (case/'frozen_config.json').write_text(json.dumps(_clean(cfg),sort_keys=True,indent=2,allow_nan=False)+'\n')
        report=run_benchmark(evaluation,cfg,case)
        records.append({'acquisition':source['name'],'subject_id':subject,'code_sha256':report['code_sha256'],
            'input_sha256':report['input_sha256'],'acquisition_metadata':meta,'runtime':report['runtime'],'results':report['results']})
        print(json.dumps({'experimental_acquisition_completed':source['name'],'methods':len(report['results']),
                          'ok':sum(r['status']=='ok' for r in report['results'])}),flush=True)
    report={'code_sha256':numerical_source_hash(),'config_sha256':hashlib.sha256(encoded.encode()).hexdigest(),
            'records':records,'parameters_frozen_before_execution':True,'scope':config['scope']}
    (out/'results.json.gz').write_bytes(gzip.compress(json.dumps(_clean(report),separators=(',',':'),allow_nan=False).encode(),mtime=0))
    return report
