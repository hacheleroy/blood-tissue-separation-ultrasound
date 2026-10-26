"""Retrospective eSDIva angle/receive-aperture information test.

Usage (with pinned eSDIva installed):
  PYTHONPATH=src:examples:/path/to/eSDIva/src:/path/to/eSDIva/examples/example22_flow_doppler \
    python examples/multidim_aperture_study.py --split development --out OUT
  ... --split evaluation --out OUT

The split called evaluation was previously seen in phase 11 and is not a new
independent RF holdout. All oracle source fitting remains inside this scorer.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from functools import lru_cache

import numpy as np
from doppler_tools import rf2iq
from esdiva.utilities import add_noise,rf_rms

from btbench.esdiva_aperture import delayed_groups,read_rf_checkpoint
from btbench.multidim import (intersection_filter,minimum_gain,oracle_joint,
                               paper_tensor,svd_unfolding)
from btbench.stap import estimate_tissue_velocity,local_stap_views


ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/multidim_aperture_v1.json'
SOURCES=[Path(__file__),ROOT/'src/btbench/multidim.py',
         ROOT/'src/btbench/esdiva_aperture.py',ROOT/'src/btbench/stap.py',CONFIG]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arrhash(array):return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()
def rel(a,b):return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-30))


@lru_cache(maxsize=2)
def aperture_case(directory,groups):
    """RF-delayed group IQ and exact existing beamformed component controls."""
    import step1_define_flow_phantom as ex
    directory=Path(directory)
    info=json.loads((directory/'report.json').read_text())
    if info['upstream_commit']!='ab7eaef6d23428762ab5b9079ab03a5aa99af085':
        raise ValueError('Changed eSDIva upstream source')
    with np.load(directory/'broad_iq.npz') as f:
        saved={k:f[k] for k in f.files}
    rf_t,t0,dt=read_rf_checkpoint(directory/'rf_tissue',info['rf_shape'][-1])
    rf_b,t0_b,dt_b=read_rf_checkpoint(directory/'rf_blood',info['rf_shape'][-1])
    if (rf_t.shape!=tuple(info['rf_shape']) or dt!=dt_b or
            not np.array_equal(t0,t0_b) or arrhash(rf_t)!=info['rf_tissue_sha256'] or
            arrhash(rf_b)!=info['rf_blood_sha256']):
        raise ValueError('RF checkpoint identity differs from report')
    rf_clean=rf_t+rf_b
    rf_mix=add_noise(rf_clean,6.,reference=rf_rms(rf_clean),
                     rng=np.random.default_rng(info['seed']+101))
    rf_n=rf_mix-rf_clean
    if arrhash(rf_n)!=info['rf_noise_sha256']:
        raise ValueError('Original common receiver-noise realization failed replay')
    _,probe,_=ex.build_simulator()
    _,nine=ex.compound_events(probe,1)
    events=[nine[j] for j in info['view_indices']]
    grid=info['grid_mm']
    z_all=np.arange(grid['z_extent'][0],grid['z_extent'][1],grid['dz'])
    x_all=np.arange(grid['x_extent'][0],grid['x_extent'][1],grid['dx'])
    idx=np.array([np.argmin(abs(z_all-z)) for z in saved['z_mm']])
    if (not np.allclose(z_all[idx],saved['z_mm']) or
            not np.allclose(x_all,saved['x_mm'])):
        raise ValueError('Saved beam grid differs from RF delay grid')
    stop=info['calibration_frames']*3
    n=info['evaluation_frames']*3
    def convert(source):
        wave=np.stack([delayed_groups(source[j],t0[j],dt,events[j%3],
                       probe.element_centers,z_all,x_all,c=ex.C,groups=groups)
                       for j in range(stop,stop+n)])
        # The paper applies HOSVD after delay and before beam sum. RF->IQ is
        # linear and uses the exact same calibrated carrier as the stored IQ.
        image=rf2iq(wave.reshape(info['evaluation_frames'],3,groups,
                                    len(z_all),len(x_all)),
                    fs=ex.C/(2*grid['dz']*1e-3),fc=info['carrier_hz'],axis=-2)
        return image[:,:,:,idx,:].transpose(0,1,3,4,2).reshape(
            info['evaluation_frames'],3,-1,groups)
    t,b,noise=(convert(source) for source in (rf_t,rf_b,rf_n))
    core=saved['core_rows']
    ref={k:saved[k][:,:,core].reshape(info['evaluation_frames'],3,-1)
         for k in ('tissue','blood','noise','views')}
    # Select 8x16 core *after* RF-to-IQ on the full axial grid.
    p=len(x_all)
    chan={"tissue":t.reshape(len(t),3,len(idx),p,groups)[:,:,core].reshape(len(t),3,-1,groups),
          "blood":b.reshape(len(b),3,len(idx),p,groups)[:,:,core].reshape(len(b),3,-1,groups),
          "noise":noise.reshape(len(noise),3,len(idx),p,groups)[:,:,core].reshape(len(noise),3,-1,groups)}
    errors={name:rel(value.sum(-1),ref[name]) for name,value in chan.items()}
    if max(errors.values())>.02:
        raise ValueError(f'Delayed receive aperture does not reproduce stored IQ: {errors}')
    return chan,ref,saved,info,errors


def score_case(directory,config,phase_error=False):
    chan,ref,saved,info,errors=aperture_case(directory,config['receive_groups'])
    t,b,n=(chan[k].copy() for k in ('tissue','blood','noise'))
    if phase_error:
        time=np.linspace(0.,1.,len(t))
        angle_phase=np.array([0.,.12,-.12])[None,:,None,None]*time[:,None,None,None]
        aperture_phase=np.linspace(-.12,.12,t.shape[-1])[None,None,None,:]
        phase=np.exp(1j*(angle_phase+aperture_phase))
        t*=phase;b*=phase;n*=phase
    mix=t+b+n
    background=t+n
    # All operands come from the same acoustic RF events. The beam sum is
    # the one fixed target for every method, including the compounded SVD.
    def representations(y):
        angle=y.sum(-1).transpose(0,2,1)   # time,pixel,transmit angle
        receive=y.mean(1)                    # time,pixel,receive group
        compounded=receive.sum(-1)
        full=y.transpose(0,2,1,3)
        return {'compound':compounded,'angle':angle,'receive':receive,'full':full}
    m,ts,bs,ns,bg=(representations(y) for y in (mix,t,b,n,background))
    last=config['score_last_frames']
    def beamsum(y,rep):
        if rep=='compound':return y
        if rep=='angle':return y.mean(axis=-1)
        if rep=='receive':return y.sum(axis=-1)
        if rep=='full':return y.sum(axis=-1).mean(axis=-1)
        raise ValueError(rep)
    target=bs['compound'][-last:]
    den_b=max(float(np.vdot(target,target).real),1e-30)
    den_t=max(float(np.vdot(ts['compound'][-last:],ts['compound'][-last:]).real),1e-30)
    def metrics(output_b,output_bg,rep):
        ob=beamsum(output_b,rep)[-last:]
        ot=beamsum(output_bg,rep)[-last:]
        return {'blood_iq_error':float(np.linalg.norm(ob-target)/np.sqrt(den_b)),
                'background_power_remaining':float(np.vdot(ot,ot).real/den_t)}
    paper_ops={
        'compounded_svd_rank1':('compound',svd_unfolding(m['compound'])),
        'jiang_ga_svd_rank1':('angle',svd_unfolding(m['angle'])),
        'ozgun_svd_a_rank1':('receive',svd_unfolding(m['receive'])),
        'jiang_ac_hosvd_adaptation':('angle',paper_tensor(m['angle'],family='jiang')),
        'ozgun_hosvd_adaptation':('receive',paper_tensor(m['receive'],family='ozgun')),
        'angle_core_intersection':('angle',intersection_filter(m['angle'])),
        'receive_core_intersection':('receive',intersection_filter(m['receive']))}
    methods={}
    for name,(rep,op) in paper_ops.items():
        out=op.apply(m[rep]); out_b=op.apply(bs[rep]);out_bg=op.apply(bg[rep])
        if rel(out,out_b+out_bg)>1e-8:raise AssertionError('Frozen filter lost source additivity')
        methods[name]={**metrics(out_b,out_bg,rep),'settings':op.settings}
    try:
        geom=info['geometry']
        vel=estimate_tissue_velocity(saved['calibration_views'][:,:,saved['core_rows']],
            saved.get('training_roi',np.ones((8,16),bool)) if 'training_roi' in saved else
            np.ones((8,16),bool),geom['tx_rad'],geom['rx_rad'],geom['prf_hz'],
            geom['f0_hz'],geom['sound_speed_mps'],geom['phase_sign'])
        stap=local_stap_views(m['angle'].transpose(0,2,1).reshape(len(t),3,8,16),
            geom['prf_hz'],geom['f0_hz'],geom['tx_rad'],geom['rx_rad'],[vel],
            sound_speed_mps=geom['sound_speed_mps'],phase_sign=geom['phase_sign'],
            view_offsets_s=geom['view_offsets_s'],taps=config['stap_taps'],adaptive=False)
        def stap_apply(z):
            return stap.apply_frozen(z.transpose(0,2,1).reshape(len(t),3,8,16)).reshape(len(t),3,-1).transpose(0,2,1)
        methods['calibrated_stap_orthogonal']={**metrics(stap_apply(bs['angle']),stap_apply(bg['angle']),'angle'),
                                             'settings':{'calibration_frames':16,'taps':config['stap_taps'],
                                                         'tissue_velocity_mps':np.asarray(vel).tolist()}}
    except (ValueError,np.linalg.LinAlgError) as e:
        methods['calibrated_stap_orthogonal']={'failure':str(e)}
    oracle={}
    for rep in ('compound','angle','receive','full'):
        for rank in config['oracle_ranks']:
            remove,energy=oracle_joint(ts[rep],rank)
            rt=remove(ts[rep]);rb=remove(bs[rep]);rbg=remove(bg[rep])
            at=beamsum(rt,rep)[-last:]
            ab=beamsum(rb,rep)[-last:]
            abg=beamsum(rbg,rep)[-last:]
            score_bg=beamsum(bg[rep],rep)[-last:]
            score_t=beamsum(ts[rep],rep)[-last:]
            gain=minimum_gain(score_bg,abg,score_t,config['tissue_power_limit'])
            name=f'{rep}_rank{rank}'
            oracle[name]={'tissue_fitted_energy_fraction':energy,
                          'blood_iq_error_full':float(np.linalg.norm(ab)/np.sqrt(den_b)),
                          'background_power_full':float(np.linalg.norm(score_bg-abg)**2/den_t),
                          'gain_to_tissue_limit':gain,
                          'blood_error_at_tissue_limit':None if gain is None else
                          float(gain*np.linalg.norm(ab)/np.sqrt(den_b))}
    return {'seed':info['seed'],'scene':Path(directory).name,'phase_error':phase_error,
            'transmissions':info['evaluation_frames']*3,'duration_ms':1000*info['evaluation_frames']/info['geometry']['prf_hz'],
            'score_frames':last,'receive_groups':config['receive_groups'],
            'source_replay_errors':errors,'component_sha256':{k:arrhash(v) for k,v in chan.items()},
            'methods':methods,'oracle':oracle,
            'normalization':{'compounded_blood_power':den_b,'compounded_tissue_power':den_t}}


def summary(cases,ranks):
    names=cases[0]['methods']
    results={}
    for name in names:
        for phase in (False,True):
            vals=[c['methods'][name] for c in cases if c['phase_error']==phase]
            ok=[v for v in vals if 'failure' not in v]
            results[f'{name}:phase_error={phase}']={'count':len(ok),
               'median_blood_error':float(np.median([v['blood_iq_error'] for v in ok])) if ok else None,
               'median_background_power':float(np.median([v['background_power_remaining'] for v in ok])) if ok else None,
               'joint_target_count':sum(v['blood_iq_error']<=.25 and v['background_power_remaining']<=.1 for v in ok)}
    for rep in ('compound','angle','receive','full'):
        for rank in ranks:
            vals=[c['oracle'][f'{rep}_rank{rank}'] for c in cases if not c['phase_error']]
            ok=[v['blood_error_at_tissue_limit'] for v in vals if v['blood_error_at_tissue_limit'] is not None]
            results[f'oracle_{rep}_rank{rank}']={'tissue_limit_feasible':len(ok),
                  'median_blood_error_if_feasible':float(np.median(ok)) if ok else None,
                  'joint_target_count':sum(v<=.25 for v in ok)}
    return results


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--split',choices=('development','evaluation'),required=True)
    cli.add_argument('--out',type=Path,required=True)
    cli.add_argument('--rf-base',type=Path,default=ROOT.parent/'esdiva_phase11_work/acquisitions')
    a=cli.parse_args()
    cfg=json.loads(CONFIG.read_text())
    sources={str(p.relative_to(ROOT)):sha(p) for p in SOURCES}
    a.out.mkdir(parents=True,exist_ok=True)
    target=a.out/f'{a.split}.json'
    if target.exists():cli.error('Refusing to overwrite a scored split')
    if a.split=='evaluation':
        freeze=json.loads((a.out/'freeze.json').read_text())
        if sources!=freeze['sources'] or sha(a.out/'development.json')!=freeze['development_sha256']:
            cli.error('Freeze/source identity changed')
    cases=[]
    for seed in cfg[f'{a.split}_seeds']:
        for scene in cfg['scenes']:
            path=a.rf_base/str(seed)/scene
            for phase in (False,True):
                result=score_case(path,cfg,phase_error=phase)
                cases.append(result)
                print(seed,scene,'phase error',phase,flush=True)
    result={'scope':'retrospective eSDIva acoustic RF, previously seen phase-11 clouds; not new independent holdout',
            'paper_status':'correct published tensor/unfolding architectures; fixed cutoff adaptations, not exact author implementations',
            'sources':sources,'cases':cases,'summary':summary(cases,cfg['oracle_ranks'])}
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if a.split=='development':
        (a.out/'freeze.json').write_text(json.dumps({'sources':sources,
          'development_sha256':sha(target),'evaluation_seeds':cfg['evaluation_seeds'],
          'note':'Prior phase-11 evaluation seed; new methods not previously scored but RF is not independent'},indent=2)+'\n')
    print(json.dumps({'split':a.split,'cases':len(cases),'sha256':sha(target),'summary':result['summary']},indent=2))


if __name__=='__main__':main()
