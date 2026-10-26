"""Source/mixture information test in tissue-frame displacement–lag space.

The source-separated arrays and tube coordinates are used only in the scorer.
The preprocessing estimates B-mode translation and a global complex-IQ phase
from the *observed mixture*. There is no source-derived registration, tuning,
component selection, or fitted blood/tissue label. The fixed score is a
diagnostic feature, not a deployable clutter filter.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import shift
from scipy.optimize import minimize
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/transport_lag_information_v1.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arrhash(x):return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
def power(x):return float(np.vdot(x,x).real)


def auc(positive,negative):
    a,b=np.asarray(positive),np.asarray(negative)
    ranks=rankdata(np.r_[a,b],method='average')
    return float((ranks[:len(a)].sum()-len(a)*(len(a)+1)/2)/(len(a)*len(b)))


def align(mix,sources):
    """Fit each frame to the first from B-mode, phase-lock on its complex IQ.

    Entire imaged field participates; there is no known tube mask or tissue
    source in the estimator. Small translations are fitted on magnitude IQ
    to avoid conflating carrier phase with geometry.
    """
    ref=mix[0]
    ref_mag=abs(ref)
    ref_weight=ref_mag/max(np.linalg.norm(ref_mag),1e-30)
    shifts=[];phasors=[];out={k:[] for k in sources}
    for t in range(len(mix)):
        if t==0:pos=np.zeros(2)
        else:
            im=abs(mix[t]);im=im/max(np.linalg.norm(im),1e-30)
            def cost(v):
                z=shift(im,v,order=1,mode='nearest',prefilter=False)
                return float(np.mean((z-ref_weight)**2))
            fit=minimize(cost,np.zeros(2),method='Powell',bounds=((-1,1),(-1,1)),
                         options={'xtol':1e-3,'ftol':1e-6,'maxiter':50})
            pos=np.asarray(fit.x)
        shifted=shift(mix[t],pos,order=1,mode='nearest',prefilter=False)
        phase=np.exp(-1j*np.angle(np.vdot(ref,shifted)))
        shifts.append(pos);phasors.append(phase)
        for name,x in sources.items():
            out[name].append(shift(x[t],pos,order=1,mode='nearest',prefilter=False)*phase)
    return {k:np.stack(v) for k,v in out.items()},np.asarray(shifts),np.asarray(phasors)


def covariance(x,lag,disp_pixels=(0.,0.),pair_mask=None):
    a=x[:-lag]
    b=x[lag:]
    if pair_mask is not None:
        a=a[pair_mask];b=b[pair_mask]
    if np.linalg.norm(disp_pixels)>0:
        b=shift(b,(0.,-disp_pixels[0],-disp_pixels[1]),order=1,mode='nearest',prefilter=False)
    numerator=np.sum(a.conj()*b,axis=0)
    den=np.sqrt(np.sum(abs(a)**2,axis=0)*np.sum(abs(b)**2,axis=0))
    c=numerator/np.maximum(den,1e-30)
    assert np.max(abs(c))<1.0000001
    return c


def masks(z,x,core):
    # Only evaluator uses the vessel radius and inclination.
    zz,xx=np.meshgrid(z,x,indexing='ij')
    theta=np.deg2rad(60)
    r=abs(xx*np.cos(theta)-(zz-16)*np.sin(theta))/2
    wall=(r>=.75)&(r<1)
    outer=(r>=1)&(r<1.25)
    assert int(wall[core].sum())==16
    return wall,outer,r


def case(directory,config,legacy=False):
    directory=Path(directory)
    report=json.loads((directory/'report.json').read_text())
    iq_file=directory/('broad_iq.npz' if legacy else 'iq_sources.npz')
    with np.load(iq_file) as data:
        saved={k:data[k] for k in data.files}
    if legacy:
        for label,key in [('views','mixture'),('tissue','tissue'),('blood','blood'),('noise','noise')]:
            assert arrhash(saved[label][:,:,saved['core_rows']])==report[key+'_sha256']
        saved['cycle_times_s']=np.arange(80)/3700.
    else:
        for k,v in report['iq_sha256'].items():assert arrhash(saved[k])==v,k
        assert arrhash(saved['cycle_times_s'])==report['cycle_times_sha256']
        assert arrhash(saved['event_times_s'])==report['event_times_sha256']
        assert report['transmissions']==240 and report['cycles']==80
        assert report['script_sha256']==digest(ROOT/'examples/transport_lag_acquire.py')
        assert report['noise_reference_rule'].startswith('first 24 identical events')
    sources={k:saved[k].mean(1) for k in ('tissue','blood','noise')}
    sources['background']=sources['tissue']+sources['noise']
    sources['mixture']=saved['views'].mean(1)
    # Evaluation-only IQ counterfactual: same first-emission blood scatterers
    # repeated at the same locations. This is a model-derived static-blood
    # control, not a measured physical flow-off acquisition.
    static_mix=sources['background']+np.repeat(sources['blood'][:1],80,axis=0)
    assert np.linalg.norm(sources['mixture']-sources['background']-sources['blood']) < 2e-4*max(np.linalg.norm(sources['mixture']),1e-30)
    aligned,positions,phasors=align(sources['mixture'],sources)
    static_aligned,_,_=align(static_mix,{'static_mixture':static_mix})
    aligned.update(static_aligned)
    assert np.linalg.norm(aligned['mixture']-aligned['background']-aligned['blood'])<2e-4*max(np.linalg.norm(aligned['mixture']),1e-30)
    core=saved['core_rows'];wall,outer,radial=masks(saved['z_mm'],saved['x_mm'],core)
    # Evaluator-only, matched-across-schedules weak subset: the first eight
    # cycles coincide in emission times and scatterer positions, so rank their
    # blood/tissue power and keep the weaker half of the 16 wall pixels.
    # The thermal noise draws differ between schedules, so selection excludes
    # noise to ensure the anatomical pixel set is exactly paired.
    early_b=np.sum(abs(sources['blood'][:8,core])**2,axis=0)
    early_t=np.sum(abs(sources['tissue'][:8,core])**2,axis=0)
    wall_f=wall[core]
    weak_flat=np.flatnonzero(wall_f.ravel())[
        np.argsort((early_b/np.maximum(early_t,1e-30)).ravel()[wall_f.ravel()],kind='stable')[:8]]
    weak=np.zeros(wall_f.size,bool);weak[weak_flat]=True;weak=weak.reshape(wall_f.shape)
    physical=saved['cycle_times_s']
    lag_s={str(i):float(np.median(physical[i:]-physical[:-i])) for i in config['lags']}
    # Exclude the nine gaps between bursts from the short lag-one feature.
    close=(np.diff(physical)<2/3700.)
    cov={name:{str(lag):covariance(x,lag,pair_mask=close if lag==1 else None) for lag in config['lags']}
         for name,x in aligned.items()}
    short,long=map(str,config['feature_lags'])
    # A high short-lag magnitude rejects incoherent noise; a reduced real
    # long-lag covariance signals motion relative to the tissue frame.
    score={k:(abs(c[short])-c[long].real) for k,c in cov.items()}
    summary={}
    for name in aligned:
        p=aligned[name][:,core][:,:, :]
        wall_p=p[:,wall_f]
        summary[name]={
            'near_wall_iq_power':power(wall_p),
            'score_median':float(np.median(score[name][core][wall_f])),
            'score_p10':float(np.quantile(score[name][core][wall_f],.1)),
            'score_p90':float(np.quantile(score[name][core][wall_f],.9)),
            'short_coherence_median':float(np.median(abs(cov[name][short][core][wall_f]))),
            'long_coherence_median':float(np.median(abs(cov[name][long][core][wall_f]))),
            'long_real_median':float(np.median(cov[name][long].real[core][wall_f])),
            'long_phase_median_rad':float(np.median(np.angle(cov[name][long][core][wall_f])))}
    blood,tissue,noise=[score[k][core][wall_f] for k in ('blood','tissue','noise')]
    mixture=score['mixture'][core][wall_f]
    background=score['background'][core][wall_f]
    geometry={'source_auc_blood_vs_tissue':auc(blood,tissue),
              'source_auc_blood_vs_noise':auc(blood,noise),
              'source_auc_blood_vs_background':auc(blood,background),
              'mixture_minus_background_score_median':float(np.median(mixture-background)),
              'mixture_over_background_score_positive_fraction':float(np.mean(mixture>background)),
              'mixture_minus_static_blood_score_median':float(np.median(mixture-score['static_mixture'][core][wall_f])),
              'mixture_over_static_blood_score_positive_fraction':float(np.mean(mixture>score['static_mixture'][core][wall_f])),
              'mixture_near_wall_vs_outer_auc':auc(mixture,score['mixture'][core][outer[core]]),
              'near_wall_pixels':int(wall_f.sum()),'outer_pixels':int(outer[core].sum()),
              'near_wall_blood_background_db':float(10*np.log10(power(aligned['blood'][:,core][:,wall_f])/
                        max(power(aligned['background'][:,core][:,wall_f]),1e-30)))}
    geometry['weak_wall']={
        'pixels':int(weak.sum()),'pixel_indices':weak_flat.tolist(),
        'selection':'lowest first-eight-cycle source blood/tissue power among fixed 16 wall pixels; evaluation only',
        'source_auc_blood_vs_tissue':auc(score['blood'][core][weak],score['tissue'][core][weak]),
        'source_auc_blood_vs_noise':auc(score['blood'][core][weak],score['noise'][core][weak]),
        'mixture_minus_background_score_median':float(np.median((score['mixture']-score['background'])[core][weak])),
        'mixture_over_background_score_positive_fraction':float(np.mean(score['mixture'][core][weak]>score['background'][core][weak])),
        'mixture_minus_static_blood_score_median':float(np.median((score['mixture']-score['static_mixture'])[core][weak])),
        'mixture_over_static_blood_score_positive_fraction':float(np.mean(score['mixture'][core][weak]>score['static_mixture'][core][weak])),
        'aligned_blood_background_db':float(10*np.log10(power(aligned['blood'][:,core][:,weak])/
                           max(power(aligned['background'][:,core][:,weak]),1e-30))),
        'blood_score_median':float(np.median(score['blood'][core][weak])),
        'tissue_score_median':float(np.median(score['tissue'][core][weak])),
        'mixture_score_median':float(np.median(score['mixture'][core][weak]))}
    # Physical-time lag-8 phase, never the ordinary lag-one estimator across
    # irregular burst gaps. Compare with blood-only IQ on the same schedule;
    # geometric velocity is a distinct reference owing to PSF and speckle.
    velocity_lag='8'
    v_ref=-np.angle(cov['blood'][velocity_lag])*1540/(4*np.pi*2555519.761864708*lag_s[velocity_lag])
    v_mix=-np.angle(cov['mixture'][velocity_lag])*1540/(4*np.pi*2555519.761864708*lag_s[velocity_lag])
    geometry['phase_velocity_lag8']={
        'blood_only_weak_median_mm_s':float(1e3*np.median(v_ref[core][weak])),
        'observed_weak_median_mm_s':float(1e3*np.median(v_mix[core][weak])),
        'weak_pixel_error_vs_blood_only_median_mm_s':float(1e3*np.median(abs((v_mix-v_ref)[core][weak]))),
        'wall_pixel_error_vs_blood_only_median_mm_s':float(1e3*np.median(abs((v_mix-v_ref)[core][wall_f]))),
        'nominal_weak_median_mm_s':float(1e3*np.median((.45*.02*np.cos(np.deg2rad(60))*(1-radial**2))[core][weak]))}
    # Optimistic evaluator-only directional test: does shifting the delayed
    # image along the *known vessel axis* expose extra blood correlation?
    # Direction and true speed never enter registration or the fixed score.
    long_lag=config['feature_lags'][1]
    dt=lag_s[str(long_lag)]
    direction=np.array([np.cos(np.deg2rad(60)),np.sin(np.deg2rad(60))])  # z,x
    shifted={}
    for name in ('tissue','blood','noise','mixture'):
        base=abs(cov[name][long][core][wall_f])
        candidates=[]
        for speed_mm_s in config['oracle_transport_speeds_mm_s']:
            dzdx=speed_mm_s*dt*direction/np.array([.1,1.])
            c=covariance(aligned[name],long_lag,dzdx)[core][wall_f]
            candidates.append(abs(c))
        best=np.max(candidates,axis=0)
        shifted[name]={'mean_gain_over_zero_shift':float(np.mean(best-base)),
                       'median_gain_over_zero_shift':float(np.median(best-base)),
                       'fraction_nonzero_speed_wins':float(np.mean(np.argmax(candidates,axis=0)>0))}
    return {'seed':report['seed'],'motion':report['motion'],'schedule':'retrospective_dense' if legacy else report['schedule'],
            'transmissions':240,'duration_ms':1000*80/3700 if legacy else report['duration_ms'],
            'source_report_sha256':digest(directory/'report.json'),
            'source_iq_npz_sha256':digest(iq_file),
            'common_noise_reference_rf_rms':None if legacy else report['noise_reference_rf_rms'],
            'lags_s':lag_s,'registration':{'estimated_shift_z_pixels_first_last':positions[[0,-1],0].tolist(),
                'estimated_shift_x_pixels_first_last':positions[[0,-1],1].tolist(),
                'phase_first_last_rad':np.angle(phasors[[0,-1]]).tolist(),
                'max_translation_pixels':float(np.max(np.linalg.norm(positions,axis=1)))},
            'sources':summary,'separation':geometry,'oracle_directional_shift':shifted,
            'feature_values':{name:{'near_wall':score[name][core][wall_f].tolist(),
                                    'weak_wall':score[name][core][weak].tolist()}
                              for name in aligned}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,default=ROOT.parent/'transport_lag_work')
    p.add_argument('--split',choices=('retrospective','development','evaluation'),required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    config=json.loads(CONFIG.read_text())
    archive=ROOT/'runs/transport_lag_information_v1'
    if a.split=='evaluation':
        frozen=json.loads((archive/'freeze.json').read_text())
        for label,path in [('config',CONFIG),('analysis',Path(__file__)),
                           ('acquisition',ROOT/'examples/transport_lag_acquire.py'),
                           ('development',archive/'development.json'),
                           ('retrospective',archive/'retrospective.json')]:
            if digest(path)!=frozen[label+'_sha256']:
                raise RuntimeError(f'Frozen {label} differs before independent evaluation')
    out=[]
    for seed in ([3907,3919] if a.split=='retrospective' else config['seeds'][a.split]):
        for motion in config['motions']:
            for schedule in (['rigid_slow' if motion=='rigid' else 'smooth_slow'] if a.split=='retrospective' else config['schedules']):
                directory=(ROOT.parent/'esdiva_phase11_work/acquisitions'/str(seed)/schedule
                           if a.split=='retrospective' else a.base/str(seed)/motion/schedule)
                item=case(directory,config,legacy=a.split=='retrospective')
                out.append(item)
                print(seed,motion,schedule,'AUC',item['separation']['source_auc_blood_vs_tissue'],flush=True)
    payload={'status':a.split,'config_sha256':digest(CONFIG),'analysis_sha256':digest(__file__),
             'acquisition_sha256':digest(ROOT/'examples/transport_lag_acquire.py'),
             'records':out}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if a.out.exists():raise FileExistsError(a.out)
    a.out.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
