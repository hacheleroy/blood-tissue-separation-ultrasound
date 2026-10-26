"""Retrospective source-separated receive-core and slow-flow profile diagnostic.

Run with pinned eSDIva and example22 on PYTHONPATH, e.g.
PYTHONPATH=src:examples:/path/to/eSDIva/src:/path/to/eSDIva/examples/example22_flow_doppler \
  python scripts/receive_core_source_profile.py --out runs/receive_core_source_profile_v1/results.json

All operators fit the observed mixture. Source truth and tube geometry enter
only the evaluator. Both clouds have been seen previously; this is not a new
holdout or a velocity-estimator development experiment.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from multidim_aperture_study import aperture_case, arrhash
from btbench.filters import svd_filter
from btbench.measurements import doppler_velocity
from btbench.multidim import factor_tensor, paper_tensor, svd_unfolding, tensor_gate
import step1_define_flow_phantom as phantom


ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/receive_core_source_profile_v1.json'
REFERENCE=ROOT/'runs/flow_background_audit_v1/multidim_exact.json'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def power(x):return float(np.vdot(x,x).real)
def db(a,b):return float(10*np.log10(max(a,1e-300)/max(b,1e-300)))


def source_powers(b,t,n):
    pb,pt,pn=map(power,(b,t,n))
    cross=float(2*np.vdot(t,n).real)
    bg=power(t+n)
    assert np.isclose(bg,pt+pn+cross,rtol=1e-12,atol=1e-9)
    return dict(blood=pb,tissue=pt,noise=pn,tissue_noise_cross=cross,background=bg,
                incoherent_noise_fraction=pn/max(pt+pn,1e-30),blood_background_db=db(pb,bg))


def geometry(saved,flow_scale,radial_bins):
    z,x=np.meshgrid(saved['z_mm'][saved['core_rows']],saved['x_mm'],indexing='ij')
    theta=phantom.THETA
    center=phantom.VESSEL_CENTER_MM
    signed=(x-center[0])*np.cos(theta)-(z-center[2])*np.sin(theta)
    along=(x-center[0])*np.sin(theta)+(z-center[2])*np.cos(theta)
    radial=np.abs(signed)/phantom.VESSEL_RADIUS_MM
    assert np.max(np.abs(along))<phantom.VESSEL_LENGTH_MM/2
    true_velocity=phantom.V_PEAK*flow_scale*np.cos(theta)*np.maximum(0,1-radial**2)
    regions={}
    for lo,hi,label in radial_bins:
        mask=(radial>=lo)&(radial<hi)
        assert np.count_nonzero(mask)>2,label
        regions[label]=mask.ravel()
    regions['whole_tile']=np.ones(radial.size,bool)
    regions['lumen']=(radial<1).ravel()
    return regions,true_velocity.ravel(),radial.ravel()


def velocity(iq,info):
    g=info['geometry']
    return doppler_velocity(iq,g['prf_hz'],g['f0_hz'],g['sound_speed_mps'],
                            g['phase_sign'],geometry_factor=2.)


def region_scores(images,source_images,regions,nominal_velocity,info,windows):
    b,t,n=images['blood'],images['tissue'],images['noise']
    raw_b=source_images['blood']
    scores={}
    for window in windows:
        bi,ti,ni=(y[-window:] for y in (b,t,n))
        ref_b=raw_b[-window:]
        input_t=source_images['tissue'][-window:]
        input_b=raw_b
        v_ref,_=velocity(ref_b,info)
        v_out,coherence=velocity((b+t+n)[-window:],info)
        v_bout,_=velocity(bi,info)
        out={}
        for label,mask in regions.items():
            p=source_powers(bi[:,mask],ti[:,mask],ni[:,mask])
            ib=power(input_b[-window:,mask]); it=power(input_t[:,mask])
            pixel_b=np.sum(abs(bi[:,mask])**2,axis=0)
            pixel_bg=np.sum(abs((ti+ni)[:,mask])**2,axis=0)
            local_db=10*np.log10(np.maximum(pixel_b,1e-300)/np.maximum(pixel_bg,1e-300))
            out[label]={**p,'pixels':int(np.count_nonzero(mask)),
                'blood_power_retention':p['blood']/max(ib,1e-30),
                'background_remaining_over_input_tissue':p['background']/max(it,1e-30),
                'blood_iq_error':float(np.linalg.norm(bi[:,mask]-ref_b[:,mask])/
                                         max(np.linalg.norm(ref_b[:,mask]),1e-30)),
                'local_sbr_median_db':float(np.median(local_db)),
                'local_sbr_p10_db':float(np.quantile(local_db,.1)),
                'fraction_pixels_blood_above_background':float(np.mean(local_db>=0)),
                'nominal_velocity_mps':float(np.median(nominal_velocity[mask])),
                'blood_only_velocity_mps':float(np.median(v_ref[mask])),
                'filtered_blood_velocity_mps':float(np.median(v_bout[mask])),
                'observed_velocity_mps':float(np.median(v_out[mask])),
                'velocity_error_vs_blood_only_mps':float(np.median(abs(v_out[mask]-v_ref[mask]))),
                'filtered_blood_velocity_error_mps':float(np.median(abs(v_bout[mask]-v_ref[mask]))),
                'velocity_error_vs_nominal_mps':float(np.median(abs(v_out[mask]-nominal_velocity[mask]))),
                'velocity_coherence_median':float(np.median(coherence[mask])),
                'positive_velocity_fraction':float(np.mean(v_out[mask]>0))}
        scores[str(window)]=out
    return scores


def mode_scores(mixture,sources,bins):
    """Orthogonal full-ensemble source energies; no source used to fit a basis.

    Both sets of energies are before final receive summation. The separate
    post-sum image powers in region_scores include all coherent cross terms.
    """
    u_time,u_space,u_recv=factor_tensor(mixture)
    tensor={}
    for name,y in sources.items():
        c=np.einsum('ti,pj,gk,tpg->ijk',u_time.conj(),u_space.conj(),
                    u_recv.conj(),y,optimize=True)
        e=abs(c)**2
        total=power(y)
        sectors={
            'first_time_first_receive':float(e[0,:,0].sum()),
            'first_time_other_receive':float(e[0,:,1:].sum()),
            'other_time_first_receive':float(e[1:,:,0].sum()),
            'other_time_other_receive':float(e[1:,:,1:].sum())}
        assert np.isclose(sum(sectors.values()),total,rtol=1e-10)
        tensor[name]={k:v/total for k,v in sectors.items()}
    mat=mixture.sum(-1)
    temporal,_,_=np.linalg.svd(mat,full_matrices=False)
    modes={}
    for name,y in sources.items():
        coeff=temporal.conj().T@y.sum(-1)
        e=np.sum(abs(coeff)**2,axis=1)
        assert np.isclose(e.sum(),power(y.sum(-1)),rtol=1e-10)
        modes[name]=[float(v) for v in e]
    bands=[]
    for lo,hi in bins:
        assert 0<=lo<hi<=len(temporal)
        energies={k:float(sum(e[lo:hi])) for k,e in modes.items()}
        bands.append({'rank_from_one':lo+1,'rank_to_one':hi,'power':energies,
                      'blood_to_noise_db':db(energies['blood'],energies['noise'])})
    return (u_time,u_space,u_recv),tensor,dict(per_mode=modes,bands=bands)


def process(case,config,reference,rf_base):
    seed,scene=case['seed'],case['scene']
    directory=rf_base/str(seed)/scene
    ch,_,saved,info,closure=aperture_case(directory,config['receive_groups'])
    assert not config['phase_error'] and info['seed']==seed and info['flow_scale']==.02
    assert {name:arrhash(value) for name,value in ch.items()}==case['component_sha256']
    assert list(ch.values())[0].shape==(80,3,128,8)
    receive={name:y.mean(1) for name,y in ch.items()}
    compound={name:y.sum(-1) for name,y in receive.items()}
    rx_mix=sum(receive.values()); cmp_mix=sum(compound.values())
    basis,tensor,svd_modes=mode_scores(rx_mix,receive,config['svd_mode_bins'])
    keep=np.ones(rx_mix.shape,bool);keep[0,:,0]=False
    receive_core=tensor_gate(rx_mix,keep,basis)
    bar=svd_filter(cmp_mix,selector='baranger')
    cut=bar.diagnostics['windows'][0]
    assert cut['upper']<=len(cmp_mix)
    bar_low=svd_filter(cmp_mix,selector='fixed',discard=cut['discard'],upper=len(cmp_mix))
    methods={
        'unfiltered':('compound',None),
        'compounded_svd_rank1':('compound',svd_unfolding(cmp_mix)),
        'receive_svd_a_rank1':('receive',svd_unfolding(rx_mix)),
        'receive_core_intersection':('receive',receive_core),
        'receive_hosvd_adaptation':('receive',paper_tensor(rx_mix,family='ozgun')),
        'baranger_lower_only':('compound',bar_low),
        'baranger_double':('compound',bar)}
    regions,nominal,radial=geometry(saved,info['flow_scale'],config['radial_bins'])
    output=[]
    for name in config['methods']:
        domain,op=methods[name]
        input_source=receive if domain=='receive' else compound
        apply=(op.apply if hasattr(op,'apply') else op.apply_frozen) if op is not None else lambda v:v
        image={k:apply(v) for k,v in input_source.items()}
        mixture_output=apply(sum(input_source.values()))
        assert np.linalg.norm(mixture_output-sum(image.values())) <= (
            1e-9*max(np.linalg.norm(mixture_output),1e-30)),name
        if domain=='receive':image={k:v.sum(-1) for k,v in image.items()}
        scored=region_scores(image,compound,regions,nominal,info,config['score_windows'])
        if name in reference:
            a=reference[name];m=scored['16']['whole_tile']
            for archived,now in [('blood_power_retention','blood_power_retention'),
                                 ('background_power_remaining','background_remaining_over_input_tissue'),
                                 ('blood_iq_error','blood_iq_error')]:
                assert np.isclose(m[now],a[archived],atol=2e-8),(seed,scene,name,archived)
        output.append({'method':name,'domain':domain,'scores':scored})
    return {'seed':seed,'scene':scene,'phase_error':False,'flow_scale':info['flow_scale'],
            'source_sha256':{'report':sha(directory/'report.json'),
                             'components':case['component_sha256']},
            'iq_channel_closure':closure,
            'radial_bin_counts':{k:int(sum(v)) for k,v in regions.items()},
            'nominal_axial_peak_mps':float(phantom.V_PEAK*info['flow_scale']*np.cos(phantom.THETA)),
            'baranger_boundaries_zero_based':{'lower':cut['discard'],'upper':cut['upper'],
                                              'modes':len(cmp_mix)},
            'tensor_mode_fraction':tensor,'compounded_svd_modes':svd_modes,
            'methods':output}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rf-base',type=Path,default=ROOT.parent/'esdiva_phase11_work/acquisitions')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    config=json.loads(CONFIG.read_text())
    baseline=json.loads((ROOT/'configs/multidim_aperture_v1.json').read_text())
    assert config['receive_groups']==baseline['receive_groups']==8
    assert config['score_windows']==[baseline['score_last_frames'],80]
    previous=json.loads(REFERENCE.read_text())
    ref={(r['seed'],r['scene'],r['phase_error']):{} for r in previous['rows']}
    for r in previous['rows']:
        ref[(r['seed'],r['scene'],r['phase_error'])][r['method']]=r
    cases=[];archive_hash={}
    for split,seed in [('development',3907),('evaluation',3919)]:
        path=ROOT/'runs/multidim_aperture_v1'/f'{split}.json'
        archive_hash[split]=sha(path)
        d=json.loads(path.read_text())
        for relative,digest in d['sources'].items():assert sha(ROOT/relative)==digest
        for scene in config['scenes']:
            case=next(x for x in d['cases'] if x['seed']==seed and x['scene']==scene and not x['phase_error'])
            cases.append(process(case,config,ref[(seed,scene,False)],a.rf_base))
            print(seed,scene,'verified',flush=True)
            aperture_case.cache_clear()
    result={'status':config['source_status'],'config_sha256':sha(CONFIG),
            'script_sha256':sha(__file__),'reference_sha256':sha(REFERENCE),
            'multidim_source_archive_sha256':archive_hash,
            'upstream_commit':'ab7eaef6d23428762ab5b9079ab03a5aa99af085',
            'cases':cases}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('cases',len(cases),'sha256',sha(a.out))


if __name__=='__main__':main()
