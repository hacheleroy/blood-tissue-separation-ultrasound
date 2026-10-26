"""Read-only source replay for the frozen multidimensional RF contrast archive.

Run with the pinned eSDIva source and its dependencies on PYTHONPATH. This
adds powers missing from the original score without changing its frozen code.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from multidim_aperture_study import aperture_case
from btbench.multidim import intersection_filter, oracle_joint, paper_tensor, svd_unfolding
from btbench.stap import local_stap_views


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(x, y):
    return float(np.linalg.norm((x-y).ravel())/max(np.linalg.norm(y.ravel()), 1e-30))


def power(x):
    return float(np.vdot(x, x).real)


def replay_case(case, directory, config):
    channels, _, _, info, _ = aperture_case(directory, config['receive_groups'])
    assert case['seed'] == info['seed'] and case['scene'] == Path(directory).name
    from multidim_aperture_study import arrhash
    assert {key:arrhash(val) for key,val in channels.items()} == case['component_sha256']
    t, b, n = (channels[k].copy() for k in ('tissue', 'blood', 'noise'))
    if case['phase_error']:
        time = np.linspace(0., 1., len(t))
        angle = np.array([0., .12, -.12])[None, :, None, None]*time[:, None, None, None]
        aperture = np.linspace(-.12, .12, t.shape[-1])[None, None, None, :]
        phase = np.exp(1j*(angle+aperture))
        t *= phase; b *= phase; n *= phase

    def reps(y):
        angle = y.sum(-1).transpose(0, 2, 1)
        receive = y.mean(1)
        return {'compound':receive.sum(-1), 'angle':angle, 'receive':receive,
                'full':y.transpose(0, 2, 1, 3)}

    def image(y, rep):
        if rep == 'compound': return y
        if rep == 'angle': return y.mean(-1)
        if rep == 'receive': return y.sum(-1)
        if rep == 'full': return y.sum(-1).mean(-1)
        raise ValueError(rep)

    mixture, blood, tissue, background = (reps(y) for y in (t+b+n, b, t, t+n))
    last = config['score_last_frames']
    target = blood['compound'][-last:]
    bden, tden = power(target), power(tissue['compound'][-last:])
    assert np.isclose(bden, case['normalization']['compounded_blood_power'], rtol=1e-10)
    assert np.isclose(tden, case['normalization']['compounded_tissue_power'], rtol=1e-10)
    out = []

    def append(name, rep, ob, og, archived, status='source_blind', gain=None):
        ob, og = image(ob, rep)[-last:], image(og, rep)[-last:]
        bp = power(ob)/bden
        bg = power(og)/tden
        err = rel(ob, target)
        target_map=np.sum(abs(target)**2,axis=0)
        active=target_map>=.1*target_map.max()
        blood_map=np.sum(abs(ob)**2,axis=0)[active]
        background_map=np.sum(abs(og)**2,axis=0)[active]
        local_db=10*np.log10(np.maximum(blood_map,1e-300)/np.maximum(background_map,1e-300))
        if 'blood_iq_error' in archived:
            assert np.isclose(err, archived['blood_iq_error'], atol=2e-8), (name, err, archived)
            assert np.isclose(bg, archived['background_power_remaining'], atol=2e-8), (name, bg, archived)
        if 'blood_error_at_tissue_limit' in archived and name.endswith('_at_limit'):
            assert np.isclose(err, archived['blood_error_at_tissue_limit'], atol=2e-8), (name, err, archived)
        if 'blood_iq_error_full' in archived and name.endswith('_full'):
            assert np.isclose(err, archived['blood_iq_error_full'], atol=2e-8), (name, err, archived)
            assert np.isclose(bg, archived['background_power_full'], atol=2e-8), (name, bg, archived)
        out.append(dict(seed=case['seed'], scene=case['scene'], phase_error=case['phase_error'],
                        method=name, status=status, signal='frozen_blood' if status=='source_blind'
                        else 'truth_guided_blood', blood_iq_error=err,
                        blood_power_retention=bp, background_power_remaining=bg,
                        input_blood_to_tissue_power=bden/tden, gain=gain,
                        local_blood_support_pixels=int(active.sum()),
                        local_sbr_median_db=float(np.median(local_db)),
                        local_sbr_p10_db=float(np.quantile(local_db,.1)),
                        local_fraction_blood_above_background=float(np.mean(local_db>=0))))

    ops = {
        'compounded_svd_rank1':('compound', svd_unfolding(mixture['compound'])),
        'jiang_ga_svd_rank1':('angle', svd_unfolding(mixture['angle'])),
        'ozgun_svd_a_rank1':('receive', svd_unfolding(mixture['receive'])),
        'jiang_ac_hosvd_adaptation':('angle', paper_tensor(mixture['angle'],family='jiang')),
        'ozgun_hosvd_adaptation':('receive', paper_tensor(mixture['receive'],family='ozgun')),
        'angle_core_intersection':('angle', intersection_filter(mixture['angle'])),
        'receive_core_intersection':('receive', intersection_filter(mixture['receive']))}
    append('unfiltered','compound',blood['compound'],background['compound'],{})
    for name,(rep,op) in ops.items():
        append(name, rep, op.apply(blood[rep]), op.apply(background[rep]), case['methods'][name])

    stap_record = case['methods']['calibrated_stap_orthogonal']
    if 'failure' not in stap_record:
        geo = info['geometry']
        velocity = stap_record['settings']['tissue_velocity_mps']
        stap = local_stap_views(mixture['angle'].transpose(0,2,1).reshape(len(t),3,8,16),
            geo['prf_hz'],geo['f0_hz'],geo['tx_rad'],geo['rx_rad'],[velocity],
            sound_speed_mps=geo['sound_speed_mps'],phase_sign=geo['phase_sign'],
            view_offsets_s=geo['view_offsets_s'],taps=config['stap_taps'],adaptive=False)
        def apply(y):
            return stap.apply_frozen(y.transpose(0,2,1).reshape(len(t),3,8,16)).reshape(len(t),3,-1).transpose(0,2,1)
        append('calibrated_stap_orthogonal','angle',apply(blood['angle']),
               apply(background['angle']),stap_record)
        out[-1]['status'] = 'calibrated'

    for rep in ('compound','angle','receive','full'):
        for rank in config['oracle_ranks']:
            label=f'{rep}_rank{rank}'
            record=case['oracle'][label]
            removed, _ = oracle_joint(tissue[rep], rank)
            rb, rg = removed(blood[rep]), removed(background[rep])
            # Source-truth projector, both at full subtraction and its
            # truth-guided minimum gain that attains the tissue target.
            full={**record, 'background_power_full':record['background_power_full']}
            append('oracle_'+label+'_full',rep,blood[rep]-rb,background[rep]-rg,
                   full,'oracle',1.)
            gain=record['gain_to_tissue_limit']
            if gain is not None:
                append('oracle_'+label+'_at_limit',rep,blood[rep]-gain*rb,
                       background[rep]-gain*rg,record,'oracle',gain)
    return out


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--rf-base',type=Path,default=ROOT.parent/'esdiva_phase11_work/acquisitions')
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args()
    config=json.loads((ROOT/'configs/multidim_aperture_v1.json').read_text())
    archive=ROOT/'runs/multidim_aperture_v1'
    inputs={name:sha(archive/f'{name}.json') for name in ('development','evaluation','freeze')}
    rows=[]
    for split in ('development','evaluation'):
        data=json.loads((archive/f'{split}.json').read_text())
        for relative, digest in data['sources'].items():
            assert sha(ROOT/relative)==digest,(relative,'frozen numerical code changed')
        for case in data['cases']:
            result=replay_case(case,a.rf_base/str(case['seed'])/case['scene'],config)
            rows.extend(dict(split=split,**r) for r in result)
            print(split,case['seed'],case['scene'],'phase',case['phase_error'],len(result),flush=True)
    output={'status':'retrospective source-only metric; frozen filters and source hashes verified',
            'audit_script_sha256':sha(__file__),'input_sha256':inputs,'rows':rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print('rows',len(rows),'sha256',sha(a.out))


if __name__=='__main__':main()
