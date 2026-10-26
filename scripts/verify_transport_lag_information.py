"""Check the frozen transport-lag information archives without raw RF.

The acoustic RF/IQ checkpoints remain local. Their report and IQ hashes are
recorded in each row and may be checked on machines that have those files.
Clean CI checks source/code/freeze identity, matched emission budgets, all
declared cases, and feature-score arithmetic from stored pixel-level values.
"""
import hashlib
import json
from pathlib import Path
import math
import statistics
import importlib.util

ROOT=Path(__file__).resolve().parents[1]
RUNS=ROOT/'runs/transport_lag_information_v1'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def auc(pos,neg):
    # Direct pair counting, including ties, avoids a numerical dependency.
    return sum((1. if a>b else .5 if a==b else 0.) for a in pos for b in neg)/(len(pos)*len(neg))


def check_feature_physics():
    """A persistent echo has zero lag contrast; phase-transport is distinct."""
    import numpy as np
    path=ROOT/'scripts/transport_lag_information.py'
    spec=importlib.util.spec_from_file_location('transport_feature',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rng=np.random.default_rng(1307)
    base=rng.normal(size=(12,16))+1j*rng.normal(size=(12,16))
    still=np.broadcast_to(base,(80,12,16))
    phase=still*np.exp(1j*.08*np.arange(80))[:,None,None]
    noise=(rng.normal(size=phase.shape)+1j*rng.normal(size=phase.shape))
    def score(x):
        return np.abs(module.covariance(x,1))-module.covariance(x,48).real
    assert np.max(np.abs(score(still)))<1e-12
    assert np.median(score(phase))>1.5
    assert np.median(np.abs(score(noise)))<.3


def check(local=None):
    check_feature_physics()
    cfg=json.loads((ROOT/'configs/transport_lag_information_v1.json').read_text())
    frozen=json.loads((RUNS/'freeze.json').read_text())
    manifest=json.loads((RUNS/'manifest.json').read_text())
    assert manifest['status']=='transport-lag information archive published sources'
    for path,expected in manifest['sha256'].items():
        assert sha(ROOT/path)==expected,path
    assert sha(RUNS/'freeze_pre_noise_amendment.json')==frozen['initial_freeze_sha256']
    assert frozen['independent_RF_results_not_inspected'] is True
    paths={'config':ROOT/'configs/transport_lag_information_v1.json',
           'analysis':ROOT/'scripts/transport_lag_information.py',
           'acquisition':ROOT/'examples/transport_lag_acquire.py',
           'development':RUNS/'development.json','retrospective':RUNS/'retrospective.json'}
    for k,p in paths.items():assert sha(p)==frozen[k+'_sha256'],k
    all_cases=0
    for split in ('retrospective','development','evaluation'):
        d=json.loads((RUNS/f'{split}.json').read_text())
        assert d['status']==split
        for k in ('config','analysis','acquisition'):
            assert d[k+'_sha256']==frozen[k+'_sha256'],(split,k)
        expected={(seed,m,s) for seed in ([3907,3919] if split=='retrospective' else cfg['seeds'][split])
                  for m in cfg['motions'] for s in (['retrospective_dense'] if split=='retrospective' else cfg['schedules'])}
        got={(r['seed'],r['motion'],r['schedule']) for r in d['records']}
        assert got==expected and len(got)==len(d['records'])
        if split!='retrospective':
            for seed in cfg['seeds'][split]:
                for motion in cfg['motions']:
                    a,b=(next(r for r in d['records'] if r['seed']==seed and r['motion']==motion and
                              r['schedule']==s) for s in ('dense','burst'))
                    assert abs(a['common_noise_reference_rf_rms']-b['common_noise_reference_rf_rms'])<1e-15
                    assert a['separation']['weak_wall']['pixel_indices']==b['separation']['weak_wall']['pixel_indices']
        for r in d['records']:
            all_cases+=1
            assert r['transmissions']==240
            assert (19<r['duration_ms']<23 if r['schedule']!='burst' else 198<r['duration_ms']<202)
            assert all(len(v)==64 for v in (r['source_report_sha256'],r['source_iq_npz_sha256']))
            assert 0<r['lags_s']['1']<r['lags_s']['8']<r['lags_s']['48']
            assert r['separation']['near_wall_pixels']==16
            w=r['separation']['weak_wall']
            assert w['pixels']==8 and len(set(w['pixel_indices']))==8
            values=r['feature_values']
            for name,axes in values.items():
                assert len(axes['near_wall'])==16 and len(axes['weak_wall'])==8
                assert all(math.isfinite(x) for v in axes.values() for x in v)
            for suffix,axis in (('', 'near_wall'),('_weak','weak_wall')):
                metrics=w if suffix else r['separation']
                a=auc(values['blood'][axis],values['tissue'][axis])
                b=auc(values['blood'][axis],values['noise'][axis])
                assert abs(a-metrics['source_auc_blood_vs_tissue'])<1e-12
                assert abs(b-metrics['source_auc_blood_vs_noise'])<1e-12
                mix=values['mixture'][axis]
                bg=values['background'][axis]
                static=values['static_mixture'][axis]
                assert abs(statistics.median([a-b for a,b in zip(mix,bg)])-
                           metrics['mixture_minus_background_score_median'])<1e-12
                assert abs(statistics.median([a-b for a,b in zip(mix,static)])-
                           metrics['mixture_minus_static_blood_score_median'])<1e-12
            if local and split!='retrospective':
                root=Path(local)/str(r['seed'])/r['motion']/r['schedule']
                assert sha(root/'report.json')==r['source_report_sha256']
                assert sha(root/'iq_sources.npz')==r['source_iq_npz_sha256']
    print('Transport-lag information verified:',all_cases,'records; local IQ' if local else 'records; clean archive')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--local',type=Path)
    a=p.parse_args()
    check(a.local)
