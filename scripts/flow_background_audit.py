"""Retrospective blood-to-background power ratios from frozen project archives.

Never infer surviving blood power from an IQ error. No filter is refitted here.
Rows distinguish exact paired/source-frozen calculations from tissue-only,
incoherent noise proxies, and rigorous error-based intervals. Physical tube
contrast is exported separately because it has no isolated blood IQ truth.
"""
import argparse
import csv
import gzip
import hashlib
import json
import lzma
import math
from pathlib import Path

import numpy as np

from btbench.shared import array_hash, shared_scene


ROOT=Path(__file__).resolve().parents[1]
RUNS=ROOT/'runs'
FIELDS=['archive','source','split','seed','condition','method','metric_kind','signal',
        'blood_iq_error','blood_power_retention','background_power_remaining',
        'input_blood_to_tissue_power','sbr_db','sbr_lower_db','sbr_upper_db',
        'secondary_mixed_fit_proxy_db','local_sbr_median_db','local_sbr_p10_db',
        'local_fraction_blood_above_background','local_blood_support_pixels','note']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    path=Path(path)
    if not path.exists() and path.name=='results.json.gz':
        parts=sorted(path.parent.glob(path.name+'.part*'))
        if parts:return json.loads(gzip.decompress(b''.join(p.read_bytes() for p in parts)))
    if path.suffix=='.xz':return json.loads(lzma.decompress(path.read_bytes()))
    if path.suffix=='.gz':return json.loads(gzip.decompress(path.read_bytes()))
    return json.loads(path.read_text())


def ratio_db(numerator,denominator):
    if numerator is None or denominator is None or numerator<0 or denominator<0:
        return None
    if numerator==0 and denominator==0:return None
    if numerator==0:return float('-inf')
    if denominator==0:return float('inf')
    return 10.*math.log10(numerator/denominator)


def err_bounds(epsilon,rho,background):
    """Triangle inequality on ||B+error||, with no phase assumption."""
    if epsilon is None or rho is None or background is None:return None,None
    low=abs(1.-epsilon)**2
    high=(1.+epsilon)**2
    return ratio_db(low*rho,background),ratio_db(high*rho,background)


def add(rows,archive,source,split,seed,condition,method,kind,signal,
        error=None,blood=None,bg=None,rho=None,note='',lower=None,upper=None,proxy=None,
        local=None):
    sbr=ratio_db(blood*rho,bg) if blood is not None and rho is not None else None
    if sbr is None and lower is None and upper is None and error is not None:
        lower,upper=err_bounds(error,rho,bg)
    local=local or {}
    rows.append(dict(zip(FIELDS,[archive,str(source),split,seed,condition,method,kind,signal,
               error,blood,bg,rho,sbr,lower,upper,proxy,
               local.get('local_sbr_median_db'),local.get('local_sbr_p10_db'),
               local.get('local_fraction_blood_above_background'),
               local.get('local_blood_support_pixels'),note])))


def legacy(rows):
    sources=[*sorted(RUNS.glob('verified_seed*/results.json')),
        RUNS/'integration_v02_verified/results.json',
        RUNS/'integration_fidelity_v03/results.json.gz',
        RUNS/'oracle_objectives_check/results.json',
        *sorted(RUNS.glob('phase_diagram_validated/case_*/results.json')),
        RUNS/'svd_literature_v02_verified/results.json',
        RUNS/'stap_v02_verified/results.json',
        RUNS/'svd_fidelity_v03/results.json.gz',
        RUNS/'stap_fidelity_v03/results.json.gz',
        RUNS/'separation_fidelity_v03/results.json.gz']
    for source in sources:
        data=load(source)
        records=data.get('records')
        if records is None:
            records=[dict(seed=data.get('acquisition',{}).get('seed'),
                          condition=data.get('acquisition',{}).get('simulation_parameters',{}),
                          results=data['results'])]
        for case in records:
            condition=case.get('condition',case.get('case',case.get('scene_parameters',{})))
            if isinstance(condition,dict):
                condition=json.dumps(condition,sort_keys=True,separators=(',',':'))
            for record in case['results']:
                m=record.get('metrics',{})
                if record.get('status','ok')!='ok' or not m.get('frozen_component_metrics_available',True):
                    continue
                tissue=m.get('frozen_tissue_to_blood_power_ratio')
                noise=m.get('frozen_noise_to_blood_power_ratio')
                if tissue is None:continue
                domain=m.get('frozen_component_input_domain','archived native IQ/view representation')
                note=f'Component domain: {domain}. Isolated sources use one mixture-fitted operator.'
                if noise is None:
                    kind='frozen_tissue_only'
                    note+=' Filtered noise was not archived: ratio excludes noise.'
                    sbr=ratio_db(1.,tissue);lower=upper=None
                else:
                    kind='frozen_incoherent_proxy'
                    note+=' The tissue-noise cross term was not archived; bounds include all possible phases.'
                    sbr=ratio_db(1.,tissue+noise)
                    low_bg=(math.sqrt(tissue)-math.sqrt(noise))**2
                    high_bg=(math.sqrt(tissue)+math.sqrt(noise))**2
                    lower,upper=ratio_db(1.,high_bg),ratio_db(1.,low_bg)
                add(rows,source.parent.name,source.relative_to(ROOT),case.get('split',''),
                    case.get('seed'),condition,record.get('id',record.get('method')),
                    kind,'frozen_blood_in_'+domain,m.get('frozen_blood_iq_nrmse'),
                    m.get('frozen_blood_power_retention'),
                    m.get('frozen_tissue_power_retention'),None,note,lower,upper)
                rows[-1]['sbr_db']=sbr


def followup(rows):
    source=RUNS/'fidelity_followup_v04/results.json.xz'
    data=load(source)
    for case in data['records']:
        baseline=next(r for r in case['results'] if r['id']=='unfiltered')['metrics']
        rho=1./baseline['frozen_tissue_to_blood_power_ratio']
        for record in case['results']:
            if record['status']!='ok':continue
            m=record['metrics']
            for alpha,tag in [(1.,'alpha1'),(.1,'alpha0.1')]:
                suffix=f'_alpha{alpha:g}'
                q=m.get('refitted_blood_response_power_retention'+suffix)
                bg_to_blood=m.get('refitted_background_to_blood_power')
                if q is None or bg_to_blood is None:continue
                add(rows,'fidelity_followup_v04',source.relative_to(ROOT),case['split'],
                    case['seed'],case['condition'],record['id'],'exact_paired_'+tag,
                    'F(T+N+alpha B)-F(T+N), divided by alpha',
                    m.get('refitted_blood_response_nrmse'+suffix),q,
                    m.get('refitted_background_to_tissue_power'),rho,
                    'Counterfactual incremental response; background independently refitted.')
                direct=ratio_db(q,bg_to_blood)
                if direct is None:
                    assert rows[-1]['sbr_db'] is None
                else:
                    assert np.isclose(rows[-1]['sbr_db'],direct,atol=1e-7)


def shared(rows,name):
    source=RUNS/name/'results.json.gz';data=load(source)
    config=load(RUNS/name/'frozen_config.json')
    for case in data['records']:
        condition=next(v for v in config['conditions'] if v['id']==case['condition'])
        scene=shared_scene(**{**config['generator'],**{k:v for k,v in condition.items() if k!='id'},
                              'seed':case['seed']})
        assert array_hash(scene.blood)==case['blood_truth_sha256']
        assert array_hash(scene.tissue)==case['tissue_truth_sha256']
        start=config['score_start']
        rho=float(np.vdot(scene.blood[start:],scene.blood[start:]).real/
                  np.vdot(scene.tissue[start:],scene.tissue[start:]).real)
        for record in case['results']:
            if record['status']!='ok':continue
            m=record['metrics']
            for tag,qkey,ekey in [('alpha1','incremental_power_ratio_alpha1',
                                   'refitted_blood_response_nrmse_alpha1'),
                                  ('alpha0.1','incremental_power_ratio_alpha0_1',
                                   'refitted_blood_response_nrmse_alpha0_1')]:
                add(rows,name,source.relative_to(ROOT),case['split'],case['seed'],case['condition'],
                    record['id'],'exact_paired_'+tag,'F(T+N+alpha B)-F(T+N), divided by alpha',
                    m[ekey],m[qkey],m['refitted_background_to_tissue_power'],rho,
                    'Input blood/tissue powers regenerated with both archived component hashes verified.')


def axis1(rows,name):
    for split in ('seen_diagnostic','development','heldout'):
        source=RUNS/name/(split+'.json')
        if not source.exists():continue
        data=load(source)
        for case in data['cases']:
            rho=case['blood_input_power']/case['tissue_input_power']
            results=case['results']
            if isinstance(results,list):results={r['method']:r for r in results}
            for method,m in results.items():
                add(rows,name,source.relative_to(ROOT),split,case['seed'],case['condition'],
                    method,'exact_paired_alpha1','F(T+N+B)-F(T+N)',m['added_blood_iq_error'],
                    m['added_blood_power_retention'],m['refitted_background_power_remaining'],rho,
                    'Same scored frames, crop and coherent angle sum.')
                if 'same_motion_background_power_remaining' in m:
                    add(rows,name,source.relative_to(ROOT),split,case['seed'],case['condition'],
                        method,'bounded_same_motion','Only response error and background power archived',
                        m['same_motion_added_blood_iq_error'],None,
                        m['same_motion_background_power_remaining'],rho,
                        'Power of incremental response for same-motion fit not archived.')


def oracle_bounds(rows):
    name='axis1_oracle_limit_v1'
    for split in ('seen_diagnostic','development','heldout'):
        source=RUNS/name/(split+'.json')
        if not source.exists():continue
        data=load(source)
        for case in data['cases']:
            rho=case['blood_input_power']/case['tissue_input_power']
            for length in case['lengths']:
                condition=f"{case['condition']}:{case['flow_mps']:.5g}mps:{length['ensemble_frames']}frames"
                for method,m in length['controls'].items():
                    add(rows,name,source.relative_to(ROOT),split,case['seed'],condition,
                        method,'bounded_paired','Added blood IQ error only',m['added_blood_iq_error'],
                        None,m['background_power_remaining'],rho,
                        'Exact output blood power not archived; interval is from triangle inequality.')
                for rank,m in length['oracles'].items():
                    add(rows,name,source.relative_to(ROOT),split,case['seed'],condition,
                        'truth_oracle_rank'+rank+'_full','bounded_oracle','Truth-guided projector',
                        m['blood_iq_error_at_full_subtraction'],None,
                        m['background_power_at_full_subtraction'],rho,
                        'Diagnostic only. Exact output blood power not archived.')
                    gain=m['minimum_gain_for_tissue_limit']
                    if gain is not None:
                        add(rows,name,source.relative_to(ROOT),split,case['seed'],condition,
                            'truth_oracle_rank'+rank+'_at_limit','bounded_oracle','Truth-guided gain',
                            m['minimum_blood_error_at_tissue_limit'],None,
                            .1 if gain>0 else None,rho,
                            'At positive minimum gain background equals declared 0.10 target; no archived response power.')


def esdiva(rows):
    names=('esdiva_phase10_prospective','esdiva_phase10_reserve',
           'esdiva_phase10_stress','esdiva_phase11_prospective')
    exact_source=RUNS/'flow_background_audit_v1/esdiva_paired_exact.json'
    exact_doc=load(exact_source)
    exact={(r['bank'],r['seed'],r['scene'],r['variant'],r['method']):r
           for r in exact_doc['rows']}
    for name in names:
        for source in sorted((RUNS/name).glob('**/score.json')):
            score=load(source)
            if 'variants' not in score:continue
            seed,scene=score['seed'],score['scene']
            candidates=[ROOT.parent/'esdiva_phase11_work/acquisitions'/str(seed)/scene,
                        ROOT.parent/'esdiva_stress_work'/str(seed)/scene]
            matching=next((p for p in candidates if (p/'iq_sources.npz').exists()),None)
            rho=None
            if matching is not None:
                with np.load(matching/'iq_sources.npz') as f:
                    acquisition=load(matching/'report.json')
                    assert sha(matching/'report.json')==score['acquisition_sha256']
                    start=len(f['blood'])-score['score_frames']
                    assert start>=0 and len(f['blood'])==acquisition['evaluation_frames']
                    b=f['blood'][start:]
                    t=f['tissue'][start:]
                    rho=float(np.vdot(b,b).real/np.vdot(t,t).real)
            for variant,group in score['variants'].items():
                for method,m in group['methods'].items():
                    matched=exact.get((name,seed,scene,variant,method))
                    if matched:
                        assert np.isclose(matched['blood_iq_error'],m['paired_added_blood_iq_nrmse'])
                        assert np.isclose(matched['background_power_remaining'],
                                          m['refitted_background_output_to_input_tissue_power'])
                        if rho is not None:assert np.isclose(rho,matched['input_blood_to_tissue_power'])
                        add(rows,name,exact_source.relative_to(ROOT),source.parts[-4],seed,scene,
                            variant+'/'+method,'exact_paired_alpha1',
                            'F(T+N+B)-F(T+N)',matched['blood_iq_error'],
                            matched['blood_power_retention'],matched['background_power_remaining'],
                            matched['input_blood_to_tissue_power'],
                            'Replay of archived paired score, original source and score checked.')
                        continue
                    error=m['paired_added_blood_iq_nrmse'];bg=m['refitted_background_output_to_input_tissue_power']
                    proxy=ratio_db(m['frozen_blood_power_ratio']*rho,bg) if rho is not None else None
                    add(rows,name,source.relative_to(ROOT),source.parts[-4],seed,scene,
                        variant+'/'+method,'bounded_paired','F(T+N+B)-F(T+N)',
                        error,None,bg,rho,
                        'Pair-response power not archived; mixed-fit proxy uses frozen B numerator and refitted background denominator.'
                        if rho is not None else 'Score-domain source power unavailable in current workspace.',proxy=proxy)


def multidim(rows):
    source=RUNS/'flow_background_audit_v1/multidim_exact.json'
    data=load(source)
    for name,digest in data['input_sha256'].items():
        assert sha(RUNS/'multidim_aperture_v1'/f'{name}.json')==digest
    for m in data['rows']:
        add(rows,'multidim_aperture_v1',source.relative_to(ROOT),m['split'],m['seed'],
            m['scene']+(':phase_error' if m['phase_error'] else ':clean'),m['method'],
            'exact_frozen' if m['status'] in ('source_blind','calibrated') else 'exact_oracle',
            m['signal'],m['blood_iq_error'],m['blood_power_retention'],
            m['background_power_remaining'],m['input_blood_to_tissue_power'],
            'Source RF replay with frozen operator; same coherent output image.'
            if m['status'] in ('source_blind','calibrated') else
            'Evaluator-only tissue-truth projector/rank/gain; not deployable.',
            local=m)


def physical():
    values=[]
    for name in ('resolvestroke_single_view','resolvestroke_aperture_pilot_v1'):
        source=RUNS/name/'report.json';data=load(source)
        clips=data.get('clips',data.get('clip_scores',[]))
        for clip in clips:
            for method,scores in clip['scores'].items():
                for tube,entry in scores.items():
                    values.append({'archive':name,'clip':clip['clip'],'method':method,
                                   'tube':tube,'tube_to_outside_db':entry['mean_power_contrast_db'],
                                   'tube_auc':entry['auc'],'tube_pixels':entry['n_tube_pixels'],
                                   'note':'Anatomical tube/background contrast, not source-separated blood-to-clutter SBR.'})
    return values


def summary(rows):
    from collections import Counter,defaultdict
    counts=Counter(r['metric_kind'] for r in rows)
    cohorts=defaultdict(list)
    for row in rows:
        if row['sbr_db'] is not None and math.isfinite(row['sbr_db']):
            cohorts[row['archive']].append(row)
    return {'status_counts':dict(sorted(counts.items())),
            'row_count':len(rows),'by_archive':{k:{'rows_with_point_ratio':len(v)}
             for k,v in sorted(cohorts.items())}}


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--out',type=Path,default=RUNS/'flow_background_audit_v1')
    a=cli.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    legacy(rows);followup(rows)
    shared(rows,'shared_axes_v05');shared(rows,'motion_repair_v06')
    for name in ('axis1_motion_order_v1','axis1_native_subtraction_v1','axis1_stable_tissue_v1'):
        axis1(rows,name)
    oracle_bounds(rows);esdiva(rows);multidim(rows)
    with (a.out/'rows.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,FIELDS);writer.writeheader();writer.writerows(rows)
    (a.out/'rows.csv.xz').write_bytes(lzma.compress((a.out/'rows.csv').read_bytes(),preset=9))
    physical_rows=physical()
    with (a.out/'physical_tube_proxy.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,physical_rows[0].keys());writer.writeheader();writer.writerows(physical_rows)
    overview=summary(rows)
    sources=sorted({r['source'] for r in rows})
    def source_digest(relative):
        path=ROOT/relative
        if path.exists():return sha(path)
        parts=sorted(path.parent.glob(path.name+'.part*'))
        if not parts:raise FileNotFoundError(path)
        return hashlib.sha256(b''.join(p.read_bytes() for p in parts)).hexdigest()
    overview.update({'sources':{p:source_digest(p) for p in sources},
                     'analysis_script_sha256':sha(__file__),
                     'rows_csv_sha256':sha(a.out/'rows.csv'),
                     'rows_xz_sha256':sha(a.out/'rows.csv.xz'),
                     'physical_tube_proxy_sha256':sha(a.out/'physical_tube_proxy.csv'),
                     'physical_tube_proxy_rows':len(physical_rows)})
    (a.out/'summary.json').write_text(json.dumps(overview,indent=2,allow_nan=False)+'\n')
    print(json.dumps(overview,indent=2,allow_nan=False))


if __name__=='__main__':main()
