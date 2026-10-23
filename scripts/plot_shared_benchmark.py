"""Regenerate shared-study scientific figures and transparent tabular summaries."""
import argparse
import csv
import gzip
import json
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


LABELS={'svd':'Joint-view SVD','weighted_svd':'Weighted-history SVD','axis1':'Axis 1',
        'axis2':'Axis 2','stap_orthogonal':'Axis 3 · orthogonal','stap_loaded':'Axis 3 · loaded',
        'lee_temporal_raw':'Lee AER · temporal','lee_spatial_raw':'Lee AER · spatial',
        'lee_unit_extension':'AER · unit extension','baranger':'Baranger SSM'}
FAMILIES=list(LABELS)
CONDITION_LABELS=['Slow\noverlap','Separated\nlow rank','Subpixel\ntranslation','Smooth\ndeformation',
                  'Unmodelled\ndeformation','Calibration\ndrift','Opposed\nflow','Identical\nstatistics']


def load_report(directory):
    directory=Path(directory)
    if (directory/'results.json.gz').exists():return json.loads(gzip.decompress((directory/'results.json.gz').read_bytes()))
    path=directory/'results.json.xz'
    data=path.read_bytes() if path.exists() else b''.join(p.read_bytes() for p in sorted(directory.glob('results.json.xz.part*')))
    return json.loads(lzma.decompress(data))


def summarize(directory, figure_directory=None):
    directory=Path(directory);report=load_report(directory);config=json.loads((directory/'frozen_config.json').read_text())
    figure_directory=Path(figure_directory or directory/'figures');figure_directory.mkdir(parents=True,exist_ok=True)
    records=report['records'];evaluation=[c for c in records if c['split']=='evaluation'];conditions=[c['id'] for c in config['conditions']]
    lookup={(c['condition'],c['seed'],r['id']):r for c in evaluation for r in c['results']}
    def values(mid,key,condition):return np.array([lookup[condition,s,mid]['metrics'][key] for s in config['evaluation_seeds']])
    def means(mid,key):return np.array([values(mid,key,c).mean() for c in conditions])
    def save(fig,name):
        fig.savefig(figure_directory/(name+'.png'),dpi=180,bbox_inches='tight')
        fig.savefig(figure_directory/(name+'.svg'),bbox_inches='tight',metadata={'Date':None})
        plt.close(fig)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'shared-benchmark-v05','axes.titleweight':'semibold'})
    selections=report['operating_points']
    points={ (s['family'],s['motion'],s['horizon_frames'],s['background_target']):s for s in selections }
    fig,axes=plt.subplots(1,2,figsize=(15.5,6.4),layout='constrained',sharey=True)
    for ax,target in zip(axes,[None,.01]):
        matrix=np.full((len(FAMILIES),len(conditions)),np.nan)
        for i,family in enumerate(FAMILIES):
            selection=points[(family,'none',64,target)]
            if selection['selected_id']:matrix[i]=means(selection['selected_id'],'blood_iq_nrmse')
        cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('#e3e6e8')
        im=ax.imshow(matrix,norm=LogNorm(vmin=.03,vmax=30),cmap=cmap,aspect='auto')
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                v=matrix[i,j]
                rgb=cmap(LogNorm(vmin=.03,vmax=30)(v))[:3] if np.isfinite(v) else (1,1,1)
                luminance=np.dot(rgb,[.299,.587,.114])
                ax.text(j,i,'—' if not np.isfinite(v) else f'{v:.2g}',ha='center',va='center',
                        color='white' if luminance<.45 else '#151515',fontsize=8)
        ax.set_xticks(range(len(conditions)),CONDITION_LABELS,fontsize=8)
        ax.set_yticks(range(len(FAMILIES)),[LABELS[f] for f in FAMILIES])
        ax.set_title('Training-best actual error' if target is None else 'Training-best subject to mean background ≤ 10⁻²')
    fig.colorbar(im,ax=axes,shrink=.78,label='Held-out actual complex IQ NRMSE · lower is better')
    fig.suptitle('Common three-view input · 64-frame history · no spatial compensation\n'
                 'Each cell: mean over three held-out seeds; one globally trained setting per family\n'
                 '— = training-infeasible · fixed-cutoff SVD selected empty output · colours saturate above error 30',fontsize=11)
    save(fig,'shared_01_operating_points')

    fig,axes=plt.subplots(1,3,figsize=(16,5.2),layout='constrained')
    for mid,label,style in [('svd_none_h16_k1','16-frame SVD','o-'),('svd_none_h64_k1','64-frame uniform SVD','s-'),
                            ('weighted_svd_none_h64_k1','64-frame weighted SVD','^-'),('axis1_none_h64_k1','64-frame Axis 1','x--')]:
        axes[0].plot(means(mid,'blood_iq_nrmse'),style,label=label,ms=5)
    axes[0].set_yscale('log');axes[0].set_title('A · Memory and its weighted-SVD control')
    axes[0].set_ylabel('Held-out actual IQ NRMSE');axes[0].legend(fontsize=8)
    for mode,label in [('none','None'),('rigid','Estimated rigid'),('smooth','Estimated smooth'),('known','Known warp · diagnostic')]:
        y=[np.mean([c['contexts'][mode]['actual']['evaluation_only_displacement_rmse_pixels'] for c in evaluation if c['condition']==cond]) for cond in conditions]
        axes[1].plot(y,'o-',label=label,ms=4)
        axes[2].plot(means(f'svd_{mode}_h64_k1','blood_iq_nrmse'),'o--' if mode=='known' else 'o-',label=label,ms=4)
    axes[1].set_title('B · Fitted displacement versus prescribed warp');axes[1].set_ylabel('Displacement RMSE · pixels');axes[1].legend(fontsize=8)
    axes[2].set_title('C · Fixed cutoff-1 SVD under the same motion arms');axes[2].set_ylabel('Held-out actual IQ NRMSE');axes[2].set_yscale('log')
    for ax in axes:
        ax.set_xticks(range(len(conditions)),[x.replace('\n',' ') for x in CONDITION_LABELS],rotation=55,ha='right',fontsize=8)
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('Motion and memory controls · means over three held-out seeds\nKnown coordinates retain interpolation and boundary losses',fontsize=12)
    save(fig,'shared_02_memory_motion')

    fig,axes=plt.subplots(1,2,figsize=(14.8,5.8),layout='constrained',sharey=True)
    for ax,motion in zip(axes,['none','smooth']):
        mid=f'axis2_{motion}_h64_r4_p0.1'
        for key,label,marker in [('blood_iq_nrmse','Actual capped heuristic','o-'),
                                 ('heuristic_uncapped_nrmse','Actual uncapped heuristic','s--'),
                                 ('best_diagonal_uncapped_nrmse','Truth-fitted diagonal weights · diagnostic','^-')]:
            samples=np.array([values(mid,key,c) if key=='blood_iq_nrmse' else
                              [lookup[c,s,mid]['metrics']['axis2_basis_diagnostic'][key] for s in config['evaluation_seeds']]
                              for c in conditions])
            y=samples.mean(1);ax.plot(y,marker,label=label,ms=5);ax.fill_between(range(len(y)),samples.min(1),samples.max(1),alpha=.12)
        ax.set_yscale('log');ax.set_title('No spatial compensation' if motion=='none' else 'Estimated smooth compensation')
        ax.set_xticks(range(len(conditions)),CONDITION_LABELS,fontsize=8);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Actual-output IQ NRMSE');axes[0].legend(fontsize=8)
    fig.suptitle('Axis 2: Motion-demodulated multi-lag separation\nRank 4, phase width 0.1; means and full ranges of three held-out seeds',fontsize=12)
    save(fig,'shared_03_axis2_basis_identification')

    fixed=[('svd_none_h64_k1','Joint-view SVD · k1'),('axis1_none_h64_k1','Axis 1 · k1'),
           ('stap_none_h64_orthogonal','Axis 3 · orthogonal'),('stap_none_h64_loaded0.01','Axis 3 · loading .01'),
           ('lee_aer_none_h64_temporal_abs26','Lee temporal AER · ε26'),('lee_aer_none_h64_spatial_abs26','Lee spatial AER · ε26')]
    fig,axes=plt.subplots(1,3,figsize=(15.2,4.8),layout='constrained',sharey=True)
    for ax,cond in zip(axes,['slow_overlap','separated_low_rank','calibration_drift']):
        for i,(mid,label) in enumerate(fixed):
            xv=values(mid,'refitted_background_to_tissue_power',cond);yv=values(mid,'blood_iq_nrmse',cond)
            ax.scatter(xv,yv,alpha=.3,s=18,color=f'C{i}')
            ax.scatter(xv.mean(),yv.mean(),marker=['o','x','s','^','D','v'][i],s=65,color=f'C{i}',label=label)
        ax.axhline(1.,color='.55',lw=.8,ls=':');ax.set_xscale('log');ax.set_yscale('log')
        ax.set_title(cond.replace('_',' ').capitalize());ax.set_xlabel('Refitted background / input tissue power');ax.grid(alpha=.15)
    axes[0].set_ylabel('Actual complex IQ NRMSE');axes[0].legend(fontsize=7,loc='best')
    fig.suptitle('Local ultrasound STAP and Stephen Lee AER · matched 64-frame support and calibration\n'
                 'Fixed declared settings, no held-out tuning; small dots = seeds, large symbols = means',fontsize=12)
    save(fig,'shared_04_stap_lee')

    # CSV exports contain every setting/case and every selected global recipe.
    flat=[]
    for case in records:
        for row in case['results']:
            item=dict(condition=case['condition'],split=case['split'],seed=case['seed'],id=row['id'],family=row['family'],evaluation_only=row['evaluation_only'])
            item.update({k:v for k,v in row['metrics'].items() if not isinstance(v,dict)})
            flat.append(item)
    def csv_write(path,rows):
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=sorted({k for r in rows for k in r}));writer.writeheader();writer.writerows(rows)
    metric_directory=directory/'metrics';metric_directory.mkdir(exist_ok=True)
    for condition in conditions:
        csv_write(metric_directory/(condition+'.csv'),[r for r in flat if r['condition']==condition])
    selected=[]
    for s in selections:
        rows=s['evaluation'];m=[r['metrics'] for r in rows]
        selected.append(dict(family=s['family'],motion=s['motion'],horizon_frames=s['horizon_frames'],
            background_target='unconstrained' if s['background_target'] is None else s['background_target'],
            status=s['status'],selected_id=s['selected_id'],held_out_cases=len(rows),
            held_out_mean_iq_nrmse=float(np.mean([r['blood_iq_nrmse'] for r in m])) if m else None,
            held_out_mean_background=float(np.mean([r['refitted_background_to_tissue_power'] for r in m])) if m else None,
            held_out_target_pass_fraction=float(np.mean([r['meets_target'] for r in rows])) if rows and s['background_target'] is not None else None))
    csv_write(directory/'operating_points.csv',selected);csv_write(directory/'method_ledger.csv',config['methods'])
    summaries=[]
    for motion in ('none','smooth'):
        for family in FAMILIES:
            for target in (None,.01,.0001):
                s=points[family,motion,64,target]
                for cond in conditions:
                    rows=[r for r in s['evaluation'] if r['condition']==cond]
                    summaries.append(dict(family=family,motion=motion,target=target,condition=cond,selected_id=s['selected_id'],status=s['status'],
                        iq_nrmse_mean=float(np.mean([r['metrics']['blood_iq_nrmse'] for r in rows])) if rows else None,
                        background_mean=float(np.mean([r['metrics']['refitted_background_to_tissue_power'] for r in rows])) if rows else None,
                        incremental_nrmse_alpha01_mean=float(np.mean([r['metrics']['refitted_blood_response_nrmse_alpha0_1'] for r in rows])) if rows else None))
    equivalence=[]
    for case in records:
        rr={r['id']:r for r in case['results']}
        for mid,r in rr.items():
            if mid.startswith('axis1_') and mid.replace('axis1_','weighted_svd_',1) in rr:
                other=rr[mid.replace('axis1_','weighted_svd_',1)]
                equivalence.append(max(abs(r['metrics'][key]-other['metrics'][key]) for key in
                    ['blood_iq_nrmse','refitted_background_to_tissue_power','refitted_blood_response_nrmse_alpha1','refitted_blood_response_nrmse_alpha0_1']))
    output_summary=dict(study=report['study'],code_sha256=report['code_sha256'],cases=len(records),
        outputs=len(flat),training_cases=sum(c['split']=='training' for c in records),evaluation_cases=len(evaluation),
        max_frozen_closure=max(r['metrics']['frozen_component_closure_nrmse'] for c in records for r in c['results']),
        max_axis1_weighted_svd_primary_metric_difference=max(equivalence),
        infeasible_operating_points=sum(s['status']!='selected' for s in selections),
        selected=selected,condition_summaries=summaries)
    (directory/'summary.json').write_text(json.dumps(output_summary,indent=2)+'\n')
    lines=['# Shared three-axis benchmark v05','',
        f"{len(records)} synthetic cases, {len(flat)} case-setting evaluations, three fits per setting. {len(evaluation)} held-out cases.",
        f"Frozen source: `{report['code_sha256']}`.",'',
        'Global operating points pool all training conditions. Means below pool the 24 held-out cases; a mean training target is not a per-case guarantee. Known-motion arms and truth-fitted Axis 2 weights are excluded from selection.','',
        '| Family | Motion | History | Training target | Selected setting | Held-out mean IQ error | Mean background | Pass fraction |',
        '|---|---|---:|---:|---|---:|---:|---:|']
    for r in selected:
        value=lambda k:'—' if r[k] is None else f'{r[k]:.4g}'
        lines.append(f"| {r['family']} | {r['motion']} | {r['horizon_frames']} | {r['background_target']} | {r['selected_id'] or 'infeasible'} | {value('held_out_mean_iq_nrmse')} | {value('held_out_mean_background')} | {value('held_out_target_pass_fraction')} |")
    lines+=['',f"Maximum frozen closure error: {output_summary['max_frozen_closure']:.3g}.",
            f"Maximum Axis 1 / matched weighted-SVD primary-metric discrepancy: {max(equivalence):.3g}.",
            '', 'The accompanying figures and CSVs are regenerated by `scripts/plot_shared_benchmark.py`. Statistical-IQ software results do not establish experimental or acoustic performance.']
    (directory/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in output_summary.items() if k not in ('selected','condition_summaries')}))
    return output_summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('--figures')
    args=p.parse_args();summarize(args.directory,args.figures)
