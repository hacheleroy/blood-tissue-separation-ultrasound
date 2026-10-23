"""Analyze every frozen condition; plot fixed calibration target and equivalence checks."""
import argparse
import gzip
import hashlib
import json
import lzma
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from threadpoolctl import threadpool_limits
from btbench.followup_experiment import make_followup_scene,numerical_source_hash
from btbench.benchmark import run_filter


FAMILIES=[('svd_short','Short SVD'),('svd_long','Long SVD'),('local_svd_short','Local short SVD'),
    ('axis1_memory','Axis 1, memory'),('axis1_rigid_memory','Axis 1, rigid + memory'),
    ('lee_temporal_raw_w256','Lee temporal ZCA'),('lee_spatial_raw_w256','Lee spatial ZCA'),
    ('axis2_original','Axis 2, original'),('axis2_ridge_capped','Axis 2, ridge + cap')]


def save(fig,path):
    fig.savefig(path.with_suffix('.png'),dpi=180,bbox_inches='tight')
    fig.savefig(path.with_suffix('.svg'),bbox_inches='tight')
    p=path.with_suffix('.svg');p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('study');p.add_argument('figures');args=p.parse_args()
    root=Path(args.study);figdir=Path(args.figures);figdir.mkdir(parents=True,exist_ok=True)
    archive=root/'results.json.xz'
    report=json.loads(lzma.decompress(archive.read_bytes()) if archive.exists()
                      else gzip.decompress((root/'results.json.gz').read_bytes()))
    config=json.loads((root/'frozen_config.json').read_text())
    if report['code_sha256']!=numerical_source_hash():raise ValueError('Analyze with the frozen numerical code')
    cases=[r for r in report['records'] if r['split']=='evaluation']
    conditions=[c['label'] for c in config['conditions']]
    titles={'slow_overlap':'Slow spectral overlap','accelerating':'Accelerating tissue','texture_jump':'Texture jump',
            'opposed_flow':'Opposed flow','translated_texture':'Padded image translation','psf_deforming':'PSF: deforming tissue',
            'equal_static_failure_control':'Both sources static (failure control)'}
    colors=plt.get_cmap('tab10').colors
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    selected=[s for s in report['matched_operating_points'] if s['target']==1e-4]
    fig,axes=plt.subplots(3,3,figsize=(13,10.5),layout='constrained')
    for ax,condition in zip(axes.flat,conditions):
        for j,(family,label) in enumerate(FAMILIES):
            point=next(s for s in selected if s['condition']==condition and s['family']==family)
            for row in point['evaluation']:
                if row['status']!='ok':continue
                m=row['metrics'];x=m['refitted_background_to_tissue_power'];y=m['blood_iq_nrmse']
                ax.scatter(max(x,1e-16),max(y,1e-6),s=34,color=colors[j],marker='o',edgecolor='white',linewidth=.4)
        ax.set(xscale='log',yscale='log',title=titles[condition],xlabel='Refitted background / input-tissue power',
               ylabel='Output-to-blood IQ error')
        ax.axvline(1e-4,color='.6',ls='--',lw=.8);ax.grid(True,which='major',alpha=.2)
    for ax in axes.flat[len(conditions):]:ax.set_visible(False)
    handles=[Line2D([],[],color=colors[j],marker='o',ls='',label=label) for j,(_,label) in enumerate(FAMILIES)]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Held-out reconstruction and refitted background response\nSettings selected on two calibration seeds at background/tissue power ≤ 10⁻⁴; each point is one held-out seed',fontsize=13)
    fig.get_layout_engine().set(rect=(0,.11,1,.85))
    save(fig,figdir/'followup_v04_matched')

    pairs=[('Global SVD','svd_short_c0_k1','svd_spatial_equivalence_w32'),
           ('AER, common normalization','aer_common_equivalence_temporal','aer_common_equivalence_spatial'),
           ('Local SVD / no-memory Axis 1','local_svd_short_c0_k1','axis1_no_memory_c0_k1')]
    specs={m['id']:m for m in config['methods']};definitions={c['label']:c for c in config['conditions']};checks=[]
    with threadpool_limits(limits=1):
        for case in cases:
            a=make_followup_scene(config,definitions[case['condition']],case['seed'])
            if hashlib.sha256(a.iq.tobytes()).hexdigest()!=case['input_sha256']:raise ValueError('Regenerated input changed')
            measured={r['id']:r for r in case['results']}
            for label,mid1,mid2 in pairs:
                left=run_filter(a.iq,specs[mid1],prf_hz=a.prf_hz);right=run_filter(a.iq,specs[mid2],prf_hz=a.prf_hz)
                x,y=measured[mid1]['metrics'],measured[mid2]['metrics']
                checks.append({'condition':case['condition'],'seed':case['seed'],'pair':label,
                    'relative_mixture_difference':float(np.linalg.norm(left.iq-right.iq)/max(np.linalg.norm(left.iq),1e-30)),
                    'relative_frozen_blood_attribution_difference':float(np.linalg.norm(left.apply_frozen(a.blood_true)-right.apply_frozen(a.blood_true))/np.linalg.norm(a.blood_true)),
                    'blood_distortion_first':x['frozen_blood_iq_nrmse'],'blood_distortion_second':y['frozen_blood_iq_nrmse'],
                    'refitted_response_first':x['refitted_blood_response_nrmse_alpha1'],'refitted_response_second':y['refitted_blood_response_nrmse_alpha1']})
    (root/'equivalence_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.4),layout='constrained')
    for j,(label,_,_) in enumerate(pairs):
        group=[r for r in checks if r['pair']==label]
        axes[0].scatter(np.full(len(group),j),[max(r['relative_mixture_difference'],1e-16) for r in group],s=24,color=colors[j],alpha=.65)
        axes[1].scatter([max(r['blood_distortion_first'],1e-6) for r in group],[max(r['blood_distortion_second'],1e-6) for r in group],s=27,color=colors[j],alpha=.7,label=label)
        axes[2].scatter([max(r['refitted_response_first'],1e-6) for r in group],[max(r['refitted_response_second'],1e-6) for r in group],s=27,color=colors[j],alpha=.7)
    axes[0].set(yscale='log',title='Equivalent mixture outputs',ylabel='Relative output difference',xticks=[0,1,2],xticklabels=['SVD','AER','Local / Axis 1'])
    for ax,title,xlabel,ylabel in [(axes[1],'Frozen attribution depends on side','First implementation: frozen blood error','Second implementation: frozen blood error'),
                                  (axes[2],'Refitted response agrees','First implementation: injection error','Second implementation: injection error')]:
        ax.set(xscale='log',yscale='log',title=title,xlabel=xlabel,ylabel=ylabel)
        lo=min(ax.get_xlim()[0],ax.get_ylim()[0]);hi=max(ax.get_xlim()[1],ax.get_ylim()[1])
        ax.plot([lo,hi],[lo,hi],ls='--',lw=.8,color='.5');ax.set_xlim(lo,hi);ax.set_ylim(lo,hi)
    for ax in axes:ax.grid(True,which='major',alpha=.2)
    fig.legend(handles=[Line2D([],[],marker='o',ls='',color=colors[j],label=label) for j,(label,_,_) in enumerate(pairs)],
               loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('A correction to the earlier preservation interpretation · all seven conditions and three held-out seeds',fontsize=13)
    fig.get_layout_engine().set(rect=(0,.12,1,.82))
    save(fig,figdir/'followup_v04_attribution')

    gains=[]
    for case in cases:
        rows={r['id']:r for r in case['results']}
        for centered in (0,1):
            raw=rows[f'axis2_original_c{centered}_m1'];cap=rows[f'axis2_cap_c{centered}']
            def gain(row):return max(p.get('operator_gain_max',1.) for p in row['diagnostics']['packets'])
            gains.append({'condition':case['condition'],'seed':case['seed'],'centered':bool(centered),
                          'original_gain':gain(raw),'capped_gain':gain(cap),
                          'original_output_error':raw['metrics']['blood_iq_nrmse'],'capped_output_error':cap['metrics']['blood_iq_nrmse']})
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3),layout='constrained')
    for j,condition in enumerate(conditions):
        group=[r for r in gains if r['condition']==condition]
        axes[0].scatter([max(r['original_gain'],1e-8) for r in group],[max(r['capped_gain'],1e-8) for r in group],s=31,color=colors[j],alpha=.8,label=titles[condition])
        axes[1].scatter([max(r['original_output_error'],1e-6) for r in group],[max(r['capped_output_error'],1e-6) for r in group],s=31,color=colors[j],alpha=.8)
    axes[0].set(xscale='log',yscale='log',xlabel='Original maximum local-operator gain',ylabel='Capped maximum local-operator gain',title='Conditioning control')
    axes[0].axhline(1,color='.4',lw=.8,ls='--')
    axes[1].set(xscale='log',yscale='log',xlabel='Original output-to-blood IQ error',ylabel='Capped output-to-blood IQ error',title='Actual reconstruction response')
    lo=min(axes[1].get_xlim()[0],axes[1].get_ylim()[0]);hi=max(axes[1].get_xlim()[1],axes[1].get_ylim()[1])
    axes[1].plot([lo,hi],[lo,hi],ls='--',color='.5',lw=.8);axes[1].set_xlim(lo,hi);axes[1].set_ylim(lo,hi)
    for ax in axes:ax.grid(True,which='major',alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Axis 2 · cap-only paired comparison, centered and uncentered · three held-out seeds\nA bounded local map does not by itself establish correct source identification',fontsize=12)
    fig.get_layout_engine().set(rect=(0,.21,1,.69))
    save(fig,figdir/'followup_v04_axis2')

    summary={'cases':len(report['records']),'outputs':sum(len(c['results']) for c in report['records']),
             'errors':[{'condition':c['condition'],'seed':c['seed'],'id':r['id'],'error':r.get('error')} for c in report['records'] for r in c['results'] if r['status']!='ok'],
             'max_equivalent_output_difference':max(r['relative_mixture_difference'] for r in checks),
             'max_frozen_attribution_difference':max(r['relative_frozen_blood_attribution_difference'] for r in checks),
             'max_original_axis2_local_gain':max(r['original_gain'] for r in gains),
             'max_capped_axis2_local_gain':max(r['capped_gain'] for r in gains),'matched_target':1e-4,
             'matched_operating_points_at_target':[
                 {'condition':s['condition'],'family':s['family'],'status':s['status'],'selected_id':s['selected_id'],
                  'evaluation':[{'seed':r['seed'],'meets_target':r['meets_target'],
                                 'blood_iq_nrmse':r['metrics']['blood_iq_nrmse'],
                                 'refitted_background_to_tissue_power':r['metrics']['refitted_background_to_tissue_power'],
                                 'refitted_blood_response_nrmse_alpha1':r['metrics']['refitted_blood_response_nrmse_alpha1']} for r in s['evaluation']]}
                 for s in selected],'axis2_cap_pairs':gains}
    (root/'ANALYSIS.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    lines=['# Frozen follow-up execution','',f"{len(report['records'])} cases × {len(config['methods'])} settings; {summary['outputs']} returned records.",
        '','Two calibration seeds; three held-out seeds per condition. No evaluation truth selects settings.',
        'Primary comparison uses refitted background response and actual output error. Frozen attribution is secondary.',
        '',f"Maximum relative difference of equivalent mixture outputs: {summary['max_equivalent_output_difference']:.6g}.",
        f"Maximum normalized difference in frozen blood attribution: {summary['max_frozen_attribution_difference']:.6g}.",
        '', '| Condition | Family | Selected setting at calibration target 1e-4 | Held-out median output IQ error | Held-out median background/tissue power | Held-out target passes |',
        '|---|---|---|---:|---:|---:|']
    for condition in conditions:
        for family,label in FAMILIES:
            point=next(r for r in selected if r['condition']==condition and r['family']==family)
            values=point['evaluation']
            if not values:lines.append(f'| {condition} | {label} | infeasible on calibration | n/a | n/a | 0/0 |');continue
            loss=np.median([r['metrics']['blood_iq_nrmse'] for r in values]);bg=np.median([r['metrics']['refitted_background_to_tissue_power'] for r in values])
            lines.append(f"| {condition} | {label} | {point['selected_id']} | {loss:.6g} | {bg:.6g} | {sum(r['meets_target'] for r in values)}/{len(values)} |")
    lines+=['','Infeasible settings are not interpolated. Passing the background-only target does not prove tissue leakage inside F(X).',
             'Absolute AER units and amplitude changes remain in the scores; no gain is fitted to source truth.',
             'All thresholds, held-out failures, source distortions and conditional operator diagnostics are retained in the archive.']
    (root/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['matched_operating_points_at_target','axis2_cap_pairs']}),flush=True)


if __name__=='__main__':main()
