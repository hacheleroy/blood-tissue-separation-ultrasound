"""Summarize the frozen v06 held-out bank without tuning its gate."""
import csv
import gzip
import json
import lzma
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
NEW=ROOT/'runs/motion_repair_v06'


def load_v05():
    parts=sorted((ROOT/'runs/shared_axes_v05').glob('results.json.xz.part*'))
    return json.loads(lzma.decompress(b''.join(p.read_bytes() for p in parts)))


def main():
    report=json.loads(gzip.decompress((NEW/'results.json.gz').read_bytes()))
    old=load_v05();config=json.loads((NEW/'frozen_config.json').read_text())
    cases=[c for c in report['records'] if c['split']=='evaluation']
    assert len(cases)==40 and len(report['records'])==64
    oldcases=[c for c in old['records'] if c['split']=='evaluation']
    lookup=lambda c,mid:next(r for r in c['results'] if r['id']==mid)
    modes=['none','rigid','smooth','known'];conditions=[c['id'] for c in config['conditions']]
    methods=['svd','axis1','axis2','stap','lee_aer']
    mids={'svd':'svd_{mode}_h64_k1','axis1':'axis1_{mode}_h64_k1',
          'axis2':'axis2_{mode}_h64_r4_p0.1','stap':'stap_{mode}_h64_loaded0.01',
          'lee_aer':'lee_aer_{mode}_h64_temporal_abs26'}
    rows=[]
    for c in cases:
        for family in methods:
            for mode in modes:
                record=lookup(c,mids[family].format(mode=mode))
                metrics=record['metrics'];motion=c['contexts'][mode]['actual']
                audit=motion['motion_fit'];counts=audit.get('reason_counts',{})
                rows.append(dict(condition=c['condition'],seed=c['seed'],family=family,motion=mode,
                    blood_iq_nrmse=metrics['blood_iq_nrmse'],
                    background_to_tissue_power=metrics['refitted_background_to_tissue_power'],
                    incremental_blood_nrmse_alpha1=metrics['refitted_blood_response_nrmse_alpha1'],
                    incremental_blood_nrmse_alpha0_1=metrics['refitted_blood_response_nrmse_alpha0_1'],
                    displacement_rmse_pixels=motion['evaluation_only_displacement_rmse_pixels'],
                    accepted_frames=audit.get('accepted_frames',0),validation_rejected=counts.get('validation',0),
                    support_rejected=counts.get('support',0),invalid_map_rejected=counts.get('invalid_map',0)))
    with (NEW/'metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    def mean(condition,family,mode,key):
        return float(np.mean([r[key] for r in rows if r['condition']==condition and r['family']==family and r['motion']==mode]))
    def pooled(family,mode,key):
        return float(np.mean([r[key] for r in rows if r['family']==family and r['motion']==mode]))
    motion_conditions=['subpixel_translation','smooth_deformation'];old_rmse={}
    for condition in motion_conditions:
        old_rmse[condition]=float(np.mean([c['contexts']['smooth']['actual']['evaluation_only_displacement_rmse_pixels']
            for c in oldcases if c['condition']==condition]))
    geometry={cond:{mode:mean(cond,'svd',mode,'displacement_rmse_pixels') for mode in modes}
              for cond in conditions}
    paired={family:{mode:{'blood_iq_nrmse':pooled(family,mode,'blood_iq_nrmse'),
        'background_to_tissue_power':pooled(family,mode,'background_to_tissue_power'),
        'incremental_blood_nrmse_alpha1':pooled(family,mode,'incremental_blood_nrmse_alpha1')}
        for mode in modes} for family in methods}
    static=next(c for c in conditions if c=='slow_overlap')
    # Static gate is checked at the primitive level, as some conditions contain
    # phase drift or imposed motion. All other criteria are evaluated on v06.
    checks={
        'translation_geometry_better_than_v05_smooth':geometry['subpixel_translation']['smooth']<old_rmse['subpixel_translation'],
        'smooth_geometry_better_than_v05_smooth':geometry['smooth_deformation']['smooth']<old_rmse['smooth_deformation'],
        'svd_blood_error_no_more_than_5_percent_above_none':paired['svd']['smooth']['blood_iq_nrmse']<=1.05*paired['svd']['none']['blood_iq_nrmse'],
        'axis1_blood_error_no_more_than_5_percent_above_none':paired['axis1']['smooth']['blood_iq_nrmse']<=1.05*paired['axis1']['none']['blood_iq_nrmse'],
    }
    summary=dict(study=report['study'],code_sha256=report['code_sha256'],config_sha256=report['config_sha256'],
        training_cases=24,heldout_cases=40,settings_per_case=36,heldout_seeds=config['evaluation_seeds'],
        old_smooth_displacement_rmse=old_rmse,displacement_rmse=geometry,
        paired_heldout_means=paired,criteria=checks,software_gate_passed=all(checks.values()),
        interpretation='Prescribed warp evaluates geometry only; selection/gating never receives it. Synthetic angular IQ, not measured RF or phantom validation.')
    (NEW/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    for mode in modes:
        axes[0].plot(range(len(conditions)),[geometry[c][mode] for c in conditions],'o-',label=mode)
    axes[0].set(ylabel='Displacement RMSE (pixels; evaluation-only)',title='Motion recovery')
    display={'svd':'Joint-view SVD','axis1':'Axis 1: Motion-compensated virtual long ensemble',
             'axis2':'Axis 2: Motion-demodulated multi-lag separation',
             'stap':'Axis 3: Local ultrasound STAP','lee_aer':'Stephen Lee AER'}
    for family in methods:
        axes[1].plot(range(len(modes)),[paired[family][m]['blood_iq_nrmse'] for m in modes],
                     'o-',label=display[family])
    axes[1].set(xticks=range(len(modes)),xticklabels=modes,ylabel='Held-out actual blood IQ NRMSE',
                title='Paired filtering outcomes (all conditions)')
    axes[0].set_xticks(range(len(conditions)),conditions,rotation=65,ha='right',fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8,loc='upper left',bbox_to_anchor=(1.02,1))
    fig.suptitle('Motion repair v06 · five reserved held-out seeds · no truth-based tuning')
    fig.savefig(ROOT/'docs/figures/motion_repair_v06.png',dpi=180,bbox_inches='tight')
    fig.savefig(ROOT/'docs/figures/motion_repair_v06.svg',bbox_inches='tight',metadata={'Date':None})
    plt.close(fig)
    svg=ROOT/'docs/figures/motion_repair_v06.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
