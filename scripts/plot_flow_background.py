"""Plot matched slow-flow signal/background ratios against blood IQ error."""
import csv
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import numpy as np


ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'runs/flow_background_audit_v1/rows.csv'
OUTPUT=ROOT/'docs/figures/flow_background_tradeoff.png'
NAMES={
    'unfiltered':'Input',
    'compounded_svd_rank1':'SVD',
    'jiang_ga_svd_rank1':'GA-SVD',
    'ozgun_svd_a_rank1':'SVD-a',
    'jiang_ac_hosvd_adaptation':'AC-HOSVD*',
    'ozgun_hosvd_adaptation':'Rx-HOSVD*',
    'angle_core_intersection':'Angle core',
    'receive_core_intersection':'Rx core',
    'calibrated_stap_orthogonal':'STAP',
    'none/roundtrip':'Input',
    'none/svd_cut1':'SVD',
    'none/axis1_cut1':'Axis 1',
    'none/axis2_default':'Axis 2',
    'none/axis2_slow_r4_cap040':'Axis 2 capped',
    'none/axis3_orthogonal':'Axis 3',
    'none/baranger_ssm':'Baranger',
    'none/lee_aer':'AER'}


def plot():
    with INPUT.open() as f:rows=list(csv.DictReader(f))
    panels=[('multidim_aperture_v1',{'rigid_slow:clean','smooth_slow:clean'},
             'Acoustic RF · coherent image',{'exact_frozen'}),
            ('esdiva_phase11_prospective',{'rigid_slow','smooth_slow'},
             'Phase 11 · three IQ views',{'exact_paired_alpha1'})]
    fig,axes=plt.subplots(1,2,figsize=(13,7.2))
    for ax,(archive,conditions,title,kinds) in zip(axes,panels):
        group={}
        for r in rows:
            if (r['archive']!=archive or r['condition'] not in conditions or
                r['metric_kind'] not in kinds or r['method'] not in NAMES or not r['sbr_db']):
                continue
            group.setdefault(r['method'],[]).append(r)
        for j,(method,items) in enumerate(group.items()):
            xx=np.array([float(r['blood_iq_error']) for r in items])
            yy=np.array([float(r['sbr_db']) for r in items])
            retained=np.array([float(r['blood_power_retention']) for r in items])
            x,y=median(xx),median(yy)
            size=55+65*np.sqrt(max(float(median(retained)),0))
            color=plt.cm.tab20(j/19)
            ax.scatter(x,y,s=size,color=color,edgecolors='black',linewidths=.4,
                       zorder=3,label=NAMES[method])
            ax.vlines(x,min(yy),max(yy),color=color,lw=1,alpha=.8)
        ax.axhline(0,color='#666',linestyle=':',lw=1)
        ax.set(title=title,xlabel='Added blood IQ error (lower is better)',
               ylabel='Blood / residual background (dB, higher is better)')
        ax.grid(alpha=.18)
        ax.legend(loc='upper center',bbox_to_anchor=(.5,-.20),ncol=3,fontsize=8,frameon=False)
    fig.subplots_adjust(left=.08,right=.98,bottom=.36,top=.86,wspace=.25)
    fig.suptitle('Each point is the median of four slow-flow RF scenes; bars span all four',fontsize=12,y=.965)
    fig.text(.5,.07,'* Fixed-rule structural adaptation of the published HOSVD method. '
             'Point size indicates median surviving blood power; axes differ between panels.',
             ha='center',fontsize=9)
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUTPUT,dpi=190,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':plot()
