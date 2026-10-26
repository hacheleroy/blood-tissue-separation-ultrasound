"""Paired matched-emission information-test figure from frozen archives."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUNS=ROOT/'runs/transport_lag_information_v1'
OUT=ROOT/'docs/figures/transport_lag_information_v1.png'


def main():
    rows=[]
    for split in ('development','evaluation'):
        rows.extend((split,r) for r in json.loads((RUNS/f'{split}.json').read_text())['records'])
    pairs={}
    for split,r in rows:pairs.setdefault((split,r['seed'],r['motion']),{})[r['schedule']]=r
    fig,axes=plt.subplots(2,2,figsize=(10.4,7.3),constrained_layout=True)
    descriptors=[
        ('Isolated-source AUC: blood vs tissue','source_auc_blood_vs_tissue',(0,1),lambda r:r['separation']['weak_wall']),
        ('Mixture minus static-blood IQ-control score','mixture_minus_static_blood_score_median',None,lambda r:r['separation']['weak_wall']),
        ('Input blood / background (dB)','aligned_blood_background_db',None,lambda r:r['separation']['weak_wall']),
        ('Lag-8 velocity error (mm/s)','weak_pixel_error_vs_blood_only_median_mm_s',None,lambda r:r['separation']['phase_velocity_lag8'])]
    colors={('development','rigid'):'#8a929a',('development','smooth'):'#b8bec4',
            ('evaluation','rigid'):'#155e91',('evaluation','smooth'):'#cd641b'}
    for ax,(title,key,limits,part) in zip(axes.ravel(),descriptors):
        for (split,seed,motion),pair in sorted(pairs.items()):
            x=[part(pair[s])[key] for s in ('dense','burst')]
            ax.plot((0,1),x,'o-',color=colors[(split,motion)],markersize=6,
                    label=f'{split}: {motion}' if ax is axes[0,0] else None)
        ax.set(xlim=(-.12,1.12),xticks=[0,1],xticklabels=['21.5 ms dense','200.1 ms bursts'],title=title)
        ax.grid(axis='y',alpha=.25)
        if limits:ax.set_ylim(*limits)
        if key=='mixture_minus_static_blood_score_median':ax.axhline(0,color='black',lw=.9,ls=':')
    axes[0,0].legend(frameon=False,fontsize=8,loc='lower right')
    fig.suptitle('Transport features at the eight weakest near-wall pixels\n'
                 'Same 240 transmissions; source truth selects pixels only for evaluation',fontsize=12)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT,dpi=170)
    plt.close(fig)
    print(OUT)


if __name__=='__main__':main()
