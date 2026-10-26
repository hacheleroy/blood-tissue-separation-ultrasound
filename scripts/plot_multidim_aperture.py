"""Rebuild the slow-flow IQ-versus-background scientific figure."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
run=ROOT/'runs/multidim_aperture_v1'
cases=[]
for name in ('development','evaluation'):
    cases += [c for c in json.loads((run/f'{name}.json').read_text())['cases']
              if not c['phase_error'] and c['scene'] in ('rigid_slow','smooth_slow')]
names={
    'compounded_svd_rank1':'Compounded SVD',
    'jiang_ga_svd_rank1':'GA-SVD (angle stacking)',
    'ozgun_svd_a_rank1':'SVD-a (receive stacking)',
    'jiang_ac_hosvd_adaptation':'Angle HOSVD adaptation',
    'ozgun_hosvd_adaptation':'Receive HOSVD adaptation',
    'angle_core_intersection':'Angle core intersection',
    'receive_core_intersection':'Receive core intersection',
    'calibrated_stap_orthogonal':'Calibrated STAP'}
fig,ax=plt.subplots(figsize=(8.5,5.8))
palette=plt.cm.tab10(np.linspace(0,1,len(names)))
for (key,label),color in zip(names.items(),palette):
    x=np.array([c['methods'][key]['blood_iq_error'] for c in cases])
    y=np.array([c['methods'][key]['background_power_remaining'] for c in cases])
    ax.scatter(x,y,color=color,alpha=.25,s=35)
    ax.scatter(np.median(x),np.median(y),color=color,s=95,edgecolors='black',
               linewidths=.5,label=label,zorder=4)
ax.axvline(.25,color='0.4',linestyle='--',linewidth=1)
ax.axhline(.1,color='0.4',linestyle='--',linewidth=1)
ax.set(xlim=(0,1.06),ylim=(.005,.8),yscale='log',xlabel='Added blood IQ error (lower is better)',
       ylabel='Remaining background power / input tissue power (lower is better)',
       title='Four slow-flow RF scenes: fixed-rule filters on the same 240 emissions')
ax.legend(fontsize=8,ncol=2,loc='upper right')
ax.grid(True,which='both',alpha=.2)
fig.tight_layout()
target=ROOT/'docs/figures/multidim_slow_tradeoff.png'
fig.savefig(target,dpi=180)
print(target)
