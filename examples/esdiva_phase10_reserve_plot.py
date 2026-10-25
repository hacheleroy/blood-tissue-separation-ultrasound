"""Plot the frozen phase-10 reserved-bank paired gate without changing its decision."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['svg.fonttype'] = 'none'


def plot(decision_path, output_path):
    result = json.loads(Path(decision_path).read_text())
    if result['reserved_test_seeds'] != [2401, 2411]:
        raise ValueError('Changed phase-10 reserved bank')
    colors = {2401:'#4c78a8', 2411:'#f58518'}
    markers = {'svd_cut1':'o','axis1_cut1':'s'}
    fig,axes = plt.subplots(1,2,figsize=(10,4.7),sharey=True)
    fig.subplots_adjust(left=.09,right=.98,top=.84,bottom=.27,wspace=.14)
    for axis,scene,title in zip(axes,('rigid_fast','smooth_slow'),
                                ('Rigid tissue, fast flow','Smooth tissue, slow flow')):
        for row in result['rows']:
            if row['scene'] != scene:
                continue
            axis.scatter(row['tissue_leakage_ratio'],row['blood_nrmse_increase'],
                         s=70,marker=markers[row['method']],c=colors[row['seed']],
                         edgecolors='black',linewidths=.5,zorder=3)
        axis.axvline(.9,color='#555555',linestyle='--',linewidth=1)
        axis.axhline(.02,color='#555555',linestyle='--',linewidth=1)
        axis.set_title(title)
        axis.set_xlabel('Tissue leakage ratio to no correction (≤ 0.900)')
        axis.grid(alpha=.2,zorder=0)
        axis.margins(.15)
    axes[0].set_ylabel('Paired blood NRMSE increase (≤ 0.020)')
    legend = [Line2D([0],[0],marker='o',color='none',markerfacecolor=color,
                     markeredgecolor='black',markersize=8,label=f'Seed {seed}')
              for seed,color in colors.items()]
    legend += [Line2D([0],[0],marker=marker,color='none',markerfacecolor='white',
                      markeredgecolor='black',markersize=8,label=label)
               for marker,label in [('o','SVD cut 1'),('s','Axis 1 cut 1')]]
    fig.legend(handles=legend,ncol=4,loc='lower center',bbox_to_anchor=(.5,.08),
               frameon=False)
    fig.suptitle('Phase 10 · frozen reserved scatterer clouds',fontsize=13)
    if result['failures'] or result['missing']:
        note = ', '.join(result['failures']+result['missing'])
        fig.text(.5,.015,'Incomplete/failed: '+note[:145],ha='center',fontsize=8)
    output_path=Path(output_path)
    output_path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output_path,bbox_inches='tight',
                metadata={'Creator':'btbench phase-10 reserved fixed decision'})
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--decision',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    plot(args.decision,args.output)
