"""Plot the archived retrospective source/velocity diagnostic."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def med(cases,method,region,key,window=16):
    return np.median([next(m for m in c['methods'] if m['method']==method)
                      ['scores'][str(window)][region][key] for c in cases])


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--input',type=Path,required=True)
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args()
    cases=json.loads(a.input.read_text())['cases']
    methods=['unfiltered','compounded_svd_rank1','receive_core_intersection']
    names=['Raw','SVD','Receive core']
    colors=['#777777','#3569a8','#cb6a32']
    regions=['center','middle','near_wall']
    fig,axs=plt.subplots(2,2,figsize=(10.6,7.5),layout='constrained')
    fig.suptitle('Source-separated slow-flow RF: four previously seen scenes',fontsize=13)

    ax=axs[0,0]
    for i,(method,name,color) in enumerate(zip(methods,names,colors)):
        case_vals=[]
        for c in cases:
            s=next(m for m in c['methods'] if m['method']==method)['scores']['16']['whole_tile']
            den=next(m for m in c['methods'] if m['method']=='unfiltered')['scores']['16']['whole_tile']['tissue']
            case_vals.append((s['tissue']/den,s['noise']/den,s['background']/den))
        vals=np.median(case_vals,axis=0)
        ax.bar(i-.15,vals[0],width=.28,color=color,alpha=.95,label='Tissue' if i==0 else None)
        ax.bar(i+.15,vals[1],width=.28,color=color,alpha=.35,label='Noise' if i==0 else None)
        ax.scatter([i],[vals[2]],s=28,marker='D',facecolor='white',edgecolor='black',zorder=3,
                   label='Tissue + noise (including cross term)' if i==0 else None)
    ax.set_xticks(range(3),names)
    ax.set_yscale('log');ax.set_ylim(.006,1.8)
    ax.set_ylabel('Power / input tissue power')
    ax.set_title('A  Remaining background, final 16 cycles',loc='left',fontweight='bold')
    ax.legend(fontsize=7,loc='upper right',frameon=False)

    ax=axs[0,1]
    source_names=['Blood','Tissue','Noise']
    labels=['First time / first receive (removed)','First time / other receive',
            'Other time / first receive','Other time / other receive']
    keys=['first_time_first_receive','first_time_other_receive',
          'other_time_first_receive','other_time_other_receive']
    col=['#cc5522','#e5ad59','#4187b2','#8ec2c8']
    bottoms=np.zeros(3)
    for key,label,color in zip(keys,labels,col):
        values=[np.mean([c['tensor_mode_fraction'][s.lower()][key] for c in cases])
                for s in source_names]
        ax.bar(range(3),values,bottom=bottoms,label=label,color=color)
        bottoms+=values
    ax.set_ylim(0,1.04);ax.set_xticks(range(3),source_names)
    ax.set_ylabel('Fraction of source energy before receive summation')
    ax.set_title('B  Mixture-fitted HOSVD core',loc='left',fontweight='bold')
    ax.legend(fontsize=6.5,loc='upper left',bbox_to_anchor=(1.01,1.),
              frameon=False)

    ax=axs[1,0]
    for method,name,color in zip(methods,names,colors):
        y=np.array([med(cases,method,r,'blood_background_db') for r in regions])
        lo=np.array([min(next(m for m in c['methods'] if m['method']==method)['scores']['16'][r]
                             ['blood_background_db'] for c in cases) for r in regions])
        hi=np.array([max(next(m for m in c['methods'] if m['method']==method)['scores']['16'][r]
                             ['blood_background_db'] for c in cases) for r in regions])
        ax.errorbar(range(3),y,yerr=np.stack((y-lo,hi-y)),marker='o',capsize=3,
                    color=color,label=name,lw=1.5)
    ax.axhline(0,color='k',ls=':',lw=.8)
    ax.set_xticks(range(3),['Center','Middle','Near wall'])
    ax.set_ylabel('Blood / (tissue + noise), dB')
    ax.set_title('C  Lumen contrast, final 16 cycles',loc='left',fontweight='bold')
    ax.legend(fontsize=8,frameon=False)

    ax=axs[1,1]
    for method,name,color in zip(methods,names,colors):
        y=[med(cases,method,r,'observed_velocity_mps',80)*1e3 for r in regions]
        ax.plot(range(3),y,marker='o',color=color,label=name)
    ax.plot(range(3),[med(cases,'unfiltered',r,'nominal_velocity_mps',80)*1e3 for r in regions],
            marker='s',color='#385b3a',ls='--',label='Nominal parabola')
    ax.plot(range(3),[med(cases,'unfiltered',r,'blood_only_velocity_mps',80)*1e3 for r in regions],
            marker='^',color='#80549b',ls='--',label='Blood-only IQ')
    ax.set_xticks(range(3),['Center','Middle','Near wall'])
    ax.set_ylabel('Median axial Doppler velocity, mm/s')
    ax.set_title('D  Velocity profile, all 80 cycles',loc='left',fontweight='bold')
    ax.legend(fontsize=7,frameon=False,ncol=2)

    a.out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.out,dpi=190)
    plt.close(fig)


if __name__=='__main__':main()
