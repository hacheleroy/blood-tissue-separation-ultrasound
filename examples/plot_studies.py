"""Generate scientific SVG/PNG figures from the saved, unchanged study records."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from btbench.stap import steering,null_projector


LABELS={"unfiltered":"Unfiltered","butterworth_30Hz":"Butterworth 30 Hz",
    "svd_joint_matched":"Joint SVD, matched 8 frames","svd_per_view_matched":"Per-view SVD, matched 8 frames",
    "svd_joint_long":"Joint SVD, 64 frames","baranger_joint_long":"Baranger SSM, 64 frames",
    "calibrated_demodulated_svd":"Calibrated demodulation + SVD",
    "stap_orthogonal":"STAP orthogonal, 8 taps","stap_covariance":"STAP covariance, 8 taps",
    "stap_covariance_16taps":"STAP covariance, 16 taps",
    "stap_training_starved_control":"STAP training-starved control",
    "oracle_svd_joint_matched":"IQ-oracle SVD, matched bands"}


def save(fig,path):
    svg=path.with_suffix('.svg')
    fig.savefig(svg,bbox_inches='tight')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(path.with_suffix('.png'),bbox_inches='tight',dpi=170)
    plt.close(fig)


def tradeoff(report,out):
    methods=["butterworth_30Hz","svd_joint_matched","svd_joint_long","baranger_joint_long",
        "calibrated_demodulated_svd","stap_orthogonal","stap_covariance","stap_covariance_16taps"]
    fig,axes=plt.subplots(1,2,figsize=(12,5.6),sharey=True)
    for ax,speed in zip(axes,[.003,.02]):
        cases=[c for c in report['records'] if c['case'].get('label')=='factorial'
            and c['case'].get('umax_relative_mps')==speed and c['case'].get('clutter_db')==40
            and c['case'].get('tissue_acceleration_mps2')==0]
        for j,name in enumerate(methods):
            rows=[r['metrics'] for c in cases for r in c['results'] if r['method']==name and r['status']=='ok']
            tx=np.array([-10*np.log10(max(r['frozen_tissue_power_retention'],1e-18)) for r in rows])
            by=np.array([max(r['frozen_blood_power_retention'],1e-18) for r in rows])
            qx,qy=np.percentile(tx,[25,50,75]),np.percentile(by,[25,50,75])
            marker=['o','s','^','D','v','P','x','h'][j]
            ax.errorbar(qx[1],qy[1],xerr=[[qx[1]-qx[0]],[qx[2]-qx[1]]],
                yerr=[[qy[1]-qy[0]],[qy[2]-qy[1]]],fmt=marker,color=f'C{j}',capsize=2,
                markersize=8 if j==5 else 6,markerfacecolor='none' if j==5 else f'C{j}',
                label=LABELS[name])
        ax.set_title(f'Relative peak flow {speed*1000:g} mm/s')
        ax.set_xlabel('Frozen tissue suppression (dB)')
        ax.set_yscale('log');ax.axhline(1,color='.5',lw=.8,ls='--')
        ax.grid(alpha=.22,which='both')
    axes[0].set_ylabel('Frozen blood energy / input blood energy')
    axes[1].legend(loc='upper center',bbox_to_anchor=(.2,-.2),ncol=2,fontsize=8,frameon=False)
    fig.suptitle('Angle-time IQ surrogate: retention and rejection\n40 dB clutter, constant tissue velocity; median and interquartile range of 3 seeds',fontsize=12)
    fig.subplots_adjust(bottom=.33,top=.82,wspace=.15)
    save(fig,out/'stap_tradeoff')


def stress(report,out):
    names=["butterworth_30Hz","svd_joint_matched","baranger_joint_long","calibrated_demodulated_svd",
        "stap_orthogonal","stap_covariance","stap_covariance_16taps"]
    cases=[('view_angle_error','Angle error'),('view_phase_error','View phase error'),
        ('training_contamination','Training contamination'),('velocity_spread','Tissue velocity spread'),
        ('identical_velocity','Identical velocities'),('opposed_flow','Opposed flow')]
    values=np.full((len(names),len(cases)),np.nan)
    invalid=np.zeros_like(values,dtype=int)
    for i,name in enumerate(names):
        for j,(case,_) in enumerate(cases):
            rows=[r for c in report['records'] if c['case']['label']==case for r in c['results'] if r['method']==name]
            val=[r.get('metrics',{}).get('axial_velocity_nrmse') for r in rows]
            valid=[v for v in val if v is not None and np.isfinite(v)]
            invalid[i,j]=len(val)-len(valid)
            if valid:values[i,j]=np.median(valid)
    fig,ax=plt.subplots(figsize=(11,5.2))
    im=ax.imshow(np.maximum(values,1e-3),norm=LogNorm(vmin=.01,vmax=10),cmap='viridis')
    ax.set_xticks(range(len(cases)),[c[1] for c in cases],rotation=20,ha='right')
    ax.set_yticks(range(len(names)),[LABELS[n] for n in names])
    for i in range(len(names)):
        for j in range(len(cases)):
            label='n/a' if not np.isfinite(values[i,j]) else f'{values[i,j]:.2g}'
            if invalid[i,j]:label+=f'\n{invalid[i,j]} invalid'
            ax.text(j,i,label,ha='center',va='center',fontsize=9,color='white' if values[i,j]<.5 else 'black')
    fig.colorbar(im,ax=ax,label='Median axial velocity NRMSE (color clipped to [0.01,10])',shrink=.8)
    ax.set_title('STAP stress cases: all three seeds retained\nCell values show actual medians; missing estimates are counted explicitly',fontsize=11)
    fig.tight_layout()
    save(fig,out/'stap_stress')


def ideal_response(out):
    tx=np.deg2rad([-15,0,15]);rx=np.zeros(3)
    dv=np.linspace(-.03,.03,601)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for taps in (4,8,16,32):
        times=np.arange(taps)/1000
        st=steering([[.001,0]],tx,rx,times,5e6)
        sb=steering(np.column_stack((.001+dv,np.zeros(len(dv)))),tx,rx,times,5e6)
        y=null_projector(st)@sb
        retention=np.sum(abs(y)**2,axis=0)/np.sum(abs(sb)**2,axis=0)
        ax.plot(dv*1000,10*np.log10(np.maximum(retention,1e-15)),label=f'{taps} taps / {taps-1} ms span')
    ax.set_ylim(-70,1);ax.set_xlim(-30,30)
    ax.set_xlabel('Blood axial velocity relative to tissue (mm/s)')
    ax.set_ylabel('Blood snapshot energy response (dB)')
    ax.set_title('Ideal orthogonal tissue null: response around the clutter velocity\n3 calibrated views, 1000 Hz revisit, start time 0; no noise or source fitting',fontsize=11)
    ax.grid(alpha=.25);ax.legend(fontsize=9)
    fig.tight_layout();save(fig,out/'stap_ideal_response')


def main():
    p=argparse.ArgumentParser();p.add_argument('--stap',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none',
                         'axes.spines.right':False,'axes.spines.top':False})
    report=json.loads(Path(a.stap).read_text())
    tradeoff(report,out);stress(report,out);ideal_response(out)


if __name__=='__main__':main()
