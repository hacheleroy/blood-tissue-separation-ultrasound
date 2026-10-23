"""Reproduce selected frozen real-IQ outputs for display, without accuracy claims."""
import argparse
import gzip
import lzma
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from btbench.io import load_acquisition
from btbench.benchmark import run_filter
from btbench.followup_experiment import numerical_source_hash


CHOICES=[('unfiltered','Unfiltered'),('svd_w256_k1','SVD (256), rank 1'),
         ('baranger_w256','Baranger SSM (256)'),('lee_temporal_w768','Lee temporal ZCA (768)'),
         ('lee_spatial_w768','Lee spatial ZCA (768)'),('song_pdf_paper','Song PDF, normalized'),
         ('axis1_memory1_rigid1','Axis 1, memory + rigid'),('axis2_c0_guard1','Axis 2, ridge + cap')]


def main():
    p=argparse.ArgumentParser();p.add_argument('study');p.add_argument('datasets');p.add_argument('output');args=p.parse_args()
    root=Path(args.study);data=Path(args.datasets);output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    archive=root/'results.json.xz'
    report=json.loads(lzma.decompress(archive.read_bytes()) if archive.exists()
                      else gzip.decompress((root/'results.json.gz').read_bytes()))
    if report['code_sha256']!=numerical_source_hash():raise ValueError('Use the frozen numerical source')
    config=json.loads((root/'frozen_config.json').read_text());checks=[]
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(4,4,figsize=(13,11),layout='constrained')
    with threadpool_limits(limits=1):
        for i,record in enumerate(report['records']):
            a=load_acquisition(data/(record['acquisition']+'.npz'))
            cfg=json.loads((root/record['acquisition']/'frozen_config.json').read_text())
            specs={s['id']:s for s in cfg['methods']};n=config['calibration_frames'];stride=config['evaluation_spatial_stride']
            views=a.iq_views[:,:,::stride,::stride];cal=views[:n];iq=views[n:].mean(axis=1)
            if hashlib.sha256(iq.tobytes()).hexdigest()!=record['input_sha256']:raise ValueError('Experimental input changed')
            power=np.mean(abs(cal)**2,axis=(0,1));roi=power>=np.quantile(power,.75)
            extent=[a.x_m[0]*1e3,a.x_m[::stride][-1]*1e3,a.z_m[::stride][-1]*1e3,a.z_m[0]*1e3]
            rows={r['id']:r for r in record['results']}
            for j,(mid,title) in enumerate(CHOICES):
                result=run_filter(iq,specs[mid],tissue_roi=roi,prf_hz=a.prf_hz)
                pd=np.mean(abs(result.iq)**2,axis=0);expected=rows[mid]['metrics']['power_output']
                if not np.isclose(pd.mean(),expected,rtol=1e-10,atol=1e-30):raise ValueError('Regenerated power differs from archived result')
                reference=max(float(np.quantile(pd,.995)),1e-30)
                image=10*np.log10(np.maximum(pd/reference,1e-12))
                ax=axes[i*2+j//4,j%4]
                im=ax.imshow(image,extent=extent,aspect='auto',cmap='inferno',vmin=-40,vmax=0)
                ratio=rows[mid]['metrics'].get('tissue_roi_output_power_ratio')
                ax.set_title(f"S{record['subject_id']} · {title}\nBright-ROI ratio {ratio:.3g}",fontsize=9)
                ax.set(xlabel='Lateral (mm)',ylabel='Depth (mm)' if j%4==0 else '')
                checks.append({'acquisition':record['acquisition'],'method':mid,'mean_power':float(pd.mean()),
                               'archived_mean_power':expected,'display_reference_99_5th_percentile':reference})
    fig.colorbar(im,ax=axes,label='Power (dB relative to each panel’s own 99.5th percentile)',shrink=.65)
    fig.suptitle('Waterloo CCA · two subjects · 768 evaluation frames after 256 disjoint calibration frames\nSeparate display normalization for each panel; bright-ROI power and appearance do not establish blood preservation',fontsize=11)
    fig.savefig(output.with_suffix('.png'),dpi=180,bbox_inches='tight')
    fig.savefig(output.with_suffix('.svg'),bbox_inches='tight')
    s=output.with_suffix('.svg');s.write_text('\n'.join(line.rstrip() for line in s.read_text().splitlines())+'\n');plt.close(fig)
    (root/'PLOT_CHECKS.json').write_text(json.dumps(checks,indent=2)+'\n')
    lines=['# Longer Waterloo experimental processing pilot','',
        'Two distinct subjects; 1024 reconstructed frames each; first 256 calibrate nuisance ROI/velocity, last 768 evaluate.',
        'All parameters fixed before evaluation. No independent source/velocity/phase truth, no oracle and no verified contrast status.',
        'Display panels are normalized separately; their appearance is not a separation-accuracy score.','',
        '| Acquisition / subject | Method | Status | Bright-ROI output/input power | Median lag-one coherence | Maximum local operator gain |',
        '|---|---|---|---:|---:|---:|']
    for record in report['records']:
        for row in record['results']:
            m=row.get('metrics',{});packets=row.get('diagnostics',{}).get('packets',[])
            gains=[p.get('operator_gain_max',1.) for p in packets]
            def fmt(v):return 'n/a' if v is None else f'{v:.6g}'
            lines.append(f"| {record['acquisition']} / {record['subject_id']} | {row['id']} | {row['status']} | {fmt(m.get('tissue_roi_output_power_ratio'))} | {fmt(m.get('coherence_flow_median'))} | {fmt(max(gains) if gains else None)} |")
    lines+=['','The ROI is the brightest calibration-power quartile; its tissue-only status is unverified.',
            'A local cap does not guarantee full adaptive-pipeline stability or correct source identification.',
            'Absolute event timing and both uncompounded views are retained; each evaluation segment lasts about 0.256 s.']
    (root/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'regenerated_display_outputs':len(checks),'all_powers_match_archived_results':True}),flush=True)


if __name__=='__main__':main()
