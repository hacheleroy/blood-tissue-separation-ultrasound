"""Measured receive-subaperture pilot on one-transmit Resolve Stroke CEUS.

No blood/tissue IQ truth or independent transmit angles are available. Native
DDC channels are geometrically delayed and carrier-rotated using Zea 0.1.6;
eight sorted receive groups are summed *after* the correction. A strict
all-channel beam-sum check guards against grouping raw channels as views.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

from btbench.multidim import intersection_filter,paper_tensor,svd_unfolding
from resolvestroke_single_view import digest,score_power,check_source


ROOT=Path(__file__).resolve().parents[1]
CFG=ROOT/'configs/resolvestroke_aperture_pilot_v1.json'
PIPE=ROOT/'configs/resolvestroke_complex.yaml'


def crop(h5,cfg):
    coords=h5['custom/computed_references/coordinates']
    x=coords[0,0,:,0];y=coords[0,:,0,1];z=coords[:,0,0,2]
    ix=int(np.argmin(abs(x-cfg['target_x_m'])))
    iy=np.flatnonzero(abs(y)<=cfg['max_abs_y_m'])
    iz=np.flatnonzero((z>=cfg['depth_m'][0])&(z<=cfg['depth_m'][1]))
    lo,hi=cfg['depth_rows_within_anatomical_grid']
    rows=iz[0]+np.arange(lo,hi)
    cols=iy[::cfg['lateral_stride']]
    grid=np.asarray(coords[int(rows[0]):int(rows[-1])+1,
                           int(cols[0]):int(cols[-1])+1,ix,:])[:,::cfg['lateral_stride']]
    mask=np.asarray(h5['custom/computed_references/tube_mask'][int(rows[0]):int(rows[-1])+1,
                 int(cols[0]):int(cols[-1])+1,ix])[:,::cfg['lateral_stride']]
    if grid.shape!=(10,39,3) or min((mask==1).sum(),(mask==2).sum())<5:
        raise ValueError('Anatomical fixed crop has changed')
    return grid,mask,{'x_index':ix,'rows':rows.tolist(),'columns':cols.tolist(),
                     'tube_pixels':{str(k):int((mask==k).sum()) for k in (1,2)}}


def rois(grid,mask,cfg,label):
    tube=mask==label
    depth=grid[...,2]
    within=(depth>=depth[tube].min()-cfg['tube_depth_margin_m'])&(depth<=depth[tube].max()+cfg['tube_depth_margin_m'])
    outside=(mask==0)&(distance_transform_edt(mask==0)>cfg['outside_guard_pixels'])&within
    if outside.sum()<40:raise ValueError('Anatomical background support too small')
    return tube,outside


def run(source,out):
    import h5py,hdf5plugin,zea
    from zea import Config,File,Pipeline
    from zea.beamform.beamformer import tof_correction

    cfg=json.loads(CFG.read_text())
    source=Path(source);out=Path(out)
    if (source.stat().st_size!=cfg['source_size_bytes'] or
            digest(source)!=cfg['source_sha256'] or zea.__version__!=cfg['zea_version']):
        raise ValueError('Dataset/beamformer identity differs')
    out.mkdir(parents=True,exist_ok=False)
    config=Config.from_path(str(PIPE))
    zea.init_device()
    report={'study':cfg['study'],'scope':cfg['status'],
            'source_sha256':cfg['source_sha256'],
            'sources':{str(p.relative_to(ROOT)):digest(p) for p in
                       (CFG,PIPE,Path(__file__),ROOT/'src/btbench/multidim.py',
                        ROOT/'scripts/resolvestroke_single_view.py')},'clips':[]}
    with h5py.File(source,'r',locking=False) as h,File(str(source)) as zf:
        legacy=json.loads((ROOT/'configs/resolvestroke_single_view.json').read_text())
        raw=check_source(h,legacy)
        grid,mask,audit=crop(h,cfg)
        params=zf.load_parameters(**config.parameters)
        pipeline=Pipeline.from_config(config)
        prepared=pipeline.prepare_parameters(params,grid=grid,flatgrid=grid.reshape(-1,3))
        args={k:v for k,v in prepared.items() if k!='grid'}
        geom=np.asarray(prepared['probe_geometry'])
        sort=np.lexsort((geom[:,0],geom[:,1]))
        if len(sort)!=256:raise ValueError('Not 256 receive elements')
        k=cfg['receive_groups']
        tgc=np.asarray(params.tgc_gain_curve,dtype=np.float32)[None,:,None,None]
        scores_roi={str(label):rois(grid,mask,cfg,label) for label in (1,2)}
        for clip in cfg['clips']:
            frames=[]
            for j in range(clip['start'],clip['start']+cfg['frames']):
                channel=np.asarray(raw[j],np.float32)/tgc
                aligned=np.asarray(tof_correction(channel,**args))[0]
                grouped=(aligned[:,sort,:].reshape(len(grid)*len(grid[0]),k,256//k,2)
                         .sum(axis=2,dtype=np.float32))
                frames.append(grouped[...,0]+1j*grouped[...,1])
            rx=np.stack(frames)
            first=np.asarray(raw[clip['start']:clip['start']+1],np.float32)/tgc[None]
            direct=np.asarray(pipeline(**{pipeline.key:first},**prepared,
                                        return_numpy=True)[pipeline.output_key])[0]
            direct=direct[...,0].ravel()+1j*direct[...,1].ravel()
            closure=float(np.linalg.norm(rx[0].sum(-1)-direct)/np.linalg.norm(direct))
            if closure>1e-4:raise AssertionError(f'Receive groups do not reproduce Zea DAS: {closure}')
            methods={'raw':rx.sum(-1),
                'compounded_svd_rank1':svd_unfolding(rx.sum(-1)).apply(rx.sum(-1)),
                'ozgun_svd_a_rank1':svd_unfolding(rx).apply(rx).sum(-1),
                'ozgun_hosvd_adaptation':paper_tensor(rx,family='ozgun').apply(rx).sum(-1),
                'receive_core_intersection':intersection_filter(rx).apply(rx).sum(-1)}
            rows={}
            for name,filtered in methods.items():
                power=np.mean(abs(filtered)**2,axis=0).reshape(mask.shape)
                rows[name]={label:score_power(power,*roi) for label,roi in scores_roi.items()}
                if name!='raw':
                    for label in rows[name]:
                        rows[name][label]['outside_power_relative_to_raw_db']=float(10*np.log10(
                            max(rows[name][label]['mean_outside_power'],1e-30)/
                            max(rows['raw'][label]['mean_outside_power'],1e-30)))
            report['clips'].append({'clip':clip['name'],'start':clip['start'],
              'number_of_transmissions':cfg['frames'],'receive_group_count':k,
              'beam_sum_relative_error':closure,'geometry':audit,'scores':rows})
            print(clip['name'],'DAS closure',closure,'tube AUC',
                  [(name,round(row['1']['auc'],3),round(row['2']['auc'],3)) for name,row in rows.items()],flush=True)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    run(args.source,args.out)
