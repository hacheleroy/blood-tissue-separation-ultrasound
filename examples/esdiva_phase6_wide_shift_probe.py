"""Evaluator-only development control: match no-motion and shift on 32 rows.

The known scatterer displacement is read only *after* both data-fitted outputs
are computed. This probe is not one of the bounded candidate trials.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import relative_error
from esdiva_phase3_bounded_context import Phase3Context
from esdiva_phase5_local_phase import tissue_mask
from esdiva_phase6_unit_shift import UnitShiftContext


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(acquisition, output):
    acquisition,output=Path(acquisition).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    report_path=acquisition/'report.json'
    broad_path=acquisition/'broad_iq.npz'
    audit_path=acquisition/'broad_iq.json'
    report=json.loads(report_path.read_text())
    audit=json.loads(audit_path.read_text())
    if audit['source_report_sha256']!=sha(report_path) or not audit['rf_sources_and_noise_verified']:
        raise ValueError('Unverified source-resolved RF replay')
    with np.load(broad_path) as d,np.load(acquisition/'iq_sources.npz') as core:
        wide={k:d[k] for k in ('views','calibration_views','tissue','blood','noise','z_mm','x_mm','core_rows')}
        known=core['known_scatterer_displacement_pixels']
    rows=wide['core_rows']
    if not np.array_equal(rows,np.arange(12,20)):
        raise ValueError('Original score interval is not centered in the wide tile')
    roi=tissue_mask(wide['z_mm'],wide['x_mm'])
    observations={'mixture':wide['views'],'background':wide['tissue']+wide['noise']}
    contexts={}
    for label,views in observations.items():
        public=PublicInput(views,wide['calibration_views'],roi,report['geometry'])
        contexts[label]={
            'none':Phase3Context(public,packet=16,score_start=64,motion='none'),
            'wide_unit_shift':UnitShiftContext(public,public,np.arange(32),
                            wide['z_mm'],wide['x_mm'],packet=16,score_start=64)}
    specs={'roundtrip':{'method':'roundtrip','horizon':64},
           'svd_cut1':{'method':'svd','horizon':64,'discard':1},
           'axis1_cut1':{'method':'axis1','horizon':64,'discard':1,'rho':.75}}
    variants={}
    for variant in ('none','wide_unit_shift'):
        methods={}
        for name,spec in specs.items():
            mixed,bg=(fit_shared(contexts[label][variant],spec) for label in
                      ('mixture','background'))
            methods[name]={
                'paired_added_blood_iq_nrmse':relative_error(
                    mixed.views[:,:,rows]-bg.views[:,:,rows],wide['blood'][64:,:,rows]),
                'refitted_background_output_to_input_tissue_power':float(
                    np.mean(abs(bg.views[:,:,rows])**2)/
                    max(np.mean(abs(wide['tissue'][64:,:,rows])**2),1e-30))}
        variants[variant]=methods
    shift=contexts['mixture']['wide_unit_shift'].shifts.mean(axis=1)
    truth=-known.mean(axis=(-2,-1))
    evaluator={axis:{'correlation':float(np.corrcoef(shift[:,i],truth[:,i])[0,1]),
                     'rmse_pixels':float(np.sqrt(np.mean((shift[:,i]-truth[:,i])**2)))}
               for i,axis in enumerate(('axial','lateral'))}
    result={'status':'rejected/eligible only as development diagnostic, no candidate selection',
            'seed':report['seed'],'motion':report['motion'],'flow_scale':report['flow_scale'],
            'report_sha256':sha(report_path),'wide_audit_sha256':sha(audit_path),
            'script_sha256':sha(__file__),
            'unit_shift_context_sha256':sha(Path(__file__).with_name('esdiva_phase6_unit_shift.py')),
            'variants':variants,'known_scatterer_motion_evaluator_only':evaluator,
            'shift_fit_audit':contexts['mixture']['wide_unit_shift'].motion_summary,
            'limitation':'cyclic Fourier translation on 32 rows; independent physical motion unavailable'}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--acquisition',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();r=run(a.acquisition,a.output)
    print(json.dumps({'seed':r['seed'],'status':r['status'],
                      'evaluator_motion':r['known_scatterer_motion_evaluator_only']}))
