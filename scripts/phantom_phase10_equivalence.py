"""Verify the physical geometry port against frozen eSDIva phase-10 IQ.

Provide local acquisition folders with report.json, iq_sources.npz and
broad_iq.npz. Only a compact hash-bound audit is written; raw arrays stay local.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'examples')]
from btbench.shared import PublicInput,fit_shared,array_hash
from esdiva_phase10_support import phase10_support,choose_phase10_mode,build_phase10_context
from phantom_phase10_port import structural_support,calibration_support,choose_mode,build_context


def digest(path):return sha256(Path(path).read_bytes()).hexdigest()


def compare(name,folder):
    folder=Path(folder).resolve()
    report=json.loads((folder/'report.json').read_text())
    with np.load(folder/'iq_sources.npz',allow_pickle=False) as c,\
         np.load(folder/'broad_iq.npz',allow_pickle=False) as w:
        core={k:c[k] for k in ('views','calibration_views','tissue','blood','noise','z_mm','x_mm')}
        wide={k:w[k] for k in ('views','calibration_views','tissue','blood','noise','z_mm','x_mm','core_rows')}
    z,x=wide['z_mm'],wide['x_mm']
    axis_x=np.tan(np.deg2rad(60.))*(z-16.)
    wide_roi,core_roi=structural_support(z,x,axis_x)
    old_wide,old_core,old_gate=phase10_support(
        wide['calibration_views'],z,x,wide['core_rows'])
    if (not np.array_equal(wide_roi,old_wide) or
            not np.array_equal(core_roi,old_core) or
            calibration_support(wide['calibration_views'],wide_roi,core_roi)!=old_gate):
        raise AssertionError(f'{name}: structural support or calibration gate differs')
    ci=PublicInput(core['views'],core['calibration_views'],core_roi,report['geometry'])
    wi=PublicInput(wide['views'],wide['calibration_views'],wide_roi,report['geometry'])
    old_mode,old_audit=choose_phase10_mode(ci,wi,wide['core_rows'],z,x)
    new_mode,new_audit=choose_mode(ci,wi,wide_roi,core_roi)
    if old_mode!=new_mode or old_audit!=new_audit:
        raise AssertionError(f'{name}: selected mode or audit differs')
    background_core=core['tissue']+core['noise']
    background_wide=wide['tissue']+wide['noise']
    inputs={
        'mixture':(ci,wi),
        'background':(
            PublicInput(background_core,core['calibration_views'],core_roi,report['geometry']),
            PublicInput(background_wide,wide['calibration_views'],wide_roi,report['geometry']))}
    rows=[];outputs={}
    specs={'roundtrip':{'method':'roundtrip','horizon':64},
           'svd_cut1':{'method':'svd','horizon':64,'discard':1},
           'axis1_cut1':{'method':'axis1','horizon':64,'discard':1,'rho':.75}}
    for label,(current_core,current_wide) in inputs.items():
        for mode in ('none','hybrid','reference'):
            old=build_phase10_context(mode,current_core,current_wide,wide['core_rows'],z,x)
            new=build_context(mode,current_core,current_wide,wide_roi,core_roi)
            if (not np.allclose(old.modulation,new.modulation,atol=1e-13,rtol=0) or
                    not np.allclose(old.aligned,new.aligned,atol=1e-13,rtol=0)):
                raise AssertionError(f'{name}: {label}/{mode} context differs')
            for method,spec in specs.items():
                a=fit_shared(old,spec).views;b=fit_shared(new,spec).views
                error=float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-30))
                if error>1e-12:
                    raise AssertionError(f'{name}: {label}/{mode}/{method} error {error}')
                outputs[(label,mode,method)]=(a,b)
                rows.append({'input':label,'mode':mode,'method':method,
                             'relative_output_error':error,'output_sha256':array_hash(b)})
    paired_rows=[]
    for mode in ('none','hybrid','reference'):
        for method in ('svd_cut1','axis1_cut1'):
            m0,m1=outputs[('mixture',mode,method)]
            b0,b1=outputs[('background',mode,method)]
            blood=core['blood'][64:]
            old_error=float(np.linalg.norm((m0-b0)-blood)/max(np.linalg.norm(blood),1e-30))
            new_error=float(np.linalg.norm((m1-b1)-blood)/max(np.linalg.norm(blood),1e-30))
            if abs(old_error-new_error)>1e-12:
                raise AssertionError(f'{name}: paired blood {mode}/{method} differs')
            paired_rows.append({'mode':mode,'method':method,
                                'paired_added_blood_nrmse_original':old_error,
                                'paired_added_blood_nrmse_port':new_error,
                                'absolute_difference':abs(old_error-new_error)})
    return {'case':name,'seed':report['seed'],'motion':report['motion'],
            'flow_scale':report['flow_scale'],
            'source_report_sha256':digest(folder/'report.json'),
            'core_iq_sha256':digest(folder/'iq_sources.npz'),
            'wide_iq_sha256':digest(folder/'broad_iq.npz'),
            'selected_mode':new_mode,'wide_pixels':int(wide_roi.sum()),
            'core_pixels':int(core_roi.sum()),'calibration_gate_pass':bool(old_gate['pass']),
            'rows':rows,'paired_rows':paired_rows}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',nargs=2,action='append',metavar=('NAME','FOLDER'),required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    if len({name for name,_ in a.case})!=len(a.case):raise ValueError('Case names must be unique')
    result={'status':'physical port matches frozen phase-10 numerical source on supplied IQ',
            'script_sha256':digest(__file__),
            'port_sha256':digest(ROOT/'examples/phantom_phase10_port.py'),
            'frozen_support_sha256':digest(ROOT/'examples/esdiva_phase10_support.py'),
            'frozen_selector_sha256':digest(ROOT/'examples/esdiva_phase8_selector.py'),
            'frozen_transport_sha256':digest(ROOT/'examples/esdiva_phase7_reference_transport.py'),
            'cases':[compare(name,path) for name,path in a.case],
            'limitation':'Engineering equivalence on Example 22 IQ, not physical phantom validation.'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'cases':len(result['cases']),
                      'modes':[row['selected_mode'] for row in result['cases']]}))


if __name__=='__main__':main()
