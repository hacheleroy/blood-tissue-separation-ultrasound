"""Record phase-10 retrospective screening on eight explicitly seen clouds.

Seed 2729 uses verified RF checkpoint IQ replay; those observations are now
seen development data and cannot become prospective confirmation.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from btbench.shared import PublicInput, fit_shared
from esdiva_example22_pilot import digest, relative_error
from esdiva_phase10_support import (
    build_phase10_context, choose_phase10_mode, phase10_support,
)


def ratio(a,b):
    return float(np.mean(abs(a)**2)/max(np.mean(abs(b)**2),1e-30))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_scene(seed, scene, phase6, phase7, phase8, phase9, replay):
    folder = ((phase6 if seed in (2309,2311) else
               phase7 if seed in (2503,2521) else
               phase8 if seed==2707 else phase9)/str(seed)/scene
              if seed != 2729 else replay/scene)
    with np.load(folder/'iq_sources.npz') as source,np.load(folder/'broad_iq.npz') as broad:
        core = {k:source[k] for k in ('views','calibration_views','blood','tissue','noise')}
        wide = {k:broad[k] for k in ('views','calibration_views','blood','tissue',
                                     'noise','z_mm','x_mm','core_rows')}
    if seed == 2729:
        audit = json.loads((folder/'replay.json').read_text())
        if (audit['script_sha256']!=sha(Path(__file__).with_name(
                 'esdiva_phase9_replay_seen.py')) or
                any(digest(core[k])!=audit['core_sha256'][k] for k in core) or
                any(digest(wide[k])!=audit['wide_sha256'][k] for k in wide)):
            raise ValueError('Seen 2729 RF/IQ replay mismatch')
        geo = json.loads((phase8/'2729/smooth_zero/report.json').read_text())['geometry']
        geo['f0_hz'] = audit['carrier_hz']
        provenance = {'replay_audit_sha256':sha(folder/'replay.json'),
                      'source':'previously seen RF checkpoint diagnostic'}
    else:
        report_path = folder/'report.json'
        report = json.loads(report_path.read_text())
        geo = report['geometry']
        provenance = {'acquisition_sha256':sha(report_path),
                      'source':'previously scored source-resolved development'}
    if any(not np.array_equal(wide[k][:,:,wide['core_rows']],core[k]) for k in core):
        raise ValueError('Wider and original scored IQ differ')
    wide_roi,core_roi,calibration = phase10_support(
        wide['calibration_views'],wide['z_mm'],wide['x_mm'],wide['core_rows'])
    core_contamination = ratio(core['blood'][:,:,core_roi],core['tissue'][:,:,core_roi])
    wide_contamination = ratio(wide['blood'][:,:,wide_roi],wide['tissue'][:,:,wide_roi])
    observations = {'mixture':(core['views'],wide['views']),
                    'background':(core['tissue']+core['noise'],
                                  wide['tissue']+wide['noise'])}
    inputs = {label:(PublicInput(cv,core['calibration_views'],core_roi,geo),
                     PublicInput(wv,wide['calibration_views'],wide_roi,geo))
              for label,(cv,wv) in observations.items()}
    selected,selector = choose_phase10_mode(
        *inputs['mixture'],wide['core_rows'],wide['z_mm'],wide['x_mm'])
    specs = {'roundtrip':{'method':'roundtrip','horizon':64},
             'svd_cut1':{'method':'svd','discard':1,'horizon':64},
             'axis1_cut1':{'method':'axis1','discard':1,'horizon':64,'rho':.75}}
    outputs = {}
    for mode in ('none',selected):
        contexts = {label:build_phase10_context(
            mode,ci,wi,wide['core_rows'],wide['z_mm'],wide['x_mm'])
            for label,(ci,wi) in inputs.items()}
        methods = {}
        names = dict(specs)
        if mode=='none':
            names['axis2_slow_r4_cap040']={
                'method':'axis2','rank':4,'gain_cap':.4,'horizon':64}
        for name,spec in names.items():
            mixture,background = (fit_shared(contexts[label],spec)
                                  for label in ('mixture','background'))
            methods[name] = {
                'paired_added_blood_iq_nrmse':relative_error(
                    mixture.views-background.views,core['blood'][64:]),
                'refitted_background_output_to_input_tissue_power':ratio(
                    background.views,core['tissue'][64:]),
            }
        outputs[mode] = methods
    bound = {'blood_delta_max':.02,'moving_tissue_ratio_max':.9,
             'roundtrip_max':1e-10,'slow_axis2_blood_max':.7,
             'slow_axis2_tissue_max':.18}
    failures = []
    for mode,methods in outputs.items():
        if methods['roundtrip']['paired_added_blood_iq_nrmse']>bound['roundtrip_max']:
            failures.append('roundtrip '+mode)
    comparisons = {}
    for method in ('svd_cut1','axis1_cut1'):
        base,corrected = outputs['none'][method],outputs[selected][method]
        delta = (corrected['paired_added_blood_iq_nrmse']-
                 base['paired_added_blood_iq_nrmse'])
        leakage = (corrected['refitted_background_output_to_input_tissue_power']/
                   base['refitted_background_output_to_input_tissue_power'])
        comparisons[method]={'blood_delta':delta,'tissue_ratio':leakage}
        if delta>bound['blood_delta_max']:
            failures.append('blood '+method)
        if scene in ('rigid_fast','smooth_slow') and leakage>bound['moving_tissue_ratio_max']:
            failures.append('leakage '+method)
    if scene=='smooth_slow':
        a=outputs['none']['axis2_slow_r4_cap040']
        if (a['paired_added_blood_iq_nrmse']>bound['slow_axis2_blood_max'] or
                a['refitted_background_output_to_input_tissue_power']>
                bound['slow_axis2_tissue_max']):
            failures.append('slow Axis 2')
    if not calibration['pass'] or max(core_contamination,wide_contamination)>.1:
        failures.append('input/source')
    return {'seed':seed,'scene':scene,'provenance':provenance,
            'calibration_support':calibration,'core_blood_to_tissue':core_contamination,
            'wide_blood_to_tissue':wide_contamination,'selected_mode':selected,
            'selection':selector,'methods':outputs,'comparisons':comparisons,
            'failures':failures}


def screen(phase6,phase7,phase8,phase9,replay,output):
    phase6,phase7,phase8,phase9,replay,output=map(lambda p:Path(p).resolve(),
                                           (phase6,phase7,phase8,phase9,replay,output))
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for seed in (2309,2311,2503,2521,2707,2729,3001,3019):
        for scene in ('stationary_fast','rigid_fast','smooth_slow','smooth_zero'):
            row=run_scene(seed,scene,phase6,phase7,phase8,phase9,replay)
            rows.append(row)
            print(seed,scene,row['selected_mode'],row['failures'],flush=True)
    report={'status':'retrospective phase-10 development screening only',
            'seen_seeds':[2309,2311,2503,2521,2707,2729,3001,3019],
            'source_sha256':sha(__file__),'support_sha256':sha(Path(__file__).with_name(
                'esdiva_phase10_support.py')),
            'replay_sha256':sha(Path(__file__).with_name(
                'esdiva_phase9_replay_seen.py')),
            'no_new_seed_opened':True,'rows':rows,
            'failures':[f"{r['seed']}/{r['scene']}: {f}"
                        for r in rows for f in r['failures']]}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('phase6','phase7','phase8','phase9','replay','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    a=parser.parse_args()
    result=screen(a.phase6,a.phase7,a.phase8,a.phase9,a.replay,a.output)
    print(json.dumps({'rows':len(result['rows']),'failures':result['failures']}))
