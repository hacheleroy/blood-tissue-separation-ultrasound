"""Re-evaluate the paired phase-10/11 score arrays without changing frozen files.

The historical code hashed all then-existing btbench modules. Verify that
exact historical hash after excluding subsequently added modules, then supply
that verified value to its unchanged source guard for this read-only replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import esdiva_phase10_stress_score as stress
import esdiva_phase11_prospective_score as prospective


ROOT=Path(__file__).resolve().parents[1]
EXTRA={'motion_first.py','native_subtraction.py','stable_tissue.py',
       'multidim.py','esdiva_aperture.py'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def historical_numerical_hash():
    h=hashlib.sha256()
    for path in sorted((ROOT/'src/btbench').glob('*.py')):
        if path.name not in EXTRA:
            h.update(path.name.encode());h.update(path.read_bytes())
    return h.hexdigest()


def ratio(a,b):
    return float(np.vdot(a,a).real/np.vdot(b,b).real)


def case(module,bank,root,seed,scene):
    source=ROOT/'runs'/bank/'fresh'/str(seed)/scene/'score.json'
    original=json.loads(source.read_text())
    acq=root/str(seed)/scene
    config=ROOT/'configs'/('esdiva_phase10_stress.json' if module is stress
                            else 'esdiva_phase11_prospective.json')
    expected=json.loads(config.read_text())['numerical_source_sha256']
    if historical_numerical_hash()!=expected:
        raise ValueError('Original numerical source has changed')
    module.numerical_source_hash=lambda parent:expected
    original_fit=module.fit_shared
    captured=[]
    def capture(*args,**kwargs):
        result=original_fit(*args,**kwargs)
        captured.append(result.views)
        return result
    try:
        module.fit_shared=capture
        with TemporaryDirectory() as temp:
            replay=module.evaluate(acq,acq/'broad_iq.npz',config,Path(temp)/'score.json')
    finally:
        module.fit_shared=original_fit
    if original['status']!=replay['status'] or original['selected_mode']!=replay['selected_mode']:
        raise ValueError('Changed score/selection')
    with np.load(acq/'iq_sources.npz') as f:
        length=original['score_frames']
        b=f['blood'][-length:]
        t=f['tissue'][-length:]
    rho=ratio(b,t)
    order=[(variant,method) for variant,group in replay['variants'].items()
           for method in group['methods']]
    if len(captured)!=2*len(order):
        raise ValueError('Unpaired mixture/background replay')
    rows=[]
    for i,(variant,method) in enumerate(order):
        archived=original['variants'][variant]['methods'][method]
        result=replay['variants'][variant]['methods'][method]
        for field in ('paired_added_blood_iq_nrmse',
                      'refitted_background_output_to_input_tissue_power',
                      'frozen_blood_power_ratio'):
            assert np.isclose(archived[field],result[field],rtol=1e-8,atol=1e-9),(source,field)
        out,base=captured[2*i:2*i+2]
        increment=out-base
        q=ratio(increment,b)
        background=ratio(base,t)
        error=float(np.linalg.norm(increment-b)/np.linalg.norm(b))
        assert np.isclose(error,archived['paired_added_blood_iq_nrmse'],atol=1e-8)
        assert np.isclose(background,archived['refitted_background_output_to_input_tissue_power'],atol=1e-8)
        rows.append({'bank':bank,'seed':seed,'scene':scene,'variant':variant,
                     'method':method,'blood_iq_error':error,'blood_power_retention':q,
                     'background_power_remaining':background,'input_blood_to_tissue_power':rho,
                     'frozen_blood_power_retention_comparator':archived['frozen_blood_power_ratio']})
    return rows,{'score':sha(source),'iq_sources':sha(acq/'iq_sources.npz'),
                 'broad_iq':sha(acq/'broad_iq.npz')}


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args()
    banks=[(stress,'esdiva_phase10_stress',ROOT.parent/'esdiva_stress_work'),
           (prospective,'esdiva_phase11_prospective',ROOT.parent/'esdiva_phase11_work/acquisitions')]
    rows=[];sources={}
    for module,bank,base in banks:
        for path in sorted((ROOT/'runs'/bank/'fresh').glob('*/*/score.json')):
            seed=int(path.parent.parent.name);scene=path.parent.name
            result,digests=case(module,bank,base,seed,scene)
            rows.extend(result);sources[f'{bank}/{seed}/{scene}']=digests
            print(bank,seed,scene,len(result),flush=True)
    output={'scope':'paired, re-fitted source counterfactual on available acoustic RF IQ',
            'audit_script_sha256':sha(__file__), 'source_hashes':sources,'rows':rows}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print('rows',len(rows),'sha256',sha(a.out))


if __name__=='__main__':main()
