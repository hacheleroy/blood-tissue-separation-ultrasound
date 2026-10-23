"""Losslessly compress frozen JSON and split it for reliable transfer."""
import gzip
import hashlib
import json
import lzma
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];items=[];code_hash=None
    for directory,kind in [('fidelity_followup_v04','synthetic'),('openhrf_followup_v04','experimental')]:
        original=root/'runs'/directory/'results.json.gz';raw=gzip.decompress(original.read_bytes())
        path=original.with_suffix('.xz');data=lzma.compress(raw,preset=6)
        report=json.loads(raw);records=report['records'];rows=[r for c in records for r in c['results']]
        if code_hash is None:code_hash=report['code_sha256']
        if report['code_sha256']!=code_hash or any(c['code_sha256']!=code_hash for c in records):
            raise ValueError('Archives do not share one frozen numerical source')
        if any(r['status']!='ok' for r in rows):raise ValueError('Unresolved method errors or skips in follow-up')
        config=path.parent/'frozen_config.json';config_sha=hashlib.sha256(config.read_bytes()).hexdigest()
        if report['config_sha256']!=config_sha:raise ValueError('Frozen configuration differs from archived fingerprint')
        parts=[]
        for i,start in enumerate(range(0,len(data),450000)):
            p=path.with_name(path.name+f'.part{i:03d}');value=data[start:start+450000]
            if p.exists() and p.read_bytes()!=value:raise FileExistsError('Different archive part already exists: '+str(p))
            p.write_bytes(value)
            parts.append({'path':p.relative_to(root).as_posix(),'bytes':len(value),'sha256':hashlib.sha256(value).hexdigest()})
        item={'path':path.relative_to(root).as_posix(),'codec':'xz','sha256':hashlib.sha256(data).hexdigest(),
              'decoded_sha256':hashlib.sha256(raw).hexdigest(),'outputs':len(rows),'cases':len(records),
              'config_path':config.relative_to(root).as_posix(),'config_sha256':config_sha,'parts':parts,
              'independent_source_truth_available':kind=='synthetic'}
        if kind=='synthetic':
            item['required_metrics']=['blood_iq_nrmse','refitted_background_to_tissue_power',
                'refitted_blood_response_nrmse_alpha1','refitted_blood_response_nrmse_alpha0.1']
        items.append(item)
    manifest={'version':'0.3.0','date':'2026-10-23','code_sha256':code_hash,'archives':items,
        'interpretation':'Actual output error and refitted background/injection responses are primary; frozen attribution is representation-dependent. Experimental truth is unavailable.'}
    target=root/'runs/followup_v04_manifest.json';target.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'archives':len(items),'outputs':sum(i['outputs'] for i in items),
                      'parts':sum(len(i['parts']) for i in items),'code_sha256':code_hash}),flush=True)


if __name__=='__main__':main()
