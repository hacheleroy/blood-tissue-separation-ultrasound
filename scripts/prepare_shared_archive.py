"""Archive only generated synthetic records and their inspectable derivatives."""
import gzip
import hashlib
import json
import lzma
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];directory=root/'runs/shared_axes_v05'
    raw=gzip.decompress((directory/'results.json.gz').read_bytes());report=json.loads(raw)
    config=directory/'frozen_config.json';cfg=config.read_bytes();rows=[r for c in report['records'] for r in c['results']]
    if report['config_sha256']!=hashlib.sha256(cfg).hexdigest():raise ValueError('Config identity changed')
    if any(c['code_sha256']!=report['code_sha256'] for c in report['records']):raise ValueError('Mixed source identities')
    if any(r['status']!='ok' for r in rows):raise ValueError('Unresolved filter failures')
    path=directory/'results.json.xz';data=lzma.compress(raw,preset=6);parts=[]
    for i,start in enumerate(range(0,len(data),450000)):
        part=path.with_name(path.name+f'.part{i:03d}');value=data[start:start+450000]
        if part.exists() and part.read_bytes()!=value:raise FileExistsError('Existing archive part differs')
        part.write_bytes(value);parts.append(dict(path=part.relative_to(root).as_posix(),bytes=len(value),sha256=hashlib.sha256(value).hexdigest()))
    artifacts=[]
    for p in sorted([*directory.glob('*.csv'),*list((directory/'metrics').glob('*.csv')),directory/'summary.json',directory/'SUMMARY.md',
                      *list((root/'docs/figures').glob('shared_*.png')),*list((root/'docs/figures').glob('shared_*.svg'))]):
        b=p.read_bytes();artifacts.append(dict(path=p.relative_to(root).as_posix(),bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
    item=dict(path=path.relative_to(root).as_posix(),codec='xz',sha256=hashlib.sha256(data).hexdigest(),
        decoded_sha256=hashlib.sha256(raw).hexdigest(),outputs=len(rows),cases=len(report['records']),
        config_path=config.relative_to(root).as_posix(),config_sha256=hashlib.sha256(cfg).hexdigest(),parts=parts,
        independent_source_truth_available=True,required_metrics=['blood_iq_nrmse','refitted_background_to_tissue_power',
            'refitted_blood_response_nrmse_alpha1','refitted_blood_response_nrmse_alpha0_1'])
    manifest=dict(version='0.4.0',study='shared_axes_v05',date='2026-10-23',code_sha256=report['code_sha256'],
        archives=[item],artifacts=artifacts,provenance='Only software-generated synthetic angular-IQ metrics and plots; no empirical arrays, private documents or source PDFs.',
        interpretation='Actual output and refitted responses are primary. Known motion and truth-fitted diagonal weights are evaluation-only; global operating points use training seeds only.')
    (root/'runs/shared_v05_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(cases=item['cases'],outputs=len(rows),raw_bytes=len(raw),xz_bytes=len(data),parts=len(parts),artifacts=len(artifacts))))


if __name__=='__main__':main()
