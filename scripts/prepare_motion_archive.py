"""Create byte-verified portable parts of the frozen v06 study."""
import gzip
import hashlib
import json
import lzma
from pathlib import Path


def sha(data):return hashlib.sha256(data).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1];run=root/'runs/motion_repair_v06'
    raw=gzip.decompress((run/'results.json.gz').read_bytes());report=json.loads(raw)
    config=(run/'frozen_config.json').read_bytes()
    assert sha(config)==report['config_sha256']
    assert len(report['records'])==64 and all(len(c['results'])==36 for c in report['records'])
    compressed=lzma.compress(raw,preset=6);parts=[]
    for i,start in enumerate(range(0,len(compressed),450000)):
        path=run/f'results.json.xz.part{i:03d}';chunk=compressed[start:start+450000]
        if path.exists() and path.read_bytes()!=chunk:raise FileExistsError(path)
        path.write_bytes(chunk);parts.append(dict(path=path.relative_to(root).as_posix(),bytes=len(chunk),sha256=sha(chunk)))
    artifacts=[]
    for path in [run/'metrics.csv',run/'summary.json',root/'docs/RESULTS_MOTION_REPAIR_V06.md',
                 root/'docs/figures/motion_repair_v06.png',root/'docs/figures/motion_repair_v06.svg']:
        data=path.read_bytes();artifacts.append(dict(path=path.relative_to(root).as_posix(),bytes=len(data),sha256=sha(data)))
    item=dict(path='runs/motion_repair_v06/results.json.xz',codec='xz',sha256=sha(compressed),decoded_sha256=sha(raw),
        cases=64,outputs=2304,config_path='runs/motion_repair_v06/frozen_config.json',config_sha256=sha(config),
        independent_source_truth_available=True,parts=parts,required_metrics=['blood_iq_nrmse',
            'refitted_background_to_tissue_power','refitted_blood_response_nrmse_alpha1',
            'refitted_blood_response_nrmse_alpha0_1'])
    manifest=dict(study='motion_repair_v06',version='0.5.0',code_sha256=report['code_sha256'],
        archives=[item],artifacts=artifacts,provenance='Synthetic angular IQ only; no empirical RF or phantom data.',
        interpretation='New reserved held-out bank; truth and known warp are evaluation-only. Gate fixed before any held-out draw.')
    (root/'runs/motion_v06_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(cases=64,outputs=2304,raw_bytes=len(raw),compressed_bytes=len(compressed),parts=len(parts))))


if __name__=='__main__':main()
