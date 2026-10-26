"""Verify the published, source-linked contrast reanalysis without raw RF."""
import csv
import hashlib
import io
import json
import lzma
import math
import argparse
from collections import Counter
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/flow_background_audit_v1'


def sha(data):return hashlib.sha256(data).hexdigest()


def source_bytes(relative):
    path=ROOT/relative
    if path.exists():return path.read_bytes()
    parts=sorted(path.parent.glob(path.name+'.part*'))
    if not parts:raise FileNotFoundError(path)
    return b''.join(p.read_bytes() for p in parts)


def check(require_all_sources=False):
    overview=json.loads((OUT/'summary.json').read_text())
    assert sha((ROOT/'scripts/flow_background_audit.py').read_bytes())==overview['analysis_script_sha256']
    unpublished=set(overview['unpublished_source_paths'])
    assert unpublished <= set(overview['sources'])
    matched=0;unavailable=[]
    for path,digest in overview['sources'].items():
        try: content=source_bytes(path)
        except FileNotFoundError:
            if path not in unpublished:raise
            unavailable.append(path)
            continue
        assert sha(content)==digest,path
        matched+=1
    if require_all_sources and unavailable:
        raise FileNotFoundError('Source archives unavailable: '+', '.join(unavailable))
    compressed=(OUT/'rows.csv.xz').read_bytes()
    assert sha(compressed)==overview['rows_xz_sha256']
    plain=lzma.decompress(compressed)
    assert sha(plain)==overview['rows_csv_sha256']
    rows=list(csv.DictReader(io.StringIO(plain.decode())))
    assert len(rows)==overview['row_count']
    assert dict(sorted(Counter(r['metric_kind'] for r in rows).items()))==overview['status_counts']
    exact=0;bounded=0
    for row in rows:
        db=row['sbr_db'];low=row['sbr_lower_db'];high=row['sbr_upper_db']
        if low and high:
            assert float(low)<=float(high)+1e-8
            bounded+=1
        if row['metric_kind'].startswith('exact_') and db:
            q=float(row['blood_power_retention']);bg=float(row['background_power_remaining'])
            rho=float(row['input_blood_to_tissue_power'])
            if q>0 and bg>0:
                assert abs(float(db)-10*math.log10(q*rho/bg))<1e-8
                exact+=1
    assert exact>15000 and bounded>500
    physical=(OUT/'physical_tube_proxy.csv').read_bytes()
    assert sha(physical)==overview['physical_tube_proxy_sha256']
    assert sum(1 for _ in csv.DictReader(io.StringIO(physical.decode())))==overview['physical_tube_proxy_rows']
    for artifact,script in [('multidim_exact.json','recompute_multidim_contrast.py'),
                            ('esdiva_paired_exact.json','recompute_esdiva_paired_contrast.py')]:
        data=json.loads((OUT/artifact).read_text())
        assert sha((ROOT/'scripts'/script).read_bytes())==data['audit_script_sha256']
    return {'rows':len(rows),'exact_checked':exact,'bounded_checked':bounded,
            'source_archives_verified':matched,'source_archives_unavailable':len(unavailable),
            'physical_proxy_rows':overview['physical_tube_proxy_rows']}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-all-sources',action='store_true')
    args=parser.parse_args()
    print(json.dumps(check(require_all_sources=args.require_all_sources),indent=2))
