"""Check v06 split, controls and independent gate audit in portable archive."""
import json
import lzma
from pathlib import Path
from btbench.shared_experiment import validate_shared_config,select_global_operating_points


def verify(root):
    root=Path(root);manifest=json.loads((root/'runs/motion_v06_manifest.json').read_text())
    item=manifest['archives'][0]
    report=json.loads(lzma.decompress(b''.join((root/p['path']).read_bytes() for p in item['parts'])))
    config=validate_shared_config(json.loads((root/item['config_path']).read_text()))
    assert config['training_seeds']==[601,787,919]
    assert config['evaluation_seeds']==[1601,1907,2203,2503,2801]
    expected={(c['id'],seed,split) for c in config['conditions'] for split,key in
              [('training','training_seeds'),('evaluation','evaluation_seeds')] for seed in config[key]}
    assert {(c['condition'],c['seed'],c['split']) for c in report['records']}==expected
    assert len(report['records'])==len(expected)==64
    specs={m['id']:m for m in config['methods']};assert len(specs)==36
    for c in report['records']:
        assert c['code_sha256']==report['code_sha256'] and c['config_sha256']==report['config_sha256']
        assert len(c['results'])==36 and {r['id'] for r in c['results']}==set(specs)
        assert c['contexts']['known']['actual']['evaluation_only_displacement_rmse_pixels']==0
        for mode,contexts in c['contexts'].items():
            assert contexts['actual']['input_sha256']==c['input_sha256']
            assert len({json.dumps(v['calibration'],sort_keys=True) for v in contexts.values()})==1
            if mode in ('rigid','smooth'):
                for context in contexts.values():
                    audit=context['motion_fit'];counts=audit['reason_counts'];rows=audit['audit_rows']
                    assert audit['gate_enabled'] and sum(counts.values())==128 and len(rows)==128
                    assert audit['accepted_frames']==sum(row[0]!=0 for row in rows)
                    assert all(row[0] in (0,1,2) and row[1]>=0 for row in rows)
        for r in c['results']:
            assert r['status']=='ok' and r['evaluation_only']==specs[r['id']]['evaluation_only']
            assert r['diagnostics']['latency_max_frames']==config['packet_frames']
    assert select_global_operating_points(report['records'],config['methods'],config['background_targets'])==report['operating_points']
    return dict(cases=len(report['records']),outputs=sum(len(c['results']) for c in report['records']),gate_audit='verified')


if __name__=='__main__':print(json.dumps(verify(Path(__file__).resolve().parents[1]),indent=2))
