"""Verify frozen RF study identity, declared coverage and measured pilot."""
import argparse
import hashlib
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(rf_dir,measured=None):
    rf_dir=Path(rf_dir)
    config=json.loads((ROOT/'configs/multidim_aperture_v1.json').read_text())
    freeze=json.loads((rf_dir/'freeze.json').read_text())
    dev=json.loads((rf_dir/'development.json').read_text())
    eva=json.loads((rf_dir/'evaluation.json').read_text())
    assert sha(rf_dir/'development.json')==freeze['development_sha256']
    assert freeze['evaluation_seeds']==config['evaluation_seeds']
    assert dev['sources']==eva['sources']==freeze['sources']
    for rel,digest in freeze['sources'].items():
        assert sha(ROOT/rel)==digest,(rel,'changed frozen source')
    for split,archive in [('development',dev),('evaluation',eva)]:
        assert len(archive['cases'])==2*len(config['scenes'])*len(config[f'{split}_seeds'])
        keys=set()
        for row in archive['cases']:
            key=(row['seed'],row['scene'],row['phase_error'])
            assert key not in keys
            keys.add(key)
            assert row['transmissions']==240 and row['score_frames']==16
            assert row['receive_groups']==config['receive_groups']
            assert max(row['source_replay_errors'].values())<.02
            assert set(row['methods'])=={
                'compounded_svd_rank1','jiang_ga_svd_rank1','ozgun_svd_a_rank1',
                'jiang_ac_hosvd_adaptation','ozgun_hosvd_adaptation',
                'angle_core_intersection','receive_core_intersection','calibrated_stap_orthogonal'}
            for method in row['methods'].values():
                if 'failure' not in method:
                    assert 0<=method['blood_iq_error']<5
                    assert 0<=method['background_power_remaining']<5
            assert set(row['oracle'])=={f'{rep}_rank{rank}'
              for rep in ('compound','angle','receive','full') for rank in config['oracle_ranks']}
            for value in row['oracle'].values():
                gain=value['gain_to_tissue_limit']
                assert gain is None or 0<=gain<=1+1e-9
                assert value['blood_error_at_tissue_limit'] is None or gain is not None
        assert keys=={(seed,scene,phase) for seed in config[f'{split}_seeds']
                      for scene in config['scenes'] for phase in (False,True)}
    output={'rf_cases':len(dev['cases'])+len(eva['cases']),
            'development_sha256':sha(rf_dir/'development.json'),
            'evaluation_sha256':sha(rf_dir/'evaluation.json'),
            'freeze_sha256':sha(rf_dir/'freeze.json')}
    if measured:
        doc=json.loads(Path(measured).read_text())
        mc=json.loads((ROOT/'configs/resolvestroke_aperture_pilot_v1.json').read_text())
        assert doc['source_sha256']==mc['source_sha256'] and len(doc['clips'])==2
        assert [v['clip'] for v in doc['clips']]==[v['name'] for v in mc['clips']]
        for rel,digest in doc['sources'].items():assert sha(ROOT/rel)==digest
        for clip in doc['clips']:
            assert clip['number_of_transmissions']==128
            assert clip['receive_group_count']==8 and clip['beam_sum_relative_error']<1e-4
            for method in clip['scores'].values():
                assert set(method)=={'1','2'}
                assert all(0<=method[k]['auc']<=1 for k in ('1','2'))
        output['measured_sha256']=sha(measured)
    return output


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--rf-dir',type=Path,default=ROOT/'runs/multidim_aperture_v1')
    cli.add_argument('--measured',type=Path)
    a=cli.parse_args()
    print(json.dumps(verify(a.rf_dir,a.measured),indent=2))
