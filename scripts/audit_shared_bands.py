"""Replay the frozen adaptive SSM bands and check actual output against archive.

This adds an inspectable band ledger without changing the frozen numerical code
or tuning on an evaluation result. It uses the original seeds and parameters.
"""
import csv
import gzip
import json
import lzma
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from btbench.shared import shared_scene, SharedContext, fit_shared
from btbench.shared_experiment import validate_shared_config
from btbench.followup_experiment import numerical_source_hash
from btbench.svd_literature import baranger_band


def main(directory):
    directory=Path(directory);config=validate_shared_config(json.loads((directory/'frozen_config.json').read_text()))
    if (directory/'results.json.gz').exists():raw=gzip.decompress((directory/'results.json.gz').read_bytes())
    else:raw=lzma.decompress(b''.join(p.read_bytes() for p in sorted(directory.glob('results.json.xz.part*'))))
    report=json.loads(raw)
    if numerical_source_hash()!=report['code_sha256']:raise ValueError('Replay needs the frozen numerical source')
    conditions={c['id']:c for c in config['conditions']};rows=[];differences=[]
    with threadpool_limits(limits=1):
        for case in report['records']:
            condition=conditions[case['condition']]
            scene=shared_scene(**{**config['generator'],**{k:v for k,v in condition.items() if k!='id'},'seed':case['seed']})
            for motion in ('none','smooth'):
                mid=f'baranger_{motion}_h64';spec=next(m for m in config['methods'] if m['id']==mid)
                ctx=SharedContext(scene.public,config['packet_frames'],config['score_start'],motion)
                replay=fit_shared(ctx,spec);truth=scene.blood[config['score_start']:]
                error=float(np.linalg.norm(replay.views-truth)/np.linalg.norm(truth))
                original=next(r for r in case['results'] if r['id']==mid)['metrics']['blood_iq_nrmse']
                delta=abs(error-original);differences.append(delta)
                if delta>1e-9:raise ValueError('Replayed output differs from archived numerical result')
                for stop in ctx.stops:
                    u,_,_=ctx.decomposition(stop,64)
                    k,h,d=baranger_band(u)
                    rows.append(dict(condition=case['condition'],seed=case['seed'],split=case['split'],motion=motion,
                        start=stop-config['packet_frames'],stop=stop,lower_cutoff=k,upper_cutoff_exclusive=h,
                        retained_modes=h-k,available_modes=u.shape[1],actual_iq_nrmse=error,
                        replay_error_difference=delta,code_sha256=report['code_sha256']))
    with (directory/'baranger_bands.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(dict(replayed_actual_outputs=len(differences),bands=len(rows),max_error_difference=max(differences))))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('directory');main(p.parse_args().directory)
