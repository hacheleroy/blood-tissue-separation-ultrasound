"""Execute the published phase-5 comparison on one reserved seed and scene.

The only change in the held-out score/decision copies is seed admission and
outcome wording. This runner checks the exact source substitutions, every
file in the published development freeze, and the pinned upstream simulator
before generating RF. It never reads prescribed motion as a filtering input.
"""

import argparse
import json
from pathlib import Path
import subprocess

from esdiva_phase2_acquire import acquire
from esdiva_phase5_freeze import verify_freeze, sha
from esdiva_phase5_rebeam import rebeam
from esdiva_phase5_heldout_score import evaluate


SCORE_OLD='report["seed"] not in [*config["design_seeds"],config["new_development_seed"]]'
SCORE_NEW='report["seed"] not in config["reserved_heldout_seeds"]'
DECISION_SUBSTITUTIONS=(
    ('seed_list=[*config["design_seeds"],config["new_development_seed"]]',
     'seed_list=config["reserved_heldout_seeds"]'),
    ('item["script_sha256"]!=sha(scripts/"esdiva_phase5_score.py")',
     'item["script_sha256"]!=sha(scripts/"esdiva_phase5_heldout_score.py")'),
    ('Eligible for a pre-held-out source/config/results freeze',
     'Frozen reserved-seed gate passed; proceed to phantom validation'),
    ('Reserved held-out seeds remain unopened',
     'Frozen reserved-seed gate failed; preserve all records without retuning'),
)


def check_adapters(root):
    examples=Path(root)/'examples'
    a=(examples/'esdiva_phase5_score.py').read_text()
    b=(examples/'esdiva_phase5_heldout_score.py').read_text()
    if a.count(SCORE_OLD)!=1 or b!=a.replace(SCORE_OLD,SCORE_NEW):
        raise ValueError('Held-out score differs from frozen scorer beyond seed admission')
    a=(examples/'esdiva_phase5_decide.py').read_text()
    b=(examples/'esdiva_phase5_heldout_decide.py').read_text()
    for old,new in DECISION_SUBSTITUTIONS:
        if a.count(old)!=1:raise ValueError('Frozen decision substitution is ambiguous')
        a=a.replace(old,new)
    if a!=b:raise ValueError('Held-out decision changes a frozen criterion')


def verify_protocol(root):
    root=Path(root).resolve()
    protocol=json.loads((root/'configs/esdiva_phase5_heldout_protocol.json').read_text())
    if (protocol['status']!='phase-5 reserved-bank execution frozen before opening' or
            protocol['development_freeze_sha256']!=sha(root/'configs/esdiva_phase5_freeze.json')):
        raise ValueError('The held-out execution protocol is not bound to the development freeze')
    for name,expected in protocol['files_sha256'].items():
        if sha(root/name)!=expected:raise ValueError(f'Held-out protocol file changed: {name}')
    check_adapters(root)
    return protocol


def run(root, upstream_root, output_root, seed, scene, resume=False):
    root,upstream_root,output_root=map(lambda p:Path(p).resolve(),
                                       (root,upstream_root,output_root))
    freeze=root/'configs/esdiva_phase5_freeze.json'
    verify_freeze(root,freeze,seed,scene)
    protocol=verify_protocol(root)
    if seed not in protocol['reserved_heldout_seeds'] or scene not in protocol['scenes']:
        raise ValueError('Acquisition outside the predeclared reserved bank')
    config_path=root/'configs/esdiva_phase5_candidate.json'
    config=json.loads(config_path.read_text())
    pinned=subprocess.check_output(['git','-C',str(upstream_root),'rev-parse','HEAD'],text=True).strip()
    if pinned!=config['upstream_commit']:
        raise ValueError('Upstream simulator is not the pinned commit')
    spec=config['scenes'][scene]
    folder=output_root/str(seed)/scene
    if resume:
        if not (folder/'report.json').is_file() or not (folder/'iq_sources.npz').is_file():
            raise FileNotFoundError('A complete, saved RF acquisition is required for replay')
    else:
        if folder.exists():raise FileExistsError(f'Will not reuse an output folder: {folder}')
        folder.parent.mkdir(parents=True,exist_ok=True)
        acquire(upstream_root,folder,seed=seed,
                scatterer_fraction=config['scatterer_fraction'],
                calibration_frames=config['calibration_frames'],
                evaluation_frames=config['evaluation_frames'],
                motion=spec['motion'],flow_scale=spec['flow_scale'])
    # Acquisition serializes the view tuple as a JSON list. Verify the saved
    # source record that all later replay and scoring stages consume.
    verify_freeze(root,freeze,seed,scene,json.loads((folder/'report.json').read_text()))
    broad=folder/'broad_iq.npz'
    if not broad.exists():rebeam(folder,upstream_root,broad)
    score_path=folder/'score.json'
    if score_path.exists():raise FileExistsError(score_path)
    score=evaluate(folder,broad,config_path,score_path)
    if score['seed']!=seed or score['scene']!=scene:
        raise AssertionError('Reserved score label mismatch')
    return {"seed":seed,"scene":scene,"status":score['status'],
            "report_sha256":sha(folder/'report.json'),
            "score_sha256":sha(folder/'score.json')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--upstream-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--scene',required=True)
    p.add_argument('--resume',action='store_true',help='Score a saved complete acquisition without regenerating RF')
    a=p.parse_args();print(json.dumps(run(a.root,a.upstream_root,a.output_root,a.seed,a.scene,a.resume)))
