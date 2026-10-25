"""Apply the frozen paired gate once to the complete eight-scene reserve."""

import argparse
import json
from pathlib import Path
import shutil
import tempfile

from esdiva_phase5_freeze import verify_freeze, sha
from esdiva_phase5_heldout_decide import decide
from esdiva_phase5_run_heldout import verify_protocol


def review(root, output_root):
    root,output_root=Path(root).resolve(),Path(output_root).resolve()
    config_path=root/'configs/esdiva_phase5_candidate.json'
    config=json.loads(config_path.read_text())
    protocol=verify_protocol(root)
    target=output_root/'decision.json'
    if target.exists():raise FileExistsError(target)
    with tempfile.TemporaryDirectory(prefix='phase5_reserved_review_') as tmp:
        flat=Path(tmp)
        for seed in protocol['reserved_heldout_seeds']:
            for scene in protocol['scenes']:
                folder=output_root/str(seed)/scene
                report_path,wide_path,score_path=(folder/name for name in
                    ('report.json','broad_iq.json','score.json'))
                # An acquisition or input failure is an explicit missing-scene
                # failure in the frozen decision. Preserve any partial RF and
                # continue reviewing the other reserved scenes.
                if not score_path.is_file():
                    if report_path.is_file():
                        verify_freeze(root,root/'configs/esdiva_phase5_freeze.json',
                                      seed,scene,json.loads(report_path.read_text()))
                    continue
                report=json.loads(report_path.read_text())
                score=json.loads(score_path.read_text())
                audit=json.loads(wide_path.read_text())
                verify_freeze(root,root/'configs/esdiva_phase5_freeze.json',seed,scene,report)
                if (score['seed']!=seed or score['scene']!=scene or
                        score['acquisition_sha256']!=sha(report_path) or
                        score['broad_audit_sha256']!=sha(wide_path) or
                        score['config_sha256']!=sha(config_path) or
                        score['script_sha256']!=sha(root/'examples/esdiva_phase5_heldout_score.py') or
                        score['numerical_source_sha256']!=config['numerical_source_sha256'] or
                        audit['source_report_sha256']!=sha(report_path) or
                        not audit['rf_sources_and_noise_verified']):
                    raise ValueError(f'Reserved scene source or score mismatch: {seed}/{scene}')
                shutil.copyfile(score_path,flat/f'phase5_candidate_{scene}_s{seed}.json')
        decision=decide(config_path,flat,target)
    return {'status':decision['status'],'failures':decision['failures'],
            'score_files':len(decision['input_scores'])}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output-root',type=Path,required=True)
    a=p.parse_args();print(json.dumps(review(a.root,a.output_root)))
