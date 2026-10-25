"""Deterministically build the declared phase-6 trial configs, without RF."""

import argparse
import copy
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generate(root, output):
    root,output=Path(root).resolve(),Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    rules_path=root/'configs/esdiva_phase6_trial_rules.json'
    rules=json.loads(rules_path.read_text())
    split_path=root/'configs/esdiva_phase6_split.json'
    split=json.loads(split_path.read_text())
    if (rules['admitted_development_seeds']!=split['design_development_seeds'] or
            rules['reserved_heldout_seeds']!=split['reserved_heldout_seeds'] or
            rules['prospective_development_seed']!=split['prospective_development_seed']):
        raise ValueError('Trial rules violate the declared seed split')
    base=json.loads((root/'configs/esdiva_phase5_candidate.json').read_text())
    output.mkdir(parents=True)
    names=[]
    for trial in rules['trials']:
        if (trial['gate'] not in rules['core_gate_models'] or
                trial['fraction'] not in rules['phase_fractions'] or
                trial['sigma'] not in rules['phase_sigmas_pixels'] or
                trial['motion'] not in ('local_phase','unit_shift','phase_plus_shift')):
            raise ValueError('Trial option outside bounded rules')
        config=copy.deepcopy(base)
        config.update(status=rules['status'],
            design_seeds=split['design_development_seeds'],
            new_development_seed=split['prospective_development_seed'],
            reserved_heldout_seeds=split['reserved_heldout_seeds'],
            admitted_development_seeds=rules['admitted_development_seeds'],
            core_gate_model=trial['gate'],
            core_validation_absolute_min=rules['core_validation_absolute_min'],
            core_validation_to_train_min=rules['core_validation_to_train_min'],
            motion_variant=trial['motion'],unit_shift=rules['shift_options'],
            phase6_split_sha256=sha(split_path),phase6_trial_rules_sha256=sha(rules_path))
        config['local_phase'].update(phase_fraction=trial['fraction'],
                                     spatial_sigma_pixels=trial['sigma'])
        name=(f"{trial['gate']}_{trial['motion']}_"
              f"f{round(trial['fraction']*100):03d}_s{round(trial['sigma'])}.json")
        if name in names:raise ValueError('Duplicate development trial')
        (output/name).write_text(json.dumps(config,indent=2,sort_keys=True)+'\n')
        names.append(name)
    return names


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(generate(a.root,a.output)))
