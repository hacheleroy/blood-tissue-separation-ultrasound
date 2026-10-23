"""Split and matched-support safeguards for the predeclared v06 study."""
import json
from pathlib import Path

from btbench.shared_experiment import default_motion_repair_config,validate_shared_config


def test_frozen_motion_config_matches_declared_study():
    root=Path(__file__).resolve().parents[1]
    config=json.loads((root/'configs/motion_repair_v06.json').read_text())
    assert config==validate_shared_config(default_motion_repair_config())
    assert config['evaluation_seeds']==[1601,1907,2203,2503,2801]
    assert len(config['conditions'])==8 and len(config['methods'])==36
    assert {m['motion'] for m in config['methods']}=={'none','rigid','smooth','known'}
    assert all(m['horizon']==64 for m in config['methods'] if m['method']!='roundtrip')
    assert all(m['evaluation_only']==(m['motion']=='known') for m in config['methods'])
