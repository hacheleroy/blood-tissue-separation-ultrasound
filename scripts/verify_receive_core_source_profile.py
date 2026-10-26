"""Verify source identities and arithmetic of the archived receive-core diagnostic.

The source-separated RF checkpoints are workspace-only. --require-rf checks
their report hashes when they are available; CI verifies the frozen numerical
archive without claiming to rerun the RF beamformer from a clean checkout.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/receive_core_source_profile_v1'
CONFIG=ROOT/'configs/receive_core_source_profile_v1.json'
SOURCE=ROOT/'scripts/receive_core_source_profile.py'
REFERENCE=ROOT/'runs/flow_background_audit_v1/multidim_exact.json'
FIGURE=ROOT/'docs/figures/receive_core_source_profile_v1.png'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def close(a,b):return math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-12)


def check(require_rf=False):
    manifest=json.loads((RUN/'manifest.json').read_text())
    result=json.loads((RUN/'results.json').read_text())
    config=json.loads(CONFIG.read_text())
    assert sha(RUN/'results.json')==manifest['results_sha256']
    assert sha(FIGURE)==manifest['figure_sha256']
    assert result['config_sha256']==sha(CONFIG)
    assert result['script_sha256']==sha(SOURCE)
    assert result['reference_sha256']==sha(REFERENCE)
    expected={(seed,scene) for seed in config['seeds'] for scene in config['scenes']}
    assert len(result['cases'])==len(expected)==4
    assert {(c['seed'],c['scene']) for c in result['cases']}==expected
    reference=json.loads(REFERENCE.read_text())
    ref={(r['seed'],r['scene'],r['phase_error'],r['method']):r for r in reference['rows']}
    matched_rf=0;unavailable_rf=0;checked_regions=0
    for split,seed in [('development',3907),('evaluation',3919)]:
        path=ROOT/'runs/multidim_aperture_v1'/f'{split}.json'
        assert sha(path)==result['multidim_source_archive_sha256'][split]
        archived=json.loads(path.read_text())
        for c in [x for x in result['cases'] if x['seed']==seed]:
            original=next(x for x in archived['cases'] if x['seed']==seed and
                          x['scene']==c['scene'] and x['phase_error'] is False)
            assert c['source_sha256']['components']==original['component_sha256']
            report=ROOT.parent/'esdiva_phase11_work/acquisitions'/str(seed)/c['scene']/'report.json'
            if report.exists():
                assert sha(report)==c['source_sha256']['report']
                matched_rf+=1
            else:unavailable_rf+=1
            assert c['nominal_axial_peak_mps']>0 and not c['phase_error']
            assert list(c['radial_bin_counts'][k] for _,_,k in config['radial_bins'])==[32,16,16,16]
            assert [m['method'] for m in c['methods']]==config['methods']
            assert c['baranger_boundaries_zero_based']['lower']<c['baranger_boundaries_zero_based']['upper']<=80
            for energies in c['tensor_mode_fraction'].values():
                assert all(0<=v<=1 for v in energies.values())
                assert close(sum(energies.values()),1.)
            svd=c['compounded_svd_modes']
            for source,energies in svd['per_mode'].items():
                assert len(energies)==80 and min(energies)>=0
                assert close(sum(b['power'][source] for b in svd['bands']),sum(energies))
            for method in c['methods']:
                assert list(map(int,method['scores']))==config['score_windows']
                for window,regions in method['scores'].items():
                    assert set(regions)==set(c['radial_bin_counts'])
                    raw=next(m for m in c['methods'] if m['method']=='unfiltered')['scores'][window]
                    for region,score in regions.items():
                        assert score['pixels']==c['radial_bin_counts'][region]
                        assert close(score['background'],score['tissue']+score['noise']+
                                     score['tissue_noise_cross'])
                        assert close(score['blood_power_retention'],
                                     score['blood']/raw[region]['blood'])
                        assert close(score['background_remaining_over_input_tissue'],
                                     score['background']/raw[region]['tissue'])
                        assert close(score['blood_background_db'],
                                     10*math.log10(score['blood']/score['background']))
                        assert -1e-10<=score['velocity_coherence_median']<=1+1e-10
                        checked_regions+=1
                if (seed,c['scene'],False,method['method']) in ref:
                    archive=ref[(seed,c['scene'],False,method['method'])]
                    score=method['scores']['16']['whole_tile']
                    assert close(score['blood_power_retention'],archive['blood_power_retention'])
                    assert close(score['background_remaining_over_input_tissue'],
                                 archive['background_power_remaining'])
                    assert close(score['blood_iq_error'],archive['blood_iq_error'])
    if require_rf and unavailable_rf:raise FileNotFoundError('RF source report missing')
    assert checked_regions==4*len(config['methods'])*len(config['score_windows'])*6
    return {'cases':len(result['cases']),'region_scores_checked':checked_regions,
            'rf_reports_verified':matched_rf,'rf_reports_unavailable':unavailable_rf}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-rf',action='store_true')
    args=parser.parse_args()
    print(json.dumps(check(require_rf=args.require_rf),indent=2))
