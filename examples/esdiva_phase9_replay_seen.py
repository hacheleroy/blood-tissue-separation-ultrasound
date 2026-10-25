"""Development-only IQ replay of seed 2729's previously generated RF.

This reconstructs the three 2729 scenes that the earlier acquisition stopped
before saving IQ. It is not a new holdout or an alternate source gate.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, digest, view_rf
from esdiva_phase5_rebeam import load_rf


def rf_samples(folder):
    paths = sorted(folder.glob('rf_event_*.npz'))
    if len(paths) != 288:
        raise ValueError('Incomplete source RF checkpoints')
    length = 0
    for path in paths:
        with np.load(path) as packet:
            length = max(length, packet['rf'].shape[1])
    return length


def replay(acquisitions, upstream, output, seed=2729):
    acquisitions, upstream, output = map(lambda p: Path(p).resolve(),
                                         (acquisitions, upstream, output))
    if seed != 2729 or output.exists():
        raise ValueError('Only unseen output for the now-seen 2729 cloud is allowed')
    head = subprocess.check_output(
        ['git', '-C', str(upstream), 'rev-parse', 'HEAD'], text=True).strip()
    if head != UPSTREAM_COMMIT:
        raise ValueError('Changed Example 22 source')
    sys.path.insert(0, str(upstream/'examples/example22_flow_doppler'))
    import step1_define_flow_phantom as example
    from doppler_tools import echo_center_frequency, rf2iq
    from esdiva.utilities import add_noise, rf_rms
    _, tx, _ = example.build_simulator()
    _, nine = example.compound_events(tx, 1)
    events = [nine[i] for i in (4,3,5)]
    output.mkdir(parents=True)
    for scene in ('stationary_fast','rigid_fast','smooth_slow','smooth_zero'):
        folder = acquisitions/str(seed)/scene
        length = max(rf_samples(folder/'rf_tissue'), rf_samples(folder/'rf_blood'))
        tissue, ct = load_rf(folder/'rf_tissue', 288, 64, length)
        blood, cb = load_rf(folder/'rf_blood', 288, 64, length)
        if (not np.array_equal(ct['t0_per_event'], cb['t0_per_event']) or
                ct['dt'] != cb['dt']):
            raise ValueError('Source RF event coordinates differ')
        clean = tissue+blood
        mixed = add_noise(clean, 6., reference=rf_rms(clean),
                          rng=np.random.default_rng(seed+101))
        noise = mixed-clean
        cal, z, x = view_rf(
            tissue[:48]+noise[:48],
            {'dt':ct['dt'], 't0_per_event':ct['t0_per_event'][:48]},
            events, tx, GRID_MM, example.C)
        coordinates = {'dt':ct['dt'], 't0_per_event':ct['t0_per_event'][48:]}
        beam = {key:view_rf(value[48:], coordinates, events, tx,
                            GRID_MM, example.C)[0]
                for key,value in (('views',mixed),('tissue',tissue),('blood',blood))}
        fs = example.C/(2*GRID_MM['dz']*1e-3)
        carrier = echo_center_frequency(np.moveaxis(cal,-2,-1), 1/fs,
                                        band=(.3*example.FC,1.7*example.FC))
        iq = {key:rf2iq(value, fs=fs, fc=carrier, axis=-2)
              for key,value in beam.items()}
        iq['calibration_views'] = rf2iq(cal, fs=fs, fc=carrier, axis=-2)
        iq['noise'] = iq['views']-iq['tissue']-iq['blood']
        cs = int(np.argmin(abs(z-15.6)))
        ws = int(np.argmin(abs(z-14.4)))
        core = {key:value[:,:,cs:cs+8].copy() for key,value in iq.items()}
        wide = {key:value[:,:,ws:ws+32].copy() for key,value in iq.items()}
        core.update(z_mm=z[cs:cs+8], x_mm=x)
        wide.update(z_mm=z[ws:ws+32], x_mm=x, core_rows=np.arange(12,20))
        if (not np.array_equal(wide['z_mm'][wide['core_rows']], core['z_mm']) or
                not all(np.array_equal(wide[key][:,:,wide['core_rows']], core[key])
                        for key in iq)):
            raise AssertionError('Wider IQ did not preserve the exact core')
        if (folder/'iq_sources.npz').exists():
            with np.load(folder/'iq_sources.npz') as saved:
                if not all(np.array_equal(core[key], saved[key]) for key in iq):
                    raise AssertionError('The completed zero-flow core replay changed')
            with np.load(folder/'broad_iq.npz') as saved:
                if not all(np.array_equal(wide[key], saved[key]) for key in iq):
                    raise AssertionError('The completed zero-flow wider replay changed')
        dest = output/scene
        dest.mkdir()
        np.savez_compressed(dest/'iq_sources.npz', **core)
        np.savez_compressed(dest/'broad_iq.npz', **wide)
        audit = {'seed':seed,'scene':scene,'status':'seen RF development replay only',
                 'upstream_commit':head,'carrier_hz':float(carrier),
                 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 'rf_tissue_sha256':digest(tissue),'rf_blood_sha256':digest(blood),
                 'core_sha256':{key:digest(value) for key,value in core.items()},
                 'wide_sha256':{key:digest(value) for key,value in wide.items()},
                 'original_acquisition_report_exists':(folder/'report.json').exists(),
                 'source_truth_is_evaluator_only':True}
        (dest/'replay.json').write_text(json.dumps(audit,indent=2)+'\n')
        print(scene, 'RF source and IQ replayed', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acquisitions',type=Path,required=True)
    parser.add_argument('--upstream',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    replay(args.acquisitions,args.upstream,args.output)
