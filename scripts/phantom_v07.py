"""Preflight and score preregistered, independently referenced phantom IQ.

No phantom data ship with this repository. Optical/encoder fields are confined
to evaluation and never passed into the data-only SharedContext.
"""
import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from btbench.followup_experiment import numerical_source_hash
from btbench.shared import PublicInput,SharedContext,fit_shared,sample_frame,array_hash


def sha(data):return hashlib.sha256(data).hexdigest()


def check_protocol(root=ROOT):
    root=Path(root);freeze=json.loads((root/'docs/PHANTOM_FREEZE_V07.json').read_text())
    for relative,expected in freeze['files'].items():
        if sha((root/relative).read_bytes())!=expected:raise ValueError('Frozen protocol changed: '+relative)
    config=json.loads((root/'configs/phantom_v07.json').read_text())
    if numerical_source_hash()!=config['numerical_source_sha256']:
        raise ValueError('Current filter code differs from frozen phantom source')
    return config,freeze


def _verified_path(base,item):
    path=(base/item['path']).resolve();root=base.resolve()
    if not path.is_relative_to(root) or not path.is_file():raise ValueError('Missing or escaping data path')
    if sha(path.read_bytes())!=item['sha256']:raise ValueError('Data SHA256 differs: '+str(path))
    return path


def _masks(path,shape,config):
    with np.load(path,allow_pickle=False) as file:
        if set(file.files)!={'registration_roi','tissue_score_roi','flow_roi'}:
            raise ValueError('Three prospectively drawn masks are required')
        masks={key:np.asarray(file[key]) for key in file.files}
    for key,mask in masks.items():
        if mask.shape!=shape or mask.dtype!=bool:raise ValueError('ROI shape/type mismatch: '+key)
    reg,tissue,flow=(masks[k] for k in ('registration_roi','tissue_score_roi','flow_roi'))
    q=config['roi_quality']
    if (reg.sum()<q['minimum_registration_pixels'] or tissue.sum()<q['minimum_tissue_score_pixels']
            or flow.sum()<q['minimum_flow_pixels'] or np.any(reg&tissue) or np.any(reg&flow) or np.any(tissue&flow)):
        raise ValueError('ROI size/overlap failed')
    from scipy.ndimage import binary_dilation
    if np.any(binary_dilation(flow,iterations=q['flow_guard_pixels'])&(reg|tissue)):
        raise ValueError('Flow ROI needs the preregistered tissue guard')
    return masks


def _capture(path,config,prf,motion):
    n=config['n_frames'];cal=config['calibration_frames'];shape=tuple(config['tile_shape']);v=config['n_views']
    keys={'views','calibration_views','timestamps_s','calibration_timestamps_s',
          'reference_displacement','reference_valid','reference_sigma_pixels'}
    with np.load(path,allow_pickle=False) as file:
        if set(file.files)!=keys:raise ValueError('Capture arrays must match the contract exactly')
        data={key:np.asarray(file[key]) for key in keys}
    views=data['views'];calibration=data['calibration_views'];d=data['reference_displacement'];valid=data['reference_valid']
    if (views.shape!=(n,v,*shape) or calibration.shape!=(cal,v,*shape) or d.shape!=(n,2,*shape)
            or valid.shape!=(n,*shape) or valid.dtype!=bool or not np.iscomplexobj(views)
            or not np.iscomplexobj(calibration) or not np.isfinite(views).all()
            or not np.isfinite(calibration).all() or not np.isfinite(d).all()):
        raise ValueError('IQ/reference shape, type or finite-value check failed')
    sigma=float(data['reference_sigma_pixels'])
    if not np.isfinite(sigma) or not 0<=sigma<=config['reference_quality']['maximum_sigma_pixels']:
        raise ValueError('Independent motion-reference uncertainty too large')
    t=np.asarray(data['timestamps_s']);ct=np.asarray(data['calibration_timestamps_s'])
    if t.shape!=(n,) or ct.shape!=(cal,) or not np.isfinite(t).all() or not np.isfinite(ct).all():
        raise ValueError('Clock shape/finite check failed')
    for clock in (t,ct):
        if not np.allclose(np.diff(clock),1/prf,atol=1e-7/prf,rtol=1e-5):
            raise ValueError('Nonuniform same-view sampling')
    if not ct[-1]<t[0]:raise ValueError('Calibration must precede scored IQ')
    if not np.all(valid[config['score_start']:,1:-1,1:-1]):
        raise ValueError('Incomplete independent optical support on interior scoring pixels')
    ref=config['reference_quality'];region=d[config['score_start']:, :,1:-1,1:-1]
    rms=float(np.sqrt(np.mean(region**2)));peak=float(np.max(abs(region)))
    if motion=='static':
        if rms>ref['static_maximum_rms_pixels']:raise ValueError('Static optical control moved')
    elif rms<ref['moving_minimum_rms_pixels'] or peak>ref['maximum_displacement_pixels']:
        raise ValueError('Moving control outside registered motion range')
    if motion=='nonrigid' and np.sqrt(np.mean((region-region.mean(axis=(-2,-1),keepdims=True))**2))<ref['nonrigid_minimum_spatial_std_pixels']:
        raise ValueError('Nonrigid control lacks independently measured deformation')
    return data


def preflight(manifest_path,config):
    manifest_path=Path(manifest_path).resolve();base=manifest_path.parent
    manifest=json.loads(manifest_path.read_text());protocol=ROOT/'configs/phantom_v07.json'
    if manifest.get('protocol_sha256')!=sha(protocol.read_bytes()):raise ValueError('Manifest protocol hash differs')
    card=json.loads(_verified_path(base,manifest['acquisition_card']).read_text())
    try:
        prereg=datetime.fromisoformat(card['preregistered_utc'].replace('Z','+00:00'))
        first=datetime.fromisoformat(card['first_capture_utc'].replace('Z','+00:00'))
    except (KeyError,ValueError,AttributeError) as exc:
        raise ValueError('Acquisition card needs ISO timestamps') from exc
    required_card=('scanner','probe','transmit_waveform','receive_aperture','optical_registration',
                   'stage_program','pump_medium','flowmeter_calibration','geometry')
    if (prereg.tzinfo is None or first.tzinfo is None or prereg>=first
            or len(card.get('beamformer_sha256',''))!=64
            or any(not card.get(key) for key in required_card)):
        raise ValueError('Acquisition card must predate capture and pin the beamformer')
    sessions=manifest.get('sessions',[])
    if len(sessions)!=2 or {s['split'] for s in sessions}!=set(config['splits']):
        raise ValueError('Independent development and held-out sessions are required')
    if len({s['id'] for s in sessions})!=2:raise ValueError('Sessions need unique IDs')
    seen=set();capture_hashes=set();rf_hashes=set();inventory=[];mask_hashes=[]
    for session in sessions:
        g=session['geometry'];v=config['n_views'];shape=tuple(config['tile_shape'])
        if g!=card['geometry']:raise ValueError('Hardware/view geometry changed across sessions')
        if min(g['prf_hz'],g['f0_hz'],g['sound_speed_mps'])<=0 or g['phase_sign'] not in (-1,1):
            raise ValueError('Invalid shared acquisition geometry')
        for key in ('tx_rad','rx_rad','view_offsets_s'):
            if np.shape(g[key])!=(v,) or not np.isfinite(g[key]).all():raise ValueError('Three measured view geometries required')
        masks=_masks(_verified_path(base,session['masks']),shape,config)
        mask_hashes.append(session['masks']['sha256'])
        for item in session['triplets']:
            key=(session['split'],item['motion'],item['speed'],item['replicate'])
            if (key in seen or item['motion'] not in config['motions'] or item['speed'] not in config['flow_speeds_mps']
                    or item['replicate'] not in range(1,config['repeats_per_cell']+1)):
                raise ValueError('Duplicate or undeclared triplet')
            seen.add(key)
            if set(item['files'])!=set(config['matched_captures']):raise ValueError('Missing control capture')
            for kind in config['matched_captures']:
                source=item['files'][kind]
                if source['reference_source'] not in ('optical_DIC','optical_fiducials','stage_encoder'):
                    raise ValueError('Ultrasound-derived motion cannot be reference truth')
                if item['motion']=='nonrigid' and source['reference_source']!='optical_DIC':
                    raise ValueError('Nonrigid reference requires independent dense optical DIC')
                if not source.get('calibration_pump_off') or not source.get('calibration_motion_rest'):
                    raise ValueError('Calibration state was not documented as pump-off/rest')
                if source['sha256'] in capture_hashes:raise ValueError('Capture reused between independent repeats')
                capture_hashes.add(source['sha256'])
                for name in ('rf_sha256','beamformer_sha256'):
                    if len(source[name])!=64 or any(c not in '0123456789abcdef' for c in source[name]):
                        raise ValueError('Missing RF/beamformer provenance hash')
                if source['beamformer_sha256']!=card['beamformer_sha256']:
                    raise ValueError('Unregistered beamformer change')
                if source['rf_sha256'] in rf_hashes:raise ValueError('Raw acquisition reused between controls')
                rf_hashes.add(source['rf_sha256'])
                capture=_capture(_verified_path(base,source),config,g['prf_hz'],item['motion'])
                score_mask=masks['registration_roi']|masks['tissue_score_roi']|masks['flow_roi']
                if not np.all(capture['reference_valid'][config['score_start']:,score_mask]):
                    raise ValueError('Optical reference missing on a scored ROI')
                if kind=='flow_on':
                    speed=float(source['flowmeter_mps']);target=config['flow_speeds_mps'][item['speed']]
                    if not np.isfinite(speed) or not .8*target<=speed<=1.2*target:
                        raise ValueError('Measured flow differs from preregistered target')
                inventory.append(dict(split=session['split'],condition=item['motion'],speed=item['speed'],
                                      replicate=item['replicate'],capture=kind,sha256=source['sha256'],
                                      reference_sha256=array_hash(capture['reference_displacement'])))
        if len(session['triplets'])!=len(config['motions'])*len(config['flow_speeds_mps'])*config['repeats_per_cell']:
            raise ValueError('Incomplete factorial session')
    expected={(split,motion,speed,rep) for split in config['splits'] for motion in config['motions']
              for speed in config['flow_speeds_mps'] for rep in range(1,config['repeats_per_cell']+1)}
    if seen!=expected:raise ValueError('Missing or extra control cell')
    if len(set(mask_hashes))!=1:raise ValueError('ROI masks changed between sessions')
    return manifest,inventory


def reference_frame(views,field):
    grid=np.indices(views.shape[-2:],dtype=float)
    return np.stack([sample_frame(frame,grid+d,zero_extension=True) for frame,d in zip(views,field)])


def band_power(views,mask,band):
    if not np.any(mask):raise ValueError('Empty scoring ROI')
    n=len(views);window=np.hanning(n);f=np.abs(np.fft.fftfreq(n));indices=(f>=band[0])&(f<=band[1])
    if not indices.any():raise ValueError('No Doppler bins in preregistered band')
    x=(views[:,:,mask]-np.mean(views[:,:,mask],axis=0))*window[:,None,None]
    spectrum=np.fft.fft(x,axis=0)
    return float(np.mean(np.sum(abs(spectrum[indices])**2,axis=0)/sum(window**2)))


def _power(x,mask):return float(np.mean(abs(x[:,:,mask])**2))


def evaluate_triplet(item,session,base,config,masks):
    captures={kind:_capture(_verified_path(base,item['files'][kind]),config,session['geometry']['prf_hz'],item['motion'])
              for kind in config['matched_captures']}
    s=config['score_start'];packet=config['packet_frames'];band=config['metrics']['doppler_band_cycles_per_frame'][item['speed']]
    if np.count_nonzero((np.abs(np.fft.fftfreq(config['n_frames']-s))>=band[0])&
                        (np.abs(np.fft.fftfreq(config['n_frames']-s))<=band[1]))<config['metrics']['minimum_band_bins']:
        raise ValueError('Insufficient Doppler-band frequency resolution')
    input_aligned={k:reference_frame(c['views'][s:],c['reference_displacement'][s:]) for k,c in captures.items()}
    input_flow={k:band_power(input_aligned[k],masks['flow_roi'],band) for k in ('flow_on','flow_off')}
    excess=input_flow['flow_on']-input_flow['flow_off']
    identifiable=excess>config['metrics']['minimum_input_pump_excess_over_off']*input_flow['flow_off']
    rows=[]
    for mode in ('none','rigid','smooth'):
        specs=[x for x in config['methods'] if x['motion']==mode]
        if not specs:continue
        contexts={}
        for kind,c in captures.items():
            public=PublicInput(c['views'],c['calibration_views'],masks['registration_roi'],session['geometry'])
            contexts[kind]=SharedContext(public,packet,s,mode,motion_options=config['motion_options'])
        for spec in specs:
            outputs={kind:fit_shared(ctx,spec).views for kind,ctx in contexts.items()}
            aligned={k:reference_frame(outputs[k],captures[k]['reference_displacement'][s:]) for k in outputs}
            flow_out={k:band_power(aligned[k],masks['flow_roi'],band) for k in ('flow_on','flow_off')}
            tissue_in=_power(input_aligned['tissue_only'],masks['tissue_score_roi'])
            tissue_out=_power(aligned['tissue_only'],masks['tissue_score_roi'])
            ctx=contexts['tissue_only'];reference=captures['tissue_only']['reference_displacement'][s:]
            geometry=float(np.sqrt(np.mean((ctx.displacement[s:,:,1:-1,1:-1]-reference[:,:,1:-1,1:-1])**2)))
            audit=ctx.motion_diagnostics
            accepted=0 if mode=='none' else sum(row[0]!=0 for row in audit['audit_rows'][s:])
            rows.append(dict(id=spec['id'],motion=mode,method=spec['method'],
                tissue_leakage_ratio=tissue_out/max(tissue_in,1e-30),
                input_pump_excess=excess,input_pump_excess_over_off=excess/max(input_flow['flow_off'],1e-30),
                pump_response_retention=(flow_out['flow_on']-flow_out['flow_off'])/excess if identifiable else None,
                identifiable=bool(identifiable),geometry_rmse_pixels=geometry,
                accepted_fraction=accepted/(config['n_frames']-s),
                output_sha256={k:array_hash(v) for k,v in outputs.items()}))
    return dict(split=session['split'],session=session['id'],motion=item['motion'],speed=item['speed'],
                replicate=item['replicate'],rows=rows)


def _percentile_mean(values,seed,n=2000):
    values=np.asarray(values,float)
    rng=np.random.default_rng(seed);samples=values[rng.integers(len(values),size=(n,len(values)))].mean(axis=1)
    return [float(x) for x in np.percentile(samples,[2.5,97.5])]


def gate(report,config):
    """Fixed candidate only; five paired held-out triplets per stratum."""
    candidate=config['release_gate']['candidate_id'];limits=config['release_gate'];checks=[]
    for motion in config['motions']:
        for speed in config['flow_speeds_mps']:
            cases=[c for c in report['records'] if c['split']=='heldout' and c['motion']==motion and c['speed']==speed]
            if len(cases)!=config['repeats_per_cell']:raise ValueError('Missing held-out replicate')
            rows=[next(r for r in c['rows'] if r['id']==candidate) for c in cases]
            base=[next(r for r in c['rows'] if r['id']=='axis1_none_h64_k1') for c in cases]
            identifiable=all(r['identifiable'] for r in rows)
            retention=None if not identifiable else _percentile_mean([r['pump_response_retention'] for r in rows],43)
            leakage=_percentile_mean([r['tissue_leakage_ratio'] for r in rows],41)
            geometry=_percentile_mean([r['geometry_rmse_pixels'] for r in rows],47)
            static=_percentile_mean([r['accepted_fraction'] for r in rows],53)
            unregistered=float(np.mean([r['geometry_rmse_pixels'] for r in base]))
            motion_ok=(static[1]<=limits['maximum_static_acceptance_fraction'] if motion=='static' else
                       geometry[1]<=limits['maximum_motion_rmse_pixels'] and
                       np.mean([r['geometry_rmse_pixels'] for r in rows])<=
                           (1-limits['minimum_geometry_improvement_fraction'])*unregistered)
            passed=bool(identifiable and retention[0]>=limits['minimum_pump_response_retention']
                        and leakage[1]<=limits['maximum_tissue_leakage_ratio'] and motion_ok)
            checks.append(dict(motion=motion,speed=speed,identifiable=identifiable,
                retention_mean_ci95=retention,tissue_leakage_mean_ci95=leakage,
                displacement_rmse_mean_ci95=geometry,static_acceptance_mean_ci95=static,
                unregistered_displacement_mean=unregistered,passed=passed))
    return dict(passed=all(c['passed'] for c in checks),cells=checks,
                interpretation='Paired pump-associated Doppler excess, not independent blood IQ truth.')


def scaffold(directory,config):
    """Write the full factorial manifest and an intentionally incomplete card."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    card=dict(preregistered_utc=None,first_capture_utc=None,scanner=None,probe=None,
        transmit_waveform=None,receive_aperture=None,optical_registration=None,
        stage_program=None,pump_medium=None,flowmeter_calibration=None,
        beamformer_sha256=None,geometry=dict(prf_hz=None,f0_hz=None,sound_speed_mps=None,
            phase_sign=None,tx_rad=[None]*config['n_views'],rx_rad=[None]*config['n_views'],
            view_offsets_s=[None]*config['n_views']))
    (directory/'acquisition_card.json').write_text(json.dumps(card,indent=2)+'\n')
    sessions=[]
    for split in config['splits']:
        triplets=[]
        for motion in config['motions']:
            for speed in config['flow_speeds_mps']:
                for replicate in range(1,config['repeats_per_cell']+1):
                    prefix=f'{split}_{motion}_{speed}_{replicate:02d}'
                    files={}
                    for kind in config['matched_captures']:
                        files[kind]=dict(path=f'captures/{prefix}_{kind}.npz',sha256=None,
                            rf_sha256=None,beamformer_sha256=None,
                            reference_source='optical_DIC' if motion=='nonrigid' else 'optical_fiducials',
                            calibration_pump_off=True,calibration_motion_rest=True)
                        if kind=='flow_on':files[kind]['flowmeter_mps']=None
                    triplets.append(dict(motion=motion,speed=speed,replicate=replicate,files=files))
        sessions.append(dict(id=split,split=split,geometry=None,
                             masks=dict(path='masks.npz',sha256=None),triplets=triplets))
    manifest=dict(protocol_sha256=sha((ROOT/'configs/phantom_v07.json').read_bytes()),
        acquisition_card=dict(path='acquisition_card.json',sha256=None),sessions=sessions)
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return len(sessions)*len(sessions[0]['triplets'])*len(config['matched_captures'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=('check-protocol','scaffold','preflight','run'))
    parser.add_argument('--manifest');parser.add_argument('--output');args=parser.parse_args()
    config,freeze=check_protocol()
    if args.command=='check-protocol':print(json.dumps(dict(status='frozen',files=freeze['files']),indent=2));return
    if args.command=='scaffold':
        if not args.output:parser.error('--output directory is required for scaffold')
        print(json.dumps(dict(directory=args.output,expected_captures=scaffold(args.output,config))))
        return
    if not args.manifest:parser.error('--manifest is required for data operations')
    manifest,inventory=preflight(args.manifest,config)
    if args.command=='preflight':print(json.dumps(dict(status='eligible',captures=len(inventory),
        sessions=len(manifest['sessions']),protocol_sha256=sha((ROOT/'configs/phantom_v07.json').read_bytes())),indent=2));return
    if not args.output:parser.error('--output is required for run')
    output=Path(args.output)
    if output.exists():raise FileExistsError('Choose a new output path; frozen reports cannot be overwritten')
    base=Path(args.manifest).resolve().parent;records=[]
    for session in manifest['sessions']:
        masks=_masks(_verified_path(base,session['masks']),tuple(config['tile_shape']),config)
        for triplet in session['triplets']:
            records.append(evaluate_triplet(triplet,session,base,config,masks))
    report=dict(study=config['study'],status='phantom_observations',
        numerical_source_sha256=numerical_source_hash(),freeze=freeze,
        manifest_sha256=sha(Path(args.manifest).read_bytes()),inventory=inventory,records=records)
    report['gate']=gate(report,config)
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output),gate_passed=report['gate']['passed'],triplets=len(records))))


if __name__=='__main__':main()
