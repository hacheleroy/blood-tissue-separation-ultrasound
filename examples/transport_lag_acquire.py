"""Matched-emission, source-resolved eSDIva acquisition for transport-lag study.

Two schedules fire the same 240 three-angle emissions from the same initial
scatterer cloud. Dense has 80 contiguous 3.7 kHz cycles; burst has ten blocks
of eight such cycles, with blocks starting 22 ms apart. Blood and tissue
positions are evaluated at actual emission timestamps, not frame indices.
This is a linear, homogeneous, point-scatterer acoustic simulation: no wall
reflection, reverberation, red-cell spectral signature, or real electronics.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from esdiva_example22_pilot import GRID_MM, UPSTREAM_COMMIT, digest, view_rf
from esdiva_phase10_acquire import VIEW_INDICES

ROOT=Path(__file__).resolve().parents[1]
UPSTREAM=ROOT.parent/'esdiva_upstream'
F_CARRIER=2555519.761864708  # Earlier phase-11 carrier from the same probe.
N_CYCLES=80
N_ANGLES=3
PRF=3700.


def timestamps(schedule):
    if schedule=='dense':
        cycles=np.arange(N_CYCLES)/PRF
    elif schedule=='burst':
        cycles=np.repeat(np.arange(10)*.022,8)+np.tile(np.arange(8)/PRF,10)
    else:raise ValueError(schedule)
    events=(cycles[:,None]+np.arange(N_ANGLES)[None,:]/(N_ANGLES*PRF)).ravel()
    assert len(events)==240 and np.all(np.diff(events)>0)
    return cycles,events


def trajectories(example,seed,motion,events,fraction=.1,flow_scale=.02):
    """Same cloud for both schedules; displacement uses physical seconds.

    Tissue follows 3 Hz lateral and 2 Hz axial motion. The smooth class has
    mild position-dependent deformation. Blood has the original Poiseuille
    speeds at 2% scale, including slower near-wall scatterers.
    """
    old=example.SEED
    try:
        example.SEED=seed
        positions,amps,truth=example.build_phantom(len(events),PRF*N_ANGLES)
    finally:example.SEED=old
    n0=truth['n_tissue']
    rng=np.random.default_rng(seed+1)
    t_ids=np.sort(rng.choice(n0,max(1,round(n0*fraction)),replace=False))
    b_ids=np.sort(rng.choice(truth['n_blood'],max(1,round(truth['n_blood']*fraction)),replace=False))
    initial=positions[0,np.r_[t_ids,n0+b_ids]].copy()
    amps=amps[np.r_[t_ids,n0+b_ids]]
    n_t=len(t_ids)
    pos=np.broadcast_to(initial[None],(len(events),len(initial),3)).copy()
    base=initial[:n_t]
    if motion!='none':
        lateral=.30*np.sin(2*np.pi*3*events)
        axial=.06*np.sin(2*np.pi*2*events)
        if motion=='smooth':
            lateral=lateral[:,None]*(1+.20*(base[None,:,2]-16.)/3.)
            axial=axial[:,None]*(1+.25*base[None,:,0]/8.)
        elif motion=='rigid':
            lateral=lateral[:,None]
            axial=axial[:,None]
        else:raise ValueError(motion)
        pos[:,:n_t,0]+=lateral
        pos[:,:n_t,2]+=axial
    relative=initial[n_t:]-example.VESSEL_CENTER_MM
    s0=relative@example.FLOW_DIR
    radial=relative-s0[:,None]*example.FLOW_DIR
    speed=flow_scale*truth['v_blood'][b_ids]
    along=(s0[None]+1e3*events[:,None]*speed[None]+example.VESSEL_LENGTH_MM/2)
    along=along%example.VESSEL_LENGTH_MM-example.VESSEL_LENGTH_MM/2
    pos[:,n_t:]=(example.VESSEL_CENTER_MM[None,None]+along[:,:,None]*example.FLOW_DIR[None,None]+
                 radial[None]).astype(np.float32)
    assert np.array_equal(pos[0],initial)
    return pos,amps,n_t


def acquire(output,seed,motion,schedule,replay_noise=False):
    import subprocess
    from doppler_tools import rf2iq
    from esdiva.utilities import add_noise,rf_rms
    path=UPSTREAM/'examples/example22_flow_doppler'
    head=subprocess.check_output(['git','-C',str(UPSTREAM),'rev-parse','HEAD'],text=True).strip()
    if head!=UPSTREAM_COMMIT:raise RuntimeError('Changed acoustic upstream')
    sys.path.insert(0,str(path))
    import step1_define_flow_phantom as example
    output=Path(output)
    previous=None
    if (output/'report.json').exists():
        if not replay_noise:raise FileExistsError(output)
        previous=json.loads((output/'report.json').read_text())
        if (previous['seed'],previous['motion'],previous['schedule'])!=(seed,motion,schedule):
            raise ValueError('Replay requested on different acoustic scene')
    output.mkdir(parents=True,exist_ok=True)
    cycles,events_t=timestamps(schedule)
    sim,probe,_=example.build_simulator()
    _,nine=example.compound_events(probe,1)
    angle_events=[nine[i] for i in VIEW_INDICES]
    positions,amps,n_t=trajectories(example,seed,motion,events_t)
    print(seed,motion,schedule,'RF sources',len(events_t),'events',flush=True)
    t_amp=amps.copy();t_amp[n_t:]=0
    b_amp=amps.copy();b_amp[:n_t]=0
    t,ct=sim.sequence_rf(positions,t_amp,angle_events*N_CYCLES,out_path=output/'rf_tissue')
    b,cb=sim.sequence_rf(positions,b_amp,angle_events*N_CYCLES,out_path=output/'rf_blood')
    if previous and (digest(t)!=previous['rf_tissue_sha256'] or digest(b)!=previous['rf_blood_sha256']):
        raise ValueError('Acoustic RF identity changed during noise replay')
    assert t.shape==b.shape and ct['dt']==cb['dt'] and np.array_equal(ct['t0_per_event'],cb['t0_per_event'])
    # One combined event checks source linearity against independently simulated RF.
    k=17
    combined,cc=sim.sequence_rf(positions[k:k+1],amps,[angle_events[k%3]])
    assert np.array_equal(cc['t0_per_event'],ct['t0_per_event'][k:k+1])
    closure=float(np.linalg.norm((t[k:k+1,:,:combined.shape[-1]]+b[k:k+1,:,:combined.shape[-1]]-combined).ravel())/
                  max(np.linalg.norm(combined.ravel()),1e-30))
    assert closure<2e-4,closure
    clean=t+b
    # The first eight cycles are identical between schedules. Use their common
    # fixed 900-sample RF support for ONE absolute receiver-noise scale;
    # full-record RMS would give the schedules different amplifiers.
    noise_reference=rf_rms(clean[:24,:,:900])
    mix=add_noise(clean,6.,reference=noise_reference,rng=np.random.default_rng(seed+101))
    n=mix-clean
    coord={'dt':ct['dt'],'t0_per_event':ct['t0_per_event']}
    def convert(rf):
        wave,z,x=view_rf(rf,coord,angle_events,probe,GRID_MM,example.C)
        iq=rf2iq(wave,fs=example.C/(2*GRID_MM['dz']*1e-3),fc=F_CARRIER,axis=-2)
        assert iq.shape[:2]==(N_CYCLES,3)
        zidx=np.flatnonzero((z>=14.4-1e-8)&(z<=17.5+1e-8))
        assert len(zidx)==32 and len(x)==16
        return iq[:,:,zidx],z[zidx],x
    ti,z,x=convert(t);bi,_,_=convert(b);ni,_,_=convert(n)
    mi,_,_=convert(mix)
    iq_closure=float(np.linalg.norm((mi-ti-bi-ni).ravel())/max(np.linalg.norm(mi.ravel()),1e-30))
    assert iq_closure<2e-4,iq_closure
    core=np.arange(12,20)
    np.savez_compressed(output/'iq_sources.npz',views=mi,tissue=ti,blood=bi,noise=ni,
                        cycle_times_s=cycles,event_times_s=events_t,z_mm=z,x_mm=x,core_rows=core)
    report={'status':'source-resolved matched-emission acoustic simulation',
            'upstream_commit':head,'seed':seed,'motion':motion,'schedule':schedule,
            'flow_scale':.02,'tissue_motion':'3 Hz x; 2 Hz z; 0.30/0.06 mm amplitude; smooth has weak spatial deformation',
            'transmissions':len(events_t),'cycles':len(cycles),'angles':[example.ANGLES_DEG[i] for i in VIEW_INDICES],
            'duration_ms':float(1e3*(events_t[-1]-events_t[0])),
            'n_tissue':n_t,'n_blood':len(amps)-n_t,
            'cycle_times_sha256':digest(cycles),'event_times_sha256':digest(events_t),
            'rf_shape':list(t.shape),'rf_tissue_sha256':digest(t),'rf_blood_sha256':digest(b),
            'rf_noise_sha256':digest(n),'rf_t0_sha256':digest(ct['t0_per_event']),
            'noise_reference_rf_rms':float(noise_reference),
            'noise_reference_rule':'first 24 identical events, first 900 RF samples; common across schedules',
            'rf_additivity_relative_error':closure,'iq_additivity_relative_error':iq_closure,
            'grid_mm':GRID_MM,'carrier_hz':F_CARRIER,
            'iq_sha256':{label:digest(v) for label,v in [('views',mi),('tissue',ti),('blood',bi),('noise',ni)]},
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'noise_replayed_from_identical_rf':bool(previous),
            'limitations':'Weak point scatterers in homogeneous linear acoustics; electronic noise is synthesized; no vessel wall reflections.'}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--motion',choices=('rigid','smooth'),required=True)
    p.add_argument('--schedule',choices=('dense','burst'),required=True)
    p.add_argument('--replay-noise-from-checkpoints',action='store_true')
    args=p.parse_args()
    print(json.dumps(acquire(args.out,args.seed,args.motion,args.schedule,
                             args.replay_noise_from_checkpoints),indent=2),flush=True)
