"""Delayed receive subapertures from eSDIva source-resolved RF checkpoints.

The delay, F#1 Hann aperture, and linear interpolation mirror eSDIva's
das_volume; only the final channel accumulation is partitioned into groups.
This is acoustic *simulation* channel data, not a measured probe recording.
"""
from pathlib import Path
import numpy as np


def delayed_groups(rf,t0,dt,event,centers,z_mm,x_mm,c=1540.,groups=8):
    """(channel,fast_time) -> (group,z,x); groups sum to ordinary DAS."""
    rf=np.asarray(rf)
    el=np.asarray(centers,float)
    if rf.ndim!=2 or el.shape!=(rf.shape[0],3) or rf.shape[0]%groups:
        raise ValueError("Receiver dimensions/groups mismatch")
    if "angles_deg" not in event:
        raise ValueError("This pinned acquisition uses plane-wave events")
    ang=np.atleast_1d(np.deg2rad(event["angles_deg"]))
    n=np.array([np.sin(ang[0]),np.sin(ang[1]) if len(ang)>1 else 0.,0.])
    n[2]=np.sqrt(1-n[0]**2-n[1]**2)
    delays=np.asarray(event.get("delays",np.zeros(len(el))),float)
    apod=np.asarray(event.get("apodization",np.ones(len(el))),float)>0
    tau=delays-delays.max()
    tref=float(np.mean(tau[apod]-el[apod]@n/c))
    xx,zz=np.meshgrid(np.asarray(x_mm)*1e-3,np.asarray(z_mm)*1e-3,indexing='xy')
    pts=np.stack([xx.ravel(),np.zeros(xx.size),zz.ravel()],axis=1)
    tx=tref+pts@n/c
    delta=pts[None,:,:]-el[:,None,:]
    radial=np.linalg.norm(delta[:,:,:2],axis=-1)
    half_ap=pts[None,:,2]/2. # F-number 1
    valid=radial<=half_ap
    rx=np.linalg.norm(delta,axis=-1)/c
    samples=(tx[None]+rx-t0)/dt
    i0=np.floor(samples).astype(np.int64)
    valid &= (i0>=0)&(i0<rf.shape[1]-1)
    idx=np.clip(i0,0,rf.shape[1]-2)
    alpha=samples-i0
    val=np.take_along_axis(rf,idx,axis=1)*(1-alpha)+np.take_along_axis(rf,idx+1,axis=1)*alpha
    taper=(.5+.5*np.cos(np.pi*np.minimum(radial/half_ap,1.)))
    val*=valid*taper
    # No arbitrary normalization: summing K groups returns the DAS RF.
    return val.reshape(groups,len(el)//groups,*xx.shape).sum(axis=1).astype(np.float32)


def read_rf_checkpoint(directory, padded_length=None):
    files=sorted(Path(directory).glob('rf_event_*.npz'))
    if not files or any(p.name!=f'rf_event_{j:04d}.npz' for j,p in enumerate(files)):
        raise ValueError("Incomplete RF event schedule")
    arrays=[];times=[];dt=None
    for p in files:
        with np.load(p) as item:
            arrays.append(item['rf']);times.append(float(item['t0']))
            if dt is None:dt=float(item['dt'])
            elif dt!=float(item['dt']):raise ValueError("Event dt changed")
    length=max(x.shape[-1] for x in arrays) if padded_length is None else padded_length
    if any(x.shape[-1]>length or x.shape[0]!=arrays[0].shape[0] for x in arrays):
        raise ValueError("Incompatible RF event traces")
    padded=np.zeros((len(arrays),arrays[0].shape[0],length),dtype=arrays[0].dtype)
    for i,x in enumerate(arrays):padded[i,:,:x.shape[-1]]=x
    return padded,np.array(times),dt
