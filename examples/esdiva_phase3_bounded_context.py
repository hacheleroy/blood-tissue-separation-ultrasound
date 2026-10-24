"""Use phase-and-RMS view calibration in the existing shared filter operators.

Only the first eight pre-evaluation tissue frames set the common modulation
and tissue reference. The fit has the same complex phase as scalar regression,
but retains each view's measured RMS instead of dividing by coherence. This
limits input noise amplification from angular speckle decorrelation; whether
it preserves an added blood signal must be checked by separate source fits.
"""

import hashlib
import json

import numpy as np

from btbench.shared import (MotionGate, SharedContext, array_hash, estimate_displacement,
                            estimate_tissue_velocity, inverse_coordinates)


def bounded_calibrate(public, gain_model="phase_times_rms"):
    """Fit tissue phase and amplitude on eight early IQ frames, without truth."""
    if gain_model not in ("phase_times_rms", "regression8"):
        raise ValueError("Unknown phase-3 gain model")
    g = public.geometry
    cal = public.calibration_views
    if len(cal) != 16:
        raise ValueError("Require 16 disjoint tissue calibration frames")
    train = cal[:8]
    roi = public.training_roi
    v = estimate_tissue_velocity(train, roi, g["tx_rad"], g["rx_rad"], g["prf_hz"],
                                 g["f0_hz"], g["sound_speed_mps"], g["phase_sign"])
    tx,rx = np.asarray(g["tx_rad"]),np.asarray(g["rx_rad"])
    geometry = np.stack([np.cos(tx)+np.cos(rx), np.sin(tx)+np.sin(rx)],axis=1)
    rate = g["phase_sign"]*2*np.pi*g["f0_hz"]/g["sound_speed_mps"]*(geometry@v)
    offsets=np.asarray(g["view_offsets_s"])
    times=(np.arange(-len(cal),-len(cal)+8)[:,None]/g["prf_hz"]+offsets[None])
    aligned=train*np.exp(-1j*times*rate[None])[:,:,None,None]
    ref=aligned[:,0,roi]
    ref_power=np.vdot(ref,ref).real
    if ref_power <= 0:
        raise ValueError("Empty tissue reference power")
    gains=[]
    for j in range(train.shape[1]):
        view=aligned[:,j,roi]
        cross=np.vdot(ref,view)
        if abs(cross)<1e-12 or np.vdot(view,view).real<=0:
            raise ValueError("Unstable angular calibration")
        rms=np.sqrt(np.vdot(view,view).real/ref_power)
        gains.append((rms*np.exp(1j*np.angle(cross))) if gain_model=="phase_times_rms"
                     else cross/ref_power)
    gains=np.asarray(gains)
    evaluation_times=np.arange(len(public.views))[:,None]/g["prf_hz"]+offsets[None]
    modulation=(gains[None]*np.exp(1j*evaluation_times*rate[None]))[:,:,None,None]
    reference=np.mean(abs(aligned/gains[None,:,None,None]),axis=(0,1))
    reference_views=np.mean(aligned/gains[None,:,None,None],axis=0)
    audit={"gain_model":gain_model, "train_frames":8,
           "estimated_velocity_zx_mps":v.tolist(),
           "relative_gains_real":gains.real.tolist(),
           "relative_gains_imag":gains.imag.tolist(),
           "calibration_sha256":array_hash(cal), "training_roi_sha256":array_hash(roi)}
    return modulation,reference,audit,reference_views


class Phase3Context(SharedContext):
    """Cache the existing shared methods with bounded data-only calibration."""

    def __init__(self, public, packet=16, score_start=64, motion="none", motion_options=None,
                 gain_model="phase_times_rms"):
        n=public.views.shape[0]
        if packet<1 or n%packet or score_start%packet or not 0<=score_start<n:
            raise ValueError("Whole packets and a nonempty scoring interval are required")
        self.public=public;self.packet=packet;self.score_start=score_start
        self.modulation,reference,cal,reference_views=bounded_calibrate(public,gain_model)
        self.displacement,md=estimate_displacement(
            public.views/self.modulation,reference,public.training_roi,motion,
            reference_views=reference_views,
            options=MotionGate(**motion_options) if motion_options else None)
        self.motion=motion;self.calibration=cal;self.motion_diagnostics=md
        if "audit_rows" in md:
            md["audit_sha256"]=hashlib.sha256(json.dumps(
                md["audit_rows"],separators=(",",":"),allow_nan=False).encode()).hexdigest()
        self.motion_summary={key:value for key,value in md.items() if key!="audit_rows"}
        self.has_displacement=bool(np.any(self.displacement))
        grid=np.indices(public.views.shape[2:],dtype=float)
        self.alignment_coordinates=grid[None]+self.displacement
        self.restoration_coordinates=(np.stack([inverse_coordinates(d)
                                                for d in self.displacement])
                                      if self.has_displacement else self.alignment_coordinates)
        self.aligned=self.align(public.views);self.x=self.aligned.reshape(n,-1).T
        self.cache={};self.stops=list(range(score_start+packet,n+1,packet))
