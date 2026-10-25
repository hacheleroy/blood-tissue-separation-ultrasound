"""Phase-10 core geometry and train-only calibration eligibility."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'examples'))

from esdiva_phase10_support import phase10_support  # noqa: E402


def test_phase10_core_guard_and_disjoint_calibration_validation():
    rng=np.random.default_rng(910)
    full_z=np.arange(13.,19.01,.1)
    first=np.argmin(abs(full_z-14.4))
    z=full_z[first:first+32]
    x=np.arange(-8,8,dtype=float)
    base=rng.normal(size=(32,16))+1j*rng.normal(size=(32,16))
    zz,xx=np.indices((32,16))
    phase=.34*xx+.27*zz
    cal=np.broadcast_to(base,(16,3,32,16)).copy()
    cal[:,1]*=np.exp(1j*phase)
    cal[:,2]*=np.exp(-1j*phase)
    cal+=.025*(rng.normal(size=cal.shape)+1j*rng.normal(size=cal.shape))
    wide,core,gate=phase10_support(cal,z,x,np.arange(12,20))
    assert wide.sum()==192 and core.sum()==40
    assert gate['pass'] and gate['wide']['minimum_coherence_validation']>.95
    assert gate['fit_frames']==gate['validation_frames']==8
    perturbed=cal.copy()
    perturbed[8:,1:]*=np.exp(1j*rng.uniform(-np.pi,np.pi,
                                             size=(8,2,32,16)))
    _,_,failed=phase10_support(perturbed,z,x,np.arange(12,20))
    assert not failed['pass']
    assert failed['phase_sha256']==gate['phase_sha256']
