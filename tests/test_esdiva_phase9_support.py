"""Phase-coherent eligibility fits only prior training calibration frames."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"examples"))

from esdiva_angle_calibration import angular_diagnostics  # noqa: E402
from esdiva_phase5_local_phase import tissue_mask  # noqa: E402
from esdiva_phase9_support import phase_coherent_support  # noqa: E402


def test_spatial_phase_gate_improves_stationary_validation_and_rejects_drift():
    rng = np.random.default_rng(914)
    # Match the beamformer grid, which starts at 13 mm before cropping.
    full_z = np.arange(13., 19.01, .1)
    first = np.argmin(abs(full_z-14.4))
    z = full_z[first:first+32]
    x = np.arange(-8, 8, dtype=float)
    base = rng.normal(size=(32,16))+1j*rng.normal(size=(32,16))
    zz, xx = np.indices((32,16))
    angle = .34*xx+.27*zz
    cal = np.broadcast_to(base, (16,3,32,16)).copy()
    cal[:,1] *= np.exp(1j*angle)
    cal[:,2] *= np.exp(-1j*angle)
    cal += .025*(rng.normal(size=cal.shape)+1j*rng.normal(size=cal.shape))
    roi = tissue_mask(z,x,5.)
    raw = angular_diagnostics(cal[:8],cal[8:],roi,(4,3,5))[1:]
    assert min(row["coherence_validation"] for row in raw)<.7
    selected,core,gate = phase_coherent_support(cal,z,x,np.arange(12,20))
    assert selected.sum()==192 and core.sum()==48
    assert gate["pass"] and gate["wide"]["minimum_coherence_validation"]>.95
    assert gate["fit_frames"]==8 and gate["validation_frames"]==8
    perturbed = cal.copy()
    perturbed[8:,1:] *= np.exp(1j*rng.uniform(-np.pi,np.pi,
                                              size=(8,2,32,16)))
    _,_,failure = phase_coherent_support(perturbed,z,x,np.arange(12,20))
    assert not failure["pass"]
    assert failure["phase_sha256"]==gate["phase_sha256"]
