"""Independent spatial check and exact complex-IQ phase roundtrip."""

from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"examples"))

from btbench.shared import PublicInput, fit_shared  # noqa: E402
from esdiva_phase5_local_phase import LocalPhaseContext, fit_local_phase, tissue_mask  # noqa: E402


def example_input():
    rng=np.random.default_rng(52)
    z=14.4+.1*np.arange(32);x=np.arange(-8,8,dtype=float);roi=tissue_mask(z,x)
    base=rng.normal(size=(3,32,16))+1j*rng.normal(size=(3,32,16))
    cal=np.broadcast_to(base[None],(16,*base.shape)).copy()
    phase=np.broadcast_to(.32+.02*(z[:,None]-16),base.shape[1:])
    views=np.broadcast_to(base[None]*np.exp(1j*phase)[None,None],(32,*base.shape)).copy()
    g={"prf_hz":1000.,"f0_hz":5e6,"sound_speed_mps":1540.,"phase_sign":-1,
       "tx_rad":[0.,-.087,.087],"rx_rad":[0.,0.,0.],
       "view_offsets_s":[0.,1/3000,2/3000]}
    wide=PublicInput(views,cal,roi,g)
    core=PublicInput(views[:,:,12:20].copy(),cal[:,:,12:20].copy(),
                     np.broadcast_to(abs(x[None])>=6.,(8,16)).copy(),g)
    return wide,core,z,x


def test_phase_fit_uses_disjoint_validation_and_ignores_rejected_views():
    wide,_,z,x=example_input()
    phase,audit=fit_local_phase(wide,np.arange(12,20),z,x)
    assert audit["fit_pixels"]+audit["validation_pixels"]==wide.training_roi.sum()
    assert audit["fit_mask_sha256"]!=audit["validation_mask_sha256"]
    assert audit["accepted_view_frames"]==96
    assert np.mean(abs(phase))>.05
    zz,xx=np.indices(wide.training_roi.shape)
    checking=wide.training_roi&((zz+2*xx)%3==0)
    altered=wide.views.copy()
    altered[:,:,checking]*=np.exp(-1.1j)
    rejected=PublicInput(altered,wide.calibration_views,wide.training_roi,wide.geometry)
    _,bad=fit_local_phase(rejected,np.arange(12,20),z,x)
    assert bad["accepted_view_frames"]<audit["accepted_view_frames"]


def test_local_phase_context_roundtrips_complex_blood():
    wide,core,z,x=example_input()
    context=LocalPhaseContext(core,wide,np.arange(12,20),z,x,score_start=16)
    result=fit_shared(context,{"method":"roundtrip","horizon":16})
    blood=np.random.default_rng(3).normal(size=core.views.shape)+1j*np.random.default_rng(4).normal(size=core.views.shape)
    assert_allclose(result.views,core.views[16:],rtol=1e-13,atol=1e-13)
    assert_allclose(result.apply_frozen(blood),blood[16:],rtol=1e-13,atol=1e-13)
