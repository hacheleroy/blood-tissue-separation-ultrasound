"""Phase tracking uses disjoint tissue pixels and preserves the identity map."""

from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"examples"))

from btbench.shared import PublicInput, fit_shared, shared_scene  # noqa: E402
from esdiva_phase4_phase_context import PhaseOnlyContext, fit_phase  # noqa: E402


def test_phase_validation_cannot_change_fitted_phase():
    rng=np.random.default_rng(84)
    reference=rng.normal(size=(3,8,16))+1j*rng.normal(size=(3,8,16))
    roi=np.ones((8,16),bool)
    g={"prf_hz":1000.,"f0_hz":5e6,"sound_speed_mps":1540.,"phase_sign":1,
       "tx_rad":[-.1,0.,.1],"rx_rad":[0.,0.,0.],"view_offsets_s":[0.,.0001,.0002]}
    observed=np.broadcast_to(reference[None]*np.exp(.4j),(16,*reference.shape)).copy()
    public=PublicInput(observed, np.broadcast_to(reference[None],(16,*reference.shape)).copy(),roi,g)
    phase,first=fit_phase(public,np.ones((16,3,1,1),complex),reference)
    z,x=np.indices(roi.shape);changed=observed.copy()
    changed[:,:,((z+2*x)%3==0)]*=np.exp(.9j)
    altered=PublicInput(changed,public.calibration_views,roi,g)
    new_phase,second=fit_phase(altered,np.ones((16,3,1,1),complex),reference)
    # Changed validation data can reject a fit, never change its fitted value.
    assert first["fit_mask_sha256"]==second["fit_mask_sha256"]
    assert first["validation_mask_sha256"]==second["validation_mask_sha256"]
    assert_allclose(phase,.4,atol=1e-14)
    assert not new_phase.any()


def test_phase_only_roundtrip_is_exact_for_complex_iq():
    scene=shared_scene(seed=94,n_frames=32,calibration_frames=16,shape=(8,16),motion="translation")
    ctx=PhaseOnlyContext(scene.public,packet=16,score_start=16)
    result=fit_shared(ctx,{"method":"roundtrip","horizon":16})
    assert_allclose(result.views,scene.public.views[16:],atol=1e-13,rtol=1e-13)
    assert_allclose(result.apply_frozen(scene.blood),scene.blood[16:],atol=1e-13,rtol=1e-13)
