"""Independent objective, analytic transfer and known-subspace checks."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import linalg,signal
from btbench.core import casorati,frames
from btbench.filters import svd_filter,temporal_highpass
from btbench.benchmark import run_filter
from btbench.svd_literature import (square_fit,spatial_similarity,gram_decompose,
    randomized_svd_filter,huang_evd_filter,huang_weights,gaussian_local_evd,song_local_filter,
    temporal_statistics)


def random(shape,seed=0):
    rng=np.random.default_rng(seed)
    return rng.normal(size=shape)+1j*rng.normal(size=shape)


@pytest.mark.parametrize("r",[5,9,13])
def test_baranger_sum_table_matches_explicit_corr2(r):
    rng=np.random.default_rng(r)
    c=rng.normal(size=(r,r));c=(c+c.T)/2
    scores=[]
    for k in range(1,r):
        for h in range(k+1,r+1):
            template=np.zeros((r,r));template[:k,:k]=1;template[k:h,k:h]=1
            scores.append((np.corrcoef(c.ravel(),template.ravel())[0,1],k,h))
    best=max(scores,key=lambda v:v[0])
    k,h,d=square_fit(c)
    assert (k,h)==best[1:]
    assert_allclose(d["ssm_template_correlation"],best[0],atol=1e-13)


def test_baranger_known_boundaries_beyond_64_modes():
    c=np.zeros((90,90));c[:4,:4]=1;c[4:75,4:75]=1
    k,h,d=square_fit(c)
    assert (k,h)==(4,75)
    assert d["ssm_modes"]==90 and not d["ssm_mode_truncation"]


def test_baranger_undefined_envelope_and_flat_fallback():
    u=np.column_stack((np.ones(7),np.arange(7),np.arange(7)[::-1]))
    c,n=spatial_similarity(u)
    assert n==1 and np.isfinite(c).all()
    assert_allclose(c[1:,1:],[[1,-1],[-1,1]],atol=1e-15)
    k,h,d=square_fit(np.ones((6,6)))
    assert (k,h)==(1,6) and not d["ssm_informative"]


@pytest.mark.parametrize("centered",[False,True])
def test_gram_evd_and_direct_svd_same_filter(centered):
    a=random((31,8,9))
    for k in (0,1,7):
        direct=svd_filter(a,discard=k,upper=22,centered=centered).iq
        gram=svd_filter(a,discard=k,upper=22,centered=centered,backend="gram").iq
        assert_allclose(gram,direct,atol=8e-13)


def test_rsvd_full_range_exact_and_deterministic():
    a=random((24,8,7))
    r=randomized_svd_filter(a,discard=3,oversampling=24,power_iterations=0,seed=7)
    assert_allclose(r.iq,svd_filter(a,discard=3).iq,atol=8e-13)
    assert_allclose(randomized_svd_filter(a,discard=3,seed=7).iq,
                    randomized_svd_filter(a,discard=3,seed=7).iq,atol=0)


def test_rsvd_rejects_clutter_and_preserves_blood_with_spatial_subsampling():
    # Spatial subsampling changes the covariance estimate. Measure component
    # rejection/preservation against known orthogonal source subspaces.
    u,_=linalg.qr(random((100,20)),mode="economic")
    v,_=linalg.qr(random((40,20),1),mode="economic")
    s=np.r_[1000.,500.,np.ones(18)]
    a=frames((u*s)@v.conj().T,(40,10,10))
    r=randomized_svd_filter(a,discard=2,spatial_fraction=.5,seed=4)
    tissue=frames((u[:,:2]*s[:2])@v[:,:2].conj().T,a.shape)
    blood=frames(u[:,2:]@v[:,2:].conj().T,a.shape)
    assert np.linalg.norm(r.apply_frozen(tissue))**2/np.linalg.norm(tissue)**2<1e-6
    assert np.linalg.norm(r.apply_frozen(blood))**2/np.linalg.norm(blood)**2>.99


def test_huang_median_uses_eigenvalues_not_singular_values():
    u=np.eye(5);s=np.array([100,9,4,1,.1])
    wt,d=huang_weights(u,s,rule="median",beta=1.43)
    assert d["discard"]==2
    assert_allclose(wt,[0,0,1,1,1])


def test_huang_ncc_known_disjoint_spatial_regions():
    u=np.zeros((40,12))
    u[:20,:5]=abs(random((20,5)))
    u[20:,5:]=abs(random((20,7),4))
    wt,d=huang_weights(u,np.arange(12,0,-1),rule="ncc",neighbors=3)
    assert d["discard"]==5 and d["minimum_image_ncc"]==0
    soft,_=huang_weights(u,np.arange(12,0,-1),rule="sigmoid",neighbors=3,beta=2)
    assert_allclose(soft[4],.5)
    assert np.all(np.diff(soft)>0)


def test_huang_covariance_weighting_before_projection():
    a=random((21,5,6));w=np.linspace(.01,1,30).reshape(5,6)
    x=casorati(a)
    u,s,vh=linalg.svd(w.reshape(-1,1)*x,full_matrices=False)
    wt,_=huang_weights(u,s,"median",beta=3)
    expected=frames(((x@vh.conj().T)*wt)@vh,a.shape)
    r=huang_evd_filter(a,rule="median",beta=3,spatial_weight=w)
    assert_allclose(r.iq,expected,atol=3e-13)


@pytest.mark.parametrize("fn",[
    lambda a:gaussian_local_evd(a,rule="median",beta=4),
    lambda a:song_local_filter(a,patch_shape=(4,4),stride=(2,2),selector_options={"prf_hz":1000.}),
    lambda a:svd_filter(a,selector="baranger"),
    lambda a:randomized_svd_filter(a,seed=8)])
def test_new_frozen_operators_do_not_refit_truth(fn):
    b=random((32,6,7));t=20*random((1,6,7),2)*np.ones((32,1,1));n=.1*random(b.shape,4)
    r=fn(b+t+n)
    assert_allclose(r.apply_frozen(b)+r.apply_frozen(t)+r.apply_frozen(n),r.iq,atol=1e-12)


def test_absolute_frequency_keeps_opposite_signs_equal():
    time=np.arange(100)/1000
    f=temporal_statistics(np.exp(2j*np.pi*np.array([50.,-50.])[:,None]*time),1000)
    assert_allclose(f,[50,50],atol=1e-12)


@pytest.mark.parametrize("phase_mode,expected",[("causal",1/np.sqrt(2)),("zero_phase",.5)])
def test_butterworth_analytic_cutoff_response(phase_mode,expected):
    # At the design cutoff one-pass amplitude is 1/sqrt(2); a two-pass
    # zero-phase application squares it. Test both Doppler signs.
    time=np.arange(8192)/1000
    a=np.exp(2j*np.pi*50*time)[:,None,None]
    r=temporal_highpass(a,50,1000,order=3,phase_mode=phase_mode)
    ratio=np.mean(abs(r.iq[1000:-1000]))
    assert_allclose(ratio,expected,rtol=1e-9)
    assert_allclose(temporal_highpass(a.conj(),50,1000,phase_mode=phase_mode).iq,r.iq.conj(),atol=1e-12)
    assert_allclose(r.apply_frozen(a),r.iq,atol=0)


def test_butterworth_dc_and_invalid_ensemble():
    a=np.full((80,2,3),3+2j)
    assert np.linalg.norm(temporal_highpass(a,30,1000).iq)<1e-12
    with pytest.raises(ValueError):temporal_highpass(a[:3],30,1000)
    with pytest.raises(ValueError):temporal_highpass(a,500,1000)


def test_frequency_filters_reject_wrong_acquisition_revisit_rate():
    a=random((32,4,5))
    with pytest.raises(ValueError,match="revisit rate"):
        run_filter(a,{"method":"highpass","params":{"cutoff_hz":30,"prf_hz":500}},prf_hz=1000)
    with pytest.raises(ValueError,match="revisit rate"):
        run_filter(a,{"method":"svd","params":{"selector":"temporal_frequency",
            "selector_options":{"cutoff_hz":30,"prf_hz":500}}},prf_hz=1000)
    r=run_filter(a,{"method":"highpass","params":{"cutoff_hz":30}},prf_hz=1000)
    assert r.diagnostics["prf_hz"]==1000
