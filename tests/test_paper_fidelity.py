import numpy as np
import pytest
from numpy.testing import assert_allclose
from btbench.core import casorati,frames,patches
from btbench.filters import svd_filter,aer_filter
from btbench.paper_filters import (lag_one_mode_frequency,song_pdf_band,
    song_pdf_local_filter,maresca_cutoff,lee_aer_filter)


def test_song_uses_lag_one_frequency_and_combines_both_lower_cutoffs():
    n=100;f=np.array([0.,4.,8.,30.,100.,200.])
    vh=np.exp(2j*np.pi*f[:,None]*np.arange(n)/1000)/np.sqrt(n)
    assert_allclose(lag_one_mode_frequency(vh,1000),f,atol=1e-12)
    k,h,d=song_pdf_band(np.array([100.,50.,49.,48.,47.,46.]),vh,1000,
                         noise_rejection=False)
    assert d['lower_slope']==2 and d['lower_frequency']==3
    assert (k,h)==(3,6)


def test_song_noise_line_is_fitted_in_log_domain():
    r=20;n=100
    log_sigma=np.r_[[4.,2.5,1.5,1.1],1.-.03*np.arange(16)]
    vh=np.exp(2j*np.pi*np.linspace(0,450,r)[:,None]*np.arange(n)/1000)
    _,_,d=song_pdf_band(10**log_sigma,vh,1000)
    assert d['noise_line_domain']=='log10 singular values'
    assert_allclose(d['noise_line_slope'],-.03,atol=1e-12)


def test_song_equation_one_normalizes_then_coherently_averages_blocks():
    rng=np.random.default_rng(12);x=rng.normal(size=(12,4,6))+1j*rng.normal(size=(12,4,6))
    options={'prf_hz':1000.,'slope_threshold':1e9,'tissue_cutoff_hz':0.,'noise_rejection':False}
    result=song_pdf_local_filter(x,patch_shape=(4,4),stride=(4,2),selector_options=options)
    expected=np.zeros_like(x);count=np.zeros(x.shape[1:])
    for sl,_ in patches(x.shape[1:],(4,4),(4,2)):
        u,s,vh=np.linalg.svd(casorati(x[(slice(None),)+sl]),full_matrices=False)
        value=(u[:,1:]*s[1:])@vh[1:]/s[1:].sum()
        expected[(slice(None),)+sl]+=frames(value,(12,4,4));count[sl]+=1
    assert_allclose(result.iq,expected/count,atol=1e-13)
    assert_allclose(result.power_doppler,np.mean(abs(result.iq)**2,axis=0))
    a=.3*x;b=.7*x
    assert_allclose(result.apply_frozen(a)+result.apply_frozen(b),result.iq,atol=1e-13)


def test_maresca_uses_normalized_singular_values_curvature_and_fixed_offset():
    s=10**(np.array([0.,-15.,-30.,-31.,-32.,-33.,-34.,-35.])/20)
    k,d=maresca_cutoff(s,offset=2)
    assert d['curvature_peak_one_based']==3
    assert k==5
    assert maresca_cutoff(s*700,offset=2)[0]==k
    assert maresca_cutoff(s,offset=26)[0]==len(s)


@pytest.mark.parametrize('domain',['temporal','spatial'])
def test_lee_algorithms_match_explicit_complex_covariance_and_regularization(domain):
    rng=np.random.default_rng(31);iq=rng.normal(size=(17,3,4))+1j*rng.normal(size=(17,3,4))
    x=casorati(iq);x=x-x.mean(axis=1,keepdims=True);eps=.37
    c=x.conj().T@x/len(x) if domain=='temporal' else x@x.conj().T/(x.shape[1]-1)
    ev,u=np.linalg.eigh(c);w=(u/np.sqrt(np.maximum(ev,0)+eps))@u.conj().T
    expected=x@w if domain=='temporal' else w@x
    r=lee_aer_filter(iq,domain=domain,epsilon_abs=eps,stream_chunk=4)
    assert_allclose(casorati(r.iq),expected,atol=2e-13)
    assert r.diagnostics['covariance_denominator']==(12 if domain=='temporal' else 16)
    assert_allclose(r.apply_frozen(.4*iq)+r.apply_frozen(.6*iq),r.iq,atol=1e-13)


def test_unit_operator_aer_is_labeled_and_does_not_fit_truth():
    rng=np.random.default_rng(41);iq=rng.normal(size=(20,3,4))+1j*rng.normal(size=(20,3,4))
    r=lee_aer_filter(iq,epsilon_rel=1e-3,gain_normalization='unit_operator')
    assert r.diagnostics['gain_normalization']=='unit_operator'
    assert r.diagnostics['operator_gain_max']==1.
    assert np.linalg.norm(r.iq)<=np.linalg.norm(iq)


@pytest.mark.parametrize('domain',['temporal','spatial'])
def test_lee_aer_treats_opposite_doppler_signs_symmetrically(domain):
    rng=np.random.default_rng(19);x=rng.normal(size=(31,4,5))+1j*rng.normal(size=(31,4,5))
    a=lee_aer_filter(x,domain=domain,epsilon_rel=1e-3)
    b=lee_aer_filter(x.conj(),domain=domain,epsilon_rel=1e-3)
    assert_allclose(a.iq.conj(),b.iq,atol=5e-13)
