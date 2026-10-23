import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import ndimage
from dataclasses import replace
from btbench.registration import (MotionGate,linear_sample,projected_residual,register_displacement)
from btbench.shared import (sample_frame,inverse_coordinates,shared_scene,SharedContext,fit_shared)


def reference(seed=601,shape=(8,16)):
    rng=np.random.default_rng(seed)
    a=ndimage.gaussian_filter(rng.normal(size=(3,*shape))+1j*rng.normal(size=(3,*shape)),(0,.65,.65))
    return a/np.sqrt(np.mean(abs(a)**2))


def displaced(ref,offset):
    d=np.broadcast_to(np.array(offset)[:,None,None],(2,*ref.shape[1:])).copy()
    # Native image uses the original generator's explicit out-of-domain zero.
    return sample_frame(ref,inverse_coordinates(d)),d


def test_continuous_sampler_preserves_edges_for_infinitesimal_motion():
    ref=np.ones((3,8,16),complex);grid=np.indices((8,16),dtype=float)
    assert_allclose(sample_frame(ref,grid+1e-9,zero_extension=True),ref,atol=3e-9)
    q=grid-.37
    assert_allclose(linear_sample(ref,q),sample_frame(ref,q,zero_extension=True),atol=2e-15)
    rng=np.random.default_rng(7);q=np.array([rng.uniform(-1,8,80),rng.uniform(-1,16,80)])
    ref=reference()
    expected=np.stack([ndimage.map_coordinates(a.real,q,order=1,mode='grid-constant')+
                       1j*ndimage.map_coordinates(a.imag,q,order=1,mode='grid-constant') for a in ref])
    assert_allclose(linear_sample(ref,q),expected,atol=2e-15)


def test_variable_projection_jacobian_matches_numerical_derivative():
    ref=reference(shape=(12,16));grid=np.indices((12,16),dtype=float)[:,2:-2,2:-2].reshape(2,-1)
    p=np.array([.23,-.31]);basis=np.ones((1,grid.shape[1]));rng=np.random.default_rng(9)
    y=linear_sample(ref,grid+np.array([.4,-.2])[:,None])*(.9+.2j)+.01*rng.normal(size=(3,grid.shape[1]))
    _,jac,_=projected_residual(p,ref,grid,y,basis)
    for k in range(2):
        dp=np.zeros(2);dp[k]=1e-6
        numerical=(projected_residual(p+dp,ref,grid,y,basis)[0]-projected_residual(p-dp,ref,grid,y,basis)[0])/2e-6
        assert_allclose(jac[:,k],numerical,atol=2e-9,rtol=1e-7)


@pytest.mark.parametrize('mode',['rigid','smooth'])
def test_static_views_with_photometric_phase_changes_return_exact_zero(mode):
    ref=reference();rng=np.random.default_rng(21);g=np.exp(1j*np.arange(12).reshape(4,3))[:,:,None,None]
    views=ref[None]*g+1e-4*(rng.normal(size=(4,*ref.shape))+1j*rng.normal(size=(4,*ref.shape)))
    d,info=register_displacement(views,np.mean(abs(ref),axis=0),np.ones(ref.shape[1:],bool),mode,reference_views=ref)
    assert np.array_equal(d,np.zeros_like(d)) and info['accepted_frames']==0
    assert info['reason_counts']['reference_consistent']==len(views)


@pytest.mark.parametrize('offset',[(.3,-.25),(-.7,.4),(.65,.7),(-.12,-.2)])
def test_forward_fit_recovers_subpixel_motion_with_boundary_loss(offset):
    ref=reference();image,truth=displaced(ref,offset);roi=np.ones(ref.shape[1:],bool)
    d,info=register_displacement(image[None],np.mean(abs(ref),axis=0),roi,'rigid',reference_views=ref)
    assert info['accepted_frames']==1,info
    assert_allclose(d[0],truth,atol=.02)


def test_gate_validation_pixels_cannot_change_the_fitted_rigid_parameters():
    ref=reference();image,_=displaced(ref,(.6,-.5));roi=np.ones(ref.shape[1:],bool)
    kwargs=dict(reference=np.mean(abs(ref),axis=0),roi=roi,mode='rigid',reference_views=ref)
    d,a=register_displacement(image[None],**kwargs)
    z,x=np.indices(roi.shape);check=roi&((z+2*x)%3==0);check[[0,-1]]=False;check[:,[0,-1]]=False
    changed=image.copy();changed[:,check]+=8*reference(787)[:,check]
    rejected,b=register_displacement(changed[None],**kwargs)
    assert a['parameters_sha256']==b['parameters_sha256']
    assert a['accepted_frames']==1 and b['accepted_frames']==0
    assert np.array_equal(rejected,np.zeros_like(rejected))


@pytest.mark.parametrize('kind',['flat','wrong_reference','small_roi'])
def test_gate_rejects_unobservable_or_inconsistent_inputs(kind):
    ref=reference();roi=np.ones(ref.shape[1:],bool)
    if kind=='flat':ref[:]=1;image=2j*ref
    elif kind=='wrong_reference':image=reference(787)
    else:image,_=displaced(ref,(.5,.4));roi[:]=False;roi[2:4,4:6]=True
    d,info=register_displacement(image[None],np.mean(abs(ref),axis=0),roi,'smooth',reference_views=ref)
    assert info['accepted_frames']==0 and not d.any()


def test_ungated_arm_is_explicit_and_preserves_diagnostic_fit():
    ref=reference();image,_=displaced(ref,(.005,.004));roi=np.ones(ref.shape[1:],bool)
    args=dict(reference=np.mean(abs(ref),axis=0),roi=roi,reference_views=ref)
    gated,info=register_displacement(image[None],mode='rigid',**args)
    raw,other=register_displacement(image[None],mode='rigid_ungated',**args)
    assert not gated.any() and raw.any() and not other['gate_enabled']
    assert info['accepted_frames']==0


def test_accepted_motion_is_packet_causal_and_frozen_map_remains_linear():
    a=shared_scene(seed=787,n_frames=64,calibration_frames=32,motion='translation')
    c=SharedContext(a.public,16,32,'rigid');assert c.motion_diagnostics['accepted_frames']>0
    spec=dict(method='weighted_svd',horizon=32,discard=1,rho=.75)
    r=fit_shared(c,spec)
    assert_allclose(r.apply_frozen(a.blood)+r.apply_frozen(a.tissue)+r.apply_frozen(a.noise),r.views,atol=1e-10)
    changed=a.public.views.copy();changed[48:]*=4*np.exp(.3j)
    other=SharedContext(replace(a.public,views=changed),16,32,'rigid')
    assert_allclose(fit_shared(other,spec).views[:16],r.views[:16],atol=2e-11)


def test_gate_fallback_is_identical_to_unregistered_filter():
    a=shared_scene(seed=601,n_frames=32,calibration_frames=16)
    c=SharedContext(a.public,16,16,'smooth');unregistered=SharedContext(a.public,16,16,'none')
    assert not c.displacement.any()
    spec=dict(method='svd',horizon=16,discard=1)
    assert_allclose(fit_shared(c,spec).views,fit_shared(unregistered,spec).views,atol=0,rtol=0)
