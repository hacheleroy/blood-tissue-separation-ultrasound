import numpy as np

from btbench.esdiva_aperture import delayed_groups
from btbench.multidim import (factor_tensor,intersection_filter,minimum_gain,
                               oracle_joint,paper_tensor,svd_unfolding,tensor_gate)


def test_joint_angle_oracle_adds_information_at_shared_temporal_mode():
    time=np.ones(16,complex)
    tissue=np.zeros((16,2,3),complex)
    blood=np.zeros_like(tissue)
    tissue[:,0,:]=time[:,None]*np.array([1.,.3,.8])
    blood[:,0,:]=time[:,None]*np.array([.4,.2,.1])
    # Orthogonalize the blood angular signature to the tissue angle signature.
    g=tissue[0,0]
    k=blood[0,0]-np.vdot(g,blood[0,0])/np.vdot(g,g)*g
    blood[:,0,:]=time[:,None]*k[None]
    blood[:,1,:]=.5*blood[:,0,:]
    joint,_=oracle_joint(tissue)
    assert np.linalg.norm(joint(tissue)-tissue)<1e-10
    assert np.linalg.norm(joint(blood))<1e-10
    temporal,_=oracle_joint(tissue.mean(axis=-1))
    assert np.linalg.norm(temporal(blood.mean(axis=-1)))>0


def test_hosvd_reconstruction_and_frozen_source_linearity():
    rng=np.random.default_rng(8)
    x=rng.normal(size=(6,7,4))+1j*rng.normal(size=(6,7,4))
    y=rng.normal(size=x.shape)+1j*rng.normal(size=x.shape)
    op=tensor_gate(x,np.ones(x.shape,bool),factor_tensor(x))
    assert np.allclose(op.apply(y),y,atol=1e-12)
    for filt in (paper_tensor(x,family='ozgun'),paper_tensor(x,family='jiang'),
                 intersection_filter(x),svd_unfolding(x)):
        assert np.allclose(filt.apply(x+y),filt.apply(x)+filt.apply(y),atol=1e-11)


def test_receive_groups_partition_same_delayed_channel_signal():
    rng=np.random.default_rng(4)
    rf=rng.normal(size=(8,130))
    centers=np.column_stack((np.linspace(-.002,.002,8),np.zeros(8),np.zeros(8)))
    event={'angles_deg':0.,'delays':np.zeros(8),'apodization':np.ones(8)}
    kw=dict(rf=rf,t0=0.,dt=2e-7,event=event,centers=centers,
            z_mm=np.array([7.,8.]),x_mm=np.array([-.2,.2]),c=1540.)
    a=delayed_groups(groups=4,**kw)
    b=delayed_groups(groups=1,**kw)
    assert a.shape==(4,2,2)
    assert np.allclose(a.sum(axis=0),b[0],rtol=1e-6,atol=1e-6)


def test_minimum_gain_matches_constraint():
    x=np.array([3.,2.],complex)
    assert minimum_gain(x,x,x,.1) is not None
    a=minimum_gain(x,x,x,.1)
    assert np.isclose(np.linalg.norm(x-a*x)**2/np.linalg.norm(x)**2,.1)
