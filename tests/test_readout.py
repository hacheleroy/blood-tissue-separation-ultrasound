import numpy as np
from btbench.synthetic import poiseuille
from btbench.measurements import evaluate_windows
from btbench.core import Acquisition


def test_shared_readout_intervals():
    a=poiseuille(n_frames=128,shape=(4,6))
    r=evaluate_windows(a.iq,a,window=32,hop=32)
    assert [(x["start"],x["stop"]) for x in r["records"]]==[(0,32),(32,64),(64,96),(96,128)]


def test_dynamic_truth_has_matching_support():
    t=64;prf=1000;f0=5e6;c=1540
    velocity=np.r_[np.full(32,.01),np.full(32,.02)][:,None,None]
    phase=np.cumsum(4*np.pi*f0/c*velocity/prf,axis=0)
    iq=np.exp(1j*phase)
    a=Acquisition(iq,prf,f0,velocity_true_mps=velocity)
    r=evaluate_windows(iq,a,32,32)
    assert all(w["metrics"]["velocity_nrmse"]<1e-12 for w in r["records"])
