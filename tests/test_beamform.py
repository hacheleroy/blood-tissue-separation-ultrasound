import numpy as np
from numpy.testing import assert_allclose
from btbench.beamform import plane_wave_das
from btbench.core import Acquisition
import pytest


def test_rf_point_target_phase_and_focus():
    # Normal incidence, one receive element: moving target shifts round-trip
    # arrival time. Away-from-probe motion yields decreasing analytic RF phase.
    fs=40e6;f0=5e6;c=1540.;depth=.008
    t=np.arange(1024)/fs
    displacement=np.array([0.,2e-6,4e-6])
    rf=[]
    for d in displacement:
        delay=2*(depth+d)/c
        pulse=np.exp(-((t-delay)/(0.3e-6))**2)*np.cos(2*np.pi*f0*(t-delay))
        rf.append(pulse)
    rf=np.asarray(rf)[:,None,:,None]
    y=plane_wave_das(rf,[[0,0]],[depth],[0],[0],fs,f0,[0],[0])[:,0,0,0]
    measured=np.angle(y[1:]*y[:-1].conj())
    expected=-4*np.pi*f0/c*np.diff(displacement)
    assert_allclose(measured,expected,atol=.006)


def test_irregular_timestamps_rejected():
    with pytest.raises(ValueError):
        Acquisition(np.ones((4,2,2),complex),1000,5e6,timestamps_s=[0,.001,.005,.006])
