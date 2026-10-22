from dataclasses import dataclass, field
from itertools import product
from typing import Any, Callable

import numpy as np


@dataclass
class Acquisition:
    """Canonical IQ layout: (slow_time, z, x), SI units, uniformly sampled.

    Reference fields are evaluation-only. No filter receives this object.
    Angles for vector Doppler are supplied separately to measurements.py.
    """
    iq: np.ndarray
    prf_hz: float
    f0_hz: float
    sound_speed_mps: float = 1540.0
    phase_sign: int = 1
    timestamps_s: np.ndarray | None = None
    iq_views: np.ndarray | None = None
    z_m: np.ndarray | None = None
    x_m: np.ndarray | None = None
    tissue_roi: np.ndarray | None = None
    flow_mask: np.ndarray | None = None
    near_wall_mask: np.ndarray | None = None
    blood_true: np.ndarray | None = None
    tissue_true: np.ndarray | None = None
    velocity_true_mps: np.ndarray | None = None
    wall: dict | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.iq = np.asarray(self.iq)
        if self.iq.ndim != 3 or not np.iscomplexobj(self.iq):
            raise ValueError("iq must be complex with shape (time,z,x); B-mode/PD is invalid")
        if self.iq.shape[0] < 3 or not np.isfinite(self.iq).all():
            raise ValueError("Need >=3 finite IQ frames")
        if min(self.prf_hz, self.f0_hz, self.sound_speed_mps) <= 0:
            raise ValueError("PRF, frequency, and sound speed must be positive")
        if self.phase_sign not in (-1, 1):
            raise ValueError("phase_sign must be explicitly +1 or -1")
        if self.timestamps_s is not None:
            t=np.asarray(self.timestamps_s)
            if t.shape!=(len(self.iq),) or not np.allclose(np.diff(t),1/self.prf_hz,rtol=1e-4,atol=1e-9):
                raise ValueError("Nonuniform/gapped slow time: split contiguous segments before filtering")
        if self.iq_views is not None:
            a=np.asarray(self.iq_views)
            if a.ndim!=4 or (a.shape[0],*a.shape[2:])!=self.iq.shape or not np.iscomplexobj(a):
                raise ValueError("iq_views must be complex (time,view,z,x), matched to iq")
        for name in ("tissue_roi", "flow_mask", "near_wall_mask"):
            a = getattr(self, name)
            if a is not None and np.shape(a) != self.iq.shape[1:]:
                raise ValueError(f"{name} must match spatial shape")
        if self.velocity_true_mps is not None and np.shape(self.velocity_true_mps) not in (self.iq.shape[1:],self.iq.shape):
            raise ValueError("velocity truth must be (z,x) or (time,z,x)")
        for name in ("blood_true", "tissue_true"):
            a = getattr(self, name)
            if a is not None and np.shape(a) != self.iq.shape:
                raise ValueError(f"{name} must match iq shape")
        for coord, n in ((self.z_m, self.iq.shape[1]), (self.x_m, self.iq.shape[2])):
            if coord is not None and (np.shape(coord) != (n,) or np.any(np.diff(coord) <= 0)):
                raise ValueError("Grid coordinates must be increasing, one-dimensional, in meters")


@dataclass
class FilterResult:
    iq: np.ndarray
    diagnostics: dict[str, Any] = field(default_factory=dict)
    apply_frozen: Callable | None = field(default=None,repr=False)


def casorati(iq):
    """(T,*space) -> (P,T); transpose, never conjugate the IQ values."""
    return np.asarray(iq).reshape(len(iq), -1).T.astype(np.complex128, copy=False)


def frames(x, shape):
    return x.T.reshape(shape)


def prepare(x, centered=False):
    mean = x.mean(axis=1, keepdims=True) if centered else np.zeros((len(x), 1), complex)
    return x - mean, mean


def chunks(n, size):
    if size < 3:
        raise ValueError("Packet/window length must be >=3")
    for start in range(0, n, size):
        stop = min(start + size, n)
        # Keep final short packet with its predecessor, avoiding two-frame fits.
        if n - stop in (1, 2):
            stop = n
        yield slice(start, stop)
        if stop == n:
            return


def patches(shape, patch_shape=(8, 8), stride=(4, 4)):
    """Gaussian overlap/add; boundaries covered and weights strictly positive."""
    ps = tuple(min(int(p), n) for p, n in zip(patch_shape, shape))
    if len(ps) != len(shape) or any(s < 1 or s > p for s, p in zip(stride, ps)):
        raise ValueError("Patch stride must lie in [1, patch_size]")
    starts = []
    for n, p, s in zip(shape, ps, stride):
        a = list(range(0, n - p + 1, s))
        if not a or a[-1] != n - p:
            a.append(n - p)
        starts.append(a)
    mesh = np.meshgrid(*[np.linspace(-1, 1, p) for p in ps], indexing="ij")
    weight = np.exp(-2 * sum(a*a for a in mesh))
    for loc in product(*starts):
        yield tuple(slice(a, a+p) for a, p in zip(loc, ps)), weight
