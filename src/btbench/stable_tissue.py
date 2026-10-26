"""Data-only, motion-stabilized tissue estimate on the original IQ grid.

The brightest spatial quartile trains one temporal tissue direction, while
every pixel is filtered by that direction. This is a deliberately narrow
candidate: a spatially restricted fit does not guarantee that slow blood is
orthogonal to the fitted temporal mode.
"""

import numpy as np
from scipy import ndimage

from .core import FilterResult
from .motion_first import MotionFit, _inverse_coordinates


def smooth_motion_path(motion, sigma_frames=1.0):
    """Smooth estimated framewise coordinates, anchoring frame zero.

    The original mixture-derived phase is retained. This regularizes a fitted
    B-mode path and does not use imposed displacement or source truth.
    """
    if not np.isfinite(sigma_frames) or sigma_frames <= 0:
        raise ValueError("Need a positive finite temporal smoothing width")
    grid = np.indices(motion.forward.shape[2:], dtype=float)
    displacement = motion.forward - grid[None]
    displacement = ndimage.gaussian_filter1d(displacement, sigma_frames,
                                               axis=0, mode="nearest")
    displacement -= displacement[0].copy()
    forward = grid[None] + displacement
    inverse = np.stack([_inverse_coordinates(q) for q in forward])
    diagnostics = {"source":"mixture-fitted envelope motion",
                   "smoothing_sigma_frames":float(sigma_frames),
                   "frame_zero_anchored":True,
                   "phase":"original per-view fitted IQ phase retained",
                   "base":motion.diagnostics}
    return MotionFit(forward, inverse, motion.phase.copy(), diagnostics,
                     motion.scalar_input)


def bright_reference_subtraction(iq, motion, fraction=0.25, phase=True):
    """Fit temporal rank one in bright pixels; subtract its native-grid image.

    The support is inferred from the observed aligned IQ median magnitude and
    shared by all views. ``fraction=1`` is the full-field control. The frozen
    replay is linear; the re-fitted filter across separate acquisitions is not.
    """
    x = np.asarray(iq)
    if (x.ndim not in (3, 4) or len(x) < 3 or not np.iscomplexobj(x)
            or not np.isfinite(x).all() or not np.isfinite(fraction)
            or not 0 < fraction <= 1):
        raise ValueError("Need finite complex IQ and support fraction in (0,1]")
    aligned = motion.align(x, phase=phase)
    a = aligned[:, None] if x.ndim == 3 else aligned
    n, views, nz, nx = a.shape
    if fraction == 1:
        support = np.ones((nz, nx), bool)
    else:
        brightness = np.median(abs(a), axis=(0, 1))
        count = max(1, int(np.ceil(fraction * brightness.size)))
        order = np.argsort(brightness.ravel(), kind="stable")
        support = np.zeros(brightness.size, bool)
        support[order[-count:]] = True
        support = support.reshape(brightness.shape)
    training = a[:, :, support].reshape(n, -1).T
    _, singular_values, vh = np.linalg.svd(training, full_matrices=False)
    mode = vh[0].conj()

    def replay(other):
        other = np.asarray(other)
        if other.shape != x.shape or not np.iscomplexobj(other) or not np.isfinite(other).all():
            raise ValueError("Frozen bright-reference IQ shape/type mismatch")
        transformed = motion.align(other, phase=phase)
        matrix = transformed.reshape(n, -1).T
        removed = (matrix @ mode)[:, None] * mode.conj()[None, :]
        return other - motion.restore(removed.T.reshape(x.shape), phase=phase)

    return FilterResult(replay(x), {
        "method":"bright_reference_native_tissue_subtraction",
        "training_fraction":float(fraction),"training_pixels":int(support.sum()),
        "rank":1,"phase_alignment":bool(phase),
        "leading_training_sigma_fraction":float(singular_values[0] /
                                                max(np.linalg.norm(singular_values), 1e-30)),
        "formula":"X - R[P_(bright-trained temporal rank one) A(X)]",
        "motion":motion.diagnostics}, replay)
