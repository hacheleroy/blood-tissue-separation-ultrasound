"""Motion-aligned tissue estimation with output on the original IQ grid.

Registration, phase, and SVD are fitted to the observed mixture. Replay uses
the frozen mapping and subspace; neither source truth nor known motion enters
the filter. This module is separate from the frozen motion-order v1 source.
"""
import numpy as np

from .core import FilterResult
from .filters import svd_filter


def native_tissue_subtraction(iq, motion, discard=1, phase=True):
    """Return X - R[P_tissue A(X)] for the mixture-fitted SVD tissue modes.

    A is the fitted B-mode geometry (with optional IQ phase), R its numerical
    restoration, and P_tissue removes exactly the first ``discard`` modes.
    This keeps the observed IQ X on its native grid. The removed tissue
    estimate alone undergoes inverse interpolation; no extra rank or spatial
    mask is fitted to the source components.
    """
    x = np.asarray(iq)
    if discard < 1 or int(discard) != discard:
        raise ValueError("Need a positive integer number of tissue modes")
    aligned = motion.align(x, phase=phase)
    fitted = svd_filter(aligned, selector="fixed", discard=int(discard))

    def replay(other):
        other = np.asarray(other)
        if other.shape != x.shape or not np.iscomplexobj(other):
            raise ValueError("Frozen tissue-subtraction shape/type mismatch")
        a = motion.align(other, phase=phase)
        removed = a - fitted.apply_frozen(a)
        return other - motion.restore(removed, phase=phase)

    return FilterResult(replay(x), {"method":"native_grid_tissue_subtraction",
        "discard":int(discard), "phase_alignment":bool(phase),
        "fitted_filter":fitted.diagnostics, "motion":motion.diagnostics,
        "formula":"X - R[P_tissue A(X)]"}, replay)
