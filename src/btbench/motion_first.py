"""Envelope-derived registration, complex-IQ replay and motion-window controls.

All estimates use the observed mixture. Source truth and imposed motion are
reserved for evaluation. Registration is offline and two-dimensional.
"""
from dataclasses import dataclass

import numpy as np
from scipy import ndimage, optimize

from .core import FilterResult
from .filters import svd_filter


def _views(iq):
    x = np.asarray(iq)
    if x.ndim not in (3, 4) or not np.iscomplexobj(x) or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("Need finite complex IQ (time,z,x) or (time,view,z,x)")
    return (x[:, None] if x.ndim == 3 else x), x.ndim == 3


def _coordinates(shape, coefficients):
    grid = np.indices(shape, dtype=float)
    z = 2 * grid[0] / (shape[0] - 1) - 1
    x = 2 * grid[1] / (shape[1] - 1) - 1
    terms = [np.ones(shape), z, x, np.sin(np.pi * z) * np.sin(np.pi * x)]
    return grid + np.einsum("ij,jzx->izx", coefficients, np.stack(terms))


def _sample(frame, coordinates):
    return (ndimage.map_coordinates(frame.real, coordinates, order=1, mode="constant", cval=0)
            + 1j * ndimage.map_coordinates(frame.imag, coordinates, order=1, mode="constant", cval=0))


def _inverse_coordinates(coordinates):
    grid = np.indices(coordinates.shape[1:], dtype=float)
    displacement = coordinates - grid
    q = grid.copy()
    for _ in range(16):
        q = grid - np.stack([ndimage.map_coordinates(d, q, order=1, mode="nearest")
                             for d in displacement])
    return q


def _envelope_sample(frame, coordinates):
    return ndimage.map_coordinates(frame, coordinates, order=1, mode="constant", cval=0)


@dataclass
class MotionFit:
    forward: np.ndarray
    inverse: np.ndarray
    phase: np.ndarray
    diagnostics: dict
    scalar_input: bool

    def align(self, other, phase=True):
        x, scalar = _views(other)
        if scalar != self.scalar_input or x.shape[:2] != self.phase.shape or x.shape[2:] != self.forward.shape[2:]:
            raise ValueError("Frozen motion shape mismatch")
        out = np.stack([np.stack([_sample(v, q) for v in frame])
                        for frame, q in zip(x, self.forward)])
        if phase:
            out *= np.exp(-1j * self.phase)[:, :, None, None]
        return out[:, 0] if scalar else out

    def restore(self, other, phase=True):
        x, scalar = _views(other)
        if scalar != self.scalar_input or x.shape[:2] != self.phase.shape or x.shape[2:] != self.forward.shape[2:]:
            raise ValueError("Frozen motion shape mismatch")
        if phase:
            x = x * np.exp(1j * self.phase)[:, :, None, None]
        out = np.stack([np.stack([_sample(v, q) for v in frame])
                        for frame, q in zip(x, self.inverse)])
        return out[:, 0] if scalar else out


def estimate_bmode_motion(iq, model="translation", max_shift_px=2.5):
    """Register envelope texture to frame zero; gated affine/smooth extensions.

    A bright, automatically determined tissue-candidate support is split into
    checkerboard fit/check pixels. Models are accepted when held-out BIC falls.
    Adjacent pixels are correlated, so this is a pragmatic gate, not a p-value.
    The same fitted coordinate map applies to every view and source component.
    """
    x, scalar = _views(iq)
    if model not in ("translation", "auto_affine", "auto_smooth") or max_shift_px <= 0:
        raise ValueError("Invalid registration model or maximum displacement")
    n, views, nz, nx = x.shape
    if min(nz, nx) < 12:
        raise ValueError("Need at least 12 pixels in each spatial dimension")
    envelope = np.mean(abs(x), axis=1)
    # Log compression removes the influence of a handful of saturated pixels.
    scale = max(float(np.median(envelope)), 1e-12)
    images = np.log1p(envelope / scale)
    ref = images[0]
    zz, xx = np.indices(ref.shape)
    margin = int(np.ceil(max_shift_px)) + 1
    interior = ((zz >= margin) & (zz < nz - margin) & (xx >= margin) & (xx < nx - margin))
    support = interior & (np.mean(images, axis=0) >= np.median(np.mean(images, axis=0)[interior]))
    train = support & ((zz + xx) % 2 == 0)
    check = support & ~train
    if min(train.sum(), check.sum()) < 12:
        raise ValueError("Insufficient fit/check pixels after boundary exclusion")
    grid = np.indices(ref.shape, dtype=float)
    coefficients = np.zeros((n, 2, 4))
    rows = []
    for t in range(n):
        if t == 0:
            rows.append({"frame":t,"model":"reference","check_ncc":1.,"failure":False,
                         "coefficients_px":coefficients[t].tolist()})
            continue
        frame = images[t]
        # Flatten order is z:(1,z,x,bend), x:(1,z,x,bend).
        # Fit the first two coefficients separately before the nested models.
        def objective_translation(d):
            c = np.zeros((2, 4));c[:, 0] = d
            warped = _envelope_sample(frame, _coordinates(ref.shape, c))
            return (warped - ref)[train]
        bounds = ([-max_shift_px] * 2, [max_shift_px] * 2)
        coarse = [np.array([a, b], float) for a in range(-int(max_shift_px), int(max_shift_px)+1)
                  for b in range(-int(max_shift_px), int(max_shift_px)+1)]
        start = min(coarse, key=lambda d: np.linalg.norm(objective_translation(d)))
        tr = optimize.least_squares(objective_translation, start, bounds=bounds, max_nfev=35).x
        candidates = [("identity", np.zeros((2, 4)), 0)]
        c = np.zeros((2, 4));c[:, 0] = tr
        candidates.append(("translation", c, 2))
        def fit_extension(previous, columns):
            slots = [(i, j) for i in range(2) for j in columns]
            def objective(p):
                trial = previous.copy()
                for (i, j), v in zip(slots, p):trial[i, j] = v
                return (_envelope_sample(frame, _coordinates(ref.shape, trial)) - ref)[train]
            values = optimize.least_squares(objective, np.zeros(len(slots)),
                bounds=(-.5*np.ones(len(slots)), .5*np.ones(len(slots))), max_nfev=35).x
            fitted = previous.copy()
            for (i, j), v in zip(slots, values):fitted[i, j] = v
            return fitted
        if model != "translation":
            candidates.append(("affine", fit_extension(c, (1, 2)), 6))
        if model == "auto_smooth":
            candidates.append(("smooth", fit_extension(candidates[-1][1], (3,)), 8))
        def stats(candidate):
            warped = _envelope_sample(frame, _coordinates(ref.shape, candidate))
            err = warped[check] - ref[check]
            rss = max(float(err @ err), 1e-20)
            a, b = warped[check], ref[check]
            ncc = float(np.dot(a-a.mean(), b-b.mean()) / max(np.linalg.norm(a-a.mean()) * np.linalg.norm(b-b.mean()), 1e-20))
            return rss, ncc
        chosen, choice, complexity = candidates[0]
        best_rss, best_ncc = stats(choice)
        best_bic = len(ref[check])*np.log(best_rss/len(ref[check]))
        for name, trial, parameters in candidates[1:]:
            rss, ncc = stats(trial)
            bic = len(ref[check])*np.log(rss/len(ref[check])) + parameters*np.log(len(ref[check]))
            if bic < best_bic:
                chosen, choice, complexity, best_rss, best_ncc, best_bic = name, trial, parameters, rss, ncc, bic
        failure = not np.isfinite(best_ncc) or best_ncc < .5
        if failure:
            choice = np.zeros((2, 4)); chosen = "identity_fallback"
        coefficients[t] = choice
        rows.append({"frame":t,"model":chosen,"check_ncc":best_ncc,"check_rss":best_rss,
                     "check_bic":float(best_bic),"free_parameters":complexity,
                     "failure":failure,"coefficients_px":choice.tolist()})
    forward = np.stack([_coordinates(ref.shape, c) for c in coefficients])
    inverse = np.stack([_inverse_coordinates(q) for q in forward])
    aligned = np.stack([np.stack([_sample(v, q) for v in frame])
                        for frame, q in zip(x, forward)])
    stable = support & (ref >= np.median(ref[support]))
    phase = np.zeros((n, views));coherence = np.ones((n, views))
    for t in range(1, n):
        for v in range(views):
            a, b = aligned[0, v, stable], aligned[t, v, stable]
            cross = np.vdot(a, b)
            phase[t, v] = np.angle(cross)
            coherence[t, v] = abs(cross)/max(np.linalg.norm(a)*np.linalg.norm(b), 1e-30)
    diag = {"motion_source":"observed log-envelope mean across views", "model_request":model,
            "fit_check":"checkerboard of automatic bright interior", "max_shift_px":max_shift_px,
            "roi_pixels":int(support.sum()), "failed_frames":sum(row["failure"] for row in rows),
            "phase_coherence_min":float(coherence.min()), "phase_source":"aligned bright IQ vs reference",
            "frames":rows}
    return MotionFit(forward, inverse, phase, diag, scalar)


def registration_order_filter(iq, motion, order="before", selector="fixed", discard=1,
                              phase=True, window=None):
    """Fit SVD/Baranger before or after identical frozen geometric maps."""
    if order not in ("before", "after"):
        raise ValueError("order must be before or after")
    target = motion.align(iq, phase=phase) if order == "after" else iq
    fitted = svd_filter(target, selector=selector, discard=discard, window=window)
    def replay(other):
        if np.shape(other) != np.shape(iq):raise ValueError("Frozen-filter shape mismatch")
        if order == "after":
            return motion.restore(fitted.apply_frozen(motion.align(other, phase=phase)), phase=phase)
        # The same geometric round trip after filtering controls interpolation loss.
        return motion.restore(motion.align(fitted.apply_frozen(other), phase=phase), phase=phase)
    return FilterResult(replay(iq), {"order":order,"selector":selector,"phase_alignment":phase,
        "svd":fitted.diagnostics,"motion":motion.diagnostics}, replay)


def motion_windows(iq, motion, prf_hz, f0_hz, dz_m, dx_m, window=48, min_window=32,
                   hop=8, adaptive=True, selector="fixed", discard=1, offset=26):
    """Offline overlap/add on windows limited by a one-wavelength motion path.

    This is a Zhang-inspired automatic-ROI adaptation, not a reproduction of
    the published manual myocardial ROI and calibrated per-ensemble offset.
    """
    if min(prf_hz, f0_hz, dz_m, dx_m) <= 0 or min_window < 3 or not 1 <= hop <= min_window or window < min_window:
        raise ValueError("Invalid acquisition constants or window geometry")
    n = len(iq)
    if n < min_window:raise ValueError("Ensemble shorter than minimum window")
    wavelength = 1540. / f0_hz
    centers = motion.forward[:, :, motion.forward.shape[2]//2, motion.forward.shape[3]//2]
    cumulative = np.r_[0., np.cumsum(np.linalg.norm(np.diff(centers, axis=0)*np.array([dz_m, dx_m]), axis=1))]
    starts = list(range(0, n-min_window+1, hop))
    if starts[-1] != n-min_window:starts.append(n-min_window)
    mass = np.zeros((n,)+(1,)*(iq.ndim-1))
    output = np.zeros_like(iq, complex);plans=[];records=[]
    for start in starts:
        end = min(n, start+window)
        if adaptive:
            allowed = np.searchsorted(cumulative, cumulative[start]+wavelength, side="right")
            end = min(end, max(start+min_window, allowed))
        sl = slice(start,end)
        r = svd_filter(iq[sl], selector=selector, discard=discard,
            selector_options={"offset":offset} if selector == "maresca" else None)
        weights = np.exp(-2*np.linspace(-1,1,end-start)**2).reshape((-1,)+(1,)*(iq.ndim-1))
        output[sl] += weights*r.iq;mass[sl] += weights
        plans.append((sl,weights,r))
        records.append({"start":start,"stop":end,"path_m":float(cumulative[end-1]-cumulative[start]),
                        "wavelength_exceeded_by_minimum":bool(cumulative[end-1]-cumulative[start]>wavelength),
                        "svd":r.diagnostics["windows"][0]})
    def replay(other):
        if np.shape(other) != np.shape(iq):raise ValueError("Frozen-window shape mismatch")
        y = np.zeros_like(other, complex)
        for sl,w,r in plans:y[sl] += w*r.apply_frozen(other[sl])
        return y/mass
    return FilterResult(output/mass,{"adaptive":adaptive,"path_budget_m":wavelength,
        "window_max":window,"window_min":min_window,"hop":hop,"selector":selector,
        "offset":offset if selector=="maresca" else None,"windows":records},replay)
