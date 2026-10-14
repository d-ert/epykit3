"""Precision-weighted local discovery at predefined genomic bandwidths.

Every scale is independently seeded and raw-effect trimmed. Greedy score
competition removes overlapping scale duplicates rather than union-merging
their boundaries. Depth/support weights are approximate precision weights,
not a claim of independent CpGs or an estimated biological covariance.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ._regional_gls_native import discover_prepared, prepare_discovery


@njit(cache=True, nogil=True)
def _smooth_tricube(pos, effect, weights, radius, starts, stops):
    smooth = np.zeros(len(pos))
    for segment in range(len(starts)):
        lo, hi = starts[segment], stops[segment]
        left, right = lo, lo
        for j in range(lo, hi):
            while left < hi and pos[left] < pos[j]-radius:
                left += 1
            while right < hi and pos[right] <= pos[j]+radius:
                right += 1
            numerator, denominator = 0.0, 0.0
            for t in range(left, right):
                distance = abs(pos[t]-pos[j])/radius
                kernel = (1.0-distance**3)**3
                weight = kernel*weights[t]
                numerator += weight*effect[t]
                denominator += weight
            if denominator > 0:
                smooth[j] = numerator/denominator
    return smooth


@njit(cache=True, nogil=True)
def _distinct_candidates(bounds, effect, tau):
    scores = np.empty(len(bounds))
    for j in range(len(bounds)):
        lo, hi = bounds[j]
        scores[j] = abs(effect[lo:hi].sum())-tau*(hi-lo)
    order = np.argsort(scores)[::-1]
    accepted = np.empty_like(bounds)
    count = 0
    for index in order:
        if scores[index] <= 0:
            continue
        lo, hi = bounds[index]
        overlap = False
        for previous in range(count):
            if lo < accepted[previous, 1] and hi > accepted[previous, 0]:
                overlap = True
                break
        if not overlap:
            accepted[count, 0], accepted[count, 1] = lo, hi
            count += 1
    result = accepted[:count]
    return result[np.argsort(result[:, 0])]


def discover_multiscale(pos, effects, precision, window_bp, min_effect, min_cpgs,
                        max_gap, min_length, boundary_min_effect=0.0):
    pos = np.asarray(pos, dtype=np.int64)
    effects, precision = np.asarray(effects, float), np.asarray(precision, float)
    if effects.shape != pos.shape or precision.shape != pos.shape or not np.isfinite(effects).all() or not np.isfinite(precision).all() or np.any(precision <= 0):
        raise ValueError("Kernel effects and positive precision must match genomic coordinates")
    if len(pos) < min_cpgs:
        return np.empty((0, 2), dtype=np.int64)
    geometry = prepare_discovery(pos, window_bp, min_cpgs, max_gap)
    weights = np.minimum(precision / np.median(precision), 4.0)
    radii = sorted({max(1, window_bp//4), max(1, window_bp//2), window_bp})
    candidates = []
    for radius in radii:
        smooth = _smooth_tricube(pos, effects, weights, radius, geometry[2], geometry[3])
        bounds = discover_prepared(pos, effects, *geometry, min_effect, min_cpgs, max_gap,
                                   min_length, boundary_min_effect, smoothed=smooth, boundary_effects=effects)
        if len(bounds):
            candidates.append(bounds)
    if not candidates:
        return np.empty((0, 2), dtype=np.int64)
    tau = .55*min_effect
    return _distinct_candidates(np.concatenate(candidates), effects, tau)
