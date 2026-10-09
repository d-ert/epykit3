"""Native batched spatial REML with analytical donor-variance profiling.

After whitening the CAR component, the sample-offset loading has one eigenvalue
1-c+c*(h'A^-1h)/median(depth); its orthogonal complement has eigenvalue 1-c.
The two residual sums of squares have (S-2) and (S-1)*(K-1) degrees of freedom.
Their variance ratio gives the exact conditional REML estimate of c. Only the
CAR length needs a bounded scalar search. No independence fallback is used.
"""

from __future__ import annotations

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads


def prepare_discovery(pos, window_bp, min_cpgs, max_gap):
    """One coordinate search per chromosome, with gaps respected at every scale."""
    pos = np.asarray(pos, dtype=np.int64)
    cuts = np.r_[0, np.flatnonzero(np.diff(pos) > max_gap) + 1, len(pos)]
    sizes = np.diff(cuts)
    left = np.maximum(np.searchsorted(pos, pos - window_bp // 2), np.repeat(cuts[:-1], sizes))
    right = np.minimum(
        np.searchsorted(pos, pos + (window_bp + 1) // 2, side="right"), np.repeat(cuts[1:], sizes)
    )
    keep = sizes >= min_cpgs
    return left, right, cuts[:-1][keep], cuts[1:][keep]


@njit(cache=True, nogil=True)
def discover_prepared(
    pos, effects, left, right, starts, stops, min_effect, min_cpgs, max_gap, min_length
):
    cumulative = np.empty(len(pos))
    regions = []
    for segment in range(len(starts)):
        lo, hi = starts[segment], stops[segment]
        total = 0.0
        for j in range(lo, hi):
            total += effects[j]
            cumulative[j] = total
        first, last, direction = -1, -1, 0
        for j in range(lo, hi):
            value = cumulative[right[j] - 1]
            if left[j] > lo:
                value -= cumulative[left[j] - 1]
            value /= right[j] - left[j]
            if abs(value) < min_effect:
                continue
            sign = 1 if value > 0 else -1
            if first >= 0 and (pos[j] - pos[last] > max_gap or sign != direction):
                if last - first + 1 >= min_cpgs and pos[last] - pos[first] + 1 >= min_length:
                    regions.append((first, last + 1))
                first = -1
            if first < 0:
                first, direction = j, sign
            last = j
        if first >= 0 and last - first + 1 >= min_cpgs and pos[last] - pos[first] + 1 >= min_length:
            regions.append((first, last + 1))
    result = np.empty((len(regions), 2), dtype=np.int64)
    for j in range(len(regions)):
        result[j, 0], result[j, 1] = regions[j]
    return result


@njit(cache=True, nogil=True)
def _prepare(yc, h, g):
    samples, k = yc.shape
    gram = 0.0
    center = np.mean(g)
    diagonal = np.zeros(k)
    adjacent = np.zeros(k - 1)
    contrast = np.zeros(k)
    for s in range(samples):
        x = g[s] - center
        gram += x * x
        for j in range(k):
            diagonal[j] += yc[s, j] * yc[s, j]
            contrast[j] += x * yc[s, j]
            if j + 1 < k:
                adjacent[j] += yc[s, j] * yc[s, j + 1]
    return diagonal, adjacent, contrast, gram, np.median(h * h)


@njit(cache=True, nogil=True)
def _components(p, yc, h, length, diagonal, adjacent, contrast, gram):
    k = len(p)
    diag = np.ones(k)
    off = np.empty(k - 1)
    logdet = 0.0
    for j in range(k - 1):
        if length == 0.0:
            a, den = 0.0, 1.0
        else:
            gap = (p[j + 1] - p[j]) / length
            a, den = np.exp(-gap), -np.expm1(-2 * gap)
        diag[j] += a * a / den
        diag[j + 1] += a * a / den
        off[j] = -a / den
        logdet += np.log(den)
    qh = diag * h
    for j in range(k - 1):
        qh[j] += off[j] * h[j + 1]
        qh[j + 1] += off[j] * h[j]
    information = np.dot(h, qh)
    total = np.dot(diag, diagonal) + 2 * np.dot(off, adjacent)
    projected = 0.0
    for s in range(yc.shape[0]):
        b = 0.0
        for j in range(k):
            b += yc[s, j] * qh[j]
        projected += b * b / information
    effect = np.dot(contrast, qh) / (gram * information)
    donor = projected - gram * information * effect * effect
    perpendicular = total - projected
    # Cancellation can leave a tiny negative sum of squares. Larger failures
    # remain unavailable and abort the entire scan in the Python boundary.
    tolerance = 1e-10 * max(total, 1.0)
    if donor < -tolerance or perpendicular < -tolerance or information <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    return information, max(perpendicular, 0.0), max(donor, 0.0), effect, logdet


@njit(cache=True, nogil=True)
def _evaluate(p, yc, h, length, c, prepared):
    diagonal, adjacent, contrast, gram, median_depth = prepared
    t, sp, sd, effect, logdet = _components(p, yc, h, length, diagonal, adjacent, contrast, gram)
    samples, k = yc.shape
    dfp, dfd = (samples - 1) * (k - 1), samples - 2
    df = dfp + dfd
    a = 1 - c
    b = a + c * t / median_depth
    rss = sp / a + sd / b
    if rss <= 0 or not np.isfinite(rss):
        return np.inf, effect, np.nan
    value = (
        df * np.log(rss / df)
        + (samples - 1) * logdet
        + dfp * np.log(a)
        + dfd * np.log(b)
        + np.log(t)
    )
    se = np.sqrt((rss / df) / (gram * t / b))
    return value, effect, se


@njit(cache=True, nogil=True)
def _profile(p, yc, h, length, prepared):
    diagonal, adjacent, contrast, gram, median_depth = prepared
    t, sp, sd, effect, logdet = _components(p, yc, h, length, diagonal, adjacent, contrast, gram)
    samples, k = yc.shape
    dfp, dfd = (samples - 1) * (k - 1), samples - 2
    if not np.isfinite(t):
        return np.inf, 0.0, effect, np.nan
    if sp <= 0:
        c = 0.98 if sd > 0 else 0.0
    else:
        ratio = (sd / dfd) / (sp / dfp)
        c = min(0.98, (ratio - 1) / (ratio - 1 + t / median_depth)) if ratio > 1 else 0.0
    a, b = 1 - c, 1 - c + c * t / median_depth
    rss = sp / a + sd / b
    df = dfp + dfd
    if rss <= 0 or not np.isfinite(rss):
        return np.inf, c, effect, np.nan
    value = (
        df * np.log(rss / df)
        + (samples - 1) * logdet
        + dfp * np.log(a)
        + dfd * np.log(b)
        + np.log(t)
    )
    se = np.sqrt((rss / df) / (gram * t / b))
    return value, c, effect, se


@njit(cache=True, nogil=True)
def _golden(p, yc, h, prepared, left, right):
    golden = 0.6180339887498949
    a = right - golden * (right - left)
    b = left + golden * (right - left)
    fa = _profile(p, yc, h, np.exp(a), prepared)[0]
    fb = _profile(p, yc, h, np.exp(b), prepared)[0]
    for _ in range(48):
        if right - left < 1e-7:
            break
        if fa < fb:
            right, b, fb = b, a, fa
            a = right - golden * (right - left)
            fa = _profile(p, yc, h, np.exp(a), prepared)[0]
        else:
            left, a, fa = a, b, fb
            b = left + golden * (right - left)
            fb = _profile(p, yc, h, np.exp(b), prepared)[0]
    return a if fa < fb else b


@njit(cache=True, nogil=True)
def _fit(p, yc, h, g):
    prepared = _prepare(yc, h, g)
    best, c, effect, se = _profile(p, yc, h, 0.0, prepared)
    best_length = 0.0
    # Search each sampled local basin; checking endpoints alone would miss
    # interior minima and optimising one start can miss a second spatial scale.
    grid = np.linspace(np.log(0.1), np.log(10000.0), 17)
    values = np.empty(len(grid))
    for j in range(len(grid)):
        length = np.exp(grid[j])
        value, cc, ee, ss = _profile(p, yc, h, length, prepared)
        values[j] = value
        if value < best:
            best, best_length, c, effect, se = value, length, cc, ee, ss
    for j in range(1, len(grid) - 1):
        if values[j] < values[j - 1] and values[j] < values[j + 1]:
            length = np.exp(_golden(p, yc, h, prepared, grid[j - 1], grid[j + 1]))
            value, cc, ee, ss = _profile(p, yc, h, length, prepared)
            if value < best:
                best, best_length, c, effect, se = value, length, cc, ee, ss
    return best_length, c, effect, se, best


@njit(cache=True, nogil=True, parallel=True)
def _batch(p, yc, h, g, bounds):
    result = np.empty((len(bounds), 5))
    for j in prange(len(bounds)):
        lo, hi = bounds[j]
        result[j, :] = _fit(p[lo:hi], yc[:, lo:hi], h[lo:hi], g)
    return result


def _validate(p, yc, h, g):
    p, yc, h, g = (np.asarray(v, dtype=np.float64) for v in (p, yc, h, g))
    if (
        yc.shape != (len(g), len(p))
        or h.shape != p.shape
        or p.ndim != 1
        or np.any(np.diff(p) <= 0)
        or np.any(h <= 0)
        or any(not np.isfinite(v).all() for v in (p, yc, h, g))
        or set(g) != {0.0, 1.0}
        or min((g == 0).sum(), (g == 1).sum()) < 2
    ):
        raise ValueError("Invalid native GLS observations or sample design")
    return p, yc, h, g


def evaluate_covariance(p, yc, h, g, length, c):
    p, yc, h, g = _validate(p, yc, h, g)
    if length < 0 or not 0 <= c < 1:
        raise ValueError("Invalid covariance parameters")
    return _evaluate(p, yc, h, length, c, _prepare(yc, h, g))


def profile_covariance(p, yc, h, g, length):
    p, yc, h, g = _validate(p, yc, h, g)
    return _profile(p, yc, h, length, _prepare(yc, h, g))


def fit_regions(p, yc, h, g, bounds, *, workers=None):
    p, yc, h, g = _validate(p, yc, h, g)
    bounds = np.asarray(bounds, dtype=np.int64).reshape(-1, 2)
    if (
        np.any(bounds[:, 0] < 0)
        or np.any(bounds[:, 1] > len(p))
        or np.any(np.diff(bounds, axis=1) < 2)
    ):
        raise ValueError("Invalid native region boundaries")
    previous = get_num_threads()
    workers = min(previous, 8) if workers is None else workers
    if not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    try:
        set_num_threads(min(workers, previous))
        result = _batch(p, yc, h, g, bounds)
    finally:
        set_num_threads(previous)
    if not np.isfinite(result).all() or np.any(result[:, 3] <= 0):
        raise ValueError(
            "Degenerate regional covariance; failed fits cannot be omitted from calibration"
        )
    return result
