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

from ._regional_gls_masked import (
    _exact_design_fit,
    _fit_masked,
    _profile_masked,
    _spectral_evaluate,
    _spectral_prepare,
)


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


@njit(inline="always")
def _trim_and_append_subregion(
    pos, effects, sub_start, sub_end, direction, tau, min_cpgs, min_length, regions
):
    n_sub = sub_end - sub_start
    if n_sub < min_cpgs:
        return

    # Each right endpoint admits a growing set of starts meeting BOTH count
    # and genomic-span constraints. The smallest admissible prefix gives the
    # globally best segment ending there; no unconstrained Kadane resets.
    prefix = np.zeros(n_sub + 1)
    for j in range(n_sub):
        prefix[j + 1] = prefix[j] + direction * effects[sub_start + j] - tau
    next_start, min_start = 0, 0
    min_prefix, best_score = np.inf, -np.inf
    best_l, best_r = -1, -1
    for end in range(min_cpgs, n_sub + 1):
        while (
            next_start <= end - min_cpgs
            and pos[sub_start + end - 1] - pos[sub_start + next_start] + 1 >= min_length
        ):
            if prefix[next_start] < min_prefix:
                min_prefix, min_start = prefix[next_start], next_start
            next_start += 1
        score = prefix[end] - min_prefix
        if score > best_score:
            best_score, best_l, best_r = score, min_start, end
    if best_l >= 0:
        regions.append((int(sub_start + best_l), int(sub_start + best_r)))


@njit(inline="always")
def _refine_and_append_region(
    pos,
    effects,
    first,
    stop_excl,
    direction,
    min_effect,
    min_cpgs,
    min_length,
    regions,
    boundary_min_effect=0.0,
):
    """Trim neutral boxcar halos and split valley bridges via penalized CUSUM."""
    tau = boundary_min_effect if boundary_min_effect > 0.0 else 0.55 * min_effect
    edge_floor = boundary_min_effect if boundary_min_effect > 0.0 else 0.35 * min_effect
    while first < stop_excl and direction * effects[first] < edge_floor:
        first += 1
    while stop_excl > first and direction * effects[stop_excl - 1] < edge_floor:
        stop_excl -= 1

    if stop_excl - first < min_cpgs:
        return

    sub_start = first
    run_low = 0
    low_start = first
    for j in range(first, stop_excl):
        if direction * effects[j] < tau:
            if run_low == 0:
                low_start = j
            run_low += 1
        else:
            if run_low >= 3 and (
                (pos[j - 1] - pos[low_start]) >= 100
                or (low_start > 0 and pos[j] - pos[low_start - 1] >= 150)
            ):
                _trim_and_append_subregion(
                    pos, effects, sub_start, low_start, direction, tau, min_cpgs, min_length, regions
                )
                sub_start = j
            run_low = 0

    if stop_excl - sub_start >= min_cpgs:
        _trim_and_append_subregion(
            pos, effects, sub_start, stop_excl, direction, tau, min_cpgs, min_length, regions
        )


@njit(cache=True, nogil=True)
def discover_prepared(
    pos,
    effects,
    left,
    right,
    starts,
    stops,
    min_effect,
    min_cpgs,
    max_gap,
    min_length,
    boundary_min_effect=0.0,
    smoothed=None,
    boundary_effects=None,
):
    cumulative = np.empty(len(pos))
    trim_effects = effects if boundary_effects is None else boundary_effects
    regions = [(np.int64(0), np.int64(0))][:0]
    threshold = boundary_min_effect if boundary_min_effect > 0.0 else min_effect

    for segment in range(len(starts)):
        lo, hi = starts[segment], stops[segment]
        total = 0.0
        for j in range(lo, hi):
            total += effects[j]
            cumulative[j] = total
        first, last, direction = -1, -1, 0
        has_seed = False
        for j in range(lo, hi):
            if smoothed is None:
                value = cumulative[right[j] - 1]
                if left[j] > lo:
                    value -= cumulative[left[j] - 1]
                value /= right[j] - left[j]
            else:
                value = smoothed[j]
            if abs(value) < threshold:
                if first >= 0:
                    if (
                        (has_seed or boundary_min_effect == 0.0)
                        and last - first + 1 >= min_cpgs
                        and pos[last] - pos[first] + 1 >= min_length
                    ):
                        _refine_and_append_region(
                            pos,
                            trim_effects,
                            first,
                            last + 1,
                            direction,
                            min_effect,
                            min_cpgs,
                            min_length,
                            regions,
                            boundary_min_effect,
                        )
                    first = -1
                    has_seed = False
                continue
            sign = 1 if value > 0 else -1
            if first >= 0 and (pos[j] - pos[last] > max_gap or sign != direction):
                if (
                    (has_seed or boundary_min_effect == 0.0)
                    and last - first + 1 >= min_cpgs
                    and pos[last] - pos[first] + 1 >= min_length
                ):
                    _refine_and_append_region(
                        pos,
                        trim_effects,
                        first,
                        last + 1,
                        direction,
                        min_effect,
                        min_cpgs,
                        min_length,
                        regions,
                        boundary_min_effect,
                    )
                first = -1
                has_seed = False
            if first < 0:
                first, direction = j, sign
                has_seed = False
            last = j
            if abs(value) >= min_effect:
                has_seed = True
        if (
            first >= 0
            and (has_seed or boundary_min_effect == 0.0)
            and last - first + 1 >= min_cpgs
            and pos[last] - pos[first] + 1 >= min_length
        ):
            _refine_and_append_region(
                pos,
                trim_effects,
                first,
                last + 1,
                direction,
                min_effect,
                min_cpgs,
                min_length,
                regions,
                boundary_min_effect,
            )

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
            y_val = yc[s, j]
            diagonal[j] += y_val * y_val
            contrast[j] += x * y_val
            if j + 1 < k:
                adjacent[j] += y_val * yc[s, j + 1]
    return diagonal, adjacent, contrast, max(gram, 1e-12), np.median(h * h)


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
    if _exact_design_fit(yc, h, g, np.ones(yc.shape, dtype=np.bool_)):
        return 0.0, 0.0, np.nan, np.nan, np.inf
    prepared = _prepare(yc, h, g)
    best, c, effect, se = _profile(p, yc, h, 0.0, prepared)
    best_length = 0.0
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
def _batch(p, yc, h, g, bounds, obs_mask):
    result = np.empty((len(bounds), 5))
    for j in prange(len(bounds)):
        lo, hi = bounds[j]
        p_sub = p[lo:hi]
        yc_sub = yc[:, lo:hi]
        h_sub = h[lo:hi]
        mask_sub = obs_mask[:, lo:hi]
        if np.all(mask_sub):
            length, c, effect, se, best = _fit(p_sub, yc_sub, h_sub, g)
        else:
            length, c, effect, se, best = _fit_masked(p_sub, yc_sub, h_sub, g, mask_sub)
        result[j, 0] = length
        result[j, 1] = c
        result[j, 2] = effect
        result[j, 3] = se
        result[j, 4] = best
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


def _observation_mask(yc, obs_mask):
    if obs_mask is None:
        return np.ones(yc.shape, dtype=np.bool_)
    mask = np.asarray(obs_mask, dtype=np.bool_)
    if mask.shape != yc.shape or np.any(mask.sum(axis=0) == 0):
        raise ValueError("Observation mask must cover each locus and match the response shape")
    return mask


def evaluate_covariance(p, yc, h, g, length, c, *, obs_mask=None):
    p, yc, h, g = _validate(p, yc, h, g)
    if _exact_design_fit(yc, h, g, _observation_mask(yc, obs_mask)):
        raise ValueError("Degenerate regional residual covariance; exact design fit")
    if length < 0 or not 0 <= c < 1:
        raise ValueError("Invalid covariance parameters")
    if obs_mask is not None:
        mask = _observation_mask(yc, obs_mask)
        return _spectral_evaluate(_spectral_prepare(p, yc, h, g, mask, length), c)
    return _evaluate(p, yc, h, length, c, _prepare(yc, h, g))


def profile_covariance(p, yc, h, g, length, *, obs_mask=None):
    p, yc, h, g = _validate(p, yc, h, g)
    if _exact_design_fit(yc, h, g, _observation_mask(yc, obs_mask)):
        raise ValueError("Degenerate regional residual covariance; exact design fit")
    if obs_mask is not None:
        mask = _observation_mask(yc, obs_mask)
        return _profile_masked(p, yc, h, g, mask, length)
    return _profile(p, yc, h, length, _prepare(yc, h, g))


def fit_regions(p, yc, h, g, bounds, *, workers=None, obs_mask=None):
    p, yc, h, g = _validate(p, yc, h, g)
    mask = _observation_mask(yc, obs_mask)
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
        result = _batch(p, yc, h, g, bounds, mask)
    finally:
        set_num_threads(previous)
    if not np.isfinite(result).all() or np.any(result[:, 3] <= 0):
        raise ValueError(
            "Degenerate regional covariance; failed fits cannot be omitted from calibration"
        )
    return result
