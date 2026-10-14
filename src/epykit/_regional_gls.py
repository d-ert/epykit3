"""Experimental region discovery and spatial GLS, retaining biological samples.

Architecture informed by dmrseq (Korthauer et al., 2019), independently
implemented. GLS t statistics are ranking scores, not parametric p-values.
Pooled full-discovery label permutations require exchangeable samples and
comparable regional null scores; their calibration must be checked empirically.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
from itertools import combinations, pairwise
from math import comb
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)


def transform_counts(methylated, coverage):
    m, n = np.asarray(methylated, float), np.asarray(coverage, float)
    if (
        m.shape != n.shape
        or m.ndim != 2
        or not np.isfinite(m).all()
        or not np.isfinite(n).all()
        or np.any(n <= 0)
        or np.any(m < 0)
        or np.any(m > n)
    ):
        raise ValueError(
            "Regional inference requires finite counts and identical positive-coverage support"
        )
    return np.arcsin(2 * (m + 0.05) / (n + 0.1) - 1)


def covariance_precision(positions, depths, length_scale, donor_fraction):
    """O(K) inverse of CAR plus a depth-aware biological sample effect.

    In weighted coordinates R=(1-c)A+c uu', u=sqrt(depth/median(depth)).
    The second term corresponds to a constant sample offset in unweighted
    transformed methylation, regardless of unequal CpG depth.
    """
    pos, depth = np.asarray(positions, float), np.asarray(depths, float)
    if (
        len(pos) < 2
        or depth.shape != pos.shape
        or np.any(np.diff(pos) <= 0)
        or np.any(depth <= 0)
        or not np.isfinite(depth).all()
        or not np.isfinite(pos).all()
        or length_scale < 0
        or not 0 <= donor_fraction < 1
    ):
        raise ValueError("Invalid regional covariance")
    a = np.exp(-np.diff(pos) / length_scale) if length_scale else np.zeros(len(pos) - 1)
    denominator = -np.expm1(-2 * np.diff(pos) / length_scale) if length_scale else np.ones(len(a))
    diag = np.ones(len(pos))
    diag[:-1] += a * a / denominator
    diag[1:] += a * a / denominator
    off = -a / denominator

    def car_apply(v):
        v = np.asarray(v)
        out = diag * v
        out[..., :-1] += off * v[..., 1:]
        out[..., 1:] += off * v[..., :-1]
        return out

    u = np.sqrt(depth / np.median(depth))
    qu = car_apply(u)
    c = donor_fraction
    factor = c / ((1 - c) ** 2 * (1 + c / (1 - c) * (u @ qu)))

    def apply(v):
        v = np.asarray(v)
        return car_apply(v) / (1 - c) - factor * np.asarray(v @ qu)[..., None] * qu

    logdet = len(pos) * np.log1p(-c) + np.log(denominator).sum() + np.log1p(c / (1 - c) * (u @ qu))
    return apply, float(logdet)


def fit_region(
    positions,
    methylated,
    coverage,
    labels,
    *,
    covariance=None,
    precomputed_yc=None,
    precomputed_depth=None,
    precomputed_h=None,
    obs_mask=None,
):
    """Fit locus intercepts and one treatment contrast across biological samples.

    Covariance is profiled by REML. A fixed (length_scale, donor_fraction)
    option supports independent dense-GLS numerical verification.
    """
    pos, g = np.asarray(positions, float), np.asarray(labels, float)
    if set(g) != {0.0, 1.0} or min((g == 0).sum(), (g == 1).sum()) < 2:
        raise ValueError("At least two independent samples per group are required")

    if len(pos) < 2 or np.any(np.diff(pos) <= 0) or not np.isfinite(pos).all():
        raise ValueError("Invalid regional covariance")
    if precomputed_yc is not None and precomputed_depth is not None and precomputed_h is not None:
        yc, h = precomputed_yc, precomputed_h
    else:
        m, n = np.asarray(methylated, float), np.asarray(coverage, float)
        if m.shape != n.shape or m.shape != (len(pos), len(g)):
            raise ValueError("Shape mismatch between counts and labels")
        if not np.isfinite(m).all() or not np.isfinite(n).all() or np.any(n < 0) or np.any(m < 0) or np.any(m > n):
            raise ValueError("Regional counts must be finite and satisfy 0 <= methylated <= coverage")
        valid = n > 0
        if np.any(valid.sum(axis=1) == 0):
            raise ValueError("Each regional locus must have at least one observed sample")
        if obs_mask is not None and not np.array_equal(np.asarray(obs_mask, bool), valid.T):
            raise ValueError("Count observation mask must match positive coverage")
        obs_mask = valid.T if not valid.all() else None
        y = np.arcsin(2*(m+.05)/(np.maximum(n, 1)+.1)-1).T
        locus_mean = np.sum(y.T*valid, axis=1) / valid.sum(axis=1)
        h = np.sqrt(np.nanmedian(np.where(valid, n, np.nan), axis=1))
        yc = np.where(valid.T, y-locus_mean, 0.0) * h
    from ._regional_gls_native import _fit, _observation_mask, _validate, evaluate_covariance

    pos, yc, h, g = _validate(pos, yc, h, g)
    mask = _observation_mask(yc, obs_mask)
    if covariance is None:
        if mask.all():
            length, fraction, effect, se, objective = _fit(pos, yc, h, g)
        else:
            from ._regional_gls_masked import _fit_masked

            length, fraction, effect, se, objective = _fit_masked(pos, yc, h, g, mask)
    else:
        length, fraction = covariance
        objective, effect, se = evaluate_covariance(pos, yc, h, g, length, fraction, obs_mask=obs_mask)
    if not np.isfinite(objective) or not np.isfinite(se) or se <= 0:
        raise ValueError("Degenerate regional residual covariance; inference cannot be calibrated")
    return dict(
        effect_transformed=float(effect),
        statistic=float(effect / se),
        residual_df=int(mask.sum() - len(pos) - 1),
        length_scale=float(length),
        donor_fraction=float(fraction),
        standard_error=float(se),
    )


def discover_regions(
    positions,
    effects,
    *,
    window_bp=500,
    min_cpgs=5,
    max_gap=500,
    min_length=50,
    min_effect=0.1,
    boundary_min_effect=None,
):
    """Effect-based discovery without a single-CpG significance prerequisite."""
    pos, effect = np.asarray(positions, np.int64), np.asarray(effects, float)
    if (
        pos.ndim != 1
        or effect.shape != pos.shape
        or np.any(np.diff(pos) <= 0)
        or np.any(pos < 0)
        or not np.isfinite(effect).all()
        or window_bp < 1
        or min_cpgs < 2
        or max_gap < 1
        or min_length < 1
        or not 0 < min_effect <= 1
        or (
            boundary_min_effect is not None
            and not (0 < boundary_min_effect <= min_effect)
        )
    ):
        raise ValueError("Invalid region discovery inputs")
    from ._regional_gls_native import discover_prepared, prepare_discovery

    geom = prepare_discovery(pos, window_bp, min_cpgs, max_gap)
    b_eff = 0.0 if boundary_min_effect is None else float(boundary_min_effect)
    res = discover_prepared(pos, effect, *geom, min_effect, min_cpgs, max_gap, min_length, b_eff)
    return [tuple(r) for r in res.tolist()]


def _precompute_discovery_geometry(pos: np.ndarray, window_bp: int, min_cpgs: int, max_gap: int):
    """Precompute chromosome-wide spatial cuts and window indices for fast discovery."""
    cuts = np.r_[0, np.flatnonzero(np.diff(pos) > max_gap) + 1, len(pos)]
    geom = []
    for lo, hi in pairwise(cuts):
        if hi - lo < min_cpgs:
            continue
        p = pos[lo:hi]
        left = np.searchsorted(p, p - window_bp // 2)
        right = np.searchsorted(p, p + (window_bp + 1) // 2, side="right")
        denom = right - left
        geom.append((lo, hi, p, left, right, denom))
    return geom


def _fast_discover_with_geom(
    geom,
    effect_arr: np.ndarray,
    min_effect: float,
    min_cpgs: int,
    max_gap: int,
    min_length: int,
    boundary_min_effect: float | None = None,
):
    """Fast candidate region discovery using precomputed chromosome geometry."""
    regions = []
    for lo, hi, p, left, right, denom in geom:
        sums = np.r_[0.0, np.cumsum(effect_arr[lo:hi])]
        smooth = (sums[right] - sums[left]) / denom
        if boundary_min_effect is not None:
            first, last, direction = -1, -1, 0
            has_seed = False
            for j in range(hi - lo):
                val = smooth[j]
                if abs(val) < boundary_min_effect:
                    if first >= 0:
                        if (
                            has_seed
                            and last - first + 1 >= min_cpgs
                            and p[last] - p[first] + 1 >= min_length
                        ):
                            regions.append((int(lo + first), int(lo + last + 1)))
                        first = -1
                        has_seed = False
                    continue
                sign = 1 if val > 0 else -1
                if first >= 0 and sign != direction:
                    if (
                        has_seed
                        and last - first + 1 >= min_cpgs
                        and p[last] - p[first] + 1 >= min_length
                    ):
                        regions.append((int(lo + first), int(lo + last + 1)))
                    first = -1
                    has_seed = False
                if first < 0:
                    first, direction = j, sign
                    has_seed = False
                last = j
                if abs(val) >= min_effect:
                    has_seed = True
            if (
                first >= 0
                and has_seed
                and last - first + 1 >= min_cpgs
                and p[last] - p[first] + 1 >= min_length
            ):
                regions.append((int(lo + first), int(lo + last + 1)))
        else:
            seeds = np.flatnonzero(np.abs(smooth) >= min_effect)
            if not len(seeds):
                continue
            splits = np.r_[
                0,
                np.flatnonzero((np.diff(p[seeds]) > max_gap) | (np.diff(np.sign(smooth[seeds])) != 0))
                + 1,
                len(seeds),
            ]
            for a, b in pairwise(splits):
                start, end = int(seeds[a]), int(seeds[b - 1]) + 1
                if end - start >= min_cpgs and p[end - 1] - p[start] + 1 >= min_length:
                    regions.append((int(lo + start), int(lo + end)))
    return regions


def _fast_decoy_statistic(
    p_slice: np.ndarray,
    yc_slice: np.ndarray,
    depth_slice: np.ndarray,
    h_slice: np.ndarray,
    g_flt: np.ndarray,
    length_scale: float,
    donor_fraction: float,
) -> float:
    """Evaluate GLS test statistic with fixed covariance in O(K) time without REML optimization."""
    x = g_flt - g_flt.mean()
    gram = float(x @ x)
    k_reg = len(p_slice)
    df = (len(g_flt) - 1) * k_reg - 1

    c = donor_fraction
    one_minus_c = 1.0 - c
    dpos = np.diff(p_slice)

    if length_scale == 0.0:
        a = np.zeros(k_reg - 1)
        den = np.ones(k_reg - 1)
    else:
        inv_len = 1.0 / length_scale
        a = np.exp(-dpos * inv_len)
        den = -np.expm1(-2.0 * dpos * inv_len)

    a_sq_den = (a * a) / den
    off = -a / den
    diag = np.ones(k_reg)
    diag[:-1] += a_sq_den
    diag[1:] += a_sq_den

    med_d = float(np.median(depth_slice))
    if med_d <= 0 or not np.isfinite(med_d):
        return np.nan
    u = np.sqrt(depth_slice / med_d)
    qu = diag * u
    qu[:-1] += off * u[1:]
    qu[1:] += off * u[:-1]

    u_qu = float(u @ qu)
    factor = c / ((one_minus_c * one_minus_c) * (1.0 + (c / one_minus_c) * u_qu))

    qh = diag * h_slice
    qh[:-1] += off * h_slice[1:]
    qh[1:] += off * h_slice[:-1]
    qh = qh / one_minus_c - factor * float(h_slice @ qu) * qu

    information = float(h_slice @ qh)
    if information <= 0 or not np.isfinite(information):
        return np.nan
    effect = float((x @ yc_slice @ qh) / (gram * information))

    qyc = diag * yc_slice
    qyc[:, :-1] += off * yc_slice[:, 1:]
    qyc[:, 1:] += off * yc_slice[:, :-1]
    qyc = qyc / one_minus_c - factor * (yc_slice @ qu)[:, None] * qu

    rss = float(np.sum(yc_slice * qyc) - gram * information * effect**2)
    if rss <= 0.0 or not np.isfinite(rss):
        return np.nan
    se = np.sqrt(rss / df / (gram * information))
    stat = float(effect / se)
    return stat if np.isfinite(stat) else np.nan


def label_assignments(labels, n_perm, seed):
    labels = np.asarray(labels, dtype=int)
    n, k = len(labels), int(labels.sum())
    if set(labels) != {0, 1} or min(k, n - k) < 2 or n_perm < 2:
        raise ValueError("Invalid independent-sample permutation design")
    total = comb(n, k)
    if n_perm >= total:
        indices = list(combinations(range(n), k))
    else:
        rng = np.random.default_rng(seed)
        indices = {tuple(np.flatnonzero(labels))}
        while len(indices) < n_perm:
            indices.add(tuple(sorted(rng.choice(n, k, replace=False))))
        indices = sorted(indices)
    result = []
    for index in indices:
        g = np.zeros(n, dtype=int)
        g[list(index)] = 1
        result.append(g)
    return result


def _permutation_fdp_qvalues(scores, null_scans):
    """Monotone expected-false-count FDP estimates over complete scan exposure.

    This is not an exact randomization test or a finite-sample FDR guarantee.
    A zero tail is a point estimate of zero false calls, not proof of zero risk.
    Score ties form one threshold and all scans, including empty scans, count.
    """
    scores = np.abs(np.asarray(scores, float))
    if not np.isfinite(scores).all():
        raise ValueError("Failed observed scores must not be dropped")
    if not null_scans:
        raise ValueError("At least one complete null scan is required")
    null = np.sort(np.abs(np.concatenate([np.asarray(s, float) for s in null_scans])))
    if not np.isfinite(null).all():
        raise ValueError("Failed permutation scores must not be dropped")
    if not len(scores):
        return np.empty(0)
    thresholds = np.unique(scores)[::-1]
    tails = len(null) - np.searchsorted(null, thresholds, side="left")
    discoveries = len(scores) - np.searchsorted(np.sort(scores), thresholds, side="left")
    fdp = np.minimum(1.0, tails / (len(null_scans) * discoveries))
    curve = np.minimum.accumulate(fdp[::-1])[::-1]
    # Thresholds are descending; search on the ascending view preserves ties.
    return curve[::-1][np.searchsorted(thresholds[::-1], scores)]


def calibrate_regions(observed, null_scans, *, correction="pooled_bh"):
    """Calibrate selected scores; pooled-BH retained for historical migrations."""
    null = (
        np.sort(np.abs(np.concatenate([np.asarray(s, float) for s in null_scans])))
        if null_scans
        else np.array([])
    )
    if not np.isfinite(null).all():
        raise ValueError("Failed permutation scores must not be dropped")
    scores = np.array([abs(r["statistic"]) for r in observed])
    p = (1 + len(null) - np.searchsorted(null, scores, side="left")) / (1 + len(null))
    order = np.argsort(p)
    q = np.empty(len(p))
    q[order] = np.minimum(
        1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    )
    if correction == "permutation_fdp":
        q = _permutation_fdp_qvalues(scores, null_scans)
    elif correction != "pooled_bh":
        raise ValueError("Unknown regional score correction")
    for r, pv, qv in zip(observed, p, q, strict=False):
        r.update(combined_pvalue=float(pv), combined_qvalue=float(qv))
    return observed, dict(
        null_scans=len(null_scans),
        null_regions=len(null),
        empty_null_scans=sum(len(s) == 0 for s in null_scans),
        correction=("full-discovery-permutation-expected-false-count-experimental"
                    if correction == "permutation_fdp" else "pooled-full-discovery-permutation-bh-experimental"),
        score_assumption="comparable selected regional null statistics; exchangeable samples",
        pvalue_semantics="descriptive pooled selected-region score tail, not an exact per-region randomization p-value",
        qvalue_semantics=("threshold-wise expected decoy discoveries per complete scan divided by observed discoveries"
                          if correction == "permutation_fdp" else "BH transform of descriptive pooled score tails"),
        finite_sample_fdr_guarantee=False,
    )


def remove_self_mirror_scores(observed, null_region_count, duplicates):
    """Recover conditional tails from a saved pool containing exact score copies.

    Used only to migrate development-v1 benchmarks. It avoids refitting unchanged
    GLS models; independently checked against directly excluding duplicate scans.
    """
    import copy

    result = copy.deepcopy(observed)
    scores = np.abs([r["statistic"] for r in result])
    tails = len(scores) - np.searchsorted(np.sort(scores), scores, side="left")
    old = np.array([r["combined_pvalue"] for r in result]) * (1 + null_region_count) - 1
    if not np.allclose(old, np.rint(old), atol=1e-6):
        raise ValueError("Saved scores do not encode integer empirical tail counts")
    remaining = null_region_count - duplicates * len(result)
    tail = np.rint(old) - duplicates * tails
    if remaining < 0 or np.any(tail < 0) or np.any(tail > remaining):
        raise ValueError("Invalid duplicate-score pool reconstruction")
    p = (1 + tail) / (1 + remaining)
    order = np.argsort(p)
    q = np.empty(len(p))
    q[order] = np.minimum(
        1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    )
    for r, pv, qv in zip(result, p, q, strict=False):
        r.update(combined_pvalue=float(pv), combined_qvalue=float(qv))
    return result


def _load_regional_cache(paths, key, q_cutoff, workers):
    table_path, meta_path = paths
    if not table_path.is_file() or not meta_path.is_file():
        return None
    try:
        meta = json.loads(meta_path.read_text())
        if (
            meta.get("cache_key") != key
            or meta.get("table_sha256") != hashlib.sha256(table_path.read_bytes()).hexdigest()
        ):
            return None
        all_regions = pl.read_parquet(table_path)
        result = (
            all_regions
            if q_cutoff is None
            else all_regions.filter(pl.col("combined_qvalue") <= q_cutoff)
        )
        meta.update(cache_hit=True, min_mean_qvalue=q_cutoff, workers=workers)
        logger.info(
            "regional_gls: reused complete cached inference (%d candidates)", len(all_regions)
        )
        return result, all_regions, meta
    except (OSError, json.JSONDecodeError, pl.exceptions.PolarsError):
        logger.warning("Invalid regional cache; recomputing the complete inference")
        return None


def _scan_group_effects(beta_clean, valid_flt, labels, min_group_obs):
    """Raw/reporting and support-shrunk/discovery effects for one label split."""
    case = np.asarray(labels, dtype=float)
    control = 1.0 - case
    n_case, n_control = valid_flt @ case, valid_flt @ control
    mean_case = (beta_clean @ case) / np.maximum(n_case, 1.0)
    mean_control = (beta_clean @ control) / np.maximum(n_control, 1.0)
    raw = mean_case - mean_control
    harmonic = 2.0 * n_case * n_control / np.maximum(n_case + n_control, 1.0)
    full_harmonic = 2.0 * case.sum() * control.sum() / len(case)
    discovery = raw * np.sqrt(np.minimum(harmonic / full_harmonic, 1.0))
    eligible = (n_case >= min_group_obs) & (n_control >= min_group_obs)
    return raw, discovery, eligible


def _prepare_label_scan(pos, yc, h, mask, beta, valid, labels, min_group_obs,
                        window_bp, min_cpgs, max_gap, min_length, min_effect, boundary):
    from ._regional_gls_native import discover_prepared, prepare_discovery

    raw, effects, eligible = _scan_group_effects(beta, valid, labels, min_group_obs)
    scan_pos, scan_h = pos[eligible], h[eligible]
    scan_yc, scan_mask = yc[:, eligible], mask[:, eligible]
    raw, effects = raw[eligible], effects[eligible]
    geometry = prepare_discovery(scan_pos, window_bp, min_cpgs, max_gap)
    bounds = discover_prepared(
        scan_pos, effects, *geometry, min_effect, min_cpgs, max_gap, min_length, boundary
    )
    return scan_pos, scan_yc, scan_h, scan_mask, raw, bounds


def run_regional_gls(
    md,
    *,
    chromosomes=None,
    canonical_only=False,
    window_bp=500,
    min_cpgs=5,
    max_gap=500,
    min_length=50,
    min_effect=0.1,
    boundary_min_effect=None,
    n_perm=100,
    seed=42,
    q_cutoff=0.05,
    pooled_decoy_covariance=False,
    calibration="permutation_fdp",
    record_scan_scores=False,
    workers=None,
    cache=True,
):
    """Stream a label-independent union; replay support, discovery and fitting."""
    if (
        window_bp < 1
        or min_cpgs < 2
        or max_gap < 1
        or min_length < 1
        or not 0 < min_effect <= 1
        or (
            boundary_min_effect is not None
            and not (0 < boundary_min_effect <= min_effect)
        )
    ):
        raise ValueError("Invalid regional discovery geometry or effect threshold")
    if q_cutoff is not None and not 0 < q_cutoff <= 1:
        raise ValueError("q cutoff must be in (0,1]")
    if calibration not in {"permutation_fdp", "pooled_bh"}:
        raise ValueError("Unknown regional calibration")
    if record_scan_scores and pooled_decoy_covariance:
        raise ValueError("Complete scan score diagnostics require per-region covariance replay")
    samples = list(md.control_ids) + list(md.treatment_ids)
    if len(samples) != len(set(samples)):
        raise ValueError("Regional GLS requires unique biological sample IDs")
    labels = np.r_[np.zeros(len(md.control_ids), int), np.ones(len(md.treatment_ids), int)]
    assignments = label_assignments(labels, n_perm, seed)
    generated_assignments = len(assignments)
    # Pooled regional calibration is conditional decoy inference, not an exact
    # per-locus randomization p-value. Self/mirror scans duplicate the observed
    # scores and impose q >= 2*m_observed/n_decoys in balanced designs. Retain
    # every other selected assignment, including scans producing no candidates.
    assignments = [
        g
        for g in assignments
        if not np.array_equal(g, labels) and not np.array_equal(g, 1 - labels)
    ]
    if not assignments:
        raise ValueError("No non-self/mirror decoy assignments; increase n_perm")
    # Complementary balanced assignments discover identical regions and have
    # exactly opposite signed GLS statistics. Fit one, retaining both scans in
    # the empirical pool, including their empty outcomes.
    representatives, scan_map, represented = [], [], {}
    for g in assignments:
        key, mirror = tuple(g), tuple(1 - g)
        if mirror in represented:
            scan_map.append((represented[mirror], -1.0))
        else:
            represented[key] = len(representatives)
            scan_map.append((len(representatives), 1.0))
            representatives.append(g)
    root = Path(md.store)
    if not root.is_dir() or any(not (root / f"sample={s}").is_dir() for s in samples):
        raise ValueError(
            "Regional count store must contain each sample partition; use md.store from ingestion"
        )
    available = {
        p.name.removeprefix("chrom=")
        for sample in samples for p in (root / f"sample={sample}").glob("chrom=*")
    }
    if not available:
        raise ValueError("Regional count store contains no chromosome partitions")
    chroms = sorted(available if chromosomes is None else set(chromosomes))
    if set(chroms) - available:
        raise ValueError("Requested chromosome is absent from count store")
    if canonical_only:
        canonical = {f"chr{x}" for x in [*range(1, 23), "X", "Y", "M", "MT"]}
        chroms = [c for c in chroms if c in canonical]
    cache_paths = None
    code_digest = hashlib.sha256(
        Path(__file__).read_bytes()
        + Path(__file__).with_name("_regional_gls_native.py").read_bytes()
        + Path(__file__).with_name("_regional_gls_masked.py").read_bytes()
    ).hexdigest()
    if cache and md.analysis_root and not pooled_decoy_covariance:
        inputs = []
        for sample in samples:
            for chrom in chroms:
                for part in sorted(
                    (root / f"sample={sample}" / f"chrom={chrom}").glob("*.parquet")
                ):
                    stat = part.stat()
                    inputs.append((str(part.resolve()), stat.st_size, stat.st_mtime_ns))
        signature = dict(
            source=code_digest,
            samples=samples,
            labels=labels.tolist(),
            inputs=inputs,
            chromosomes=chroms,
            window_bp=window_bp,
            min_cpgs=min_cpgs,
            max_gap=max_gap,
            min_length=min_length,
            min_effect=min_effect,
            boundary_min_effect=boundary_min_effect,
            n_perm=n_perm,
            seed=seed,
            calibration=calibration,
            record_scan_scores=record_scan_scores,
            numpy=np.__version__,
        )
        key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
        directory = Path(md.analysis_root) / ".cache" / "regional_gls"
        cache_paths = directory / f"{key}.parquet", directory / f"{key}.json"
        cached = _load_regional_cache(cache_paths, key, q_cutoff, workers)
        if cached is not None:
            return cached
    observed, null_scans = [], [[] for _ in assignments]
    representative_scans = [[] for _ in representatives]
    common_support, union_support = {}, {}

    for chrom in chroms:
        # Optimized Parquet scan and multi-sample common support join
        frames = []
        for sample in samples:
            sample_dir = root / f"sample={sample}" / f"chrom={chrom}"
            parts = sorted(sample_dir.glob("*.parquet"))
            if not parts:
                frames.append(pl.DataFrame(schema={"pos": pl.Int64, "N_meth": pl.Int64, "coverage": pl.Int64}))
                continue
            pattern = str(sample_dir / "*.parquet")
            frame = (
                pl.scan_parquet(pattern)
                .select(["pos", "N_meth", "coverage"])
                .filter(pl.col("coverage") > 0)
                .with_columns(pl.all().cast(pl.Int64))
                .collect()
            )
            if frame["pos"].n_unique() != len(frame):
                raise ValueError(f"Duplicate CpG coordinates in {sample}/{chrom}")
            frames.append(frame)

        # Only label-independent support is applied to the shared universe.
        renamed = [
            f.rename({"N_meth": f"m_{i}", "coverage": f"n_{i}"})
            for i, f in enumerate(frames)
        ]
        res = renamed[0]
        for f in renamed[1:]:
            res = res.join(f, on="pos", how="full", coalesce=True)

        m_cols = [f"m_{i}" for i in range(len(samples))]
        n_cols = [f"n_{i}" for i in range(len(samples))]
        res = res.with_columns(
            [pl.col(c).fill_null(0) for c in m_cols + n_cols]
        )

        # Require >= 2 valid replicates per group (or min(2, group_size)) and >= 50% total sample support
        min_group_obs = min(2, len(md.control_ids), len(md.treatment_ids))
        min_total_obs = max(min_group_obs * 2, int(np.ceil(0.5 * len(samples))))

        res = res.filter(
            pl.sum_horizontal([pl.col(c) > 0 for c in n_cols]) >= min_total_obs
        ).sort("pos")

        union_support[chrom] = len(res)
        if len(res) < min_cpgs:
            common_support[chrom] = len(res)
            continue

        pos = res["pos"].to_numpy()

        m = res.select(m_cols).to_numpy().astype(np.float64)
        n = res.select(n_cols).to_numpy().astype(np.float64)
        valid_mask = n > 0

        # Safe arcsin transformation on positive coverage
        n_safe = np.where(valid_mask, n, 1.0)
        m_safe = np.where(valid_mask, m, 0.0)
        y_raw = np.arcsin(2.0 * (m_safe + 0.05) / (n_safe + 0.1) - 1.0)  # (K, S)

        # Precomputed matrices for fast group effect computation with partial-support shrinkage
        beta_clean = np.where(valid_mask, m_safe / n_safe, 0.0)  # (K, S)
        valid_flt = valid_mask.astype(np.float64)  # (K, S)

        # Center y_raw at each CpG using ONLY the valid observed samples at that CpG.
        # Missing entries are set to 0.0 in yc_full so they are 100% exchangeable across
        # all label permutations g and never leak observed group means into decoy fits!
        obs_counts = np.maximum(np.sum(valid_mask, axis=1), 1.0)
        y_locus_mean = np.sum(y_raw * valid_mask, axis=1) / obs_counts
        y_centered = np.where(valid_mask, y_raw - y_locus_mean[:, None], 0.0).T  # (S, K)

        # Missing sample support is already accounted for in the observed-data
        # information matrix. Depth weights use only positive observed depths.
        n_sorted = np.where(valid_mask, n, np.nan)
        depth_full = np.nanmedian(n_sorted, axis=1)
        h_full = np.sqrt(np.maximum(depth_full, 1.0))
        yc_full = y_centered * h_full

        b_eff = 0.0 if boundary_min_effect is None else float(boundary_min_effect)
        scan_args = (pos, yc_full, h_full, valid_mask.T, beta_clean, valid_flt)
        discovery_args = (min_group_obs, window_bp, min_cpgs, max_gap, min_length, min_effect, b_eff)
        obs_pos, obs_yc, obs_h, obs_mask, obs_effects, obs_regions = _prepare_label_scan(
            *scan_args, labels, *discovery_args
        )
        common_support[chrom] = len(obs_pos)
        current = []
        from ._regional_gls_native import fit_regions

        native = fit_regions(obs_pos, obs_yc, obs_h, labels, obs_regions, workers=workers, obs_mask=obs_mask)
        for (lo, hi), (length, fraction, effect, se, _objective) in zip(
            obs_regions, native, strict=True
        ):
            fit = dict(
                effect_transformed=float(effect),
                statistic=float(effect / se),
                residual_df=int(obs_mask[:, lo:hi].sum() - (hi - lo) - 1),
                length_scale=float(length),
                donor_fraction=float(fraction),
                standard_error=float(se),
            )
            fit.update(
                chrom=chrom,
                start=int(obs_pos[lo]),
                end=int(obs_pos[hi - 1]) + 1,
                n_cpgs=hi - lo,
                mean_meth_diff=float(obs_effects[lo:hi].mean()),
                n_case=int(labels.sum()),
                n_control=int(len(labels) - labels.sum()),
            )
            current.append(fit)

        observed.extend(current)
        logging.getLogger(__name__).info(
            "regional_gls %s: %d common CpGs, %d observed candidates", chrom, len(pos), len(current)
        )

        if pooled_decoy_covariance:
            logger.warning(
                "Experimental pooled covariance changes the decoy fitting procedure; nominal calibration is unvalidated"
            )
            valid_lens = [f["length_scale"] for f in current if f["length_scale"] > 0]
            valid_fracs = [f["donor_fraction"] for f in current]
            chrom_len = float(np.median(valid_lens)) if valid_lens else 50.0
            chrom_donor = float(np.median(valid_fracs)) if valid_fracs else 0.05

            for j, g in enumerate(assignments):
                scan_p, scan_y, scan_h, scan_mask, _, regions = _prepare_label_scan(
                    *scan_args, g, *discovery_args
                )
                for lo, hi in regions:
                    from ._regional_gls_native import evaluate_covariance

                    _, effect, se = evaluate_covariance(
                        scan_p[lo:hi], scan_y[:, lo:hi], scan_h[lo:hi], g,
                        chrom_len, chrom_donor, obs_mask=scan_mask[:, lo:hi]
                    )
                    if not np.isfinite(effect / se):
                        raise ValueError("Failed pooled covariance score cannot be dropped")
                    null_scans[j].append(effect / se)
        else:
            unique_scores = []
            for index, g in enumerate(representatives):
                scan_p, scan_y, scan_h, scan_mask, _, regions = _prepare_label_scan(
                    *scan_args, g, *discovery_args
                )
                fit = fit_regions(scan_p, scan_y, scan_h, g, regions, workers=workers, obs_mask=scan_mask)
                unique_scores.append(fit[:, 2] / fit[:, 3])
                if record_scan_scores:
                    representative_scans[index].extend(float(v) for v in unique_scores[-1])
            for j, (index, sign) in enumerate(scan_map):
                null_scans[j].extend(sign * unique_scores[index])

    observed, meta = calibrate_regions(observed, null_scans, correction=calibration)
    if record_scan_scores:
        meta["scan_scores"] = dict(
            observed=[r["statistic"] for r in observed], unique_null=representative_scans,
            unique_labels=[g.tolist() for g in representatives],
            observed_labels=labels.tolist(), scan_map=scan_map,
        )
    schema = {
        "chrom": pl.String,
        "start": pl.Int64,
        "end": pl.Int64,
        "n_cpgs": pl.Int64,
        "mean_meth_diff": pl.Float64,
        "effect_transformed": pl.Float64,
        "statistic": pl.Float64,
        "residual_df": pl.Int64,
        "length_scale": pl.Float64,
        "donor_fraction": pl.Float64,
        "standard_error": pl.Float64,
        "n_case": pl.Int64,
        "n_control": pl.Int64,
        "combined_pvalue": pl.Float64,
        "combined_qvalue": pl.Float64,
    }
    all_regions = pl.DataFrame(observed, schema=schema) if observed else pl.DataFrame(schema=schema)
    result = (
        all_regions
        if q_cutoff is None
        else all_regions.filter(pl.col("combined_qvalue") <= q_cutoff)
    )
    meta.update(
        method="regional_gls",
        experimental=True,
        promoted=False,
        algorithm_revision="spatial-gls-observed-data-scan-count-reml-v7",
        accelerator="numba-native-analytical-donor-reml",
        masked_accelerator="numba-native-observed-CAR-residual-donor-eigensystem-REML",
        likelihood="observed-cells-with-locus-intercepts-and-shared-treatment",
        unique_null_fits=len(representatives),
        workers=workers,
        approximate_decoy_covariance=pooled_decoy_covariance,
        decoy_covariance="observed-pooled-reml" if pooled_decoy_covariance else "per-region-reml",
        excluded_self_or_mirror=generated_assignments - len(assignments),
        chromosomes=chroms,
        common_cpgs=common_support,
        union_cpgs=union_support,
        support_filter="label-independent-total-union; minimum-group-support-replayed-per-scan",
        residual_variance="unmodified-profiled-REML",
        window_bp=window_bp,
        min_cpgs=min_cpgs,
        dis_merge_bp=max_gap,
        minlen_bp=min_length,
        min_abs_meth_diff=min_effect,
        boundary_min_effect=boundary_min_effect,
        perm_seed=seed,
        requested_permutations=n_perm,
        min_mean_qvalue=q_cutoff,
        coordinate_convention="zero-based-half-open",
        covariance="depth-weighted-CAR-plus-biological-sample-offset-REML",
        inference_unit="effect-selected-region-full-discovery-replayed",
    )
    meta.update(cache_hit=False, source_sha256=code_digest)
    if cache_paths:
        table_path, meta_path = cache_paths
        table_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=table_path.parent, suffix=".parquet", delete=False
        ) as tmp:
            temporary = Path(tmp.name)
        try:
            all_regions.write_parquet(temporary)
            meta.update(
                cache_key=table_path.stem,
                table_sha256=hashlib.sha256(temporary.read_bytes()).hexdigest(),
            )
            os.replace(temporary, table_path)
            with tempfile.NamedTemporaryFile(
                mode="w", dir=table_path.parent, suffix=".json", delete=False
            ) as tmp:
                json.dump(meta, tmp)
                temporary = Path(tmp.name)
            os.replace(temporary, meta_path)
        finally:
            temporary.unlink(missing_ok=True)
    return result, all_regions, meta
