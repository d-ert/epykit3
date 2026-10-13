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
        y = transform_counts(methylated, coverage).T
        if y.shape != (len(g), len(pos)):
            raise ValueError("Shape mismatch between counts and labels")
        h = np.sqrt(np.median(coverage, axis=1))
        yc = (y - y.mean(axis=0)) * h
    from ._regional_gls_native import _fit, _validate, evaluate_covariance

    pos, yc, h, g = _validate(pos, yc, h, g)
    if covariance is None:
        length, fraction, effect, se, objective = _fit(pos, yc, h, g)
    else:
        length, fraction = covariance
        objective, effect, se = evaluate_covariance(pos, yc, h, g, length, fraction)
    if not np.isfinite(objective) or not np.isfinite(se) or se <= 0:
        raise ValueError("Degenerate regional residual covariance; inference cannot be calibrated")
    return dict(
        effect_transformed=float(effect),
        statistic=float(effect / se),
        residual_df=int((len(g) - 1) * len(pos) - 1),
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
    regions = []
    cuts = np.r_[0, np.flatnonzero(np.diff(pos) > max_gap) + 1, len(pos)]
    for lo, hi in pairwise(cuts):
        if hi - lo < min_cpgs:
            continue
        p = pos[lo:hi]
        left = np.searchsorted(p, p - window_bp // 2)
        right = np.searchsorted(p, p + (window_bp + 1) // 2, side="right")
        sums = np.r_[0.0, np.cumsum(effect[lo:hi])]
        smooth = (sums[right] - sums[left]) / (right - left)
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


def calibrate_regions(observed, null_scans):
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
    for r, pv, qv in zip(observed, p, q, strict=False):
        r.update(combined_pvalue=float(pv), combined_qvalue=float(qv))
    return observed, dict(
        null_scans=len(null_scans),
        null_regions=len(null),
        empty_null_scans=sum(len(s) == 0 for s in null_scans),
        correction="pooled-full-discovery-permutation-bh-experimental",
        score_assumption="comparable selected regional null statistics; exchangeable samples",
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
    workers=None,
    cache=True,
):
    """Stream one common-support chromosome at a time; repeat the entire search."""
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
        p.name.removeprefix("chrom=") for p in (root / f"sample={samples[0]}").glob("chrom=*")
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
            numpy=np.__version__,
        )
        key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
        directory = Path(md.analysis_root) / ".cache" / "regional_gls"
        cache_paths = directory / f"{key}.parquet", directory / f"{key}.json"
        cached = _load_regional_cache(cache_paths, key, q_cutoff, workers)
        if cached is not None:
            return cached
    observed, null_scans = [], [[] for _ in assignments]
    common_support = {}

    for chrom in chroms:
        # Optimized Parquet scan and multi-sample common support join
        frames = []
        missing_sample = False
        for sample in samples:
            sample_dir = root / f"sample={sample}" / f"chrom={chrom}"
            if not sample_dir.is_dir():
                missing_sample = True
                break
            parts = sorted(sample_dir.glob("*.parquet"))
            if not parts:
                missing_sample = True
                break
            pattern = str(sample_dir / "*.parquet")
            frame = (
                pl.scan_parquet(pattern)
                .select(["pos", "N_meth", "coverage"])
                .filter(pl.col("coverage") > 0)
                .collect()
            )
            if frame["pos"].n_unique() != len(frame):
                raise ValueError(f"Duplicate CpG coordinates in {sample}/{chrom}")
            frames.append(frame)

        if missing_sample or not frames:
            common_support[chrom] = 0
            continue

        # Union join across samples with minimum group presence filter
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
        ctrl_n_cols = [f"n_{i}" for i in np.flatnonzero(labels == 0)]
        case_n_cols = [f"n_{i}" for i in np.flatnonzero(labels == 1)]
        min_group_obs = min(2, len(ctrl_n_cols), len(case_n_cols))
        min_total_obs = max(min_group_obs * 2, int(np.ceil(0.5 * len(samples))))

        res = res.filter(
            (pl.sum_horizontal([pl.col(c) > 0 for c in ctrl_n_cols]) >= min_group_obs)
            & (pl.sum_horizontal([pl.col(c) > 0 for c in case_n_cols]) >= min_group_obs)
            & (pl.sum_horizontal([pl.col(c) > 0 for c in n_cols]) >= min_total_obs)
        ).sort("pos")

        if len(res) < min_cpgs:
            common_support[chrom] = len(res)
            continue

        pos = res["pos"].to_numpy()
        common_support[chrom] = len(pos)

        m = res.select(m_cols).to_numpy().astype(np.float64)
        n = res.select(n_cols).to_numpy().astype(np.float64)
        valid_mask = n > 0

        # Safe arcsin transformation on positive coverage; missing filled via group mean
        n_safe = np.where(valid_mask, n, 1.0)
        m_safe = np.where(valid_mask, m, 0.0)
        y_raw = np.arcsin(2.0 * (m_safe + 0.05) / (n_safe + 0.1) - 1.0)  # (K, S)

        # Precomputed matrices for fast group effect computation with partial support
        beta_clean = np.where(valid_mask, m_safe / n_safe, 0.0)  # (K, S)
        valid_flt = valid_mask.astype(np.float64)  # (K, S)

        def _calc_effects(g_arr: np.ndarray) -> np.ndarray:
            g_c = g_arr.astype(np.float64)
            g_k = 1.0 - g_c
            k_case = valid_flt @ g_c
            k_ctrl = valid_flt @ g_k
            mean_case = np.where(k_case > 0, (beta_clean @ g_c) / np.maximum(k_case, 1.0), 0.0)
            mean_ctrl = np.where(k_ctrl > 0, (beta_clean @ g_k) / np.maximum(k_ctrl, 1.0), 0.0)
            return mean_case - mean_ctrl

        obs_effects = _calc_effects(labels)

        # Fill missing y_{k, s} with group-conditional transformed mean so missing
        # entries introduce zero artificial variance and zero attenuation of treatment effect
        ctrl_mask = valid_mask[:, labels == 0]
        case_mask = valid_mask[:, labels == 1]
        y_ctrl_mean = np.sum(y_raw[:, labels == 0] * ctrl_mask, axis=1) / np.maximum(
            np.sum(ctrl_mask, axis=1), 1
        )
        y_case_mean = np.sum(y_raw[:, labels == 1] * case_mask, axis=1) / np.maximum(
            np.sum(case_mask, axis=1), 1
        )
        y_filled = y_raw.copy()
        for s_idx in range(len(samples)):
            miss = ~valid_mask[:, s_idx]
            if np.any(miss):
                y_filled[miss, s_idx] = (
                    y_case_mean[miss] if labels[s_idx] == 1 else y_ctrl_mean[miss]
                )

        y_full = y_filled.T  # (S, K)
        # Effective depth per CpG: median over observed samples, scaled by observed sample fraction
        n_sorted = np.where(valid_mask, n, np.nan)
        obs_frac = np.sum(valid_mask, axis=1) / float(len(samples))
        depth_full = np.nanmedian(n_sorted, axis=1) * obs_frac
        h_full = np.sqrt(np.maximum(depth_full, 1.0))
        yc_full = (y_full - y_full.mean(axis=0)) * h_full

        # Precompute discovery geometry for this chromosome
        from ._regional_gls_native import discover_prepared, prepare_discovery

        geom = prepare_discovery(pos, window_bp, min_cpgs, max_gap)

        b_eff = 0.0 if boundary_min_effect is None else float(boundary_min_effect)
        obs_regions = discover_prepared(
            pos, obs_effects, *geom, min_effect, min_cpgs, max_gap, min_length, b_eff
        )
        current = []
        from ._regional_gls_native import fit_regions

        native = fit_regions(pos, yc_full, h_full, labels, obs_regions, workers=workers)
        for (lo, hi), (length, fraction, effect, se, _objective) in zip(
            obs_regions, native, strict=True
        ):
            fit = dict(
                effect_transformed=float(effect),
                statistic=float(effect / se),
                residual_df=int((len(labels) - 1) * (hi - lo) - 1),
                length_scale=float(length),
                donor_fraction=float(fraction),
                standard_error=float(se),
            )
            fit.update(
                chrom=chrom,
                start=int(pos[lo]),
                end=int(pos[hi - 1]) + 1,
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
                g_flt = g.astype(float)
                effects = _calc_effects(g)
                regions = discover_prepared(
                    pos, effects, *geom, min_effect, min_cpgs, max_gap, min_length, b_eff
                )
                for lo, hi in regions:
                    stat = _fast_decoy_statistic(
                        pos[lo:hi],
                        yc_full[:, lo:hi],
                        depth_full[lo:hi],
                        h_full[lo:hi],
                        g_flt,
                        chrom_len,
                        chrom_donor,
                    )
                    null_scans[j].append(stat)
        else:
            unique_scores = []
            for g in representatives:
                effects = _calc_effects(g)
                regions = discover_prepared(
                    pos, effects, *geom, min_effect, min_cpgs, max_gap, min_length, b_eff
                )
                fit = fit_regions(pos, yc_full, h_full, g, regions, workers=workers)
                unique_scores.append(fit[:, 2] / fit[:, 3])
            for j, (index, sign) in enumerate(scan_map):
                null_scans[j].extend(sign * unique_scores[index])

    observed, meta = calibrate_regions(observed, null_scans)
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
        algorithm_revision="spatial-gls-native-reml-v4",
        accelerator="numba-native-analytical-donor-reml",
        unique_null_fits=len(representatives),
        workers=workers,
        approximate_decoy_covariance=pooled_decoy_covariance,
        decoy_covariance="observed-pooled-reml" if pooled_decoy_covariance else "per-region-reml",
        excluded_self_or_mirror=generated_assignments - len(assignments),
        chromosomes=chroms,
        common_cpgs=common_support,
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
