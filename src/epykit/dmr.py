"""Differentially Methylated Region (DMR) calling.

Two algorithms:

``call_dmr_tile_based(methylstore_path, samples_treatment, samples_control,
...)``
    Tile aggregation: sums (N_meth, coverage) across CpGs per sample
    within each fixed-size tile, then runs a full DMC test on the
    tile-level counts. Recommended path -- read-pooled tile tests have
    dramatically more power than per-CpG p-value combination at typical
    WGBS coverage.

``call_dmr_sliding_window(dmc_results, ...)``
    Operates on a precomputed per-CpG DMC table; combines p-values per
    window via signed Stouffer's Z. Faster but lower-power; sign comes
    from each CpG's meth_diff so mixed-direction windows are downweighted.
    Direction (hyper / hypo / mixed) is set by the sign of the mean
    meth_diff rather than a raw site tally.

``smooth_methylation_gaussian`` is a coverage-weighted Gaussian-kernel
smoother -- see its own docstring.
"""

from __future__ import annotations

import gc
import hashlib
import logging
import math
import tempfile
from pathlib import Path
from typing import Iterator, Union

import numpy as np
import polars as pl
from scipy import stats as sp_stats

from ._dmc_store import DMCStore

logger = logging.getLogger(__name__)


def _dmr_sliding_cache_key(
    store: DMCStore,
    window_bp: int,
    min_cpgs: int,
    min_sites_significant: int,
    alpha: float,
    min_abs_meth_diff: float,
    p_col: str,
) -> str:
    """SHA-256 fingerprint of DMR-sliding inputs that affect the result.

    Combines the DMC store's input signature with the DMR parameters so
    a rerun with identical arguments reads from cache. ``step_bp`` is
    deliberately omitted -- the new two-pointer sweep ignores it.
    """
    base_sig = store.manifest.get("input_sig", "")
    h = hashlib.sha256()
    h.update(b"|base="); h.update(base_sig.encode())
    h.update(b"|win=");  h.update(str(int(window_bp)).encode())
    h.update(b"|mincp="); h.update(str(int(min_cpgs)).encode())
    h.update(b"|minsig="); h.update(str(int(min_sites_significant)).encode())
    h.update(b"|alpha=");   h.update(f"{float(alpha):.10g}".encode())
    h.update(b"|delta=");   h.update(f"{float(min_abs_meth_diff):.10g}".encode())
    h.update(b"|pcol=");    h.update(p_col.encode())
    return h.hexdigest()


def _dmr_chain_merge_cache_key(
    store: DMCStore,
    alpha: float,
    min_abs_meth_diff: float,
    dis_merge_bp: int,
    min_cpgs: int,
    pct_sig: float,
    minlen_bp: int,
    p_col: str,
) -> str:
    """SHA-256 fingerprint of chain-merge DMR inputs (DSS callDMR semantics).

    Mirrors :func:`_dmr_sliding_cache_key` for the chain-merge caller.
    """
    base_sig = store.manifest.get("input_sig", "")
    h = hashlib.sha256()
    h.update(b"|base=");    h.update(base_sig.encode())
    h.update(b"|alpha=");   h.update(f"{float(alpha):.10g}".encode())
    h.update(b"|delta=");   h.update(f"{float(min_abs_meth_diff):.10g}".encode())
    h.update(b"|dismerge="); h.update(str(int(dis_merge_bp)).encode())
    h.update(b"|mincp=");   h.update(str(int(min_cpgs)).encode())
    h.update(b"|pctsig=");  h.update(f"{float(pct_sig):.10g}".encode())
    h.update(b"|minlen=");  h.update(str(int(minlen_bp)).encode())
    h.update(b"|pcol=");    h.update(p_col.encode())
    return h.hexdigest()

_DMR_EMPTY_SCHEMA = {
    "chrom":            pl.Utf8,
    "start":            pl.Int32,
    "end":              pl.Int32,
    "n_cpgs":           pl.Int32,
    "n_significant":    pl.Int32,
    "mean_meth_diff":   pl.Float32,
    "combined_pvalue":  pl.Float64,
    "combined_qvalue":  pl.Float64,
    "dmr_type":         pl.Utf8,
}

_DMR_TILE_SCHEMA = {
    "chrom":            pl.Utf8,
    "start":            pl.Int32,
    "end":              pl.Int32,
    "n_cpgs":           pl.Int32,
    "n_case":           pl.Int32,
    "n_control":        pl.Int32,
    "mean_beta_case":   pl.Float32,
    "mean_beta_control": pl.Float32,
    "meth_diff":        pl.Float32,
    "log2_odds_ratio":  pl.Float64,
    "pvalue":           pl.Float64,
    "qvalue":           pl.Float64,
    "dmr_type":         pl.Utf8,
}

_SMOOTH_EMPTY_SCHEMA = {
    "chrom":        pl.Utf8,
    "pos":          pl.Int32,
    "sample":       pl.Utf8,
    "beta_raw":     pl.Float32,
    "beta_smooth":  pl.Float32,
}

# cap merged DMR size to prevent biologically implausible mega-DMRs.
# Mammalian DMRs are typically 200 bp - 5 kb; 10 kb is a generous ceiling.
_MAX_DMR_BP: int = 10_000

# a window's direction is called "mixed" when the fraction of
# valid sites agreeing with the sign of the mean is below this threshold.
_MIXED_DIRECTION_THRESHOLD: float = 0.6


# Internal helpers -- DMC input streaming for sliding-window DMR

def _dmc_store_columns(store: DMCStore) -> set[str]:
    """Return the column set present in the store's per-chrom parquets.

    Reads the schema of the first chromosome's parquet. All chroms are
    assumed to share the same schema (enforced by
    ``process_chromosomes_dmc``).
    """
    chroms = store.chroms()
    if not chroms:
        return set()
    schema = pl.read_parquet_schema(str(store.path / f"chrom={chroms[0]}.parquet"))
    return set(schema.keys())


def _iter_dmc_store_chroms(
    store: DMCStore,
    p_col: str,
) -> Iterator[tuple[str, pl.DataFrame]]:
    """Yield ``(chrom, sorted per-chrom DataFrame)`` from a ``DMCStore``.

    Only the columns the sliding-window sweep needs are read; this keeps
    per-chrom IO at ~50 MB even on chr1 (1.8M CpGs).
    """
    cols = ["chrom", "pos", "meth_diff", "pvalue"]
    if p_col == "qvalue":
        cols.append("qvalue")
    for chrom in store.chroms():
        df = store.read_chrom(chrom, columns=cols)
        if len(df) == 0:
            yield chrom, df
            continue
        yield chrom, df.sort("pos")


def _iter_dataframe_chroms(
    df: pl.DataFrame,
    p_col: str,
) -> Iterator[tuple[str, pl.DataFrame]]:
    """Yield per-chrom slices of an in-memory DMC DataFrame in sorted order."""
    for chrom in sorted(df["chrom"].unique().to_list()):
        yield chrom, df.filter(pl.col("chrom") == chrom).sort("pos")


# Internal helpers -- p-value combination

def _stouffer_combine_signed(
    pvals: np.ndarray,
    meth_diffs: np.ndarray,
    weights: np.ndarray | None = None,
) -> float:
    """Combine per-CpG two-sided p-values via signed Stouffer's Z.

    Each CpG contributes a signed Z-score:
        z_i = sign(meth_diff_i) * Phi^-^1(1 - p_i / 2)
    so that hyper-methylated CpGs contribute positive Z and hypo-methylated
    CpGs contribute negative Z. The combined statistic is
        Z = Sigma w_i * z_i  /  sqrt(Sigma w_i^2)
    which is two-sided-tested. When all CpGs in a window agree in direction
    the |Z| grows as sqrtk and the combined p-value gets correspondingly
    small; when directions are mixed, contributions cancel and the
    combined p-value stays large.

    this replaces the previous Brown's method implementation, which
    required a correlation matrix of the per-CpG test statistics. The old
    code estimated it from genomic distances between CpGs -- a proxy for the
    correlation of methylation STATES, not of the test statistics -- which
    systematically over-inflated the variance correction f and weakened
    combined p-values regardless of how strong the per-CpG signal was.
    Stouffer's Z does NOT model correlation between adjacent CpGs. Because
    neighbouring WGBS CpGs are positively correlated, the true variance of
    ``Sigma w_i z_i`` exceeds ``Sigma w_i^2``, so dividing by
    ``sqrt(Sigma w_i^2)`` understates the SD and makes the combined p-value
    *anti-conservative* (too small) in dense regions -- the opposite of the
    over-smoothing the old Brown's-method proxy caused. Treat the region
    p-value as a ranking signal; for calibrated region-level inference use the
    permutation empirical FDR (``tl.dmr(..., empirical_fdr=True)``).

    Parameters
    ----------
    pvals : np.ndarray
        Two-sided p-values for each CpG in the window.
    meth_diffs : np.ndarray
        Signed per-CpG effect sizes (mean_beta_case - mean_beta_ctrl).
        Used only for direction; magnitude is ignored.
    weights : np.ndarray, optional
        Per-CpG weights (e.g. coverage). Defaults to equal weights.

    Returns
    -------
    float
        Two-sided combined p-value.
    """
    pvals      = np.asarray(pvals,      dtype=np.float64)
    meth_diffs = np.asarray(meth_diffs, dtype=np.float64)

    valid = (
        ~np.isnan(pvals) & (pvals > 0.0) & (pvals <= 1.0)
        & ~np.isnan(meth_diffs)
    )
    if not np.any(valid):
        return float("nan")

    p_valid    = np.clip(pvals[valid], np.finfo(float).tiny, 1.0 - 1e-15)
    diff_valid = meth_diffs[valid]

    # Magnitude Z from two-sided p-value: |z| = Phi^-^1(1 - p/2)
    z_mag = sp_stats.norm.isf(p_valid / 2.0)
    # Signed contribution: hyper (+) vs hypo (-). meth_diff == 0 -> no
    # contribution (sign = 0), which is correct: zero-effect CpGs neither
    # add nor subtract evidence.
    z_signed = np.sign(diff_valid) * z_mag

    if weights is None:
        w = np.ones_like(z_signed)
    else:
        w = np.asarray(weights, dtype=np.float64)[valid]
        w = np.where(np.isfinite(w) & (w >= 0), w, 0.0)

    w_sq_sum = float(np.sum(w * w))
    if w_sq_sum <= 0.0:
        return float("nan")

    z_combined = float(np.sum(w * z_signed) / np.sqrt(w_sq_sum))
    # Two-sided normal tail
    return float(2.0 * sp_stats.norm.sf(abs(z_combined)))


def _merge_intervals(starts: list[int], ends: list[int]) -> list[tuple[int, int]]:
    """Merge overlapping (start, end) integer intervals."""
    if not starts:
        return []
    pairs = sorted(zip(starts, ends), key=lambda x: x[0])
    merged: list[tuple[int, int]] = [pairs[0]]
    for s, e in pairs[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def _classify_direction(
    mean_diff: float,
    n_hyper: int,
    n_hypo: int,
) -> str:
    """Classify DMR direction from mean effect + per-site sign tally.

    previously the direction was the larger of n_hyper / n_hypo.
    A window with 6 hyper sites at +0.12 and 4 hypo sites at -0.40 was
    called "hyper" even though the mean effect was clearly negative.

    The new rule:
      - Mean direction governs (sign of mean_meth_diff). This matches
        what the tile-based path does when it pools reads.
      - If the per-site sign tally is too close (< 60 % majority), the
        window is labelled "mixed" so downstream filters can drop it.
    """
    total_signed = n_hyper + n_hypo
    if total_signed == 0:
        return "mixed"
    consensus = max(n_hyper, n_hypo) / total_signed
    if consensus < _MIXED_DIRECTION_THRESHOLD:
        return "mixed"
    if np.isnan(mean_diff):
        # Fall back to majority tally if mean is undefined for some reason.
        return "hyper" if n_hyper > n_hypo else "hypo"
    return "hyper" if mean_diff > 0 else "hypo"


def _recompute_dmr_stats(
    chrom: str,
    start: int,
    end: int,
    positions: np.ndarray,
    meth_diffs: np.ndarray,
    pvals: np.ndarray,
    is_sig: np.ndarray,
    min_cpgs: int,
    min_sites_significant: int,
    # Optional pre-computed prefix arrays for O(log n) slice
    cum_sig: np.ndarray | None = None,
) -> dict | None:
    """Recompute accurate per-site statistics over a merged interval.

    When ``cum_sig`` is supplied the significance count is computed in O(1)
    via prefix-sum lookup; otherwise falls back to a boolean mask scan.
    """
    # reject biologically implausible mega-DMRs that arise when
    # overlapping candidate windows collapse across many megabases.
    if (end - start) > _MAX_DMR_BP:
        return None

    # Use searchsorted for O(log n) slice instead of boolean mask
    lo = int(np.searchsorted(positions, start,  side="left"))
    hi = int(np.searchsorted(positions, end,    side="left"))

    n_cpgs = hi - lo
    if n_cpgs < min_cpgs:
        return None

    if cum_sig is not None:
        n_sig = int(cum_sig[hi] - cum_sig[lo])
    else:
        n_sig = int(is_sig[lo:hi].sum())

    if n_sig < min_sites_significant:
        return None

    window_diffs = meth_diffs[lo:hi]
    window_pvals = pvals[lo:hi]

    valid_diffs = window_diffs[~np.isnan(window_diffs)]
    if len(valid_diffs) == 0:
        return None

    n_hyper = int((valid_diffs > 0).sum())
    n_hypo  = int((valid_diffs < 0).sum())

    # signed Stouffer's Z. Sign comes from per-CpG meth_diff so the
    # combined statistic naturally cancels mixed-direction windows.
    combined_p = _stouffer_combine_signed(window_pvals, window_diffs)
    if np.isnan(combined_p):
        return None

    mean_diff = float(np.nanmean(window_diffs))
    dmr_type  = _classify_direction(mean_diff, n_hyper, n_hypo)

    return {
        "chrom":           chrom,
        "start":           start,
        "end":             end,
        "n_cpgs":          n_cpgs,
        "n_significant":   n_sig,
        "mean_meth_diff":  float(np.float32(mean_diff)),
        "combined_pvalue": float(combined_p),
        "dmr_type":        dmr_type,
    }


# Public API -- sliding-window DMR calling (works from a DMC table)

def call_dmr_sliding_window(
    dmc_results: Union[pl.DataFrame, DMCStore, str, Path],
    window_bp: int = 500,
    step_bp: int = 250,
    min_cpgs: int = 5,
    min_sites_significant: int = 3,
    alpha: float = 0.05,
    min_abs_meth_diff: float = 0.1,
) -> pl.DataFrame:
    """Call DMRs by aggregating DMC sites into overlapping sliding windows.

    This method takes a precomputed DMC table and combines per-CpG p-values
    region-by-region with signed Stouffer's Z. It is fast and reuses an
    existing DMC call, but has lower power than the tile-based path
    (`call_dmr_tile_based`) because it cannot pool reads -- windows whose
    individual CpGs aren't significant won't gather enough sig sites to
    pass the `min_sites_significant` gate.

    Memory scaling
    --------------
    The implementation uses a two-pointer CpG-anchored sweep -- each
    chromosome's peak memory scales with the number of CpGs on that
    chromosome (a few hundred MB for full human autosomes), independent
    of genomic span. The earlier bp-grid enumeration would have
    materialised one window position every ``step_bp`` across the whole
    chromosome (5M positions on chr1 alone) and OOM'd on full-genome
    inputs. ``step_bp`` is still accepted for backward compatibility
    but is effectively a no-op in the new implementation: every CpG
    anchors a candidate and the downstream merge collapses redundancy.

    Parameters
    ----------
    dmc_results : pl.DataFrame | DMCStore | str | Path
        DMC results to process. Accepts:

        * ``pl.DataFrame`` -- in-memory table (legacy path); held in
          memory for the whole DMR pass.
        * ``DMCStore`` -- handle to a persistent per-chrom parquet
          directory (returned by
          ``process_chromosomes_dmc(..., return_store=True)``).
          Chromosomes are streamed from disk; peak memory is
          O(largest chromosome).
        * ``str`` / ``Path`` -- path to a populated DMC store
          directory; opened via :meth:`DMCStore.open`.

        Required columns: chrom, pos, meth_diff, pvalue.
        Optional: qvalue (used in preference to pvalue when present).
    window_bp : int
        Window width in base pairs.
    step_bp : int
        Accepted but no longer meaningful (see "Memory scaling" above).
        Kept in the signature so old call sites don't break.
    min_cpgs : int
        Minimum CpG count in a *merged* DMR.
    min_sites_significant : int
        Minimum significant CpG sites in a window for it to be a candidate.
    alpha : float
        Significance threshold for qvalue / pvalue.
    min_abs_meth_diff : float
        Minimum |meth_diff| for a site to count as significant.

    Returns
    -------
    pl.DataFrame
        Columns: chrom, start, end, n_cpgs, n_significant,
                 mean_meth_diff, combined_pvalue, combined_qvalue,
                 dmr_type ("hyper" | "hypo" | "mixed").
        ``combined_qvalue`` is BH-corrected genome-wide across DMR
        candidates that passed the per-window gate.
    """
    if isinstance(dmc_results, (str, Path)):
        dmc_results = DMCStore.open(dmc_results)

    if step_bp > window_bp:
        raise ValueError(
            f"step_bp ({step_bp}) must be <= window_bp ({window_bp})"
        )

    # DMR cache: when the input is a DMCStore with a known input_sig,
    # the entire sliding-window output is a pure function of that sig
    # plus the DMR params. Stash the result inside the DMC store dir.
    dmr_cache_path: Path | None = None
    if isinstance(dmc_results, DMCStore) and dmc_results.manifest.get("input_sig"):
        available_cols = _dmc_store_columns(dmc_results)
        p_col_for_key = "qvalue" if "qvalue" in available_cols else "pvalue"
        key = _dmr_sliding_cache_key(
            dmc_results, window_bp, min_cpgs, min_sites_significant,
            alpha, min_abs_meth_diff, p_col_for_key,
        )
        dmr_cache_path = dmc_results.path / f".dmr_sliding_{key[:16]}.parquet"
        if dmr_cache_path.exists():
            cached = pl.read_parquet(str(dmr_cache_path))
            logger.info(
                "DMR sliding-window cache hit at %s (%d DMR(s)); skipping recompute.",
                dmr_cache_path.name, len(cached),
            )
            return cached

    if isinstance(dmc_results, DMCStore):
        available_cols = _dmc_store_columns(dmc_results)
        required = {"chrom", "pos", "meth_diff", "pvalue"}
        missing  = required - available_cols
        if missing:
            raise ValueError(f"DMC results missing required columns: {missing}")
        p_col = "qvalue" if "qvalue" in available_cols else "pvalue"
        chrom_iter = _iter_dmc_store_chroms(dmc_results, p_col)
    else:
        required = {"chrom", "pos", "meth_diff", "pvalue"}
        missing  = required - set(dmc_results.columns)
        if missing:
            raise ValueError(f"DMC results missing required columns: {missing}")
        p_col = "qvalue" if "qvalue" in dmc_results.columns else "pvalue"
        chrom_iter = _iter_dataframe_chroms(dmc_results, p_col)

    logger.info(
        "call_dmr_sliding_window: window=%d bp, step=%d bp, "
        "min_cpgs=%d, min_sig=%d, alpha=%.3f, min_|Deltabeta|=%.2f, p_col=%s",
        window_bp, step_bp, min_cpgs, min_sites_significant,
        alpha, min_abs_meth_diff, p_col,
    )

    all_records: list[dict] = []

    for chrom, chrom_df in chrom_iter:
        if len(chrom_df) == 0:
            continue

        positions  = chrom_df["pos"].to_numpy()
        meth_diffs = chrom_df["meth_diff"].to_numpy(allow_copy=True).astype(np.float32)
        # The significance gate uses the FDR-controlled column (qvalue if
        # present); the Stouffer combine MUST use the raw per-CpG p-values,
        # which are ~U(0,1) under the null. q-values are not uniform, so
        # combining them does not yield a valid p-value (M-DMR1). chain-merge
        # already keeps these two roles separate; sliding-window now matches.
        sig_vals   = chrom_df[p_col].to_numpy(allow_copy=True).astype(np.float64)
        raw_pvals  = chrom_df["pvalue"].to_numpy(allow_copy=True).astype(np.float64)

        is_sig = (
            (~np.isnan(sig_vals))
            & (sig_vals < alpha)
            & (~np.isnan(meth_diffs))
            & (np.abs(meth_diffs) >= min_abs_meth_diff)
        )

        # ---------------------------------------------------------------
        # Two-pointer CpG-anchored sweep.
        #
        # Replaces the prior bp-grid enumeration (np.arange over the whole
        # chromosome with step_bp) which scaled with genome span and
        # OOM'd on full-genome 22M-CpG input (chr1 alone produced 5M
        # window positions). The new pass is O(n_CpGs) per chrom and
        # bounded in memory by the size of the per-chrom arrays.
        #
        # For each CpG i, the candidate window is
        # ``[positions[i], positions[i] + window_bp)``; we maintain a
        # running right pointer j and a rolling significant-CpG count.
        # ``step_bp`` is accepted for backward compatibility but is no
        # longer meaningful -- every CpG anchors a candidate, and the
        # downstream merge collapses redundancy without any change in
        # the final DMR set.
        # ---------------------------------------------------------------
        is_sig_int = is_sig.astype(np.int8, copy=False)
        n_pos = len(positions)
        cand_starts: list[int] = []
        cand_ends: list[int] = []

        j = 0
        n_sig_running = 0
        for i in range(n_pos):
            limit = positions[i] + window_bp
            # advance j until positions[j] no longer fits in [pos_i, pos_i+window_bp)
            while j < n_pos and positions[j] < limit:
                n_sig_running += int(is_sig_int[j])
                j += 1
            # window [i, j) -- all CpGs in [positions[i], positions[i]+window_bp)
            n_cpgs_win = j - i
            if n_cpgs_win >= min_cpgs and n_sig_running >= min_sites_significant:
                cand_starts.append(int(positions[i]))
                cand_ends.append(int(positions[i]) + window_bp)
            # drop CpG i from the rolling count before advancing i
            n_sig_running -= int(is_sig_int[i])

        if not cand_starts:
            logger.info("  %s: no candidate windows", chrom)
            del chrom_df, positions, meth_diffs, sig_vals, raw_pvals, is_sig, is_sig_int
            gc.collect()
            continue

        merged_spans = _merge_intervals(cand_starts, cand_ends)
        chrom_dmrs = 0

        # _recompute_dmr_stats can still use a prefix-sum array for O(1)
        # range counts in the final pass -- that's per-merged-span (a
        # small number) so it doesn't break the O(n_CpGs) budget.
        cum_sig = np.empty(n_pos + 1, dtype=np.int32)
        cum_sig[0] = 0
        np.cumsum(is_sig_int.astype(np.int32, copy=False), out=cum_sig[1:])

        for start, end in merged_spans:
            rec = _recompute_dmr_stats(
                chrom, start, end,
                positions, meth_diffs, raw_pvals, is_sig,
                min_cpgs, min_sites_significant,
                cum_sig=cum_sig,
            )
            if rec is not None:
                all_records.append(rec)
                chrom_dmrs += 1

        logger.info(
            "  %s: %d candidate span(s) -> %d DMR(s)",
            chrom, len(merged_spans), chrom_dmrs,
        )
        # Free per-chrom buffers before next iteration. Critical when
        # streaming from a DMCStore on a 22M-site genome: without this
        # the buffer-pool references can pile up and defeat the
        # whole-table-streaming win.
        del chrom_df, positions, meth_diffs, sig_vals, raw_pvals, is_sig, is_sig_int, cum_sig
        gc.collect()

    if not all_records:
        logger.warning("No DMRs found with current filters")
        empty = pl.DataFrame(schema=_DMR_EMPTY_SCHEMA)
        if dmr_cache_path is not None:
            empty.write_parquet(str(dmr_cache_path))
        return empty

    dmr_df = (
        pl.DataFrame(all_records)
        .with_columns([
            pl.col("start").cast(pl.Int32),
            pl.col("end").cast(pl.Int32),
            pl.col("n_cpgs").cast(pl.Int32),
            pl.col("n_significant").cast(pl.Int32),
            pl.col("mean_meth_diff").cast(pl.Float32),
        ])
        .sort(["chrom", "start"])
    )

    # BH-correct DMR-level combined p-values so downstream filters
    # are operating on q-values. Without this the sliding-window output
    # was effectively un-corrected at the region level.
    from .dmc import apply_multiple_testing_correction

    dmr_df = apply_multiple_testing_correction(
        dmr_df,
        method="fdr_bh",
        pvalue_col="combined_pvalue",
        qvalue_col="combined_qvalue",
    )
    if dmr_cache_path is not None:
        dmr_df.write_parquet(str(dmr_cache_path))
        logger.info(
            "DMR sliding-window result cached at %s", dmr_cache_path.name,
        )
    return dmr_df


# Public API -- chain-and-merge DMR calling (DSS callDMR semantics)

DMR_PRESETS: dict[str, dict] = {
    # "strict": very confident DMRs only. Use when downstream uses cannot
    # tolerate false positives (e.g. follow-up validation experiments).
    # Higher alpha bar, larger effect-size floor, more CpGs required.
    "strict": dict(
        alpha=1e-6, min_abs_meth_diff=0.20, dis_merge_bp=250,
        min_cpgs=5, pct_sig=0.5, minlen_bp=100,
    ),
    # "default": balanced preset for general WGBS analyses. alpha=1e-4 is
    # one order looser than DSS callDMR's default (1e-5) -- empirically
    # this captures real-but-moderate signal that DSS's strict gate rejects
    # without crashing PPV. The 10% per-CpG effect-size floor is kept
    # (matches DSS delta=0.1) so individual measurement noise can't anchor
    # chains. Recommended starting point for most users; use 'strict' for
    # validation-ready DMRs (DSS-strict alpha) or 'permissive' for
    # exploratory / recall-oriented analyses.
    "default": dict(
        alpha=1e-4, min_abs_meth_diff=0.10, dis_merge_bp=500,
        min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    ),
    # "permissive": recall-oriented. Loosens alpha and the gap rule, drops
    # the effect-size floor halfway. Useful for exploratory analyses, gene-
    # set enrichment, or comparisons where false negatives are costlier
    # than false positives. Expect noticeably lower PPV.
    "permissive": dict(
        alpha=1e-4, min_abs_meth_diff=0.05, dis_merge_bp=1000,
        min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    ),
}


def call_dmr_chain_merge(
    dmc_results: Union[pl.DataFrame, DMCStore, str, Path],
    *,
    preset: str | None = None,
    alpha: float = 0.05,
    min_abs_meth_diff: float = 0.1,
    dis_merge_bp: int = 500,
    min_cpgs: int = 3,
    pct_sig: float = 0.5,
    minlen_bp: int = 50,
    use_q_for_sig: bool = False,
) -> pl.DataFrame:
    """Call DMRs by chaining contiguous significant CpGs and merging gaps.

    Reimplements DSS's ``callDMR`` semantics on top of an epykit DMC
    table:

      1. Mark significant CpGs:
         ``(pvalue < alpha) AND (|meth_diff| >= min_abs_meth_diff)``.
         When ``use_q_for_sig=True`` and a ``qvalue`` column is present,
         the q-value drives the significance gate instead.
      2. Walk sorted positions per chromosome. A new sig CpG joins the
         current chain when its distance to the *previous significant*
         CpG is <= ``dis_merge_bp``; otherwise it starts a new chain.
      3. A chain's span is ``[first_sig_pos, last_sig_pos + 1)``.
      4. Apply filters: span length >= ``minlen_bp``; total CpGs in span
         (sig + non-sig) >= ``min_cpgs``; fraction significant >=
         ``pct_sig``.
      5. Combine the per-CpG p-values in the surviving span via signed
         Stouffer's Z (reused from ``call_dmr_sliding_window``); classify
         direction (``hyper`` / ``hypo`` / ``mixed``) from the mean
         ``meth_diff`` and per-site sign tally.
      6. BH-correct the surviving combined p-values genome-wide.

    Geometrically this is more permissive than ``call_dmr_sliding_window``
    for sparse-cluster signal: two sig CpGs 140 bp apart cannot fit in a
    100 bp anchored window, but they chain in DSS at ``dis_merge_bp=100``
    because their gap is <= 100 bp from one sig CpG to the next.

    Tuning guidance
    ---------------
    If recall is too low for your use case, **loosen ``dis_merge_bp``
    first** (e.g. 500 -> 1000). It's the single highest-leverage knob:
    CpG-poor intergenic and intronic regions need wider gap-merging just
    to recover real DMRs whose CpGs are spread out. Loosening
    ``dis_merge_bp`` typically gains 5-10 pp of recall for only ~6 pp of
    PPV cost -- a genuine Pareto improvement, unlike loosening ``alpha``
    (which crashes PPV by 15-20 pp for only 2-4 pp of recall on most
    datasets we've benchmarked).

    The ``pct_sig`` knob is **effectively dead at strict ``alpha``** (e.g.
    1e-5 to 1e-4): the few CpGs that pass cluster tightly enough that
    >50% of any candidate's CpGs are already significant, so the 0.5
    threshold never bites. Don't bother tuning it unless you've also
    loosened ``alpha`` substantially.

    For common scenarios prefer the ``preset=`` bundles instead of
    hand-tuning: ``"strict"`` (validation-ready DMRs), ``"default"``
    (DSS-equivalent, recommended), ``"permissive"`` (exploratory / recall-
    oriented).

    Parameters
    ----------
    dmc_results
        DMC results -- same input contract as
        :func:`call_dmr_sliding_window` (DataFrame, DMCStore, or path).
    preset : {"strict", "default", "permissive"}, optional
        Apply a named parameter bundle from :data:`DMR_PRESETS`. Any
        explicit kwarg passed alongside ``preset`` overrides the bundled
        value (so you can pick a preset and tweak one knob).
    alpha
        Significance threshold for ``pvalue`` (or ``qvalue`` when
        ``use_q_for_sig=True``).
    min_abs_meth_diff
        Minimum ``|meth_diff|`` for a CpG to count as significant.
        Default ``0.1`` matches the convention in the methylation
        literature (DSS, methylKit, BSmooth all report DMRs with at
        least 10% methylation difference).
    dis_merge_bp
        Maximum distance (in bp) between consecutive significant CpGs
        for them to belong to the same chain. DSS default: 100. See the
        tuning guidance above -- this is the first knob to loosen if
        recall is too low.
    min_cpgs
        Minimum total CpGs (sig + non-sig) inside a chain's span. DSS
        default: 3.
    pct_sig
        Minimum fraction of CpGs in the span that must be significant.
        DSS default: 0.5. Note: this knob is effectively dead at strict
        ``alpha`` values.
    minlen_bp
        Minimum span length in base pairs. DSS default: 50.
    use_q_for_sig
        If True and ``qvalue`` is present, use it for the significance
        gate.

    Returns
    -------
    pl.DataFrame
        Same schema as :func:`call_dmr_sliding_window`:
        ``chrom, start, end, n_cpgs, n_significant, mean_meth_diff,
        combined_pvalue, combined_qvalue, dmr_type``.
    """
    # Resolve preset bundle. Caller-provided kwargs override bundled values
    # by tracking which params arrived at their default values vs explicit.
    # Easiest way: re-bind locals from the preset only if they match the
    # signature defaults (signaling "user didn't set this").
    if preset is not None:
        if preset not in DMR_PRESETS:
            raise ValueError(
                f"Unknown preset {preset!r}. Choose from {list(DMR_PRESETS)}."
            )
        bundle = DMR_PRESETS[preset]
        # Default sentinels match the signature defaults above.
        _SIG_DEFAULTS = dict(
            alpha=0.05, min_abs_meth_diff=0.1, dis_merge_bp=500,
            min_cpgs=3, pct_sig=0.5, minlen_bp=50,
        )
        if alpha == _SIG_DEFAULTS["alpha"]:                       alpha = bundle["alpha"]
        if min_abs_meth_diff == _SIG_DEFAULTS["min_abs_meth_diff"]: min_abs_meth_diff = bundle["min_abs_meth_diff"]
        if dis_merge_bp == _SIG_DEFAULTS["dis_merge_bp"]:         dis_merge_bp = bundle["dis_merge_bp"]
        if min_cpgs == _SIG_DEFAULTS["min_cpgs"]:                 min_cpgs = bundle["min_cpgs"]
        if pct_sig == _SIG_DEFAULTS["pct_sig"]:                   pct_sig = bundle["pct_sig"]
        if minlen_bp == _SIG_DEFAULTS["minlen_bp"]:               minlen_bp = bundle["minlen_bp"]
    if isinstance(dmc_results, (str, Path)):
        dmc_results = DMCStore.open(dmc_results)

    # Cache lookup (mirror of the sliding-window cache).
    dmr_cache_path: Path | None = None
    if isinstance(dmc_results, DMCStore) and dmc_results.manifest.get("input_sig"):
        available_cols = _dmc_store_columns(dmc_results)
        sig_col = "qvalue" if (use_q_for_sig and "qvalue" in available_cols) else "pvalue"
        key = _dmr_chain_merge_cache_key(
            dmc_results,
            alpha=alpha,
            min_abs_meth_diff=min_abs_meth_diff,
            dis_merge_bp=dis_merge_bp,
            min_cpgs=min_cpgs,
            pct_sig=pct_sig,
            minlen_bp=minlen_bp,
            p_col=sig_col,
        )
        dmr_cache_path = dmc_results.path / f".dmr_chain_merge_{key[:16]}.parquet"
        if dmr_cache_path.exists():
            cached = pl.read_parquet(str(dmr_cache_path))
            logger.info(
                "DMR chain-merge cache hit at %s (%d DMR(s)); skipping recompute.",
                dmr_cache_path.name, len(cached),
            )
            return cached

    if isinstance(dmc_results, DMCStore):
        available_cols = _dmc_store_columns(dmc_results)
        required = {"chrom", "pos", "meth_diff", "pvalue"}
        missing = required - available_cols
        if missing:
            raise ValueError(f"DMC results missing required columns: {missing}")
        sig_col = "qvalue" if (use_q_for_sig and "qvalue" in available_cols) else "pvalue"
        chrom_iter = _iter_dmc_store_chroms(dmc_results, sig_col)
    else:
        required = {"chrom", "pos", "meth_diff", "pvalue"}
        missing = required - set(dmc_results.columns)
        if missing:
            raise ValueError(f"DMC results missing required columns: {missing}")
        sig_col = "qvalue" if (use_q_for_sig and "qvalue" in dmc_results.columns) else "pvalue"
        chrom_iter = _iter_dataframe_chroms(dmc_results, sig_col)

    logger.info(
        "call_dmr_chain_merge: alpha=%.3g, min_|Deltabeta|=%.2f, dis_merge=%d bp, "
        "min_cpgs=%d, pct_sig=%.2f, minlen=%d bp, sig_col=%s",
        alpha, min_abs_meth_diff, dis_merge_bp,
        min_cpgs, pct_sig, minlen_bp, sig_col,
    )

    all_records: list[dict] = []

    for chrom, chrom_df in chrom_iter:
        if len(chrom_df) == 0:
            continue

        positions  = chrom_df["pos"].to_numpy()
        meth_diffs = chrom_df["meth_diff"].to_numpy(allow_copy=True).astype(np.float32)
        # pvals is always the raw pvalue column (used by Stouffer's Z).
        pvals = chrom_df["pvalue"].to_numpy(allow_copy=True).astype(np.float64)
        # sig_vals is what we threshold against alpha -- either pvalue or qvalue.
        sig_vals = chrom_df[sig_col].to_numpy(allow_copy=True).astype(np.float64)

        is_sig = (
            (~np.isnan(sig_vals))
            & (sig_vals < alpha)
            & (~np.isnan(meth_diffs))
            & (np.abs(meth_diffs) >= min_abs_meth_diff)
        )
        sig_idx = np.flatnonzero(is_sig)
        if sig_idx.size == 0:
            del chrom_df, positions, meth_diffs, pvals, sig_vals, is_sig, sig_idx
            gc.collect()
            continue

        # Chain consecutive sig CpGs whose pairwise gap <= dis_merge_bp.
        # The chain is a list of (first_sig_idx, last_sig_idx) pairs.
        sig_pos = positions[sig_idx]
        gaps = np.diff(sig_pos)
        # Boundary indices: positions where a new chain begins (gap > dis_merge_bp).
        break_after = np.flatnonzero(gaps > dis_merge_bp)
        # Convert to chain boundaries over the sig_idx array.
        chain_lo = np.concatenate([[0], break_after + 1])
        chain_hi = np.concatenate([break_after, [sig_idx.size - 1]])

        # Prefix-sum on is_sig for O(1) significant-count lookups.
        is_sig_int = is_sig.astype(np.int32, copy=False)
        cum_sig = np.empty(positions.size + 1, dtype=np.int32)
        cum_sig[0] = 0
        np.cumsum(is_sig_int, out=cum_sig[1:])

        chrom_dmrs = 0
        for lo, hi in zip(chain_lo, chain_hi):
            i_first = int(sig_idx[lo])
            i_last  = int(sig_idx[hi])
            start_pos = int(positions[i_first])
            # End is exclusive in epykit's schema; +1 to include i_last.
            end_pos = int(positions[i_last]) + 1
            if (end_pos - start_pos) < minlen_bp:
                continue

            # All CpGs (sig + non-sig) inside [start_pos, end_pos).
            i_span_lo = int(np.searchsorted(positions, start_pos, side="left"))
            i_span_hi = int(np.searchsorted(positions, end_pos, side="left"))
            n_cpgs_span = i_span_hi - i_span_lo
            if n_cpgs_span < min_cpgs:
                continue
            n_sig_span = int(cum_sig[i_span_hi] - cum_sig[i_span_lo])
            if n_sig_span / max(n_cpgs_span, 1) < pct_sig:
                continue

            span_pvals = pvals[i_span_lo:i_span_hi]
            span_diffs = meth_diffs[i_span_lo:i_span_hi]

            valid_diffs = span_diffs[~np.isnan(span_diffs)]
            if len(valid_diffs) == 0:
                continue
            n_hyper = int((valid_diffs > 0).sum())
            n_hypo  = int((valid_diffs < 0).sum())

            combined_p = _stouffer_combine_signed(span_pvals, span_diffs)
            if np.isnan(combined_p):
                continue

            mean_diff = float(np.nanmean(span_diffs))
            dmr_type  = _classify_direction(mean_diff, n_hyper, n_hypo)

            all_records.append({
                "chrom":           chrom,
                "start":           start_pos,
                "end":             end_pos,
                "n_cpgs":          n_cpgs_span,
                "n_significant":   n_sig_span,
                "mean_meth_diff":  float(np.float32(mean_diff)),
                "combined_pvalue": float(combined_p),
                "dmr_type":        dmr_type,
            })
            chrom_dmrs += 1

        logger.info(
            "  %s: %d sig CpG(s) -> %d chain(s) -> %d DMR(s)",
            chrom, int(sig_idx.size), len(chain_lo), chrom_dmrs,
        )
        del chrom_df, positions, meth_diffs, pvals, sig_vals
        del is_sig, is_sig_int, sig_idx, sig_pos, cum_sig
        gc.collect()

    if not all_records:
        logger.warning("No DMRs found with current filters")
        empty = pl.DataFrame(schema=_DMR_EMPTY_SCHEMA)
        if dmr_cache_path is not None:
            empty.write_parquet(str(dmr_cache_path))
        return empty

    dmr_df = (
        pl.DataFrame(all_records)
        .with_columns([
            pl.col("start").cast(pl.Int32),
            pl.col("end").cast(pl.Int32),
            pl.col("n_cpgs").cast(pl.Int32),
            pl.col("n_significant").cast(pl.Int32),
            pl.col("mean_meth_diff").cast(pl.Float32),
        ])
        .sort(["chrom", "start"])
    )

    # BH correction (mirrors call_dmr_sliding_window) so downstream
    # filters can operate on q-values consistently.
    from .dmc import apply_multiple_testing_correction

    dmr_df = apply_multiple_testing_correction(
        dmr_df,
        method="fdr_bh",
        pvalue_col="combined_pvalue",
        qvalue_col="combined_qvalue",
    )
    if dmr_cache_path is not None:
        dmr_df.write_parquet(str(dmr_cache_path))
        logger.info(
            "DMR chain-merge result cached at %s", dmr_cache_path.name,
        )
    return dmr_df


# Public API -- tile-based DMR calling

def _aggregate_sample_to_tiles(
    src_part_file: Path,
    chrom: str,
    tile_size_bp: int,
) -> pl.DataFrame | None:
    """Aggregate one sample/chromosome's per-CpG counts to per-tile sums.

    Returns a DataFrame with columns: chrom, pos (= tile start), strand,
    N_meth, N_unmeth, coverage, n_cpgs.  Or None when the source is missing.
    """
    if not src_part_file.exists():
        return None

    df = pl.read_parquet(str(src_part_file))
    if len(df) == 0:
        return None

    # Tile assignment: pos // tile * tile gives left-inclusive boundary.
    tile_col = (pl.col("pos") // tile_size_bp) * tile_size_bp
    tiled = (
        df.with_columns(tile_col.cast(pl.Int32).alias("tile_start"))
        .group_by("tile_start")
        .agg([
            pl.sum("N_meth").alias("N_meth"),
            pl.sum("coverage").alias("coverage"),
            pl.len().alias("n_cpgs"),
            # Preserve a strand value: first non-"*" if available, else "*"
            (
                pl.when(pl.col("strand") != "*").then(pl.col("strand")).otherwise(None)
                .drop_nulls().first()
            ).alias("strand_real")
            if "strand" in df.columns
            else pl.lit("*").alias("strand_real"),
        ])
        .with_columns([
            pl.lit(chrom).alias("chrom"),
            pl.col("tile_start").alias("pos"),
            pl.col("strand_real").fill_null("*").alias("strand"),
            (pl.col("coverage") - pl.col("N_meth")).alias("N_unmeth"),
        ])
        .drop("tile_start", "strand_real")
        .sort("pos")
    )
    return tiled


def _merge_adjacent_tiles(dmr_df: pl.DataFrame) -> pl.DataFrame:
    """Merge adjacent significant tiles on the same chromosome with same direction.

    The combined p-value uses signed Stouffer with the correct two-sided
    -> one-sided conversion ``z = isf(p/2)`` and a running ``(sum_z, n)``
    accumulator so chains of length > 2 combine as ``sum_z / sqrt(n)``,
    not by iterative pairwise /sqrt(2).
    """
    if dmr_df.is_empty():
        return dmr_df

    sorted_df = dmr_df.sort(["chrom", "start"])
    rows = sorted_df.to_dicts()
    merged: list[dict] = []
    current: dict | None = None
    current_sum_z = 0.0
    current_n_z = 0

    def _abs_z(p: float) -> float:
        # Two-sided p -> magnitude of one-sided z; clamp p to avoid inf.
        return float(sp_stats.norm.isf(max(p, 1e-300) / 2.0))

    def _finalise(c: dict, sum_z: float, n_z: int) -> dict:
        if n_z > 0:
            z_comb = sum_z / math.sqrt(n_z)
            c["pvalue"] = float(2.0 * sp_stats.norm.sf(abs(z_comb)))
        return c

    for row in rows:
        if current is None:
            current = dict(row)
            current["_count"] = 1
            current_sum_z = _abs_z(row["pvalue"])
            current_n_z = 1
            continue
        if (
            row["chrom"] == current["chrom"]
            and row["start"] <= current["end"]
            and row["dmr_type"] == current["dmr_type"]
        ):
            prev_n = current["_count"]
            current["end"] = row["end"]
            current["n_cpgs"] = current["n_cpgs"] + row["n_cpgs"]
            for col in ("meth_diff", "log2_odds_ratio",
                        "mean_beta_case", "mean_beta_control"):
                if col in current and current[col] is not None and row.get(col) is not None:
                    current[col] = (current[col] * prev_n + row[col]) / (prev_n + 1)
            for col in ("n_case", "n_control"):
                if col in current and current[col] is not None and row.get(col) is not None:
                    current[col] = max(current[col], row[col])
            current_sum_z += _abs_z(row["pvalue"])
            current_n_z += 1
            current["_count"] = prev_n + 1
        else:
            merged.append(_finalise(current, current_sum_z, current_n_z))
            current = dict(row)
            current["_count"] = 1
            current_sum_z = _abs_z(row["pvalue"])
            current_n_z = 1
    if current is not None:
        merged.append(_finalise(current, current_sum_z, current_n_z))

    for m in merged:
        m.pop("_count", None)

    if not merged:
        return dmr_df.clear()

    result = pl.DataFrame(merged, schema=dmr_df.schema)
    from .dmc import apply_multiple_testing_correction
    result = apply_multiple_testing_correction(result, method="fdr_bh")
    return result.sort(["chrom", "start"])


def call_dmr_tile_based(
    methylstore_path: str,
    samples_treatment: list[str] | None = None,
    samples_control: list[str] | None = None,
    tile_size_bp: int = 1000,
    test: str = "lr",
    chromosomes: list[str] | None = None,
    min_cpgs_per_tile: int = 5,
    alpha: float = 0.05,
    min_abs_meth_diff: float = 0.1,
    unite: bool = True,
    min_samples_treatment: int | None = None,
    min_samples_control: int = 0,
    dispersion: str = "site",
    reference: str = "chi2",
    design_full: np.ndarray | None = None,
    design_reduced: np.ndarray | None = None,
    coef_idx: int | None = None,
    *,
    backend: str = "sequential",
    n_workers: int | None = None,
    merge_adjacent: bool = True,
) -> pl.DataFrame:
    """Call DMRs by aggregating read counts within fixed-size tiles.

    Per-sample, per-chromosome, the engine sums N_meth and coverage across
    all CpGs in each tile, then runs a single tile-level DMC test. The
    sliding-window alternative tests each CpG individually and combines
    p-values, which has dramatically lower power at typical WGBS coverage:
    a tile with 20 CpGs at +15 % effect might have zero individually-
    significant CpGs but still trivially pass when its 600 pooled reads
    are tested.

    Implementation
    --------------
    1. For each sample and chromosome, aggregate (N_meth, coverage) per tile
       and write a "tiled methylstore" to a temp directory. Tiles with
       fewer than ``min_cpgs_per_tile`` CpGs (in that sample) are dropped
       before writing.
    2. Run ``process_chromosomes_dmc`` on the tiled store with the requested
       test. Tiles are treated as "sites" with pos = tile_start.
    3. BH-correct the tile-level p-values.
    4. Filter on qvalue and |meth_diff|.
    5. Reshape into the DMR output schema.

    Parameters
    ----------
    methylstore_path : str
        Path to the filtered partitioned Parquet methylstore.
    samples_treatment, samples_control : list[str]
        Sample IDs.
    tile_size_bp : int
        Tile width in bp (default 1000). Adjacent tiles do not overlap.
    test : str
        Statistical test for tile-level counts. Defaults to ``"lr"``
        (quasi-binomial likelihood-ratio), the recommended default
        when tile-level pooled counts are available.
    chromosomes : list[str], optional
        Chromosomes to process. Auto-detected when None.
    min_cpgs_per_tile : int
        Skip tiles with fewer than this many CpGs (per sample) during the
        per-sample aggregation step. Default 5 to reduce noise at sparse
        coverage.
    alpha : float
        q-value threshold for significance.
    min_abs_meth_diff : float
        Minimum |meth_diff| for a tile to be called significant.
    unite : bool
        If True (default), only test tiles covered in every sample.
    min_samples_treatment, min_samples_control : int
        Per-tile minimum number of samples required to be present in each
        group (only relevant when unite=False).

    Returns
    -------
    pl.DataFrame
        Columns: chrom, start, end, n_cpgs, n_case, n_control,
                 mean_beta_case, mean_beta_control, meth_diff,
                 log2_odds_ratio, pvalue, qvalue, dmr_type.
    """
    from .dmc import (
        process_chromosomes_dmc,
        apply_multiple_testing_correction,
    )

    if samples_treatment is None:
        raise TypeError("Missing required argument: samples_treatment")
    if samples_control is None:
        raise TypeError("Missing required argument: samples_control")
    if min_samples_treatment is None:
        min_samples_treatment = 0
    samples_case = samples_treatment
    min_samples_case = min_samples_treatment

    store       = Path(methylstore_path)
    all_samples = samples_case + samples_control

    if chromosomes is None:
        chromosomes = sorted({
            d.name.removeprefix("chrom=")
            for s in store.glob("sample=*")
            for d in s.glob("chrom=*")
        })

    if not chromosomes:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

    logger.info(
        "call_dmr_tile_based: tile=%d bp, test=%s, n_case=%d, n_control=%d, "
        "min_cpgs/tile=%d, alpha=%.3f, min_|Deltabeta|=%.2f, unite=%s",
        tile_size_bp, test, len(samples_case), len(samples_control),
        min_cpgs_per_tile, alpha, min_abs_meth_diff, unite,
    )

    with tempfile.TemporaryDirectory(prefix="epykit_tile_") as tmpdir:
        tile_store = Path(tmpdir) / "tiled_store"

        # ----- Phase 1: aggregate per-sample, per-chromosome counts -----
        # Track per-tile CpG counts (max across samples) for output column.
        # Stored as {(chrom, tile_start): n_cpgs}.
        tile_n_cpgs: dict[tuple[str, int], int] = {}

        for sample in all_samples:
            for chrom in chromosomes:
                src = store / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
                tiled = _aggregate_sample_to_tiles(src, chrom, tile_size_bp)
                if tiled is None or len(tiled) == 0:
                    continue
                tiled = tiled.filter(pl.col("n_cpgs") >= min_cpgs_per_tile)
                if len(tiled) == 0:
                    continue

                # Record per-tile CpG count (use max across samples so the
                # output reflects the most CpG-dense observation of the
                # tile).
                for tile_start, n_cpgs_val in zip(
                    tiled["pos"].to_list(), tiled["n_cpgs"].to_list()
                ):
                    key = (chrom, int(tile_start))
                    if n_cpgs_val > tile_n_cpgs.get(key, 0):
                        tile_n_cpgs[key] = int(n_cpgs_val)

                out_dir = tile_store / f"sample={sample}" / f"chrom={chrom}"
                out_dir.mkdir(parents=True, exist_ok=True)
                (
                    tiled
                    .select(["chrom", "pos", "strand", "N_meth", "N_unmeth", "coverage"])
                    .write_parquet(str(out_dir / "part-0.parquet"))
                )

        # ----- Phase 2: run DMC on the tiled store -----
        if not list(tile_store.glob("sample=*/chrom=*/part-0.parquet")):
            logger.warning("Tile aggregation produced no rows; returning empty DMR set")
            return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

        tile_dmc = process_chromosomes_dmc(
            methylstore_path=str(tile_store),
            samples_treatment=samples_case,
            samples_control=samples_control,
            test=test,
            chromosomes=chromosomes,
            unite=unite,
            min_samples_treatment=min_samples_case,
            min_samples_control=min_samples_control,
            dispersion=dispersion,
            reference=reference,
            design_full=design_full,
            design_reduced=design_reduced,
            coef_idx=coef_idx,
            backend=backend,
            n_workers=n_workers,
        )

    if len(tile_dmc) == 0:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

    # ----- Phase 3: BH at tile level -----
    tile_dmc = apply_multiple_testing_correction(tile_dmc, method="fdr_bh")

    # ----- Phase 4: filter and reshape -----
    # Attach n_cpgs from the per-sample aggregation.
    n_cpgs_rows = [
        {"chrom": c, "pos": p, "n_cpgs": n}
        for (c, p), n in tile_n_cpgs.items()
    ]
    n_cpgs_df = pl.DataFrame(
        n_cpgs_rows,
        schema={"chrom": pl.Utf8, "pos": pl.Int32, "n_cpgs": pl.Int32},
    )

    # P1-11: DMC output no longer has a real-valued 'log2_odds_ratio' column;
    # the backend-specific names are 'log2_odds_ratio_pooled' (lr/fisher) and
    # 'coef_treatment_log2' (glm).  Normalise to 'log2_odds_ratio' for the
    # DMR output schema (DMR schema rename is deferred to 0.8).
    _log2_src = (
        "coef_treatment_log2"
        if "coef_treatment_log2" in tile_dmc.columns
        else "log2_odds_ratio_pooled"
    )
    if _log2_src in tile_dmc.columns:
        tile_dmc = tile_dmc.with_columns(
            pl.col(_log2_src).alias("log2_odds_ratio")
        )

    dmr_df = (
        tile_dmc
        .join(n_cpgs_df, on=["chrom", "pos"], how="left")
        .with_columns(pl.col("n_cpgs").fill_null(0))
        .filter(
            (pl.col("qvalue") < alpha)
            & (pl.col("meth_diff").abs() >= min_abs_meth_diff)
            & (~pl.col("pvalue").is_nan())
        )
        .with_columns([
            pl.col("pos").alias("start"),
            (pl.col("pos") + tile_size_bp).cast(pl.Int32).alias("end"),
            pl.when(pl.col("meth_diff") > 0)
              .then(pl.lit("hyper"))
              .otherwise(pl.lit("hypo"))
              .alias("dmr_type"),
        ])
    )

    out_cols = [
        "chrom", "start", "end", "n_cpgs",
        "n_case", "n_control",
        "mean_beta_case", "mean_beta_control",
        "meth_diff", "log2_odds_ratio",
        "pvalue", "qvalue",
        "dmr_type",
    ]
    # GLM path adds adjusted log-odds effect size for the treatment coefficient.
    for extra in ("coef_treatment", "coef_se"):
        if extra in dmr_df.columns:
            out_cols.append(extra)
    dmr_df = dmr_df.select(out_cols).sort(["chrom", "start"])

    if merge_adjacent:
        dmr_df = _merge_adjacent_tiles(dmr_df)

    logger.info("Tile-based DMR: %s tiles -> %s significant DMRs",
                f"{len(tile_dmc):,}", f"{len(dmr_df):,}")

    gc.collect()
    return dmr_df


# Permutation-based empirical FDR

def empirical_fdr_for_dmr(
    methylstore_path: str,
    samples_treatment: list[str],
    samples_control: list[str],
    observed_dmr: pl.DataFrame,
    *,
    n_perm: int = 100,
    seed: int = 42,
    n_jobs: int = 1,
    empirical_strata: "dict[str, list[str]] | None" = None,
    **dmr_kwargs,
) -> pl.DataFrame:
    """Empirical (permutation) FDR for tile-based DMRs.

    Re-runs ``call_dmr_tile_based`` ``n_perm`` times with treatment / control
    labels shuffled. For each observed DMR, the empirical p-value is
    estimated from the fraction of null DMRs (across all permutations) with
    raw p-value <= observed. The result is BH-adjusted to
    ``empirical_qvalue``.

    Parameters
    ----------
    methylstore_path, samples_treatment, samples_control
        Same arguments passed to :func:`call_dmr_tile_based`.
    observed_dmr
        The DMR DataFrame returned by the observed (unpermuted) run.
        Empirical columns are appended to a copy of this frame.
    n_perm
        Number of permutations.
    seed
        Seed for the per-permutation label shuffler.
    n_jobs
        joblib parallel worker count. -1 uses all cores. Falls back to
        serial execution when joblib is not installed.
    empirical_strata : dict[str, list[str]] or None
        When supplied, a mapping from stratum label to the list of sample
        IDs belonging to that stratum.  Labels are shuffled **within** each
        stratum rather than globally.  Build this dict from ``md.obs`` in
        the caller (see :func:`epykit.tl.dmr`).  When ``None`` (default),
        the standard global shuffle is used.
    **dmr_kwargs
        Forwarded to ``call_dmr_tile_based`` for each permutation; should
        match the observed run's settings (tile_size_bp, test, alpha,
        min_abs_meth_diff, dispersion, reference, etc.).

    Returns
    -------
    pl.DataFrame
        ``observed_dmr`` with added columns ``empirical_pvalue`` and
        ``empirical_qvalue``. The full null pool (per-DMR raw pvalues from
        every permutation) is cached on ``observed_dmr.attrs`` only if the
        caller upstream wires it -- this function just returns the
        annotated table.
    """
    if len(observed_dmr) == 0:
        return observed_dmr.with_columns([
            pl.lit(None, dtype=pl.Float64).alias("empirical_pvalue"),
            pl.lit(None, dtype=pl.Float64).alias("empirical_qvalue"),
        ])

    n_treat = len(samples_treatment)
    n_ctrl = len(samples_control)
    if n_treat == 1 and n_ctrl == 1:
        raise ValueError(
            "empirical DMR FDR requires n>=2 per group; got n_treat=1, "
            "n_ctrl=1. Use Fisher-derived p-values directly via "
            "tl.dmc(test='fisher')."
        )

    pool = list(samples_treatment) + list(samples_control)
    rng = np.random.default_rng(seed)  # noqa: F841  (kept for reproducibility docs)

    def _run_one_perm(perm_idx: int) -> np.ndarray:
        # Local RNG so parallel workers stay deterministic.
        local_rng = np.random.default_rng(seed + perm_idx + 1)
        if empirical_strata is not None:
            # Shuffle within each stratum and re-assemble the pool.
            shuffled: list[str] = []
            for group_samples in empirical_strata.values():
                shuffled.extend(local_rng.permutation(group_samples).tolist())
        else:
            shuffled = pool.copy()
            local_rng.shuffle(shuffled)
        perm_treat = shuffled[:n_treat]
        perm_ctrl = shuffled[n_treat:]
        # Force test='lr' or whatever observed used; do not run annotation.
        kwargs = dict(dmr_kwargs)
        kwargs.pop("samples_case", None)
        kwargs.pop("min_samples_case", None)
        try:
            null_df = call_dmr_tile_based(
                methylstore_path=methylstore_path,
                samples_treatment=perm_treat,
                samples_control=perm_ctrl,
                **kwargs,
            )
        except Exception as exc:
            logger.warning("permutation %d failed: %s", perm_idx, exc)
            return np.array([], dtype=np.float64)
        if "pvalue" not in null_df.columns or len(null_df) == 0:
            return np.array([], dtype=np.float64)
        return null_df.get_column("pvalue").drop_nulls().to_numpy()

    null_pvals_list: list[np.ndarray]
    if n_jobs == 1:
        null_pvals_list = [_run_one_perm(i) for i in range(n_perm)]
    else:
        try:
            from joblib import Parallel, delayed
            null_pvals_list = Parallel(n_jobs=n_jobs)(
                delayed(_run_one_perm)(i) for i in range(n_perm)
            )
        except ImportError:
            logger.warning("joblib not installed; falling back to serial execution.")
            null_pvals_list = [_run_one_perm(i) for i in range(n_perm)]

    if all(len(arr) == 0 for arr in null_pvals_list):
        logger.warning(
            "All %d permutations produced zero null DMRs. Empirical "
            "p-values default to 1 / (1 + n_perm).",
            n_perm,
        )

    # Per-permutation tail count (max-T style):
    # For each observed DMR with p = p_obs, count the number of
    # permutations that produced at least one null DMR with p <= p_obs.
    # Denominator is n_perm + 1 (not |pooled null| + 1) so the empirical
    # p-value floor is 1 / (n_perm + 1), independent of how many DMRs
    # each permutation emitted. This is the standard region-statistic
    # FDR; pooling p-values across perms would inflate the denominator
    # and anti-conservatively shrink emp_p.
    obs_p = observed_dmr.get_column("pvalue").to_numpy()
    obs_finite_mask = np.isfinite(obs_p)
    obs_safe = np.where(obs_finite_mask, obs_p, 1.0)

    # min_null_p_per_perm[i] = min p across all null DMRs from perm i
    # (1.0 if perm produced no DMRs). The number of perms with at least
    # one null <= p_obs is then sum(min_null_p_per_perm <= p_obs).
    min_null_p_per_perm = np.array([
        float(arr.min()) if len(arr) > 0 else 1.0
        for arr in null_pvals_list
    ], dtype=np.float64)
    min_null_sorted = np.sort(min_null_p_per_perm)
    # For each obs p, count perms with min_null <= obs p.
    counts = np.searchsorted(min_null_sorted, obs_safe, side="right")
    emp_p = (counts + 1.0) / (n_perm + 1.0)
    emp_p = np.clip(emp_p, 0.0, 1.0)
    emp_p = np.where(obs_finite_mask, emp_p, np.nan)

    # BH-adjust to empirical q-value
    from statsmodels.stats.multitest import multipletests
    finite = np.isfinite(emp_p)
    emp_q = np.full_like(emp_p, np.nan, dtype=np.float64)
    if finite.any():
        _, q_finite, _, _ = multipletests(emp_p[finite], method="fdr_bh")
        emp_q[finite] = q_finite

    return observed_dmr.with_columns([
        pl.Series("empirical_pvalue", emp_p),
        pl.Series("empirical_qvalue", emp_q),
    ])


# BSmooth-style local-polynomial smoother (spec-faithful)

# Compiled-on-first-call helper. Numba is a core epykit dep but is otherwise
# unused; we import lazily so a numba-less debug install (uncommon) still
# falls back to a pure-numpy path.
_BSMOOTH_NJIT_FN = None


def _bsmooth_make_njit():
    """Build and cache the numba-compiled per-chrom BSmooth kernel."""
    global _BSMOOTH_NJIT_FN
    if _BSMOOTH_NJIT_FN is not None:
        return _BSMOOTH_NJIT_FN
    try:
        from numba import njit
    except ImportError:
        njit = None

    def _bsmooth_one_chrom(
        positions: np.ndarray,    # (n,) float64, sorted ascending
        n_meth:    np.ndarray,    # (n,) float64
        coverage:  np.ndarray,    # (n,) float64
        ns:        int,
        h_min:     float,
        degree:    int,
        min_cpgs_for_smooth: int,
    ) -> np.ndarray:
        """Local-polynomial smoother -- one chromosome, one sample.

        Per site i:
          * adaptive half-window h_i = max(distance to ns-th nearest CpG, h_min)
          * weights w_j = tricube(|x_j - x_i| / h_i) * coverage_j
          * weighted least squares of degree `degree` (1 or 2), centered at x_i
          * smoothed value = polynomial intercept, clipped to [0, 1]

        Sites with zero coverage anywhere in the window contribute zero weight.
        Sites with fewer than ``min_cpgs_for_smooth`` valid neighbors fall back
        to the raw beta.
        """
        n = positions.shape[0]
        out = np.full(n, np.nan, dtype=np.float64)

        # Raw beta (NaN where coverage == 0)
        beta_raw = np.empty(n, dtype=np.float64)
        for i in range(n):
            if coverage[i] > 0.0:
                beta_raw[i] = n_meth[i] / coverage[i]
            else:
                beta_raw[i] = np.nan

        for i in range(n):
            x_i = positions[i]

            # ---- 1. Find ns-th nearest CpG distance via two-pointer expand
            a = i
            b = i
            while (b - a + 1) < ns and (a > 0 or b < n - 1):
                d_left = x_i - positions[a - 1] if a > 0 else np.inf
                d_right = positions[b + 1] - x_i if b < n - 1 else np.inf
                if d_left <= d_right:
                    a -= 1
                else:
                    b += 1

            if (b - a + 1) >= ns:
                ns_dist = positions[b] - x_i
                if x_i - positions[a] > ns_dist:
                    ns_dist = x_i - positions[a]
            else:
                ns_dist = 0.0  # very small chrom; fall through to h_min

            h_i = h_min if ns_dist < h_min else ns_dist

            # ---- 2. Widen [lo, hi] to all CpGs within h_i of x_i
            lo = a
            hi = b
            while lo > 0 and (x_i - positions[lo - 1]) <= h_i:
                lo -= 1
            while hi < n - 1 and (positions[hi + 1] - x_i) <= h_i:
                hi += 1

            # ---- 3. Accumulate weighted moments
            s0 = 0.0
            s1 = 0.0
            s2 = 0.0
            s3 = 0.0
            s4 = 0.0
            t0 = 0.0   # X' W y
            t1 = 0.0
            t2 = 0.0
            n_valid = 0
            for j in range(lo, hi + 1):
                if not np.isfinite(beta_raw[j]):
                    continue
                t = positions[j] - x_i
                u = t / h_i
                if u < 0.0:
                    u = -u
                if u >= 1.0:
                    continue
                tri = 1.0 - u * u * u
                tri = tri * tri * tri
                w = tri * coverage[j]
                if w <= 0.0:
                    continue
                y = beta_raw[j]
                t2_loc = t * t
                s0 += w
                s1 += w * t
                s2 += w * t2_loc
                s3 += w * t2_loc * t
                s4 += w * t2_loc * t2_loc
                t0 += w * y
                t1 += w * t * y
                t2 += w * t2_loc * y
                n_valid += 1

            if n_valid < min_cpgs_for_smooth or s0 <= 0.0:
                out[i] = beta_raw[i]
                continue

            # ---- 4. Solve WLS for the intercept only (we don't need slope/curvature)
            if degree == 2:
                det = (
                    s0 * (s2 * s4 - s3 * s3)
                    - s1 * (s1 * s4 - s3 * s2)
                    + s2 * (s1 * s3 - s2 * s2)
                )
                if det == 0.0 or not np.isfinite(det):
                    out[i] = t0 / s0   # singular -> weighted mean fallback
                    continue
                num = (
                    t0 * (s2 * s4 - s3 * s3)
                    - s1 * (t1 * s4 - s3 * t2)
                    + s2 * (t1 * s3 - s2 * t2)
                )
                intercept = num / det
            else:  # degree == 1
                det = s0 * s2 - s1 * s1
                if det == 0.0 or not np.isfinite(det):
                    out[i] = t0 / s0
                    continue
                intercept = (t0 * s2 - s1 * t1) / det

            if intercept < 0.0:
                intercept = 0.0
            elif intercept > 1.0:
                intercept = 1.0
            out[i] = intercept

        return out

    if njit is not None:
        _BSMOOTH_NJIT_FN = njit(cache=True)(_bsmooth_one_chrom)
    else:
        _BSMOOTH_NJIT_FN = _bsmooth_one_chrom
    return _BSMOOTH_NJIT_FN


def smooth_methylation_bsmooth(
    methylstore_path: str,
    samples: list[str],
    *,
    ns: int = 70,
    h_bp: int = 1000,
    degree: int = 2,
    min_cpgs_for_smooth: int = 3,
    output_path: str | None = None,
) -> pl.DataFrame | None:
    """BSmooth-style local-polynomial smoother (Hansen et al. 2012).

    For each CpG, fits a local weighted-polynomial regression on the
    neighboring CpGs:

      * **Adaptive bandwidth**: ``h_i = max(distance to ns-th nearest CpG, h_bp)``.
        Sparse regions widen to capture ``ns`` CpGs; dense regions are
        capped below by ``h_bp`` so the kernel never collapses to
        immediate neighbors.
      * **Tricube distance weights x coverage**:
        ``w_j = (1 - (|x_j - x_i| / h_i)^3)^3 * N_j``.
      * **Polynomial degree** (``1`` or ``2``; default ``2`` matches BSmooth).
        Quadratic captures local curvature; linear is a faster fallback.
      * Smoothed value = polynomial intercept at the focal CpG, clipped
        to ``[0, 1]``.

    See :func:`smooth_methylation_gaussian` for the faster Gaussian-
    kernel approximation that previously occupied this slot.

    Performance: the per-site fit is compiled via ``numba.njit`` (numba
    is already an epykit core dep). First call incurs a one-time
    compilation latency (~1 s); subsequent calls run at native speed.

    Parameters
    ----------
    methylstore_path
        Path to the partitioned methylstore.
    samples
        Sample ids to smooth.
    ns
        Target number of CpGs per local window (BSmooth default 70).
    h_bp
        Minimum half-window in base pairs (BSmooth default 1000).
    degree
        Local-polynomial degree, ``1`` or ``2``. Default ``2``.
    min_cpgs_for_smooth
        Fall back to raw beta if fewer than this many valid neighbors
        contribute weight (default 3).
    output_path
        When set, write a sidecar parquet store at this root and return
        ``None``; otherwise concatenate everything and return a frame.

    Returns
    -------
    pl.DataFrame or None
        Long-form frame (chrom, pos, sample, beta_raw, beta_smooth) when
        ``output_path`` is ``None``, otherwise ``None`` after sidecar write.
    """
    if degree not in (1, 2):
        raise ValueError(f"degree must be 1 or 2; got {degree}")
    if ns < 2:
        raise ValueError(f"ns must be >= 2; got {ns}")
    if h_bp <= 0:
        raise ValueError(f"h_bp must be > 0; got {h_bp}")

    smoother = _bsmooth_make_njit()

    store = Path(methylstore_path)
    records: list[pl.DataFrame] = []
    out_root = Path(output_path) if output_path else None
    if out_root:
        out_root.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        sample_dir = store / f"sample={sample}"
        if not sample_dir.exists():
            logger.warning("Sample '%s' not found in %s; skipping", sample, store)
            continue

        for chrom_dir in sorted(sample_dir.glob("chrom=*")):
            chrom = chrom_dir.name.removeprefix("chrom=")
            parts = list(chrom_dir.glob("part-*.parquet"))
            if not parts:
                continue

            df = pl.concat([
                pl.read_parquet(str(p), columns=["pos", "N_meth", "coverage"])
                for p in parts
            ]).sort("pos")
            n = df.height
            if n == 0:
                continue

            pos = df["pos"].to_numpy().astype(np.float64)
            meth = df["N_meth"].to_numpy().astype(np.float64)
            cov = df["coverage"].to_numpy().astype(np.float64)

            beta_raw = np.where(cov > 0, meth / np.maximum(cov, 1.0), np.nan)
            if n < min_cpgs_for_smooth:
                logger.debug(
                    "  %s / %s: only %d sites; skipping bsmooth", sample, chrom, n,
                )
                beta_smooth = beta_raw.copy()
            else:
                beta_smooth = smoother(
                    pos, meth, cov,
                    int(ns), float(h_bp), int(degree), int(min_cpgs_for_smooth),
                )

            chunk = pl.DataFrame({
                "chrom":       pl.Series([chrom] * n, dtype=pl.Utf8),
                "pos":         df["pos"],
                "sample":      pl.Series([sample] * n, dtype=pl.Utf8),
                "beta_raw":    pl.Series(beta_raw.astype(np.float32)),
                "beta_smooth": pl.Series(beta_smooth.astype(np.float32)),
            })

            if out_root is not None:
                part_dir = out_root / f"sample={sample}" / f"chrom={chrom}"
                part_dir.mkdir(parents=True, exist_ok=True)
                chunk.write_parquet(
                    str(part_dir / "part-0.parquet"), compression="zstd",
                )
            else:
                records.append(chunk)

    if out_root is not None:
        return None

    if not records:
        return pl.DataFrame(schema=_SMOOTH_EMPTY_SCHEMA)
    return pl.concat(records)


# Public API -- fast Gaussian smoothing (replaces statsmodels LOESS)

def smooth_methylation_gaussian(
    methylstore_path: str,
    samples: list[str],
    bandwidth: int = 1000,
    grid_resolution_bp: int | None = None,
    output_path: str | None = None,
) -> pl.DataFrame | None:
    """Smooth per-sample beta values with a fast Gaussian kernel.

    .. note::
       This is a Gaussian-kernel approximation, not the local-LOESS smoother
       used in Hansen et al.'s BSmooth. A true LOESS-based smoother is on
       the roadmap.

    Within each chromosome and sample, raw beta values are smoothed along
    the genomic axis. The implementation projects raw betas onto a regular
    grid, applies a coverage-weighted Gaussian filter
    (``scipy.ndimage.gaussian_filter1d``), then interpolates back to the
    original CpG positions. This is O(G) where G is the grid size,
    versus O(n^2) for LOESS -- typically 100-500x faster on WGBS-scale data.
    """
    try:
        from scipy.ndimage import gaussian_filter1d
    except ImportError as exc:
        raise ImportError(
            "scipy is required for Gaussian smoothing. "
            "Install with: pip install scipy"
        ) from exc

    store   = Path(methylstore_path)
    records: list[pl.DataFrame] = []

    # Determine grid resolution once (same for all samples/chroms)
    _grid_res = max(1, bandwidth // 20) if grid_resolution_bp is None else grid_resolution_bp
    _out_root = Path(output_path) if output_path else None
    if _out_root:
        _out_root.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        sample_dir = store / f"sample={sample}"
        if not sample_dir.exists():
            logger.warning("Sample '%s' not found in %s; skipping", sample, store)
            continue

        for chrom_dir in sorted(sample_dir.glob("chrom=*")):
            chrom = chrom_dir.name.removeprefix("chrom=")
            parts = list(chrom_dir.glob("part-*.parquet"))
            if not parts:
                continue

            df = pl.concat([
                pl.read_parquet(str(p), columns=["pos", "N_meth", "coverage"])
                for p in parts
            ]).sort("pos")

            pos  = df["pos"].to_numpy().astype(np.float64)
            meth = df["N_meth"].to_numpy().astype(np.float64)
            cov  = df["coverage"].to_numpy().astype(np.float64)

            with np.errstate(invalid="ignore", divide="ignore"):
                beta_raw = np.where(cov > 0, meth / cov, np.nan).astype(np.float32)

            beta_smooth = beta_raw.copy()
            valid       = ~np.isnan(beta_raw)
            n_valid     = int(valid.sum())

            if n_valid >= 4:
                pos_valid   = pos[valid]
                beta_valid  = beta_raw[valid].astype(np.float64)

                # Build a regular grid spanning the valid positions.
                grid_start  = int(pos_valid[0])
                grid_end    = int(pos_valid[-1]) + _grid_res
                grid_pos    = np.arange(grid_start, grid_end, _grid_res,
                                        dtype=np.float64)

                # Coverage-weighted interpolation onto the regular grid
                cov_valid = cov[valid].astype(np.float64)
                grid_beta = np.interp(grid_pos, pos_valid, beta_valid)
                grid_weights = np.interp(grid_pos, pos_valid, cov_valid)
                grid_weights = np.maximum(grid_weights, 0.1)  # avoid exact zeros

                # Apply weighted Gaussian: smooth numerator and denominator separately
                sigma_grid = max(bandwidth / _grid_res, 0.5)
                grid_num = gaussian_filter1d(grid_beta * grid_weights, sigma=sigma_grid, mode="nearest")
                grid_den = gaussian_filter1d(grid_weights, sigma=sigma_grid, mode="nearest")
                smoothed_grid = grid_num / np.maximum(grid_den, 1e-9)
                np.clip(smoothed_grid, 0.0, 1.0, out=smoothed_grid)

                # Interpolate smoothed values back to original CpG positions
                smoothed_at_cpgs = np.interp(pos_valid, grid_pos, smoothed_grid)
                beta_smooth[valid] = smoothed_at_cpgs.astype(np.float32)
            else:
                logger.debug(
                    "  %s / %s: only %d valid sites; skipping smoothing",
                    sample, chrom, n_valid,
                )

            chunk = pl.DataFrame({
                "chrom":       pl.Series([chrom]  * len(df), dtype=pl.Utf8),
                "pos":         df["pos"],
                "sample":      pl.Series([sample] * len(df), dtype=pl.Utf8),
                "beta_raw":    pl.Series(beta_raw),
                "beta_smooth": pl.Series(beta_smooth),
            })

            if _out_root is not None:
                part_dir = _out_root / f"sample={sample}" / f"chrom={chrom}"
                part_dir.mkdir(parents=True, exist_ok=True)
                chunk.write_parquet(
                    str(part_dir / "part-0.parquet"), compression="zstd"
                )
            else:
                records.append(chunk)

    if _out_root is not None:
        return None

    if not records:
        return pl.DataFrame(schema=_SMOOTH_EMPTY_SCHEMA)

    return pl.concat(records).sort(["chrom", "pos", "sample"])


