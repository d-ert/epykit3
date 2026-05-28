"""_epykit_scoring.py -- shared epykit benchmark scoring helpers.

Pure helpers used by both ``run_epykit_study1.py`` (Phase 4 Task 2,
Piao-as-distributed) and ``run_epykit_simulator.py`` (Phase 4 Task 3,
held-out intrinsic-truth simulator). Single-sourcing this contract is what
prevents the two studies from drifting on:

* the lr+ ``pvalue_combined`` / ``qvalue_combined`` selection rule,
* the stale-engine filter on reassembly,
* the DMC threshold grid (``P_THRESHOLDS``, ``Q_THRESHOLD``),
* the meth_diff binning + DMR overlap threshold,
* the narrow exception clause used in per-engine loops (only catch what
  the engines actually raise; let ``MemoryError`` propagate),
* the AUROC all-equal-scores edge case (NaN, not 0.5).

This module is **pure**: no module-level paths, no I/O against a fixed
location, no global state. I/O paths and the cell-grid driver stay in
the consumer scripts. Threshold *constants* live here because they are
part of the scoring contract.
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Threshold constants (scoring contract)
# ---------------------------------------------------------------------------

# DMC threshold grid (matches legacy evaluate.py).
P_THRESHOLDS = (0.001, 0.005, 0.01, 0.05)
Q_THRESHOLD = 0.05
METH_DIFF_BINS = ("0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
DMR_OVERLAP_THRESHOLD = 0.8


# ---------------------------------------------------------------------------
# Engine-loop narrow exception clause
# ---------------------------------------------------------------------------

# What the surviving DMC/DMR engines actually raise. MemoryError and
# KeyboardInterrupt MUST propagate -- they are not "expected" engine
# failures. ``BLE001`` flagged the prior ``except Exception`` as too broad.
ENGINE_EXCEPTIONS: tuple[type[BaseException], ...] = (
    RuntimeError,
    ValueError,
    FileNotFoundError,
    pl.exceptions.PolarsError,
    ArithmeticError,
    OSError,
)


# ---------------------------------------------------------------------------
# lr+ power-stack kwargs (single source of truth for the contract)
# ---------------------------------------------------------------------------


def _dmc_kwargs(test: str, *, allow_n1: bool) -> tuple[str, dict]:
    """Translate a logical test name to (backend_test, ep.tl.dmc kwargs).

    Returns ``backend_test`` so the caller knows which
    ``md.varm['dmc_<test>']`` key to read after the call -- lr+ uses the lr
    backend with the README/CLAUDE.md-validated power stack on top:

        fdr_method='fdr_tsbh', neighbour_combine=True, sep_fallback=True,
        dispersion='eb', neighbour_bp=500, sep_threshold=0.9.

    Any change to this dict is a benchmark-contract change and must
    re-run the relevant ablations under ``benchmark/scripts/ab_*.py``.
    """
    if test == "lr+":
        kwargs = dict(
            test="lr",
            allow_n1=allow_n1,
            neighbour_combine=True,
            neighbour_bp=500,  # default per CLAUDE.md guidance
            sep_fallback=True,
            sep_threshold=0.9,
            fdr_method="fdr_tsbh",
            dispersion="eb",
        )
        return "lr", kwargs
    return test, dict(test=test, allow_n1=allow_n1)


# ---------------------------------------------------------------------------
# Scoring primitives
# ---------------------------------------------------------------------------


def _confusion(joined: pl.DataFrame, sig_col: str) -> dict:
    """TP/FP/TN/FN + TPR/FPR/precision/F1 from a joined-with-truth frame.

    ``joined`` must carry the bool columns ``sig_col`` (predicted positive)
    and ``is_dmc`` (truth). All four counts are returned alongside the
    derived rates so Wilson CIs downstream have integers to binomtest on.
    """
    df = joined.with_columns(pred=pl.col(sig_col), truth=pl.col("is_dmc"))
    tp = df.filter(pl.col("pred") & pl.col("truth")).height
    fp = df.filter(pl.col("pred") & ~pl.col("truth")).height
    fn = df.filter(~pl.col("pred") & pl.col("truth")).height
    tn = df.filter(~pl.col("pred") & ~pl.col("truth")).height
    n_pos = tp + fn
    n_neg = fp + tn
    n_pred = tp + fp
    return {
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "tpr": tp / n_pos if n_pos else 0.0,
        "fpr": fp / n_neg if n_neg else 0.0,
        "precision": tp / n_pred if n_pred else 0.0,
        "f1": (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) else 0.0,
    }


def _auroc(joined: pl.DataFrame) -> float:
    """AUROC from ranked ``1 - pvalue`` scores against ``is_dmc``.

    Returns NaN in three degenerate cases:
      * no positives (``n_pos == 0``),
      * no negatives (``n_neg == 0``),
      * all p-values identical (``pvalue.n_unique() <= 1``) -- the
        ranking is uninformative, so reporting 0.5 would mask the
        degeneracy. This guard was added in Phase 4 Task 3 to prevent
        silent misreporting when an engine returns a constant p-value
        (e.g. a pathological all-zero or all-one case on tiny data).
    """
    df = joined.select("is_dmc", "pvalue").with_columns(
        score=1.0 - pl.col("pvalue").fill_null(1.0),
    )
    n_pos = df.filter(pl.col("is_dmc")).height
    n_neg = df.height - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # All-equal-scores guard: an engine that emitted a constant p-value
    # produces a degenerate ranking; AUROC is undefined here. Report
    # NaN rather than the misleading 0.5 from the tied-rank formula.
    if df["pvalue"].fill_null(1.0).n_unique() <= 1:
        return float("nan")
    df = df.with_columns(rank=pl.col("score").rank(method="average"))
    sum_ranks_pos = df.filter(pl.col("is_dmc"))["rank"].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2
    return u / (n_pos * n_neg)


def _join_with_truth(epy_df: pl.DataFrame, truth: pl.DataFrame) -> pl.DataFrame:
    """Project pvalue/qvalue columns and inner-join to the truth table.

    For lr+ (neighbour_combine=True) the raw per-CpG p-value stays in
    ``pvalue`` by Phase-3 contract, and the neighbour-combined value
    lives in ``pvalue_combined`` / ``qvalue_combined``. lr+ is *defined*
    as the combined call, so when those columns are present we use them
    as the canonical pvalue/qvalue for scoring (otherwise downstream
    metrics ignore the power-stack improvement and lr+ collapses back
    to lr).
    """
    pcol = "pvalue_combined" if "pvalue_combined" in epy_df.columns else "pvalue"
    qcol = "qvalue_combined" if "qvalue_combined" in epy_df.columns else "qvalue"
    projected = epy_df.select(
        ["chrom", "pos"]
        + ([pcol] if pcol != "pvalue" else ["pvalue"])
        + ([qcol] if qcol != "qvalue" else ["qvalue"])
    )
    if pcol != "pvalue":
        projected = projected.rename({pcol: "pvalue"})
    if qcol != "qvalue":
        projected = projected.rename({qcol: "qvalue"})
    projected = projected.with_columns(
        pl.col("pvalue").cast(pl.Float64),
        pl.col("qvalue").cast(pl.Float64),
    )
    return truth.join(projected, on=["chrom", "pos"], how="left")


def score_dmc_parquet(
    parquet: Path, truth: pl.DataFrame,
    tool: str, scenario: str, parameter: str,
    parameter_value, test: str,
) -> list[dict]:
    """Score one per-cell DMC parquet against the truth table.

    Emits the full threshold grid (P_THRESHOLDS pvalue rows + one
    Q_THRESHOLD qvalue row, all on ``meth_diff_bin='all'``) plus per-bin
    TPR rows stratified by ``meth_diff_bin``. AUROC is computed once on
    the all-bins join and replicated on the all-bins rows; per-bin rows
    carry ``auroc=NaN`` since AUROC is defined over the full ROC.

    ``truth`` must carry the columns ``chrom``, ``pos``, ``is_dmc``, and
    ``meth_diff_bin``. The latter is optional in the per-bin loop --
    rows whose ``meth_diff_bin`` does not match any of ``METH_DIFF_BINS``
    are silently treated as not in any bin (no per-bin row emitted for
    them).
    """
    df = pl.read_parquet(parquet)
    if "pvalue" not in df.columns or "qvalue" not in df.columns:
        logger.warning("skip (no pvalue/qvalue): %s", parquet.name)
        return []
    joined = _join_with_truth(df, truth)

    rows: list[dict] = []
    auroc = _auroc(joined)

    # All-bins p-value thresholds
    for cut in P_THRESHOLDS:
        m = _confusion(joined.with_columns(sig=pl.col("pvalue") < cut), "sig")
        rows.append({
            "tool": tool, "scenario": scenario, "parameter": parameter,
            "parameter_value": parameter_value, "test": test,
            "meth_diff_bin": "all", "threshold_kind": "pvalue",
            "threshold": cut, **m, "auroc": auroc,
        })

    # All-bins q-value @ 0.05
    m = _confusion(joined.with_columns(sig=pl.col("qvalue") < Q_THRESHOLD), "sig")
    rows.append({
        "tool": tool, "scenario": scenario, "parameter": parameter,
        "parameter_value": parameter_value, "test": test,
        "meth_diff_bin": "all", "threshold_kind": "qvalue",
        "threshold": Q_THRESHOLD, **m, "auroc": auroc,
    })

    # Per-bin TPR stratified -- only meaningful when truth carries a
    # meth_diff_bin column. The simulator's truth does; some legacy
    # tables may not.
    if "meth_diff_bin" not in joined.columns:
        return rows
    for bin_label in METH_DIFF_BINS:
        sub = joined.with_columns(
            is_dmc_in_bin=pl.col("is_dmc") & (pl.col("meth_diff_bin") == bin_label),
            sig=pl.col("qvalue") < Q_THRESHOLD,
        )
        tp = sub.filter(pl.col("sig") & pl.col("is_dmc_in_bin")).height
        fn = sub.filter(~pl.col("sig") & pl.col("is_dmc_in_bin")).height
        fp_g = sub.filter(pl.col("sig") & ~pl.col("is_dmc")).height
        tn_g = sub.filter(~pl.col("sig") & ~pl.col("is_dmc")).height
        n_pos = tp + fn
        n_neg = fp_g + tn_g
        rows.append({
            "tool": tool, "scenario": scenario, "parameter": parameter,
            "parameter_value": parameter_value, "test": test,
            "meth_diff_bin": bin_label, "threshold_kind": "qvalue",
            "threshold": Q_THRESHOLD,
            "tp": tp, "fp": fp_g, "tn": tn_g, "fn": fn,
            "tpr": tp / n_pos if n_pos else 0.0,
            "fpr": fp_g / n_neg if n_neg else 0.0,
            "precision": tp / (tp + fp_g) if (tp + fp_g) else 0.0,
            "f1": (2 * tp / (2 * tp + fp_g + fn)) if (2 * tp + fp_g + fn) else 0.0,
            "auroc": float("nan"),
        })
    return rows


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge a list of half-open ``[start, end)`` intervals."""
    if not intervals:
        return []
    sorted_iv = sorted(intervals)
    merged = [sorted_iv[0]]
    for s, e in sorted_iv[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def score_dmr_parquet(
    parquet: Path, dmr_truth: pl.DataFrame,
    tool: str, scenario: str, parameter_value, method: str,
    min_overlap: float = DMR_OVERLAP_THRESHOLD,
) -> list[dict]:
    """Same overlap-fraction DMR scoring as the legacy evaluate.py.

    A truth region is "detected" when the union of called intervals
    covers at least ``min_overlap`` of its length. ``n_false`` is the
    count of called regions that fail to overlap *any* truth region
    (note: this is a slight asymmetry vs the per-truth detection rule;
    flagged for a separate follow-up but left as-is for backwards
    compat with the legacy summary).
    """
    df = pl.read_parquet(parquet)
    if "qvalue" in df.columns:
        called = df.filter(pl.col("qvalue") < Q_THRESHOLD)
    else:
        called = df

    truth_rows = dmr_truth.to_dicts()
    called_rows = called.select(["chrom", "start", "end"]).to_dicts()

    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for r in called_rows:
        by_chrom.setdefault(r["chrom"], []).append((int(r["start"]), int(r["end"])))

    n_detected = 0
    overlap_fracs: list[float] = []
    for t in truth_rows:
        t_start, t_end = int(t["start"]), int(t["end"])
        t_len = max(1, t_end - t_start)
        calls = by_chrom.get(t["chrom"], [])
        clipped = [
            (max(t_start, s), min(t_end, e))
            for s, e in calls
            if min(t_end, e) > max(t_start, s)
        ]
        merged = _merge_intervals(clipped)
        union_cov = sum(e - s for s, e in merged)
        frac = union_cov / t_len
        overlap_fracs.append(frac)
        if frac >= min_overlap:
            n_detected += 1

    n_called = called.height
    n_truth = len(truth_rows)

    n_false = 0
    if called_rows and truth_rows:
        truth_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for t in truth_rows:
            truth_by_chrom.setdefault(t["chrom"], []).append(
                (int(t["start"]), int(t["end"]))
            )
        for r in called_rows:
            s, e = int(r["start"]), int(r["end"])
            hit = False
            for ts, te in truth_by_chrom.get(r["chrom"], []):
                if min(e, te) > max(s, ts):
                    hit = True
                    break
            if not hit:
                n_false += 1
    else:
        n_false = n_called

    tpr = n_detected / n_truth if n_truth else 0.0
    precision = (n_called - n_false) / n_called if n_called else 0.0
    f1 = (2 * tpr * precision / (tpr + precision)) if (tpr + precision) else 0.0

    return [{
        "tool": tool, "scenario": scenario, "parameter": "coverage",
        "parameter_value": parameter_value, "test": method,
        "meth_diff_bin": "all", "threshold_kind": "dmr_overlap",
        "threshold": float(min_overlap),
        "tp": None, "fp": None, "tn": None, "fn": None,
        "tpr": tpr, "fpr": float("nan"),
        "precision": precision, "f1": f1, "auroc": float("nan"),
    }]


# ---------------------------------------------------------------------------
# CI helper (has/no-counts split workaround)
# ---------------------------------------------------------------------------


def split_for_ci(df: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Split eval rows into (with_counts, without_counts) for CI handling.

    ``evaluate.py --ci-only`` requires non-null integer tp/fp/tn/fn so
    Wilson ``binomtest`` has counts to operate on. DMR rows (and any
    transcribed baseline rows that lack confusion-matrix entries) carry
    nulls in those columns; running ``--ci-only`` over them poisons the
    output with NaNs and aborts the bootstrap loop.

    The pattern (originally in run_epykit_study1.py) is: split off the
    rows with counts, run ``--ci-only`` on those alone, then concatenate
    the no-counts rows back with NaN CI columns aligned to the with-CI
    schema. This helper is the split half; the concat half is small
    enough to inline in the caller.
    """
    has_counts = df.filter(
        pl.col("tp").is_not_null()
        & pl.col("fp").is_not_null()
        & pl.col("tn").is_not_null()
        & pl.col("fn").is_not_null()
    )
    no_counts = df.filter(
        pl.col("tp").is_null()
        | pl.col("fp").is_null()
        | pl.col("tn").is_null()
        | pl.col("fn").is_null()
    )
    return has_counts, no_counts


# ---------------------------------------------------------------------------
# Stale-engine filter (reassembly)
# ---------------------------------------------------------------------------


STALE_EPYKIT_TOOLS: frozenset[str] = frozenset({
    # Phase 3 removed engines; rows are dropped on reassembly.
    "epykit_bb_lr",
    # Old-runner tool labels that we now replace with fresh post-Phase-3 ones.
    "epykit_lr", "epykit_lrplus", "epykit_welch_t", "epykit_fisher",
    "epykit_dmr_tile", "epykit_dmr_merge",
    # New DMR runners (in case of partial reruns).
    "epykit_dmr_chain_merge", "epykit_dmr_sliding_window", "epykit_dmr_segment",
})


def reassemble_eval_summary(
    new_rows: list[dict],
    eval_summary_old_path: Path,
) -> pl.DataFrame:
    """Concat fresh epykit rows with the non-epykit baseline rows from
    a pre-existing eval_summary.parquet.

    Drops every row whose ``tool`` is in ``STALE_EPYKIT_TOOLS`` or
    starts with ``epykit_`` (defensive: catches any post-fix label we
    forgot to enumerate). Schemas are aligned by filling missing
    columns with nulls before vertical concat.

    ``eval_summary_old_path`` is taken as a parameter (not a module
    constant) so the same helper works for both Study 1 (Piao) and any
    other study that wants the same reassembly contract.
    """
    old = pl.read_parquet(eval_summary_old_path)
    non_epykit = old.filter(~pl.col("tool").is_in(list(STALE_EPYKIT_TOOLS)))
    # Also explicitly drop anything still starting with epykit_ to be safe.
    non_epykit = non_epykit.filter(~pl.col("tool").str.starts_with("epykit_"))

    if not new_rows:
        logger.warning("no new epykit rows -- writing baseline-only summary")
        return non_epykit

    new_df = pl.DataFrame(new_rows)
    all_cols = sorted(set(non_epykit.columns) | set(new_df.columns))
    for c in all_cols:
        if c not in non_epykit.columns:
            non_epykit = non_epykit.with_columns(pl.lit(None).alias(c))
        if c not in new_df.columns:
            new_df = new_df.with_columns(pl.lit(None).alias(c))
    combined = pl.concat(
        [non_epykit.select(all_cols), new_df.select(all_cols)],
        how="vertical_relaxed",
    )
    return combined
