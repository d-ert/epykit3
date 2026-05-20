"""Shared pytest fixtures and accuracy helpers for the epykit test suite.

The synthetic dataset is built once per pytest session (medium size: 8
samples x 5 chromosomes x ~10 000 CpGs, ~500 scattered DMCs and 10 seeded
DMRs) and reused by every test that needs read data. Per-test isolation is
preserved by giving each test its own fresh ``MethylData`` instance pointing
at the same on-disk methylstore.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl
import pytest

from tests.fixtures.synth import SimConfig, generate


# Force matplotlib's Agg backend before any test imports plotting code, so
# CI without a display still runs the pl.* smokes.
import matplotlib
matplotlib.use("Agg", force=True)


@dataclass
class SynthBundle:
    """Bundle of paths + truth table + ids passed around by fixtures."""

    samplesheet: str
    truth: pl.DataFrame
    store_root: str
    treatment_ids: list[str]
    control_ids: list[str]
    n_dmcs_true: int
    n_dmrs: int
    config: SimConfig


@pytest.fixture(scope="session")
def synth_bundle(tmp_path_factory) -> SynthBundle:
    """Generate the Bismark .cov fixture once per session."""
    cfg = SimConfig()  # default = medium fixture (8 x 5 x 2000 = 10k sites)
    out_dir = tmp_path_factory.mktemp("synth")
    result = generate(cfg, out_dir)
    truth = pl.read_parquet(result["truth"])
    return SynthBundle(
        samplesheet=result["samplesheet"],
        truth=truth,
        store_root=str(out_dir / "methyl_store"),
        treatment_ids=result["treatment_ids"],
        control_ids=result["control_ids"],
        n_dmcs_true=result["n_dmcs_true"],
        n_dmrs=result["n_dmrs"],
        config=cfg,
    )


@pytest.fixture
def synth_md(synth_bundle: SynthBundle, tmp_path):
    """Fresh MethylData pointing at the session methylstore.

    Each test gets its own object so mutations to ``md.uns`` / ``md.varm``
    don't leak across tests. The underlying Parquet methylstore is reused.
    """
    import epykit as ep
    # Each test points read_bismark at a per-test store_dir so the cache
    # write doesn't collide with other tests running in parallel.
    return ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )


@pytest.fixture
def synth_md_filtered(synth_md):
    """MethylData that has been filter_coverage'd; ready for DMC."""
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="intersect")
    return synth_md



# Accuracy metric helpers


def join_truth(dmc_df: pl.DataFrame, truth: pl.DataFrame) -> pl.DataFrame:
    """Left-join DMC results onto truth on (chrom, pos).

    epykit's DMC output stores pos as Int32; the synthetic truth uses Int64.
    Cast both to Int64 before joining so polars doesn't refuse the join.
    """
    truth_cast = truth.with_columns(pl.col("pos").cast(pl.Int64))
    dmc_cast   = dmc_df.with_columns(pl.col("pos").cast(pl.Int64))
    return truth_cast.join(dmc_cast, on=["chrom", "pos"], how="left")


def power_at_threshold(
    dmc_df: pl.DataFrame,
    truth: pl.DataFrame,
    alpha: float = 0.05,
    p_col: str = "qvalue",
) -> float:
    """Fraction of truly-DMC sites called significant.

    Power = P(reject | H1 true) = TP / (TP + FN).
    """
    joined = join_truth(dmc_df, truth)
    truly_dmc = joined.filter(pl.col("is_dmc"))
    if len(truly_dmc) == 0:
        return float("nan")
    if p_col not in truly_dmc.columns:
        # Fall back to pvalue if qvalue absent (e.g. some engines)
        p_col = "pvalue"
    called = truly_dmc.filter(pl.col(p_col) < alpha).height
    return called / len(truly_dmc)


def fdr_at_threshold(
    dmc_df: pl.DataFrame,
    truth: pl.DataFrame,
    alpha: float = 0.05,
    p_col: str = "qvalue",
) -> float:
    """Empirical false-discovery rate among significant sites.

    FDR = FP / (FP + TP). Returns 0.0 when no sites are called.
    """
    joined = join_truth(dmc_df, truth)
    if p_col not in joined.columns:
        p_col = "pvalue"
    called = joined.filter(pl.col(p_col) < alpha)
    if len(called) == 0:
        return 0.0
    fp = called.filter(~pl.col("is_dmc")).height
    return fp / len(called)


def meth_diff_bias(
    dmc_df: pl.DataFrame,
    truth: pl.DataFrame,
) -> tuple[float, float]:
    """Mean signed bias and mean absolute error of recovered Deltabeta on true DMCs.

    Returns (mean(estimated - true), mean(|estimated - true|)).
    """
    joined = join_truth(dmc_df, truth)
    dmcs = joined.filter(pl.col("is_dmc"))
    if len(dmcs) == 0 or "meth_diff" not in dmcs.columns:
        return float("nan"), float("nan")
    # epykit reports meth_diff as treatment - control; truth uses the same sign convention.
    est = dmcs["meth_diff"].to_numpy()
    tru = dmcs["true_meth_diff"].to_numpy()
    mask = np.isfinite(est) & np.isfinite(tru)
    if not mask.any():
        return float("nan"), float("nan")
    diff = est[mask] - tru[mask]
    return float(diff.mean()), float(np.abs(diff).mean())


def dmr_recovery(
    dmr_df: pl.DataFrame,
    truth: pl.DataFrame,
    cfg: SimConfig,
    alpha: float = 0.05,
    q_col: Optional[str] = None,
) -> tuple[int, int]:
    """Number of seeded DMRs recovered by at least one significant call.

    A seeded DMR is "recovered" if any DMR call in ``dmr_df`` overlaps any
    CpG of the seed and has q-value < alpha. Returns
    (n_recovered, n_seeded).
    """
    if len(dmr_df) == 0:
        return 0, cfg.n_dmrs

    # Pick q-column heuristically: tile path uses 'qvalue', sliding-window
    # uses 'combined_qvalue'.
    if q_col is None:
        if "qvalue" in dmr_df.columns:
            q_col = "qvalue"
        elif "combined_qvalue" in dmr_df.columns:
            q_col = "combined_qvalue"
        else:
            q_col = "combined_pvalue"

    sig_dmrs = dmr_df.filter(pl.col(q_col) < alpha)
    if len(sig_dmrs) == 0:
        return 0, cfg.n_dmrs

    # Build per-DMR seed intervals from truth: (chrom, min_pos, max_pos).
    seed_intervals = (
        truth.filter(pl.col("in_dmr"))
        .group_by("dmr_id")
        .agg([
            pl.col("chrom").first().alias("chrom"),
            pl.col("pos").min().alias("seed_start"),
            pl.col("pos").max().alias("seed_end"),
        ])
        .sort("dmr_id")
    )

    # For each seeded DMR, ask if any called DMR overlaps it on the same
    # chromosome. We loop in Python -- n_dmrs <= 20 in tests so this is fine.
    recovered = 0
    sig_rows = sig_dmrs.to_dicts()
    for seed in seed_intervals.to_dicts():
        s_chrom = seed["chrom"]
        s_lo, s_hi = seed["seed_start"], seed["seed_end"]
        for call in sig_rows:
            if call.get("chrom") != s_chrom:
                continue
            c_lo = call.get("start", call.get("pos"))
            c_hi = call.get("end", call.get("pos"))
            if c_lo is None or c_hi is None:
                continue
            if c_lo <= s_hi and c_hi >= s_lo:
                recovered += 1
                break
    return recovered, cfg.n_dmrs
