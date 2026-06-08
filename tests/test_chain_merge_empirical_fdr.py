"""Integration tests for empirical_fdr_for_chain_merge.

The per-permutation engine `_chain_merge_perm_survivors` (which recomputes the
full per-CpG DMC then chain-merges -- too heavy for a unit test) is monkeypatched
to return fixed null survivor pools, so the harness + count-ratio aggregation is
exercised end-to-end. Mirrors test_dmr_region_fdr_mode.py for the tile path.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from epykit import dmr as dmr_mod

_SCHEMA = {
    "chrom": pl.Utf8, "start": pl.Int32, "end": pl.Int32,
    "n_cpgs": pl.Int32, "mean_meth_diff": pl.Float64,
    "combined_pvalue": pl.Float64,
}


def _observed(n=40, p=1e-6):
    return pl.DataFrame({
        "chrom": ["chr1"] * n,
        "start": list(range(0, n * 1000, 1000)),
        "end": list(range(500, n * 1000 + 500, 1000)),
        "n_cpgs": [8] * n,
        "mean_meth_diff": [0.3] * n,
        "combined_pvalue": [p] * n,
    }, schema=_SCHEMA)


def test_chain_merge_region_mode_significant(monkeypatch):
    """40 observed regions at combined_pvalue=1e-6; each shuffle yields only
    2 extreme null regions -> count-ratio q=2/40=0.05-ish, significant-ish."""
    monkeypatch.setattr(
        dmr_mod, "_chain_merge_perm_survivors",
        lambda **kw: np.array([1e-6, 1e-6], dtype=np.float64),
    )
    out = dmr_mod.empirical_fdr_for_chain_merge(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2", "t3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmr=_observed(40, 1e-6),
        dmc_kwargs={"test": "lr"},
        chain_merge_kwargs={"preset": "default"},
        min_mean_qvalue=0.05,
        n_perm=20, seed=0, n_jobs=1, fdr_method="region",
    )
    q = out.get_column("empirical_qvalue").to_numpy()
    assert "empirical_pvalue" in out.columns
    assert "empirical_fdr_set" in out.columns
    # V(1e-6)=2 per shuffle, R=40 -> set FDR = 2/40 = 0.05; gradient <= that.
    assert out.get_column("empirical_fdr_set")[0] == pytest.approx(0.05, abs=1e-6)
    assert np.all(q <= 0.05 + 1e-9)


def test_chain_merge_uses_combined_pvalue_not_pvalue(monkeypatch):
    """The harness must read the region statistic from `combined_pvalue`."""
    monkeypatch.setattr(
        dmr_mod, "_chain_merge_perm_survivors",
        lambda **kw: np.array([0.5, 0.5], dtype=np.float64),
    )
    out = dmr_mod.empirical_fdr_for_chain_merge(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2", "t3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmr=_observed(10, 1e-8),   # very extreme observed
        dmc_kwargs={"test": "lr"},
        chain_merge_kwargs={"preset": "default"},
        min_mean_qvalue=0.05,
        n_perm=10, seed=0, n_jobs=1,
    )
    # observed p=1e-8 << all null (0.5) -> strongest region q ~ 0
    assert np.nanmin(out.get_column("empirical_qvalue").to_numpy()) < 0.01


def test_chain_merge_runtime_and_small_n_warn(monkeypatch):
    """chain_merge permutation recomputes the whole DMC per shuffle -> must
    warn about runtime; <4/group also warns about underpowered permutation."""
    monkeypatch.setattr(
        dmr_mod, "_chain_merge_perm_survivors",
        lambda **kw: np.array([1e-6], dtype=np.float64),
    )
    with pytest.warns(UserWarning):
        dmr_mod.empirical_fdr_for_chain_merge(
            methylstore_path="/dev/null",
            samples_treatment=["t1", "t2", "t3"],
            samples_control=["c1", "c2", "c3"],
            observed_dmr=_observed(10, 1e-6),
            dmc_kwargs={"test": "lr"},
            chain_merge_kwargs={"preset": "default"},
            min_mean_qvalue=0.05,
            n_perm=5, seed=0, n_jobs=1,
        )


def test_chain_merge_empty_observed_returns_empty_cols(monkeypatch):
    monkeypatch.setattr(
        dmr_mod, "_chain_merge_perm_survivors",
        lambda **kw: np.array([1e-6], dtype=np.float64),
    )
    out = dmr_mod.empirical_fdr_for_chain_merge(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmr=pl.DataFrame(schema=_SCHEMA),
        dmc_kwargs={"test": "lr"},
        chain_merge_kwargs={"preset": "default"},
        min_mean_qvalue=0.05,
        n_perm=5, seed=0, n_jobs=1,
    )
    assert "empirical_qvalue" in out.columns
    assert out.height == 0
