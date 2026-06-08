"""Integration tests for empirical_fdr_for_dmr fdr_method dispatch.

Uses the monkeypatch pattern (real call_dmr_tile_based is too heavy): each
"permutation" returns a fixed null survivor set, so the estimator math is
exercised end-to-end with hand-checkable expectations.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from epykit import dmr as dmr_mod

_SCHEMA = {
    "chrom": pl.Utf8, "start": pl.Int32, "end": pl.Int32,
    "n_cpgs": pl.Int32, "meth_diff": pl.Float64, "pvalue": pl.Float64,
}


def _observed(n=50, p=1e-6):
    return pl.DataFrame({
        "chrom": ["chr1"] * n,
        "start": list(range(0, n * 1000, 1000)),
        "end": list(range(1000, n * 1000 + 1000, 1000)),
        "n_cpgs": [10] * n,
        "meth_diff": [0.3] * n,
        "pvalue": [p] * n,
    }, schema=_SCHEMA)


def test_region_mode_significant_where_maxt_saturates(monkeypatch):
    """50 real tiles at p=1e-6; every shuffle yields only 2 extreme nulls.
    Region count-ratio: V(1e-6)=2, R=50 -> q=0.04 (significant).
    max-T would saturate every tile at emp_p=1.0."""
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based",
                        lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]}))
    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2", "t3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmr=_observed(50, 1e-6),
        n_perm=20, seed=0, n_jobs=1, fdr_method="region",
    )
    q = out.get_column("empirical_qvalue").to_numpy()
    assert np.all(q < 0.05), f"region q not significant: max={q.max()}"
    assert "empirical_fdr_set" in out.columns
    assert out.get_column("empirical_fdr_set")[0] == pytest.approx(0.04, abs=1e-6)


def test_region_is_the_default(monkeypatch):
    """Omitting fdr_method uses region (count-ratio), not max-T."""
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based",
                        lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]}))
    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2", "t3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmr=_observed(50, 1e-6),
        n_perm=20, seed=0, n_jobs=1,
    )
    assert np.all(out.get_column("empirical_qvalue").to_numpy() < 0.05)


def test_maxt_mode_still_saturates(monkeypatch):
    """fdr_method='max_t' preserves the FWER min-P behaviour."""
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based",
                        lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]}))
    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2", "t3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmr=_observed(50, 1e-6),
        n_perm=20, seed=0, n_jobs=1, fdr_method="max_t",
    )
    emp_p = out.get_column("empirical_pvalue").to_numpy()
    assert np.all(emp_p > 0.5), f"max-T should saturate; got max emp_p={emp_p.max()}"


def test_region_small_n_warns(monkeypatch):
    """Permutation inference at <4/group is underpowered -> UserWarning."""
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based",
                        lambda **kw: pl.DataFrame({"pvalue": [1e-6]}))
    with pytest.warns(UserWarning, match="permutation"):
        dmr_mod.empirical_fdr_for_dmr(
            methylstore_path="/dev/null",
            samples_treatment=["t1", "t2", "t3"],
            samples_control=["c1", "c2", "c3"],
            observed_dmr=_observed(10, 1e-6),
            n_perm=10, seed=0, n_jobs=1, fdr_method="region",
        )


def test_unknown_fdr_method_raises(monkeypatch):
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based",
                        lambda **kw: pl.DataFrame({"pvalue": [1e-6]}))
    with pytest.raises(ValueError, match="fdr_method"):
        dmr_mod.empirical_fdr_for_dmr(
            methylstore_path="/dev/null",
            samples_treatment=["t1", "t2"],
            samples_control=["c1", "c2"],
            observed_dmr=_observed(5, 1e-6),
            n_perm=5, seed=0, n_jobs=1, fdr_method="bogus",
        )
