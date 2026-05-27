"""P0-2b: empirical_fdr_for_dmc must use per-permutation tail counts,
not pooled-null counts. Same fix shape as the P0-2 DMR fix.

Strategy: mock process_chromosomes_dmc (the real heavy call inside
_run_one_perm) to return a lightweight stub DMCStore yielding prebuilt
null p-value arrays. This avoids touching the real Parquet store and
lets us control null distributions precisely.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from epykit import dmc as dmc_mod


class _MockDMCStore:
    """Minimal stand-in for DMCStore used inside empirical_fdr_for_dmc."""

    def __init__(self, pvals: np.ndarray):
        self._pvals = pvals

    @property
    def total_sites(self) -> int:
        return len(self._pvals)

    def iter_chroms(self, columns=None):
        # Yield a single fake chromosome with the requested p-value array.
        df = pl.DataFrame({"pvalue": self._pvals})
        yield "chr1", df

    def cleanup(self) -> None:
        pass  # nothing to clean up in the mock


def test_dmc_emp_p_uses_n_perm_denominator_not_pooled(monkeypatch):
    """5 permutations, each producing 100 null sites at p ~ U(0,1).
    For an observed p_obs = 0.02, every perm has ~2 null sites with
    p <= 0.02 (uniform), so 'at least one' fires in ~5/5 perms.
    Per-permutation formula => emp_p ~ 1.0.
    Old pooled formula => emp_p ~ 0.02 (total_null = 500).
    """
    rng = np.random.default_rng(42)
    n_perm = 5
    perm_pvals = [np.sort(rng.uniform(0, 1, size=100)) for _ in range(n_perm)]
    call_count = {"i": 0}

    def fake_process_chromosomes_dmc(**kwargs):
        i = call_count["i"] % n_perm
        call_count["i"] += 1
        return _MockDMCStore(perm_pvals[i])

    monkeypatch.setattr(dmc_mod, "process_chromosomes_dmc", fake_process_chromosomes_dmc)

    observed = pl.DataFrame({
        "chrom": ["chr1"],
        "pos": [100],
        "pvalue": [0.02],
    })

    out = dmc_mod.empirical_fdr_for_dmc(
        methylstore_path="/fake/path",
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmc=observed,
        n_perm=n_perm,
        seed=0,
        n_jobs=1,
    )
    emp_p = out["empirical_pvalue"][0]
    # Every perm has several null sites with p <= 0.02 (uniform from 100
    # draws), so 'at least one' fires in all 5 perms.  The per-perm formula
    # yields emp_p = (5 + 1) / (5 + 1) = 1.0. The old pooled formula would
    # have given ~ 0.02 * 500 / 501 ≈ 0.02. Require > 0.5 as a loose bound.
    assert emp_p > 0.5, (
        f"emp_p={emp_p:.4f} -- looks like the pooled-null formula is still "
        f"in use. Expected ~1.0 from per-perm tail counts."
    )


def test_dmc_emp_p_correct_on_pure_null(monkeypatch):
    """If observed p is smaller than every null p in every perm, emp_p
    should hit the floor 1 / (n_perm + 1)."""
    n_perm = 20
    perm_pvals_queue = [np.full(50, 0.5) for _ in range(n_perm)]
    call_count = {"i": 0}

    def fake_process_chromosomes_dmc(**kwargs):
        i = call_count["i"]
        call_count["i"] += 1
        return _MockDMCStore(perm_pvals_queue[i])

    monkeypatch.setattr(dmc_mod, "process_chromosomes_dmc", fake_process_chromosomes_dmc)

    observed = pl.DataFrame({
        "chrom": ["chr1"],
        "pos": [100],
        "pvalue": [1e-9],  # smaller than all null p's (0.5)
    })

    out = dmc_mod.empirical_fdr_for_dmc(
        methylstore_path="/fake/path",
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmc=observed,
        n_perm=n_perm,
        seed=0,
        n_jobs=1,
    )
    emp_p = out["empirical_pvalue"][0]
    expected_floor = 1.0 / (n_perm + 1)
    assert abs(emp_p - expected_floor) < 1e-6, (
        f"emp_p={emp_p:.6f}; expected floor 1/(n_perm+1) = {expected_floor:.6f}"
    )
