"""max_t mode: empirical_fdr_for_dmr(fdr_method="max_t") must use
per-permutation tail counts (Westfall-Young min-P / FWER), not pooled-null
counts. These guard the opt-in ``max_t`` mode; the DEFAULT ``region`` mode
deliberately uses the pooled count-ratio target-decoy FDR instead (see
tests/test_region_count_ratio_fdr.py and test_dmr_region_fdr_mode.py).
We test the math via a small monkeypatch that fakes the permutation output --
the real call_dmr_tile_based is too heavy here."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from epykit import dmr as dmr_mod


_DMR_SCHEMA = {
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "n_cpgs": pl.Int32,
    "meth_diff": pl.Float64,
    "pvalue": pl.Float64,
}


def _mk_dmr(rows: list[dict]) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=_DMR_SCHEMA)


def test_emp_p_uses_n_perm_denominator_not_pooled(monkeypatch):
    """5 permutations, each producing 100 null DMRs at p ~ U(0, 1).
    For an observed p_obs = 0.02, 'at least one null <= 0.02' should fire
    in roughly 5/5 perms (perm-wise tail), giving emp_p ~ 1.0 -- *not*
    ~ 0.02 (which is what the old pooled formula gives at total_null=500)."""
    rng = np.random.default_rng(0)
    perm_pvals = [np.sort(rng.uniform(0, 1, size=100)) for _ in range(5)]
    # Replace _run_one_perm with one that yields the precomputed null p's.
    call_count = {"i": 0}

    def fake_inner(perm_idx: int) -> np.ndarray:
        i = call_count["i"]
        call_count["i"] += 1
        return perm_pvals[i % len(perm_pvals)]

    # Patch call_dmr_tile_based inside the module so each "perm" run
    # returns one of perm_pvals as the DMR pvalue column.
    def fake_call_dmr_tile_based(**kwargs):
        # Return a small DataFrame whose 'pvalue' column is one perm's array.
        arr = fake_inner(0)
        return pl.DataFrame({"pvalue": arr})

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake_call_dmr_tile_based)

    observed = _mk_dmr([
        {"chrom": "chr1", "start": 0, "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": 0.02},
    ])

    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",  # unused with monkeypatch
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmr=observed,
        n_perm=5,
        seed=0,
        n_jobs=1,
        fdr_method="max_t",
    )
    emp_p = out["empirical_pvalue"][0]
    # Every perm has ~2 null DMRs with p <= 0.02 (uniform draws of 100),
    # so 'at least one' fires in ~5/5 perms -- emp_p should be ~ 1.0,
    # not ~ 0.02 (pooled-null formula).
    assert emp_p > 0.5, (
        f"emp_p={emp_p:.4f} -- looks like the pooled-null formula is "
        f"still in use. Expected ~1.0 from per-perm tail counts."
    )


def test_emp_p_correct_on_pure_null():
    """If the observed p is smaller than every single null p in every
    perm, emp_p should hit the floor 1 / (n_perm + 1)."""
    n_perm = 20
    # Each perm yields nulls bounded below at 0.5.
    perm_pvals = [np.full(50, 0.5) for _ in range(n_perm)]

    def fake_call_dmr_tile_based(**kwargs):
        return pl.DataFrame({"pvalue": perm_pvals.pop(0)})

    import importlib
    importlib.reload(dmr_mod)  # reset module state between tests

    from epykit import dmr as dmr2
    dmr2.call_dmr_tile_based = fake_call_dmr_tile_based

    observed = _mk_dmr([
        {"chrom": "chr1", "start": 0, "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.5, "pvalue": 1e-9},
    ])
    out = dmr2.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmr=observed,
        n_perm=n_perm,
        seed=0,
        n_jobs=1,
        fdr_method="max_t",
    )
    emp_p = out["empirical_pvalue"][0]
    assert abs(emp_p - 1.0 / (n_perm + 1)) < 1e-6, (
        f"emp_p={emp_p:.6f}; expected 1/(n_perm+1) = "
        f"{1.0 / (n_perm + 1):.6f}"
    )
