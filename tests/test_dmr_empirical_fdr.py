"""Permutation empirical FDR for DMRs.

* Paired-design awareness (labels shuffle within strata) and the n=1,1
  refusal.
* The empirical p-value denominator: per-permutation tail counts, not a
  pooled null, and failed permutations excluded from the count.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit import dmr as dmr_mod
from epykit.dmr import _empirical_pvalues_from_null_pool


def test_n_one_each_raises(synth_md_filtered):
    """With one treatment and one control sample, empirical FDR must
    refuse with a ValueError mentioning n>=2."""
    md = synth_md_filtered
    # Trim obs to 1 treatment + 1 control by filtering the obs frame.
    one_treat = md.treatment_ids[:1]
    one_ctrl = md.control_ids[:1]
    kept = one_treat + one_ctrl
    md.obs = md.obs.filter(pl.col("sample_id").is_in(kept))

    # allow_n1=True lets dmc proceed (Fisher fallback); the refusal we are
    # testing happens later, inside empirical_fdr_for_dmr.
    ep.tl.dmc(md, test="fisher", allow_n1=True)
    with pytest.raises(ValueError, match="n.*2"):
        ep.tl.dmr(md, method="tile", empirical_fdr=True, n_perm=10, allow_n1=True)


@pytest.mark.slow
def test_paired_design_shuffles_within_strata(synth_md_filtered):
    """When empirical_strata= is supplied, the shuffle must permute
    within strata. Smoke test: run completes, empirical_qvalue populated."""
    md = synth_md_filtered
    # Build a paired covariate: each treatment sample paired with one ctrl.
    n_pair = min(len(md.treatment_ids), len(md.control_ids))
    pair_ids = [f"Pair{i}" for i in range(n_pair)] * 2
    all_ids = list(md.treatment_ids[:n_pair]) + list(md.control_ids[:n_pair])
    pair_series = pl.Series("subject_id", pair_ids[: len(all_ids)])
    md.obs = md.obs.with_columns(pair_series)

    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(
        md,
        method="tile",
        empirical_fdr=True,
        n_perm=20,
        empirical_strata="subject_id",
    )
    dmrs = md.uns["dmr"]
    assert "empirical_qvalue" in dmrs.columns, (
        f"empirical_qvalue missing; got {dmrs.columns}"
    )


# --- Empirical p-value denominator: per-permutation tail counts --------------
# empirical_fdr_for_dmr must count, per permutation, whether at least one
# null region is as extreme as the observed one -- not pool every null
# region across permutations. The math is checked by faking the permutation
# output; the real call_dmr_tile_based is too heavy here.

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
    call_count = {"i": 0}

    def fake_call_dmr_tile_based(**kwargs):
        i = call_count["i"]
        call_count["i"] += 1
        return pl.DataFrame({"pvalue": perm_pvals[i % len(perm_pvals)]})

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
    )
    emp_p = out["empirical_pvalue"][0]
    assert emp_p > 0.5, (
        f"emp_p={emp_p:.4f} -- looks like the pooled-null formula is "
        f"still in use. Expected ~1.0 from per-perm tail counts."
    )


def test_emp_p_correct_on_pure_null(monkeypatch):
    """If the observed p is smaller than every single null p in every
    perm, emp_p should hit the floor 1 / (n_perm + 1)."""
    n_perm = 20
    perm_pvals = [np.full(50, 0.5) for _ in range(n_perm)]

    def fake_call_dmr_tile_based(**kwargs):
        return pl.DataFrame({"pvalue": perm_pvals.pop(0)})

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake_call_dmr_tile_based)

    observed = _mk_dmr([
        {"chrom": "chr1", "start": 0, "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.5, "pvalue": 1e-9},
    ])
    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=["t1", "t2"],
        samples_control=["c1", "c2"],
        observed_dmr=observed,
        n_perm=n_perm,
        seed=0,
        n_jobs=1,
    )
    emp_p = out["empirical_pvalue"][0]
    assert abs(emp_p - 1.0 / (n_perm + 1)) < 1e-6, (
        f"emp_p={emp_p:.6f}; expected 1/(n_perm+1) = "
        f"{1.0 / (n_perm + 1):.6f}"
    )


# --- Failed permutations must not count in the denominator -------------------
# Pre-fix, the denominator was always n_perm + 1, which biased empirical p
# downward by the failure rate. Under small-n DMR calling some permuted
# splits produce zero candidate regions, so the failure rate is not trivial.


def test_failed_permutations_excluded_from_denominator():
    observed = np.array([1e-6, 1e-3, 0.05])
    null_pool = np.array([1e-4, 1e-2, 0.5])  # only 3 successful perms
    n_perm_successful = 3

    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=n_perm_successful,
    )
    # For observed[0]=1e-6: 0 nulls <= it -> emp = (0+1)/(3+1) = 0.25.
    # The pre-fix code used (10+1) in the denominator -> 0.0909, biased low.
    assert abs(emp[0] - 0.25) < 1e-9, (
        f"Denominator must be n_perm_successful+1=4 -> emp=0.25; got {emp[0]}"
    )


def test_zero_successful_perms_returns_ones():
    """Defensive: if every permutation failed, empirical p must be 1.0
    (the most conservative possible value), not div-by-zero or undefined."""
    observed = np.array([0.01, 0.5])
    null_pool = np.array([], dtype=np.float64)
    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=0,
    )
    np.testing.assert_array_equal(emp, np.array([1.0, 1.0]))


def test_searchsorted_counts_match_expected_p():
    """Sanity: every observed p gets the right (count<=obs + 1)/(n+1)."""
    observed = np.array([0.001, 0.05, 0.1, 0.5])
    null_pool = np.array([0.0005, 0.01, 0.02, 0.06, 0.07, 0.2, 0.3, 0.6])
    n_perm_successful = 4
    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=n_perm_successful,
    )
    # counts: obs=0.001 -> 1 (only 0.0005); 0.05 -> 3; 0.1 -> 5; 0.5 -> 7.
    # denominator = 4+1 = 5
    expected = np.array([(1 + 1) / 5, (3 + 1) / 5, (5 + 1) / 5, (7 + 1) / 5])
    np.testing.assert_allclose(emp, expected, atol=1e-12)
