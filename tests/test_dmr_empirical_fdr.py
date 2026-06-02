"""P1-6: DMR empirical FDR paired-design awareness + n=1,1 refusal."""
from __future__ import annotations

import pytest
import polars as pl
import epykit as ep


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
