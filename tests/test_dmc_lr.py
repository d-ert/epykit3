"""P1-11: log2_odds_ratio renamed per backend.

lr / fisher: log2_odds_ratio_pooled (semantics unchanged, name clearer).
glm:         coef_treatment_log2  (was always the logit coefficient, not log2(OR)).

Both backends emit a transitional log2_odds_ratio column NaN-filled
with a FutureWarning once per tl.dmc call.
"""
from __future__ import annotations
import warnings
import numpy as np
import pytest
import epykit as ep


def test_lr_emits_log2_odds_ratio_pooled(synth_md_filtered):
    md = synth_md_filtered
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="lr")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns, (
        f"missing log2_odds_ratio_pooled; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns, (
        "transitional log2_odds_ratio column must be present (NaN-filled)"
    )
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all(), (
        f"transitional column must be all-NaN; got {legacy[:5]}"
    )
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message), (
        f"expected FutureWarning mentioning log2_odds_ratio; got {[str(w.message) for w in fut]}"
    )


def test_fisher_emits_log2_odds_ratio_pooled(synth_md_filtered):
    md = synth_md_filtered
    ep.tl.dmc(md, test="fisher")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns
    assert "log2_odds_ratio" in df.columns
    assert np.isnan(df["log2_odds_ratio"].to_numpy()).all()
