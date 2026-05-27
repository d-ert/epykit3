"""P1-11 GLM half: coef_treatment_log2."""
from __future__ import annotations
import warnings
import numpy as np
import polars as pl
import pytest
import epykit as ep


def test_glm_emits_coef_treatment_log2(synth_md_filtered):
    md = synth_md_filtered
    md.obs = md.obs.with_columns(
        (pl.col("group") == "treatment").cast(int).alias("treatment")
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="glm", formula="~ treatment")
    df = md.dmc
    assert "coef_treatment_log2" in df.columns, (
        f"missing coef_treatment_log2; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all()
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message)
