"""Sliding-window DMR combine must use raw per-CpG p-values, not q-values (M-DMR1).

The region statistic is a signed Stouffer combination, which assumes ~U(0,1)
inputs under the null. BH q-values are not uniform, so combining them does not
yield a valid p-value. The significance *gate* still uses q-values (FDR
control); only the combine must use raw p.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from epykit.dmr import _stouffer_combine_signed, call_dmr_sliding_window


def test_sliding_window_combines_raw_p_not_q():
    pos = [1000, 1100, 1200, 1300, 1400]
    pvalue = [1e-4, 2e-4, 1e-3, 5e-4, 3e-4]
    # q-values inflated relative to p, but still < alpha so the gate passes.
    qvalue = [4e-2, 4e-2, 4e-2, 4e-2, 4e-2]
    meth_diff = [0.30, 0.35, 0.25, 0.40, 0.30]
    df = pl.DataFrame({
        "chrom": ["chr1"] * 5,
        "pos": pos,
        "meth_diff": meth_diff,
        "pvalue": pvalue,
        "qvalue": qvalue,
    })

    out = call_dmr_sliding_window(
        df, window_bp=500, min_cpgs=5, min_sites_significant=5,
        alpha=0.05, min_abs_meth_diff=0.1,
    )
    assert out.height >= 1, "expected one DMR over the tight hyper cluster"

    got = float(out["combined_pvalue"][0])
    exp_raw = _stouffer_combine_signed(np.array(pvalue), np.array(meth_diff))
    exp_q = _stouffer_combine_signed(np.array(qvalue), np.array(meth_diff))

    assert abs(got - exp_raw) < 1e-9, (
        f"combine should use raw p: got {got}, raw-Stouffer {exp_raw}, "
        f"q-Stouffer {exp_q}"
    )
    # And it must genuinely differ from the (wrong) q-value combine.
    assert abs(got - exp_q) > 1e-6, "combine still appears to use q-values"
