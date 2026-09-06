"""P1-8 + P2-4: dmr_hmm renamed to dmr_segment; per-segment p-values
are Stouffer-combined from constituent CpG p-values (not NaN as in the
pre-Phase-3 implementation)."""
from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

import epykit as ep


def test_call_dmr_rule_segment_emits_finite_pvalues(synth_md_filtered):
    """The renamed engine must emit finite p/q-values per segment, not NaN."""
    from epykit.dmr_segment import call_dmr_rule_segment

    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    dmrs = call_dmr_rule_segment(md.dmc, min_cpgs=3, min_abs_meth_diff=0.05)
    assert dmrs.height > 0, "expected at least one called segment on synth data"
    pvals = dmrs["pvalue"].to_numpy()
    qvals = dmrs["qvalue"].to_numpy()
    finite_p = np.isfinite(pvals)
    finite_q = np.isfinite(qvals)
    assert finite_p.all(), (
        f"P2-4 fix should produce finite per-segment pvalues; "
        f"got {(~finite_p).sum()}/{len(pvals)} NaN"
    )
    assert finite_q.all(), (
        f"P2-4 fix should produce finite per-segment qvalues; "
        f"got {(~finite_q).sum()}/{len(qvals)} NaN"
    )
    assert ((pvals >= 0) & (pvals <= 1)).all(), "p-values must lie in [0, 1]"
    assert ((qvals >= 0) & (qvals <= 1)).all(), "q-values must lie in [0, 1]"


def test_segment_uses_signed_stouffer_cancelling_mixed_directions():
    """D1: the segment caller must combine per-CpG p-values with the SIGNED
    Stouffer Z (shared with the tile/sliding-window callers), not an unsigned
    two-sided combine. Unsigned addition of |z| makes a region's p shrink
    toward 0 as it grows even when directions are mixed -- anti-conservative.
    The signed combine cancels opposing directions."""
    from scipy.stats import norm

    from epykit.dmr import _stouffer_combine_signed

    pvals = np.full(6, 1e-4)                       # all individually strong
    mixed = np.array([0.3, 0.3, 0.3, -0.3, -0.3, -0.3])
    coherent = np.full(6, 0.3)

    # Reconstruct the OLD unsigned two-sided Stouffer that was removed.
    z = norm.isf(pvals / 2.0)
    p_unsigned = float(2.0 * norm.sf(abs(z.sum() / np.sqrt(6))))

    p_signed_mixed = _stouffer_combine_signed(pvals, mixed)
    p_signed_coherent = _stouffer_combine_signed(pvals, coherent)

    # Unsigned: anti-conservatively tiny regardless of direction.
    assert p_unsigned < 1e-6
    # Signed: mixed directions cancel -> large p; coherent stays significant.
    assert p_signed_mixed > 0.5
    assert p_signed_coherent < 1e-6


def test_dmr_segment_has_no_unsigned_stouffer_combine():
    """The unsigned `_stouffer_combine` helper was removed in the D1 fix;
    guard against it silently coming back."""
    import epykit.dmr_segment as seg
    assert not hasattr(seg, "_stouffer_combine"), (
        "dmr_segment._stouffer_combine (unsigned) must stay removed; use the "
        "shared _stouffer_combine_signed."
    )


def test_dmr_hmm_shim_warns_on_import_and_re_exports():
    """Old `epykit.dmr_hmm.call_dmr_hmm` import path must keep working
    with a DeprecationWarning, and must re-export the new function."""
    import importlib
    import importlib.util

    # Unload if cached to force re-import (triggers shim warning)
    import sys
    for mod in list(sys.modules.keys()):
        if "dmr_hmm" in mod:
            del sys.modules[mod]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        import epykit.dmr_hmm as legacy
        assert hasattr(legacy, "call_dmr_hmm")
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert dep, "expected DeprecationWarning on importing epykit.dmr_hmm"
        assert "dmr_segment" in str(dep[0].message).lower(), (
            f"shim warning must point users to dmr_segment; got: {dep[0].message}"
        )
