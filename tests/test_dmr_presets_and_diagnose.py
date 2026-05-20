"""DMR_PRESETS bundle + ep.tl.diagnose_dmr_calling.

Covers the contract:

* PRESETS:
    - All three presets resolve to known kwargs
    - Calling with preset+explicit kwarg has the explicit kwarg win
    - Unknown preset raises ValueError
    - Preset applied via tl.dmr produces same DMRs as the equivalent
      explicit kwarg call (round-trip equivalence)

* diagnose_dmr_calling:
    - Returns the expected dict shape
    - Buckets sum to the reference DMR count
    - A reference DMR overlapping an existing called DMR -> SUCCESS_OVERLAP
    - A reference DMR in a chromosome with zero CpGs -> H1_NO_CPGS
    - A reference DMR whose CpGs are not significant -> H2_NO_SIG_CPGS
    - Raises if no DMC table is present
    - Raises if no DMR table is present
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from epykit.dmr import DMR_PRESETS, call_dmr_chain_merge


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

def _fake_dmc_df(n: int = 200, seed: int = 0) -> pl.DataFrame:
    """Tiny synthetic DMC table for preset round-trip tests."""
    rng = np.random.default_rng(seed)
    positions = np.arange(1000, 1000 + n * 50, 50, dtype=np.int64)[:n]
    # ~10% significant CpGs with strong effect; rest are noise
    is_sig = rng.random(n) < 0.10
    pvalue = np.where(is_sig, rng.uniform(1e-9, 1e-7, n), rng.uniform(0.1, 1.0, n))
    qvalue = np.where(is_sig, rng.uniform(1e-7, 1e-5, n), rng.uniform(0.1, 1.0, n))
    meth_diff = np.where(is_sig, rng.uniform(0.15, 0.45, n) * np.sign(rng.standard_normal(n)),
                         rng.uniform(-0.02, 0.02, n))
    return pl.DataFrame({
        "chrom":     ["chr1"] * n,
        "pos":       positions,
        "pvalue":    pvalue,
        "qvalue":    qvalue,
        "meth_diff": meth_diff,
    })


def test_presets_dict_has_expected_keys():
    assert set(DMR_PRESETS) == {"strict", "default", "permissive"}
    required = {"alpha", "min_abs_meth_diff", "dis_merge_bp",
                "min_cpgs", "pct_sig", "minlen_bp"}
    for name, bundle in DMR_PRESETS.items():
        assert required <= set(bundle), f"{name} missing keys"


def test_presets_are_ordered_strict_to_permissive():
    """The bundles should be monotonic across strict -> default -> permissive
    on the knobs where direction is unambiguous."""
    s, d, p = DMR_PRESETS["strict"], DMR_PRESETS["default"], DMR_PRESETS["permissive"]
    # Stricter alpha for stricter bundle
    assert s["alpha"] <= d["alpha"] <= p["alpha"]
    # Stricter effect-size floor for stricter bundle
    assert s["min_abs_meth_diff"] >= d["min_abs_meth_diff"] >= p["min_abs_meth_diff"]
    # More CpGs required for stricter bundle
    assert s["min_cpgs"] >= d["min_cpgs"]
    # Longer minimum length for stricter
    assert s["minlen_bp"] >= d["minlen_bp"]


def test_unknown_preset_raises():
    dmc = _fake_dmc_df(n=50)
    with pytest.raises(ValueError, match="Unknown preset"):
        call_dmr_chain_merge(dmc, preset="overly_permissive")


def test_explicit_kwarg_overrides_preset():
    """Caller-provided kwargs must beat the bundled values. Otherwise
    ``preset='default'`` would be a foot-cannon for users who pass one knob
    explicitly and expect the rest to follow the preset."""
    dmc = _fake_dmc_df(n=200, seed=42)

    # preset="default" sets alpha=1e-5; override to 1e-3
    out_overridden = call_dmr_chain_merge(dmc, preset="default", alpha=1e-3)
    # Bare call with the explicit kwargs that should result
    out_explicit = call_dmr_chain_merge(
        dmc, alpha=1e-3, min_abs_meth_diff=0.1,
        dis_merge_bp=100, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    # Same number of rows + same DMR boundaries (sort first for stability)
    assert out_overridden.height == out_explicit.height
    if out_overridden.height > 0:
        a = out_overridden.sort(["chrom", "start"]).select(["chrom", "start", "end"])
        b = out_explicit.sort(["chrom", "start"]).select(["chrom", "start", "end"])
        assert a.equals(b)


def test_preset_default_matches_explicit_bundle_call():
    """preset='default' with no extra kwargs must produce identical output
    to passing all the bundled values explicitly."""
    dmc = _fake_dmc_df(n=200, seed=7)
    out_preset = call_dmr_chain_merge(dmc, preset="default")
    out_explicit = call_dmr_chain_merge(dmc, **DMR_PRESETS["default"])
    assert out_preset.height == out_explicit.height
    if out_preset.height > 0:
        a = out_preset.sort(["chrom", "start"])
        b = out_explicit.sort(["chrom", "start"])
        # Compare column-by-column; some computed fields (qvalue) recompute
        # identically so equality should hold.
        assert a.select(["chrom", "start", "end"]).equals(
            b.select(["chrom", "start", "end"])
        )


# ---------------------------------------------------------------------------
# diagnose_dmr_calling
# ---------------------------------------------------------------------------

class _StubMD:
    """Minimal duck-typed stand-in for MethylData -- avoids the full
    methylstore fixture. ``diagnose_dmr_calling`` only touches ``uns`` and
    ``varm``, so this is sufficient."""
    def __init__(self, dmc_df: pl.DataFrame, dmr_df: pl.DataFrame, dmc_key: str = "dmc_lr"):
        self.varm = {dmc_key: dmc_df}
        self.uns = {
            "dmc": {"last_key": dmc_key},
            "dmr": dmr_df,
        }


def _ref_dmr(chrom: str, start: int, end: int) -> dict:
    return {"chrom": chrom, "start": start, "end": end}


def test_diagnose_returns_expected_shape():
    from epykit.tl import diagnose_dmr_calling

    dmc = pl.DataFrame({
        "chrom":  ["chr1"] * 10,
        "pos":    list(range(1000, 1100, 10)),
        "pvalue": [0.5] * 10,
        "qvalue": [0.5] * 10,
    })
    dmr = pl.DataFrame({
        "chrom": ["chr1"], "start": [2000], "end": [2100],
        "dmr_type": ["hyper"],
    })
    md = _StubMD(dmc, dmr)
    ref = pl.DataFrame([_ref_dmr("chr1", 1000, 1100)])

    out = diagnose_dmr_calling(md, ref)
    assert set(out) == {"counts", "bucket_indices", "n_reference",
                        "summary", "alpha_threshold"}
    assert out["n_reference"] == 1
    assert sum(out["counts"].values()) == 1
    assert isinstance(out["summary"], str)


def test_diagnose_buckets_sum_to_total():
    from epykit.tl import diagnose_dmr_calling

    # 8 reference DMRs; build a DMC + DMR set that exercises each bucket
    dmc = pl.DataFrame({
        "chrom":  ["chr1"] * 5 + ["chr2"] * 5,
        "pos":    list(range(1000, 1500, 100)) + list(range(5000, 5500, 100)),
        # chr1: all non-sig; chr2: first two strong, last three middle-strength
        "pvalue": [0.5] * 5 + [1e-9, 1e-9, 1e-3, 1e-3, 1e-3],
        "qvalue": [0.5] * 5 + [1e-7, 1e-7, 1e-3, 1e-3, 1e-3],
    })
    dmr = pl.DataFrame({
        "chrom": ["chr2"], "start": [5000], "end": [5050],
        "dmr_type": ["hyper"],
    })
    md = _StubMD(dmc, dmr)
    ref = pl.DataFrame([
        _ref_dmr("chr2", 5000, 5050),     # SUCCESS_OVERLAP (our DMR overlaps)
        _ref_dmr("chr3", 1000, 2000),     # H1 (no CpGs on chr3)
        _ref_dmr("chr1", 1000, 1500),     # H2 (CpGs but q>=0.05)
        _ref_dmr("chr2", 5200, 5500),     # H3a (q in [1e-5, 0.05) range)
    ])
    out = diagnose_dmr_calling(md, ref)
    assert sum(out["counts"].values()) == out["n_reference"] == 4


def test_diagnose_bucket_classification_is_correct():
    """Each bucket is exercised by a hand-crafted reference DMR."""
    from epykit.tl import diagnose_dmr_calling

    dmc = pl.DataFrame({
        "chrom":  ["chr1"] * 5 + ["chr2"] * 5,
        "pos":    list(range(1000, 1500, 100)) + list(range(5000, 5500, 100)),
        "pvalue": [0.5] * 5 + [1e-9, 1e-9, 1e-3, 1e-3, 1e-3],
        "qvalue": [0.5] * 5 + [1e-7, 1e-7, 1e-3, 1e-3, 1e-3],
    })
    dmr = pl.DataFrame({
        "chrom": ["chr2"], "start": [5000], "end": [5050],
        "dmr_type": ["hyper"],
    })
    md = _StubMD(dmc, dmr)
    ref = pl.DataFrame([
        _ref_dmr("chr2", 5000, 5050),     # 0: SUCCESS_OVERLAP
        _ref_dmr("chr3", 1000, 2000),     # 1: H1_NO_CPGS
        _ref_dmr("chr1", 1000, 1500),     # 2: H2_NO_SIG_CPGS
        _ref_dmr("chr2", 5200, 5500),     # 3: H3a (q=1e-3 in range)
    ])
    out = diagnose_dmr_calling(md, ref, alpha_threshold=1e-5)
    assert 0 in out["bucket_indices"]["SUCCESS_OVERLAP"]
    assert 1 in out["bucket_indices"]["H1_NO_CPGS"]
    assert 2 in out["bucket_indices"]["H2_NO_SIG_CPGS"]
    assert 3 in out["bucket_indices"]["H3a_WEAK_ALPHA"]


def test_diagnose_h3b_when_q_below_alpha_threshold():
    """A reference DMR with CpGs at q < alpha_threshold but no DMR overlap
    should land in H3b_STRUCTURE (chain-merge dropped on structural filter)."""
    from epykit.tl import diagnose_dmr_calling

    dmc = pl.DataFrame({
        "chrom":  ["chr1"] * 3,
        "pos":    [1000, 1010, 1020],
        "pvalue": [1e-9, 1e-9, 1e-9],
        "qvalue": [1e-7, 1e-7, 1e-7],
    })
    # No DMR was called at chr1:1000-1020 (e.g. failed minlen or min_cpgs)
    dmr = pl.DataFrame({
        "chrom": ["chr99"], "start": [1], "end": [2],   # something elsewhere
        "dmr_type": ["hyper"],
    })
    md = _StubMD(dmc, dmr)
    ref = pl.DataFrame([_ref_dmr("chr1", 1000, 1030)])
    out = diagnose_dmr_calling(md, ref, alpha_threshold=1e-5)
    assert out["bucket_indices"]["H3b_STRUCTURE"] == [0]


def test_diagnose_raises_when_no_dmc():
    from epykit.tl import diagnose_dmr_calling

    class _NoDMC:
        varm = {}
        uns = {"dmr": pl.DataFrame({"chrom": ["chr1"], "start": [1], "end": [2], "dmr_type": ["hyper"]})}

    with pytest.raises(ValueError, match="No DMC table"):
        diagnose_dmr_calling(_NoDMC(), pl.DataFrame([_ref_dmr("chr1", 1, 2)]))


def test_diagnose_raises_when_no_dmr():
    from epykit.tl import diagnose_dmr_calling

    class _NoDMR:
        varm = {"dmc_lr": pl.DataFrame({
            "chrom": ["chr1"], "pos": [10], "pvalue": [0.1], "qvalue": [0.1],
        })}
        uns = {"dmc": {"last_key": "dmc_lr"}}

    with pytest.raises(ValueError, match="No DMR table"):
        diagnose_dmr_calling(_NoDMR(), pl.DataFrame([_ref_dmr("chr1", 1, 100)]))


def test_diagnose_chromosomes_filter_restricts_analysis():
    """Passing chromosomes= should drop reference DMRs on other chromosomes."""
    from epykit.tl import diagnose_dmr_calling

    dmc = pl.DataFrame({
        "chrom":  ["chr1", "chr2"],
        "pos":    [1000, 5000],
        "pvalue": [0.5, 0.5],
        "qvalue": [0.5, 0.5],
    })
    dmr = pl.DataFrame({
        "chrom": ["chr99"], "start": [1], "end": [2], "dmr_type": ["hyper"],
    })
    md = _StubMD(dmc, dmr)
    ref = pl.DataFrame([
        _ref_dmr("chr1", 900, 1100),
        _ref_dmr("chr2", 4900, 5100),
        _ref_dmr("chr3", 100, 200),
    ])
    out = diagnose_dmr_calling(md, ref, chromosomes=["chr1"])
    assert out["n_reference"] == 1   # chr2 + chr3 filtered out


def test_diagnose_public_api_exposed_on_tl_namespace():
    """Smoke check that ep.tl.diagnose_dmr_calling resolves."""
    import epykit as ep
    assert callable(ep.tl.diagnose_dmr_calling)


def test_dmr_presets_exposed_on_top_level_and_tl():
    import epykit as ep
    from epykit.tl import DMR_PRESETS as tl_presets
    assert "default" in ep.DMR_PRESETS
    assert ep.DMR_PRESETS is tl_presets  # same object, no shadow
