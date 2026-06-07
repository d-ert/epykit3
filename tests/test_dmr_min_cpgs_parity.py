"""M10: --min-cpgs parity for the chain_merge DMR caller.

Pre-fix bugs this pins:

  1. The CLI ``epykit dmr --method chain_merge`` branch never forwarded
     ``--min-cpgs`` to ``call_dmr_chain_merge`` -- the engine default
     applied regardless, so ``--min-cpgs 10 --method chain_merge`` was
     silently dropped (different DMR set from the API on the DEFAULT
     caller).

  2. ``call_dmr_chain_merge`` used ``min_cpgs == 3`` as the
     "user-didn't-set-this" sentinel for preset overrides. Any caller
     passing a non-3 value (e.g. ``tl.dmr`` passing its default 5)
     suppressed a preset's ``min_cpgs``; and a user genuinely asking for
     ``min_cpgs=3`` *with* a preset had it overridden.

The fix uses a ``None`` sentinel end-to-end. Resolution order for
``call_dmr_chain_merge``: explicit value wins; else preset's value (if a
preset is active); else 3 (DSS engine default). ``tl.dmr`` resolves None
to the preset's value (preset active) or 5 (no preset -- its documented
chain_merge default, which preserves the benchmark paper's numbers).
"""

from __future__ import annotations

import argparse

import polars as pl
import pytest

import epykit as ep
from epykit.dmr import DMR_PRESETS, call_dmr_chain_merge


def _mk_dmc_frame(rows: list[dict]) -> pl.DataFrame:
    """Build a minimal DMC-shaped DataFrame for call_dmr_chain_merge."""
    return pl.DataFrame(
        rows,
        schema={
            "chrom":     pl.Utf8,
            "pos":       pl.Int32,
            "meth_diff": pl.Float32,
            "pvalue":    pl.Float64,
        },
    )


def _strict_region(n: int) -> pl.DataFrame:
    """A single chain of ``n`` strongly-significant CpGs that satisfies
    every ``preset='strict'`` filter *except* min_cpgs.

    strict = alpha=1e-6, min_abs_meth_diff=0.20, dis_merge_bp=250,
    pct_sig=0.5, minlen_bp=100. The CpGs are 80 bp apart (gap <= 250),
    span >= 100 bp for n >= 3, |meth_diff|=0.30 (>= 0.20), pvalue=1e-8
    (< 1e-6), and 100% significant (>= 0.5). So only min_cpgs decides
    whether the region survives.
    """
    rows = [
        {"chrom": "chr1", "pos": 1000 + 80 * i, "meth_diff": 0.30, "pvalue": 1e-8}
        for i in range(n)
    ]
    return _mk_dmc_frame(rows)


# 1. Engine: a preset's min_cpgs is actually applied.


def test_chain_merge_engine_preset_min_cpgs_applied():
    """``preset='strict'`` (min_cpgs=5) drops a 4-CpG region; the same
    region survives once min_cpgs is explicitly lowered to 4."""
    assert DMR_PRESETS["strict"]["min_cpgs"] == 5  # guards the test premise
    df = _strict_region(4)

    dropped = call_dmr_chain_merge(df, preset="strict")
    assert len(dropped) == 0, "4-CpG region must fail the preset's min_cpgs=5"

    kept = call_dmr_chain_merge(df, preset="strict", min_cpgs=4)
    assert len(kept) == 1, "explicit min_cpgs=4 must override the preset's 5"


# 2. Engine: an explicit min_cpgs always overrides the preset --
#    including the genuine min_cpgs=3 case that the old sentinel ate.


def test_chain_merge_explicit_min_cpgs_overrides_preset():
    """A user passing ``min_cpgs=3`` alongside ``preset='strict'`` must
    get 3 (the old code treated 3 as "unset" and used the preset's 5)."""
    df = _strict_region(3)

    kept = call_dmr_chain_merge(df, preset="strict", min_cpgs=3)
    assert len(kept) == 1, (
        "explicit min_cpgs=3 must be honored, not swallowed by the preset"
    )

    dropped = call_dmr_chain_merge(df, preset="strict", min_cpgs=7)
    assert len(dropped) == 0, "explicit min_cpgs=7 must override the preset"


# 3. Engine: with no preset and no explicit value, the DSS default (3)
#    applies (preserves direct-engine-caller backward compat).


def test_chain_merge_engine_no_preset_default_is_three():
    df = _strict_region(3)  # 3 CpGs, strict geometry but no preset gates
    out = call_dmr_chain_merge(
        df, alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=250, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 1, "engine default min_cpgs must be 3 (3 CpGs kept)"


# 4. API: tl.dmr must not suppress a preset's min_cpgs.


def test_tl_dmr_chain_merge_preset_min_cpgs_not_suppressed(synth_md_filtered):
    """``tl.dmr(method='chain_merge', preset='default')`` must resolve
    min_cpgs to the preset's value (3), not tl.dmr's own default of 5."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    assert DMR_PRESETS["default"]["min_cpgs"] == 3  # premise guard
    ep.tl.dmr(md, method="chain_merge", preset="default", min_mean_qvalue=None)
    assert md.uns["dmr_params"]["min_cpgs"] == 3, (
        "preset min_cpgs (3) was suppressed by tl.dmr's default 5"
    )


# 5. API: the bare chain_merge default stays 5 (paper preservation).


def test_tl_dmr_chain_merge_bare_default_preserved(synth_md_filtered):
    """``tl.dmr(method='chain_merge')`` with no preset and no min_cpgs must
    still behave as min_cpgs=5 -- the benchmark paper's default caller."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="chain_merge", min_mean_qvalue=None)
    assert md.uns["dmr_params"]["min_cpgs"] == 5, (
        "bare tl.dmr(method='chain_merge') effective min_cpgs must stay 5"
    )


# 6. CLI: --min-cpgs is forwarded to chain_merge (the drop bug).


def test_cli_chain_merge_forwards_min_cpgs(tmp_path, monkeypatch):
    """``_cmd_dmr`` with method=chain_merge must pass --min-cpgs through to
    call_dmr_chain_merge. Pre-fix the branch dropped it entirely."""
    import epykit.dmr as ep_dmr
    from epykit.cli import _cmd_dmr

    dmc = _strict_region(5)
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(dmc_path)

    captured: dict = {}
    real = ep_dmr.call_dmr_chain_merge

    def _capturing(dmc_results, **kwargs):
        captured.update(kwargs)
        return real(dmc_results, **kwargs)

    monkeypatch.setattr(ep_dmr, "call_dmr_chain_merge", _capturing)

    args = argparse.Namespace(
        method="chain_merge", empirical_fdr=False,
        dmc_results=str(dmc_path), preset=None,
        alpha=0.05, min_abs_meth_diff=0.1, dis_merge_bp=500,
        pct_sig=0.5, minlen_bp=50, use_q_for_sig=False,
        min_cpgs=10, output=str(tmp_path / "out.parquet"), no_tsv=True,
    )
    _cmd_dmr(args)
    assert captured.get("min_cpgs") == 10, (
        "CLI --min-cpgs 10 must reach call_dmr_chain_merge"
    )


def test_cli_chain_merge_unset_min_cpgs_flows_none(tmp_path, monkeypatch):
    """With --min-cpgs unset (argparse default None), the CLI forwards None
    so the engine resolves via preset/3 rather than a hard-coded value."""
    import epykit.dmr as ep_dmr
    from epykit.cli import _cmd_dmr

    dmc = _strict_region(5)
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(dmc_path)

    captured: dict = {}
    real = ep_dmr.call_dmr_chain_merge

    def _capturing(dmc_results, **kwargs):
        captured.update(kwargs)
        return real(dmc_results, **kwargs)

    monkeypatch.setattr(ep_dmr, "call_dmr_chain_merge", _capturing)

    args = argparse.Namespace(
        method="chain_merge", empirical_fdr=False,
        dmc_results=str(dmc_path), preset="strict",
        alpha=0.05, min_abs_meth_diff=0.1, dis_merge_bp=500,
        pct_sig=0.5, minlen_bp=50, use_q_for_sig=False,
        min_cpgs=None, output=str(tmp_path / "out.parquet"), no_tsv=True,
    )
    _cmd_dmr(args)
    assert "min_cpgs" in captured, "chain_merge branch must forward min_cpgs"
    assert captured["min_cpgs"] is None, (
        "unset --min-cpgs must flow None so the preset's value wins"
    )
