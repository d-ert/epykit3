"""D11: ``epykit dmr --min-mean-qvalue`` region-level q-filter parity.

``tl.dmr`` applies a post-hoc region q-value filter (default
``min_mean_qvalue=0.05``) to the chain_merge and sliding_window callers, but the
CLI ``_cmd_dmr`` never did -- so CLI chain_merge/sliding_window output was less
filtered than the equivalent API call. This pins the new ``--min-mean-qvalue``
flag and that ``_cmd_dmr`` mirrors ``tl.dmr``'s exact filter (same column
fallback: ``combined_qvalue`` else ``combined_pvalue``).
"""

from __future__ import annotations

import argparse

import polars as pl

from epykit.cli import _cmd_dmr, build_parser
from epykit.dmr import call_dmr_chain_merge, call_dmr_sliding_window


def _two_region_dmc() -> pl.DataFrame:
    """Two well-separated DMR chains: chr1 (extremely significant -> tiny
    combined_qvalue) and chr2 (significant but much weaker)."""
    rows = []
    for i in range(5):
        rows.append({"chrom": "chr1", "pos": 1000 + 80 * i, "meth_diff": 0.4, "pvalue": 1e-10})
    for i in range(5):
        rows.append({"chrom": "chr2", "pos": 1000 + 80 * i, "meth_diff": 0.15, "pvalue": 2e-3})
    return pl.DataFrame(
        rows,
        schema={"chrom": pl.Utf8, "pos": pl.Int32, "meth_diff": pl.Float32, "pvalue": pl.Float64},
    )


def _chain_merge_namespace(tmp_path, dmc_path, *, min_mean_qvalue):
    return argparse.Namespace(
        method="chain_merge", empirical_fdr=False,
        dmc_results=str(dmc_path), preset=None,
        alpha=0.05, min_abs_meth_diff=0.1, dis_merge_bp=500,
        pct_sig=0.5, minlen_bp=50, use_q_for_sig=False,
        min_cpgs=5, min_mean_qvalue=min_mean_qvalue,
        output=str(tmp_path / "out.parquet"), no_tsv=True,
    )


# argparse wiring


def test_dmr_parser_has_min_mean_qvalue_default_005():
    parser = build_parser()
    args = parser.parse_args(["dmr", "--method", "chain_merge", "--output", "x.parquet"])
    assert hasattr(args, "min_mean_qvalue")
    assert args.min_mean_qvalue == 0.05, "must match tl.dmr's default of 0.05"


# chain_merge filtering parity


def test_cli_chain_merge_applies_q_filter(tmp_path):
    """A threshold strictly between the two regions' combined_qvalues must keep
    only the stronger region -- and match tl.dmr's exact filter."""
    dmc = _two_region_dmc()
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(dmc_path)

    unfiltered = call_dmr_chain_merge(
        dmc, alpha=0.05, min_abs_meth_diff=0.1, dis_merge_bp=500,
        pct_sig=0.5, minlen_bp=50, min_cpgs=5,
    )
    assert len(unfiltered) == 2, "premise: both regions called pre-filter"
    qs = sorted(unfiltered["combined_qvalue"].to_list())
    threshold = (qs[0] * qs[1]) ** 0.5  # geometric midpoint, splits the two

    # Mirror tl.dmr's exact post-filter to compute the expected survivor set.
    expected = unfiltered.filter(pl.col("combined_qvalue") < threshold)
    assert len(expected) == 1

    args = _chain_merge_namespace(tmp_path, dmc_path, min_mean_qvalue=threshold)
    _cmd_dmr(args)
    got = pl.read_parquet(args.output)
    assert got.sort("chrom").equals(expected.sort("chrom")), (
        "CLI chain_merge output must equal tl.dmr's q-filtered set exactly"
    )


def test_cli_chain_merge_none_disables_q_filter(tmp_path):
    """``min_mean_qvalue=None`` keeps every called region (filter disabled)."""
    dmc = _two_region_dmc()
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(dmc_path)

    args = _chain_merge_namespace(tmp_path, dmc_path, min_mean_qvalue=None)
    _cmd_dmr(args)
    got = pl.read_parquet(args.output)
    assert len(got) == 2, "None must disable the q-filter (both regions kept)"


# sliding_window filtering parity


def test_cli_sliding_window_applies_q_filter(tmp_path):
    dmc = _two_region_dmc()
    dmc_path = tmp_path / "dmc.parquet"
    dmc.write_parquet(dmc_path)

    unfiltered = call_dmr_sliding_window(
        dmc, window_bp=500, step_bp=250, min_cpgs=5,
        min_sites_significant=3, alpha=0.05, min_abs_meth_diff=0.1,
    )
    assert len(unfiltered) >= 1
    q_col = "combined_qvalue" if "combined_qvalue" in unfiltered.columns else "combined_pvalue"
    qs = sorted(unfiltered[q_col].to_list())
    # Threshold below the largest q so at least one region is dropped.
    threshold = qs[-1]
    expected = unfiltered.filter(pl.col(q_col) < threshold)

    args = argparse.Namespace(
        method="sliding_window", empirical_fdr=False,
        dmc_results=str(dmc_path),
        window_bp=500, step_bp=250, min_cpgs=5, min_sites_significant=3,
        alpha=0.05, min_abs_meth_diff=0.1,
        min_mean_qvalue=threshold,
        output=str(tmp_path / "sw.parquet"), no_tsv=True,
    )
    _cmd_dmr(args)
    got = pl.read_parquet(args.output)
    assert got.sort("chrom").equals(expected.sort("chrom")), (
        "CLI sliding_window must apply the same q-filter as tl.dmr exactly "
        "(not merely the same row count)"
    )


# tile filtering parity (D11 follow-up: the CLI tile branch was missing the
# qvalue post-filter that tl.dmr's tile path applies, so `epykit dmr
# --method tile` diverged from tl.dmr(method='tile') when --min-mean-qvalue
# was tighter than alpha).


def test_cli_tile_applies_qvalue_filter(tmp_path, monkeypatch):
    """CLI tile output must equal tl.dmr tile's ``qvalue < min_mean_qvalue``
    post-filtered set. The tile engine is stubbed so the test pins the CLI's
    post-filter behavior (the part that was missing) rather than the tile
    statistics. With a threshold tighter than alpha the filter must drop the
    marginal tile -- proving the parity gap is closed."""
    import epykit.cli as ep_cli
    import epykit.dmr as ep_dmr

    # Two tiles that both already passed the engine's alpha=0.05 filter: one
    # strongly significant, one marginal. tl.dmr tile filters on ``qvalue``.
    tile_frame = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [100, 100],
            "end": [200, 200],
            "pvalue": [1e-6, 0.04],
            "qvalue": [1e-6, 0.04],
            "meth_diff": [0.4, 0.2],
        }
    )
    monkeypatch.setattr(ep_dmr, "call_dmr_tile_based", lambda *a, **k: tile_frame)
    monkeypatch.setattr(
        ep_cli, "_read_samplesheet_groups", lambda *a, **k: (["t0", "t1"], ["c0", "c1"])
    )
    monkeypatch.setattr(ep_cli, "_cli_n1_and_footgun_checks", lambda *a, **k: None)

    threshold = 0.01  # tighter than alpha=0.05 -> must drop the q=0.04 tile
    # The expected survivor set is exactly tl.dmr tile's filter: qvalue-only.
    expected = tile_frame.filter(pl.col("qvalue") < threshold)
    assert len(expected) == 1, "premise: threshold drops the marginal tile"

    args = argparse.Namespace(
        method="tile", empirical_fdr=False,
        methylstore="store", samplesheet="ss.csv",
        treatment_group="A", control_group="B",
        tile_size_bp=1000, test="lr", min_cpgs_per_tile=3,
        alpha=0.05, min_abs_meth_diff=0.1, unite=True,
        min_samples_treatment=2, min_samples_control=2,
        min_mean_qvalue=threshold,
        output=str(tmp_path / "tile.parquet"), no_tsv=True,
    )
    _cmd_dmr(args)
    got = pl.read_parquet(args.output)
    assert got.sort("chrom").equals(expected.sort("chrom")), (
        "CLI tile output must equal tl.dmr tile's qvalue-filtered set exactly"
    )


def test_cli_tile_none_disables_q_filter(tmp_path, monkeypatch):
    """``min_mean_qvalue=None`` keeps every tile (filter disabled)."""
    import epykit.cli as ep_cli
    import epykit.dmr as ep_dmr

    tile_frame = pl.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [100, 100],
            "end": [200, 200],
            "pvalue": [1e-6, 0.04],
            "qvalue": [1e-6, 0.04],
            "meth_diff": [0.4, 0.2],
        }
    )
    monkeypatch.setattr(ep_dmr, "call_dmr_tile_based", lambda *a, **k: tile_frame)
    monkeypatch.setattr(
        ep_cli, "_read_samplesheet_groups", lambda *a, **k: (["t0", "t1"], ["c0", "c1"])
    )
    monkeypatch.setattr(ep_cli, "_cli_n1_and_footgun_checks", lambda *a, **k: None)

    args = argparse.Namespace(
        method="tile", empirical_fdr=False,
        methylstore="store", samplesheet="ss.csv",
        treatment_group="A", control_group="B",
        tile_size_bp=1000, test="lr", min_cpgs_per_tile=3,
        alpha=0.05, min_abs_meth_diff=0.1, unite=True,
        min_samples_treatment=2, min_samples_control=2,
        min_mean_qvalue=None,
        output=str(tmp_path / "tile.parquet"), no_tsv=True,
    )
    _cmd_dmr(args)
    got = pl.read_parquet(args.output)
    assert len(got) == 2, "None must disable the q-filter (both tiles kept)"
