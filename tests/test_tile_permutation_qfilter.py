"""Observed and permuted tiles must pass the same q-value post-filter."""

import polars as pl
import pytest

from epykit import MethylData, tl
from epykit import dmr as dmr_mod
from epykit.cli import build_parser


@pytest.mark.parametrize(("cutoff", "expected_fdr"), [(None, 1.0), (0.01, 0.0), (0.05, 1.0)])
def test_tile_region_fdr_uses_observed_postfilter(tmp_path, monkeypatch, cutoff, expected_fdr):
    md = MethylData(
        obs=pl.DataFrame({"sample_id": ["t1", "t2", "c1", "c2"], "treatment": [1, 1, 0, 0]}),
        store=str(tmp_path),
    )
    observed = pl.DataFrame({"pvalue": [0.001], "qvalue": [0.005]})
    null = pl.DataFrame({"pvalue": [0.001], "qvalue": [0.04]})
    monkeypatch.setattr(tl, "call_dmr_tile_based", lambda **kw: observed)
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", lambda **kw: null)

    tl.dmr(
        md,
        method="tile",
        empirical_fdr=True,
        fdr_method="region",
        min_mean_qvalue=cutoff,
        n_perm=20,
        tsv=False,
    )

    assert md.uns["dmr"]["empirical_qvalue"][0] == expected_fdr
    assert md.uns["dmr"]["empirical_fdr_set"][0] == expected_fdr


def test_cli_passes_tile_postfilter_to_permutations(tmp_path, monkeypatch):
    sheet = tmp_path / "samples.csv"
    sheet.write_text("sample_id,group\nt1,T\nt2,T\nc1,C\nc2,C\n")
    tiles = pl.DataFrame({"pvalue": [0.001, 0.002], "qvalue": [0.005, 0.04]})
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", lambda **kw: tiles)
    received = {}

    def empirical(**kwargs):
        received.update(kwargs)
        return kwargs["observed_dmr"]

    monkeypatch.setattr(dmr_mod, "empirical_fdr_for_dmr", empirical)
    args = build_parser().parse_args(
        [
            "dmr",
            "--method",
            "tile",
            "--methylstore",
            str(tmp_path),
            "--samplesheet",
            str(sheet),
            "--treatment-group",
            "T",
            "--control-group",
            "C",
            "--output",
            str(tmp_path / "dmr.parquet"),
            "--no-tsv",
            "--empirical-fdr",
            "--min-mean-qvalue",
            "0.01",
        ]
    )
    args.func(args)

    assert received["min_mean_qvalue"] == 0.01
    assert received["observed_dmr"].equals(tiles.head(1))
