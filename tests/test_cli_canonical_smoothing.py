"""The ``--canonical-only`` and ``--smoothing`` CLI flags (1.2).

``epykit convert``, ``epykit dmc`` and ``epykit dmr --method tile`` take
``--canonical-only``; ``epykit dmc`` takes ``--smoothing`` and
``--smoothing-span-bp``. Each handler runs in-process through the real parser
and its result is compared with the API call the flag forwards to: the
scaffold cohort of ``test_canonical_chrom_filter`` (chr1 plus one unplaced
contig) proves the filter, the synthetic truth-table cohort proves the
smoothing. Combinations the engines do not consume exit with a usage error
instead of being dropped.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit import dmr as dmr_mod
from epykit.cli import build_parser
from epykit.convert import convert_sample
from epykit.dmr import call_dmr_tile_based
from tests.fixtures.synth import SimConfig, generate
from tests.test_canonical_chrom_filter import SCAFFOLD, _partitions, _read, _write_cohort

_CHROMS_LOGGER = "epykit._chroms"


def _run_cli(*argv: object) -> argparse.Namespace:
    """Parse ``argv`` with the real parser and run the handler in-process."""
    args = build_parser().parse_args([str(a) for a in argv])
    args.func(args)
    return args


def _sites(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.select(["chrom", "pos", "pvalue", "qvalue"])
        .with_columns(pl.col("chrom").cast(pl.Utf8))
        .sort(["chrom", "pos"])
    )


def _chroms(df: pl.DataFrame) -> set[str]:
    return set(df["chrom"].cast(pl.Utf8).unique().to_list()) if df.height else set()


def _assert_same_sites(cli: pl.DataFrame, api: pl.DataFrame) -> None:
    assert cli.height == api.height, (cli.height, api.height)
    assert cli.select(["chrom", "pos"]).equals(api.select(["chrom", "pos"]))
    for col in ("pvalue", "qvalue"):
        assert np.allclose(cli[col].to_numpy(), api[col].to_numpy(), atol=1e-12, equal_nan=True), col


_GROUPS = ("--treatment-group", "treatment", "--control-group", "control")


@pytest.fixture
def cohort(tmp_path) -> str:
    return _write_cohort(tmp_path)


@pytest.fixture(scope="module")
def synth_sheet(tmp_path_factory) -> str:
    """Varying counts along chr1 (50-400 bp apart), so box smoothing changes
    the test; the constant-count scaffold cohort cannot show that."""
    cfg = SimConfig(
        n_per_group=3,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    return generate(cfg, tmp_path_factory.mktemp("cli_smoothing_synth"))["samplesheet"]


# Parser wiring


def test_flags_default_off_and_parse():
    parser = build_parser()
    conv = ["convert", "--input", "x.cov", "--sample-id", "s", "--output-dir", "o"]
    dmc = ["dmc", "--methylstore", "m", "--samplesheet", "s", *_GROUPS, "--output", "o"]
    dmr = ["dmr", "--output", "o"]

    assert parser.parse_args(conv).canonical_only is False
    bare = parser.parse_args(dmc)
    assert (bare.canonical_only, bare.smoothing, bare.smoothing_span_bp) == (False, False, 500)
    assert parser.parse_args(dmr).canonical_only is False

    assert parser.parse_args([*conv, "--canonical-only"]).canonical_only is True
    on = parser.parse_args([*dmc, "--canonical-only", "--smoothing", "--smoothing-span-bp", "200"])
    assert (on.canonical_only, on.smoothing, on.smoothing_span_bp) == (True, True, 200)
    assert parser.parse_args([*dmr, "--method", "tile", "--canonical-only"]).canonical_only is True


# convert


def test_convert_canonical_only_matches_convert_sample(cohort, tmp_path):
    cov = Path(cohort).read_text().splitlines()[1].split(",")[2]  # sample t0

    _run_cli("convert", "--input", cov, "--sample-id", "t0", "--output-dir", tmp_path / "cli")
    assert _partitions(str(tmp_path / "cli")) == {"t0": {"chr1", SCAFFOLD}}

    _run_cli(
        "convert",
        "--input",
        cov,
        "--sample-id",
        "t0",
        "--output-dir",
        tmp_path / "cli_canon",
        "--canonical-only",
    )
    convert_sample(cov, "t0", str(tmp_path / "api_canon"), canonical_only=True)
    assert _partitions(str(tmp_path / "cli_canon")) == {"t0": {"chr1"}}
    assert _partitions(str(tmp_path / "api_canon")) == {"t0": {"chr1"}}
    part = Path("sample=t0") / "chrom=chr1" / "part-0.parquet"
    assert pl.read_parquet(tmp_path / "cli_canon" / part).equals(
        pl.read_parquet(tmp_path / "api_canon" / part)
    )


# dmc, binary path


def test_dmc_canonical_only_matches_tl_dmc(cohort, tmp_path):
    md_api = _read(cohort, tmp_path / "api")
    ep.tl.dmc(
        md_api,
        test="lr",
        canonical_only=True,
        min_samples_treatment=0,
        min_samples_control=0,
        tsv=False,
    )
    api = _sites(md_api.dmc)
    assert _chroms(api) == {"chr1"}

    md_cli = _read(cohort, tmp_path / "cli")
    out = tmp_path / "cli_dmc.parquet"
    _run_cli(
        "dmc",
        "--methylstore",
        md_cli.store,
        "--samplesheet",
        cohort,
        *_GROUPS,
        "--output",
        out,
        "--test",
        "lr",
        "--canonical-only",
        "--no-tsv",
    )
    cli = _sites(pl.read_parquet(out))
    assert _chroms(cli) == {"chr1"}
    _assert_same_sites(cli, api)

    # Without the flag the same store tests both contigs.
    out_all = tmp_path / "cli_dmc_all.parquet"
    _run_cli(
        "dmc", "--methylstore", md_cli.store, "--samplesheet", cohort, *_GROUPS,
        "--output", out_all, "--test", "lr", "--no-tsv",
    )  # fmt: skip
    assert _chroms(pl.read_parquet(out_all)) == {"chr1", SCAFFOLD}


def test_dmc_smoothing_matches_tl_dmc(synth_sheet, tmp_path):
    span = 600  # non-default, so the span flag is proven as well as the switch
    md_api = ep.read_bismark(
        synth_sheet,
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "api"),
    )
    common = dict(test="lr", min_samples_treatment=0, min_samples_control=0, tsv=False)
    ep.tl.dmc(md_api, smoothing=True, smoothing_span_bp=span, **common)
    smoothed = _sites(md_api.dmc)
    ep.tl.dmc(md_api, **common)
    raw = _sites(md_api.dmc)
    assert not np.allclose(smoothed["pvalue"].to_numpy(), raw["pvalue"].to_numpy(), equal_nan=True), (
        "premise: smoothing must change the test on this fixture"
    )

    md_cli = ep.read_bismark(
        synth_sheet,
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "cli"),
    )
    out = tmp_path / "cli_dmc.parquet"
    _run_cli(
        "dmc", "--methylstore", md_cli.store, "--samplesheet", synth_sheet, *_GROUPS,
        "--output", out, "--test", "lr", "--min-samples-treatment", "0",
        "--min-samples-control", "0", "--smoothing", "--smoothing-span-bp", span, "--no-tsv",
    )  # fmt: skip
    _assert_same_sites(_sites(pl.read_parquet(out)), smoothed)


# dmc, formula / contrast path


def test_dmc_contrast_path_canonical_only_matches_tl_dmc(cohort, tmp_path):
    md_api = _read(cohort, tmp_path / "api")
    ep.tl.dmc(md_api, formula="~ group", contrast="group", canonical_only=True, tsv=False)
    api = _sites(md_api.dmc)
    assert _chroms(api) == {"chr1"}

    out = tmp_path / "cli_contrast.parquet"
    _run_cli(
        "dmc", "--methylstore", tmp_path / "cli", "--samplesheet", cohort, *_GROUPS,
        "--formula", "~ group", "--contrast", "group", "--canonical-only", "--output", out,
    )  # fmt: skip
    cli = _sites(pl.read_parquet(out))
    assert _chroms(cli) == {"chr1"}
    _assert_same_sites(cli, api)


# dmr, tile


def test_dmr_tile_canonical_only_matches_tile_caller(cohort, tmp_path):
    md = _read(cohort, tmp_path / "store")
    tile = dict(tile_size_bp=200, min_cpgs_per_tile=2)
    api = call_dmr_tile_based(
        methylstore_path=md.store,
        samples_treatment=md.treatment_ids,
        samples_control=md.control_ids,
        unite=False,
        min_samples_treatment=0,
        min_samples_control=0,
        canonical_only=True,
        **tile,
    )
    assert api.height > 0 and _chroms(api) == {"chr1"}

    out = tmp_path / "tile.parquet"
    _run_cli(
        "dmr", "--method", "tile", "--methylstore", md.store, "--samplesheet", cohort, *_GROUPS,
        "--tile-size-bp", tile["tile_size_bp"], "--min-cpgs-per-tile", tile["min_cpgs_per_tile"],
        "--min-mean-qvalue", "1.1", "--canonical-only", "--output", out, "--no-tsv",
    )  # fmt: skip
    cli = pl.read_parquet(out)
    assert _chroms(cli) == {"chr1"}
    assert cli.equals(api)


def test_dmr_tile_permutations_receive_the_canonical_list(cohort, tmp_path, monkeypatch, caplog):
    """The observed tile run and every ``--empirical-fdr`` permutation get the
    same resolved list, and the audit line fires once for the whole command."""
    md = _read(cohort, tmp_path / "store")
    seen: list[list[str] | None] = []
    real = dmr_mod.call_dmr_tile_based

    def recording(*args, **kwargs):
        seen.append(kwargs.get("chromosomes"))
        return real(*args, **kwargs)

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", recording)

    n_perm = 3
    out = tmp_path / "tile.parquet"
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        _run_cli(
            "dmr", "--method", "tile", "--methylstore", md.store, "--samplesheet", cohort,
            *_GROUPS, "--tile-size-bp", "200", "--min-cpgs-per-tile", "2", "--canonical-only",
            "--empirical-fdr", "--n-perm", n_perm, "--output", out, "--no-tsv",
        )  # fmt: skip

    result = pl.read_parquet(out)
    assert result.height > 0 and "empirical_qvalue" in result.columns
    assert seen == [["chr1"]] * (1 + n_perm)
    assert len([r for r in caplog.records if r.name == _CHROMS_LOGGER]) == 1, caplog.text


# Unsupported combinations


@pytest.mark.parametrize("method", ["chain_merge", "sliding_window", "segment"])
def test_dmr_non_tile_methods_reject_canonical_only(method, tmp_path):
    with pytest.raises(SystemExit, match="--method tile only") as excinfo:
        _run_cli(
            "dmr", "--method", method, "--dmc-results", "missing.parquet",
            "--canonical-only", "--output", tmp_path / "o.parquet",
        )  # fmt: skip
    assert "epykit dmc --canonical-only" in str(excinfo.value)


def _dmc_argv(cohort, tmp_path, *extra: object) -> list[object]:
    # The store never exists: every rejection fires before it is opened.
    return [
        "dmc", "--methylstore", tmp_path / "missing", "--samplesheet", cohort, *_GROUPS,
        "--output", tmp_path / "o.parquet", *extra,
    ]  # fmt: skip


def test_dmc_smoothing_rejected_with_a_non_lr_engine(cohort, tmp_path):
    with pytest.raises(SystemExit, match="option of the lr engine") as excinfo:
        _run_cli(*_dmc_argv(cohort, tmp_path, "--test", "welch_t", "--smoothing"))
    assert "--test welch_t does not consume it" in str(excinfo.value)


def test_dmc_smoothing_rejected_after_allow_n1_resolves_to_fisher(cohort, tmp_path):
    """``--allow-n1`` turns ``--test lr`` into fisher at n=1; the smoothing
    check runs after that resolution."""
    lines = Path(cohort).read_text().splitlines()
    one_each = tmp_path / "one_each.csv"
    one_each.write_text("\n".join([lines[0], lines[1], lines[3]]) + "\n")
    with pytest.raises(SystemExit, match="--test fisher does not consume it"):
        _run_cli(*_dmc_argv(one_each, tmp_path, "--test", "lr", "--allow-n1", "--smoothing"))


def test_dmc_smoothing_rejected_on_the_contrast_path(cohort, tmp_path):
    with pytest.raises(SystemExit, match="--formula / --contrast path does not consume it"):
        _run_cli(
            *_dmc_argv(cohort, tmp_path, "--formula", "~ group", "--contrast", "group", "--smoothing")
        )


def test_dmc_smoothing_requires_a_positive_span(cohort, tmp_path):
    with pytest.raises(SystemExit, match="positive number of base pairs"):
        _run_cli(*_dmc_argv(cohort, tmp_path, "--smoothing", "--smoothing-span-bp", "0"))
