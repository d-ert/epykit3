"""Opt-in canonical chromosome filtering at ingestion, in DMC and in tile DMR.

The session ``synth`` fixture is canonical-only (chr1..chr5), so it cannot
prove the filter. These tests hand-build a tiny 2 + 2 cohort with CpGs on
``chr1`` plus one real GRCh38 unplaced contig (``chrUn_KI270742v1``, the
contig that topped the GSE263850 DMR list before the filter existed) and pin
the ``canonical_only`` contract of ``read_*`` / ``convert``, of ``tl.dmc`` /
``process_chromosomes_dmc`` on the binary and the contrast path, and of
``tl.dmr(method="tile")``. The CLI flags are covered in
``test_cli_canonical_smoothing.py`` on the same cohort.
"""

from __future__ import annotations

import gzip
import inspect
import json
import logging
from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit import _dmc_stages as stages_mod
from epykit import dmc as dmc_mod
from epykit import dmr as dmr_mod
from epykit import tl as tl_mod
from epykit.convert import _can_reuse_sample, _manifest_path, ensure_converted_sample
from epykit.dmc import process_chromosomes_dmc
from epykit.dmr import _DMR_TILE_SCHEMA, call_dmr_tile_based
from epykit.pl._compute import compute_manhattan_data

SCAFFOLD = "chrUn_KI270742v1"
_CHROMS = ("chr1", SCAFFOLD)
_CHROMS_LOGGER = "epykit._chroms"

# Two treatment + two control samples with a clear hyper-vs-hypo difference
# and per-sample jitter so the LR engine sees within-group variance.
_TREAT_COUNTS = [(18, 2), (16, 4)]  # (N_meth, N_unmeth) per treatment sample
_CTRL_COUNTS = [(3, 17), (5, 15)]  # per control sample
_POSITIONS = list(range(1001, 1001 + 60 * 20, 20))  # 60 CpGs/chrom, 20 bp apart
_TILE_KW = {"tile_size_bp": 200, "min_cpgs_per_tile": 2}


def _cov_lines(
    counts: tuple[int, int], *, one_based: bool, chroms: tuple[str, ...] = _CHROMS
) -> list[str]:
    m, u = counts
    pct = 100.0 * m / (m + u)
    return [
        f"{chrom}\t{pos if one_based else pos - 1}\t{pos}\t{pct:.2f}\t{m}\t{u}"
        for chrom in chroms
        for pos in _POSITIONS
    ]


def _write_bismark(path: Path, counts: tuple[int, int], chroms: tuple[str, ...] = _CHROMS) -> None:
    path.write_text("\n".join(_cov_lines(counts, one_based=True, chroms=chroms)) + "\n")


def _write_methyldackel(
    path: Path, counts: tuple[int, int], chroms: tuple[str, ...] = _CHROMS
) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        fh.write('track type="bedGraph" description="CpG methylation levels"\n')
        fh.write("\n".join(_cov_lines(counts, one_based=False, chroms=chroms)) + "\n")


def _write_combined_bed(
    path: Path, counts: tuple[int, int], chroms: tuple[str, ...] = _CHROMS
) -> None:
    m, u = counts
    t = m + u
    pct = 100.0 * m / t
    rows = [
        f"{chrom}\t{pos - 1}\t{pos}\t{m}\t{t}\t{pct:.2f}\t0\t0\t0.00\t{m}\t{t}\t{pct:.2f}"
        for chrom in chroms
        for pos in _POSITIONS
    ]
    path.write_text("\n".join(rows) + "\n")


_FORMATS = {
    "bismark": (_write_bismark, ".cov", ep.read_bismark),
    "methyldackel": (_write_methyldackel, ".bedGraph.gz", ep.read_methyldackel),
    "combined_strand_bed": (_write_combined_bed, ".bed", ep.read_combined_strand_bed),
}


def _write_cohort(tmp_path: Path, fmt: str = "bismark", chroms: tuple[str, ...] = _CHROMS) -> str:
    """Write the 2 + 2 cohort in ``fmt`` and return the samplesheet path."""
    writer, suffix, _ = _FORMATS[fmt]
    src_dir = tmp_path / f"src_{fmt}"
    src_dir.mkdir(exist_ok=True)
    sheet = ["sample_id,group,path"]
    for group, counts_list in (("treatment", _TREAT_COUNTS), ("control", _CTRL_COUNTS)):
        for i, counts in enumerate(counts_list):
            f = src_dir / f"{group[0]}{i}{suffix}"
            writer(f, counts, chroms)
            sheet.append(f"{group[0]}{i},{group},{f}")
    sheet_path = tmp_path / f"samplesheet_{fmt}.csv"
    sheet_path.write_text("\n".join(sheet) + "\n")
    return str(sheet_path)


def _read(sheet: str, store_dir: Path, fmt: str = "bismark", **kwargs):
    reader = _FORMATS[fmt][2]
    return reader(
        sheet,
        treatment_group="treatment",
        control_group="control",
        assembly="hg38",
        store_dir=str(store_dir),
        **kwargs,
    )


def _partitions(store: str) -> dict[str, set[str]]:
    """{sample_id: set of chromosome partitions on disk}."""
    out: dict[str, set[str]] = {}
    for sample_dir in Path(store).glob("sample=*"):
        out[sample_dir.name.removeprefix("sample=")] = {
            d.name.removeprefix("chrom=") for d in sample_dir.glob("chrom=*")
        }
    return out


def _manifest_flag(store: str, sample_id: str) -> object:
    payload = json.loads(_manifest_path(Path(store) / f"sample={sample_id}").read_text())
    return payload.get("canonical_only", "missing")


def _chroms_records(caplog) -> list[logging.LogRecord]:
    return [r for r in caplog.records if r.name == _CHROMS_LOGGER]


def _dmr_chroms(md) -> set[str]:
    dmr = md.uns["dmr"]
    return set(dmr["chrom"].unique().to_list()) if dmr.height else set()


def _dmc_chroms(md) -> set[str]:
    dmc = md.dmc
    return set(dmc["chrom"].cast(pl.Utf8).unique().to_list()) if dmc.height else set()


def _frame_chroms(df: pl.DataFrame) -> set[str]:
    return set(df["chrom"].cast(pl.Utf8).unique().to_list()) if df.height else set()


@pytest.fixture
def cohort_sheet(tmp_path) -> str:
    return _write_cohort(tmp_path)


@pytest.fixture
def scaffold_md(cohort_sheet, tmp_path):
    """A MethylData whose raw store holds chr1 and the scaffold for every sample."""
    md = _read(cohort_sheet, tmp_path / "store")
    ep.pp.set_unite_type(md, type="intersect")
    return md


# Ingestion


def test_ingestion_default_keeps_every_contig(cohort_sheet, tmp_path, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        md = _read(cohort_sheet, tmp_path / "store")

    parts = _partitions(md.store)
    assert set(parts) == {"t0", "t1", "c0", "c1"}
    assert all(chroms == {"chr1", SCAFFOLD} for chroms in parts.values()), parts
    assert _manifest_flag(md.store, "t0") is False
    assert _chroms_records(caplog) == []


@pytest.mark.parametrize("fmt", sorted(_FORMATS))
def test_ingestion_canonical_only_drops_only_the_scaffold(fmt, tmp_path, caplog):
    """Every reader forwards the option; only the non-canonical contig goes."""
    sheet = _write_cohort(tmp_path, fmt)
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        md = _read(sheet, tmp_path / f"store_{fmt}", fmt, canonical_only=True)

    parts = _partitions(md.store)
    assert set(parts) == {"t0", "t1", "c0", "c1"}
    assert all(chroms == {"chr1"} for chroms in parts.values()), parts
    assert _manifest_flag(md.store, "t0") is True

    # One summary line per converted sample, naming the dropped contig.
    records = _chroms_records(caplog)
    assert len(records) == 4, caplog.text
    assert all(SCAFFOLD in r.getMessage() and "dropping 1 contig(s)" in r.getMessage() for r in records)


def test_ingestion_false_true_false_rebuilds_the_partition_set(cohort_sheet, tmp_path):
    """Flipping the option on a cached store regenerates every sample, and the
    replaced sample directory leaves no stale partition behind for a glob."""
    store_dir = tmp_path / "store"

    md_default = _read(cohort_sheet, store_dir)
    assert all(c == {"chr1", SCAFFOLD} for c in _partitions(md_default.store).values())
    n_sites_all = md_default.uns["n_sites_raw"]

    md_canon = _read(cohort_sheet, store_dir, canonical_only=True)
    assert md_canon.store == md_default.store
    assert all(c == {"chr1"} for c in _partitions(md_canon.store).values())
    assert not (Path(md_canon.store) / "sample=t0" / f"chrom={SCAFFOLD}").exists()
    assert _manifest_flag(md_canon.store, "t0") is True
    assert md_canon.uns["n_sites_raw"] * 2 == n_sites_all

    md_back = _read(cohort_sheet, store_dir, canonical_only=False)
    assert all(c == {"chr1", SCAFFOLD} for c in _partitions(md_back.store).values())
    assert _manifest_flag(md_back.store, "t0") is False
    assert md_back.uns["n_sites_raw"] == n_sites_all


def test_legacy_manifest_without_key_is_reusable_for_false_only(tmp_path):
    cov = tmp_path / "s.cov"
    _write_bismark(cov, _TREAT_COUNTS[0])
    store = tmp_path / "raw"
    assert ensure_converted_sample(str(cov), "s", str(store)) is True

    # Simulate a manifest written before the key existed.
    manifest_path = _manifest_path(store / "sample=s")
    payload = json.loads(manifest_path.read_text())
    del payload["canonical_only"]
    manifest_path.write_text(json.dumps(payload))

    reuse = dict(input_path=cov, sample_dir=store / "sample=s", row_group_size=1_000_000)
    assert _can_reuse_sample(**reuse, canonical_only=False) is True
    assert _can_reuse_sample(**reuse, canonical_only=True) is False

    assert ensure_converted_sample(str(cov), "s", str(store)) is False  # cache hit
    assert ensure_converted_sample(str(cov), "s", str(store), canonical_only=True) is True
    assert _partitions(str(store))["s"] == {"chr1"}
    assert _manifest_flag(str(store), "s") is True


# Tile DMR


def test_canonical_only_is_keyword_only():
    for fn in (call_dmr_tile_based, ep.tl.dmr):
        param = inspect.signature(fn).parameters["canonical_only"]
        assert param.kind is inspect.Parameter.KEYWORD_ONLY, fn.__name__
        assert param.default is False, fn.__name__


def test_tile_default_keeps_scaffold(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmr(scaffold_md, method="tile", **_TILE_KW)

    assert _dmr_chroms(scaffold_md) == {"chr1", SCAFFOLD}
    assert scaffold_md.uns["dmr_params"]["canonical_only"] is False
    assert _chroms_records(caplog) == []


def test_tile_canonical_only_drops_scaffold_and_logs_once(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmr(scaffold_md, method="tile", canonical_only=True, **_TILE_KW)

    assert _dmr_chroms(scaffold_md) == {"chr1"}
    assert scaffold_md.uns["dmr_params"]["canonical_only"] is True
    records = _chroms_records(caplog)
    assert len(records) == 1, caplog.text
    msg = records[0].getMessage()
    assert "[dmr/tile]" in msg and SCAFFOLD in msg


def test_tile_explicit_chromosomes_take_precedence(scaffold_md, caplog):
    """An explicit list, even one naming only the scaffold, is used verbatim."""
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmr(
            scaffold_md, method="tile", chromosomes=[SCAFFOLD], canonical_only=True, **_TILE_KW
        )
    assert _dmr_chroms(scaffold_md) == {SCAFFOLD}
    assert _chroms_records(caplog) == []


def test_tile_explicit_empty_list_takes_precedence(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmr(scaffold_md, method="tile", chromosomes=[], canonical_only=True, **_TILE_KW)
    dmr = scaffold_md.uns["dmr"]
    assert dmr.height == 0
    assert dmr.schema == pl.Schema(_DMR_TILE_SCHEMA)
    assert _chroms_records(caplog) == []


def test_engine_filters_auto_detected_chromosomes_only(scaffold_md):
    kwargs = dict(
        methylstore_path=scaffold_md.store,
        samples_treatment=scaffold_md.treatment_ids,
        samples_control=scaffold_md.control_ids,
        **_TILE_KW,
    )
    auto = call_dmr_tile_based(**kwargs, canonical_only=True)
    assert set(auto["chrom"].unique().to_list()) == {"chr1"}
    explicit = call_dmr_tile_based(**kwargs, chromosomes=[SCAFFOLD], canonical_only=True)
    assert set(explicit["chrom"].unique().to_list()) == {SCAFFOLD}


@pytest.mark.parametrize("method", ["chain_merge", "sliding_window", "segment"])
def test_non_tile_methods_reject_canonical_only(scaffold_md, method):
    """DMC-derived callers refuse the option before any DMC lookup and point
    at upstream filtering; nothing is stored on md."""
    with pytest.raises(ValueError, match=r"method='tile' only") as excinfo:
        ep.tl.dmr(scaffold_md, method=method, canonical_only=True)
    assert "ep.tl.dmc" in str(excinfo.value)
    assert "dmr" not in scaffold_md.uns


def test_tile_permutations_use_the_observed_universe(scaffold_md, monkeypatch, caplog):
    """The observed run and every permutation receive the same resolved
    chromosome list, and the audit line is emitted once for the whole call."""
    seen: list[list[str] | None] = []
    real = dmr_mod.call_dmr_tile_based

    def recording(*args, **kwargs):
        seen.append(kwargs.get("chromosomes"))
        return real(*args, **kwargs)

    # tl.py binds the engine name at import; the harness calls the dmr module global.
    monkeypatch.setattr(tl_mod, "call_dmr_tile_based", recording)
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", recording)

    n_perm = 3
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmr(
            scaffold_md,
            method="tile",
            canonical_only=True,
            empirical_fdr=True,
            n_perm=n_perm,
            **_TILE_KW,
        )

    dmr = scaffold_md.uns["dmr"]
    assert dmr.height > 0, "the fixture must yield observed tiles for permutations to run"
    assert "empirical_qvalue" in dmr.columns
    assert len(seen) == 1 + n_perm
    assert seen == [["chr1"]] * (1 + n_perm)
    assert len(_chroms_records(caplog)) == 1, caplog.text


# DMC


def test_dmc_canonical_only_is_keyword_only():
    for fn in (process_chromosomes_dmc, ep.tl.dmc):
        param = inspect.signature(fn).parameters["canonical_only"]
        assert param.kind is inspect.Parameter.KEYWORD_ONLY, fn.__name__
        assert param.default is False, fn.__name__


def test_dmc_default_keeps_scaffold(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(scaffold_md, test="lr", tsv=False)

    assert _dmc_chroms(scaffold_md) == {"chr1", SCAFFOLD}
    assert scaffold_md.uns["dmc"]["canonical_only"] is False
    assert _chroms_records(caplog) == []


def test_dmc_canonical_only_drops_scaffold_and_logs_once(scaffold_md, caplog):
    """The filter shapes the store the engine writes, so the q-values are
    corrected over chr1 alone; one audit line names the dropped contig."""
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(scaffold_md, test="lr", canonical_only=True, tsv=False)

    assert _dmc_chroms(scaffold_md) == {"chr1"}
    assert scaffold_md.uns["dmc"]["canonical_only"] is True
    assert set(scaffold_md.dmc_store.chroms()) == {"chr1"}
    records = _chroms_records(caplog)
    assert len(records) == 1, caplog.text
    msg = records[0].getMessage()
    assert "[dmc]" in msg and SCAFFOLD in msg


def test_dmc_explicit_chromosomes_take_precedence(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(scaffold_md, test="lr", chromosomes=[SCAFFOLD], canonical_only=True, tsv=False)
    assert _dmc_chroms(scaffold_md) == {SCAFFOLD}
    assert _chroms_records(caplog) == []


def test_dmc_explicit_empty_list_takes_precedence(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(scaffold_md, test="lr", chromosomes=[], canonical_only=True, tsv=False)
    assert scaffold_md.uns["dmc"]["n_sites"] == 0
    assert scaffold_md.dmc.height == 0
    assert _chroms_records(caplog) == []


def test_dmc_canonical_only_with_no_canonical_contig_is_a_valid_empty_run(tmp_path):
    """A store holding only the scaffold filters down to nothing: the run
    completes with an empty result rather than failing."""
    sheet = _write_cohort(tmp_path, chroms=(SCAFFOLD,))
    md = _read(sheet, tmp_path / "store")
    ep.pp.set_unite_type(md, type="intersect")

    ep.tl.dmc(md, test="lr", canonical_only=True, tsv=False)

    assert md.uns["dmc"]["n_sites"] == 0
    assert md.dmc.height == 0
    assert md.dmc_store.chroms() == []


def test_dmc_contrast_path_honours_canonical_only(scaffold_md, caplog):
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(scaffold_md, formula="~ group", contrast="group", canonical_only=True, tsv=False)

    assert scaffold_md.uns["dmc"]["last_key"] == "dmc_glm_contrast"
    assert _dmc_chroms(scaffold_md) == {"chr1"}
    assert scaffold_md.uns["dmc"]["canonical_only"] is True
    assert len(_chroms_records(caplog)) == 1


def test_dmc_engine_filters_auto_detected_chromosomes_only(scaffold_md):
    """``process_chromosomes_dmc`` honours the option for direct callers, and
    the resolved list is part of the low-level cache identity: the default
    run after a filtered one recomputes instead of serving the chr1-only
    store."""
    kwargs = dict(
        methylstore_path=scaffold_md.store,
        samples_treatment=scaffold_md.treatment_ids,
        samples_control=scaffold_md.control_ids,
        test="lr",
    )
    auto = process_chromosomes_dmc(**kwargs, canonical_only=True)
    assert _frame_chroms(auto) == {"chr1"}
    explicit = process_chromosomes_dmc(**kwargs, chromosomes=[SCAFFOLD], canonical_only=True)
    assert _frame_chroms(explicit) == {SCAFFOLD}
    default = process_chromosomes_dmc(**kwargs)
    assert _frame_chroms(default) == {"chr1", SCAFFOLD}


def test_dmc_resume_false_true_false_invalidates_and_identical_call_resumes(scaffold_md):
    md = scaffold_md
    ep.tl.dmc(md, test="lr", resumable=True, tsv=False)
    assert md.uns["dmc"]["resumed"] is False
    assert _dmc_chroms(md) == {"chr1", SCAFFOLD}

    ep.tl.dmc(md, test="lr", resumable=True, canonical_only=True, tsv=False)
    assert md.uns["dmc"]["resumed"] is False
    assert _dmc_chroms(md) == {"chr1"}

    ep.tl.dmc(md, test="lr", resumable=True, canonical_only=False, tsv=False)
    assert md.uns["dmc"]["resumed"] is False
    assert _dmc_chroms(md) == {"chr1", SCAFFOLD}

    ep.tl.dmc(md, test="lr", resumable=True, canonical_only=False, tsv=False)
    assert md.uns["dmc"]["resumed"] is True
    assert md.uns["dmc"]["canonical_only"] is False
    assert _dmc_chroms(md) == {"chr1", SCAFFOLD}


def test_dmc_permutations_use_the_observed_universe(scaffold_md, monkeypatch, caplog):
    """The observed engine run and every ``empirical_fdr`` permutation
    receive the same resolved chromosome list; the audit line fires once."""
    seen: list[list[str] | None] = []
    real = dmc_mod.process_chromosomes_dmc

    def recording(*args, **kwargs):
        seen.append(kwargs.get("chromosomes"))
        return real(*args, **kwargs)

    # The stage binds the engine name at import; the permutation harness in
    # dmc.py calls the module global.
    monkeypatch.setattr(stages_mod, "process_chromosomes_dmc", recording)
    monkeypatch.setattr(dmc_mod, "process_chromosomes_dmc", recording)

    n_perm = 3
    with caplog.at_level(logging.INFO, logger=_CHROMS_LOGGER):
        ep.tl.dmc(
            scaffold_md,
            test="lr",
            canonical_only=True,
            empirical_fdr=True,
            n_perm=n_perm,
            tsv=False,
        )

    assert "empirical_qvalue" in scaffold_md.dmc.columns
    assert len(seen) == 1 + n_perm
    assert seen == [["chr1"]] * (1 + n_perm)
    assert len(_chroms_records(caplog)) == 1, caplog.text


# Plot order


def test_manhattan_order_follows_the_shared_ucsc_list():
    dmc = pl.DataFrame(
        {
            "chrom": [SCAFFOLD, "chrX", "chr2", "chr1", "chr10"],
            "pos": [10, 20, 30, 40, 50],
            "pvalue": [0.5, 0.01, 0.2, 0.001, 0.3],
        }
    )
    order = [b["chrom"] for b in compute_manhattan_data(None, dmc=dmc).chrom_blocks]
    assert order == ["chr1", "chr2", "chr10", "chrX"]
    order_all = [
        b["chrom"] for b in compute_manhattan_data(None, dmc=dmc, canonical_only=False).chrom_blocks
    ]
    assert order_all == ["chr1", "chr2", "chr10", "chrX", SCAFFOLD]
