"""Integration tests for the canonical_only DMC/DMR default + ingestion opt-in.

The session ``synth`` fixture is canonical-only (chr1..chr5), so it cannot
exercise scaffold filtering. These tests hand-build a tiny Bismark ``.cov``
cohort with CpGs on ``chr1`` plus one real GRCh38 unplaced contig
(``chrUn_KI270742v1`` -- the contig that was the top DMR on the GSE263850 real
cohort before this filter existed) and assert the canonical_only contract.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

SCAFFOLD = "chrUn_KI270742v1"

# Two treatment + two control samples; slight per-sample jitter so the LR
# engine sees non-zero within-group variance. Clear hyper-vs-hypo difference.
_TREAT_COUNTS = [(18, 2), (16, 4)]   # (N_meth, N_unmeth) per treatment sample
_CTRL_COUNTS = [(3, 17), (5, 15)]    # per control sample
_POSITIONS = list(range(1001, 1001 + 60 * 20, 20))  # 60 CpGs/chrom, 20 bp apart


def _write_cov(path: Path, counts: tuple[int, int]) -> None:
    """Write a tiny 1-based Bismark .cov with chr1 + scaffold CpGs."""
    m, u = counts
    cov = m + u
    pct = 100.0 * m / cov
    lines = [
        f"{chrom}\t{pos}\t{pos}\t{pct:.2f}\t{m}\t{u}"
        for chrom in ("chr1", SCAFFOLD)
        for pos in _POSITIONS
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_cohort(tmp_path: Path) -> str:
    """Write 2+2 .cov files and a samplesheet; return the samplesheet path."""
    cov_dir = tmp_path / "cov"
    cov_dir.mkdir(exist_ok=True)
    sheet = ["sample_id,group,path"]
    for i, counts in enumerate(_TREAT_COUNTS):
        f = cov_dir / f"t{i}.cov"
        _write_cov(f, counts)
        sheet.append(f"t{i},treatment,{f}")
    for i, counts in enumerate(_CTRL_COUNTS):
        f = cov_dir / f"c{i}.cov"
        _write_cov(f, counts)
        sheet.append(f"c{i},control,{f}")
    sheet_path = tmp_path / "samplesheet.csv"
    sheet_path.write_text("\n".join(sheet) + "\n")
    return str(sheet_path)


@pytest.fixture
def scaffold_md(tmp_path):
    """A MethylData whose store has CpGs on chr1 AND a scaffold contig."""
    import epykit as ep

    sheet = _write_cohort(tmp_path)
    md = ep.read_bismark(
        sheet,
        treatment_group="treatment",
        control_group="control",
        assembly="hg38",
        store_dir=str(tmp_path / "store"),
    )
    ep.pp.set_unite_type(md, type="intersect")
    return md


def _dmc_chroms(md) -> set[str]:
    return set(md.dmc["chrom"].unique().to_list())


def test_dmc_default_drops_scaffold(scaffold_md):
    """The new default (canonical_only=True) excludes unplaced contigs."""
    import epykit as ep

    ep.tl.dmc(scaffold_md, test="lr")
    chroms = _dmc_chroms(scaffold_md)
    assert "chr1" in chroms
    assert SCAFFOLD not in chroms


def test_dmc_all_contigs_keeps_scaffold(scaffold_md):
    """canonical_only=False restores the pre-change behaviour."""
    import epykit as ep

    ep.tl.dmc(scaffold_md, test="lr", canonical_only=False)
    chroms = _dmc_chroms(scaffold_md)
    assert "chr1" in chroms
    assert SCAFFOLD in chroms


def test_dmc_explicit_chromosomes_overrides_filter(scaffold_md):
    """An explicit chromosomes= list wins over the canonical filter."""
    import epykit as ep

    # Even with canonical_only at its True default, an explicit request for the
    # scaffold must be honoured verbatim (and chr1 excluded).
    ep.tl.dmc(scaffold_md, test="lr", chromosomes=[SCAFFOLD])
    chroms = _dmc_chroms(scaffold_md)
    assert SCAFFOLD in chroms
    assert "chr1" not in chroms


def test_dmc_audit_log_names_dropped_contig(scaffold_md, caplog):
    """Dropping a contig emits one INFO line naming it (never silent)."""
    import epykit as ep

    with caplog.at_level(logging.INFO, logger="epykit._chroms"):
        ep.tl.dmc(scaffold_md, test="lr")
    msgs = "\n".join(r.getMessage() for r in caplog.records)
    assert "canonical_only" in msgs
    assert SCAFFOLD in msgs


def test_tile_dmr_default_drops_scaffold(scaffold_md):
    """The tile DMR path filters its own chromosome enumeration by default."""
    import epykit as ep

    ep.tl.dmr(scaffold_md, method="tile", tile_size_bp=200, min_cpgs_per_tile=2)
    dmr = scaffold_md.uns["dmr"]
    chroms = set(dmr["chrom"].unique().to_list()) if dmr.height else set()
    assert SCAFFOLD not in chroms


def test_tile_dmr_all_contigs_keeps_scaffold(scaffold_md):
    """method=tile with canonical_only=False still tests scaffolds."""
    import epykit as ep

    ep.tl.dmr(
        scaffold_md, method="tile", tile_size_bp=200,
        min_cpgs_per_tile=2, canonical_only=False,
    )
    dmr = scaffold_md.uns["dmr"]
    chroms = set(dmr["chrom"].unique().to_list()) if dmr.height else set()
    assert SCAFFOLD in chroms


def test_ingestion_canonical_only_omits_scaffold_partition(tmp_path):
    """read_bismark(canonical_only=True) never writes the scaffold partition."""
    import epykit as ep

    sheet = _write_cohort(tmp_path)

    md_default = ep.read_bismark(
        sheet, treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / "store_default"),
    )
    raw_default = Path(md_default.store)
    assert list(raw_default.glob("sample=*/chrom=chr1"))
    assert list(raw_default.glob(f"sample=*/chrom={SCAFFOLD}")), (
        "default ingestion should keep all contigs"
    )

    md_canon = ep.read_bismark(
        sheet, treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / "store_canon"),
        canonical_only=True,
    )
    raw_canon = Path(md_canon.store)
    assert list(raw_canon.glob("sample=*/chrom=chr1"))
    assert not list(raw_canon.glob(f"sample=*/chrom={SCAFFOLD}")), (
        "canonical_only=True must not write the scaffold partition"
    )
