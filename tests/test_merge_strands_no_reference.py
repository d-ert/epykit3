"""Regression tests for M5: merge_strands=True must merge CpG dyads without a
reference FASTA (strand-free / position-based merging).

Pre-fix: merge_strands=True was a silent no-op when reference_fasta was not
given, so each CpG dyad survived as two rows at half coverage.  These tests
verify the corrected behaviour end-to-end through read_bismark / convert_sample
and also exercise the _merge_cpg_pairs_by_position helper directly.

Coordinate note: Bismark .cov is 1-based (start == end).  We write start = end
= (0-based pos + 1) so auto-detection maps back to the correct 0-based pos.

Dyad pairing rule (strand-free):
  + strand cytosine is at 0-based pos P  -> Bismark start = P+1
  - strand cytosine is at 0-based pos P+1 -> Bismark start = P+2
  After coordinate conversion pos(+) = P, pos(-) = P+1.
  The merger pairs P with P+1 and emits a single site at pos P.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path

import polars as pl

import epykit as ep
from epykit.convert import _merge_cpg_pairs_by_position, convert_sample

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_cov(path: Path, rows) -> None:
    """Write a plain (uncompressed) Bismark .cov file.

    rows: iterable of (chrom, start, end, pct, N_meth, N_unmeth)
    start == end for 1-based Bismark .cov.
    """
    with open(path, "w", newline="") as fh:
        for chrom, start, end, pct, m, u in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{pct:.2f}\t{m}\t{u}\n")


def _write_sheet(path: Path, sample_paths: dict) -> None:
    """Write a minimal samplesheet CSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "group", "path"])
        writer.writeheader()
        for sid, (grp, p) in sample_paths.items():
            writer.writerow({"sample_id": sid, "group": grp, "path": str(p)})


def _read_store(md, sample: str, chrom: str) -> pl.DataFrame:
    part = (
        Path(md.store)
        / f"sample={sample}"
        / f"chrom={chrom}"
        / "part-0.parquet"
    )
    return pl.read_parquet(part).sort("pos")


# ---------------------------------------------------------------------------
# Unit tests for the helper directly
# ---------------------------------------------------------------------------

class TestMergeCpgPairsByPosition:
    """Direct unit tests for the strand-free merge helper."""

    def _make_df(self, rows: list[tuple]) -> pl.DataFrame:
        """rows: (chrom, pos, N_meth, N_unmeth, coverage, strand, context, sample)"""
        if not rows:
            return pl.DataFrame({
                "chrom":    pl.Series([], dtype=pl.Utf8),
                "pos":      pl.Series([], dtype=pl.Int32),
                "strand":   pl.Series([], dtype=pl.Utf8),
                "context":  pl.Series([], dtype=pl.Utf8),
                "N_meth":   pl.Series([], dtype=pl.Int32),
                "N_unmeth": pl.Series([], dtype=pl.Int32),
                "coverage": pl.Series([], dtype=pl.Int32),
                "sample":   pl.Series([], dtype=pl.Utf8),
            })
        chroms, poses, nm, nu, cov, strand, ctx, samp = zip(*rows, strict=False)
        return pl.DataFrame({
            "chrom":    list(chroms),
            "pos":      pl.Series(list(poses), dtype=pl.Int32),
            "strand":   list(strand),
            "context":  list(ctx),
            "N_meth":   pl.Series(list(nm), dtype=pl.Int32),
            "N_unmeth": pl.Series(list(nu), dtype=pl.Int32),
            "coverage": pl.Series(list(cov), dtype=pl.Int32),
            "sample":   list(samp),
        })

    def test_basic_dyad_pair_merged(self):
        """Two rows for the same CpG dyad (pos P and P+1) collapse to one."""
        df = self._make_df([
            ("chr1", 100, 8, 2, 10, "*", "CpG", "s1"),  # + strand C
            ("chr1", 101, 6, 4, 10, "*", "CpG", "s1"),  # - strand C
        ])
        out = _merge_cpg_pairs_by_position(df)
        assert out.height == 1
        row = out.row(0, named=True)
        assert row["pos"] == 100
        assert row["N_meth"] == 14
        assert row["N_unmeth"] == 6
        assert row["coverage"] == 20
        assert row["strand"] == "*"

    def test_multiple_clean_dyads(self):
        """Three dyads far apart all collapse to 3 sites."""
        df = self._make_df([
            ("chr1", 200, 5, 5, 10, "*", "CpG", "s1"),
            ("chr1", 201, 3, 7, 10, "*", "CpG", "s1"),
            ("chr1", 500, 9, 1, 10, "*", "CpG", "s1"),
            ("chr1", 501, 7, 3, 10, "*", "CpG", "s1"),
            ("chr1", 800, 4, 6, 10, "*", "CpG", "s1"),
            ("chr1", 801, 2, 8, 10, "*", "CpG", "s1"),
        ])
        out = _merge_cpg_pairs_by_position(df).sort("pos")
        assert out.height == 3
        assert out["pos"].to_list() == [200, 500, 800]
        # Each merged dyad has coverage 20
        assert out["coverage"].to_list() == [20, 20, 20]

    def test_cgcg_adjacency_pairs_correctly(self):
        """CGCG sequence: positions {P, P+1, P+2, P+3}.

        Greedy left-to-right must pair (P, P+1) and (P+2, P+3),
        not collapse all four into one.
        """
        P = 300
        df = self._make_df([
            ("chr1", P,     5, 5, 10, "*", "CpG", "s1"),
            ("chr1", P + 1, 5, 5, 10, "*", "CpG", "s1"),
            ("chr1", P + 2, 3, 7, 10, "*", "CpG", "s1"),
            ("chr1", P + 3, 3, 7, 10, "*", "CpG", "s1"),
        ])
        out = _merge_cpg_pairs_by_position(df).sort("pos")
        assert out.height == 2
        assert out["pos"].to_list() == [P, P + 2]
        assert out["coverage"].to_list() == [20, 20]

    def test_singleton_passthrough(self):
        """A site with no ±1 neighbour passes through unmodified."""
        df = self._make_df([
            ("chr1", 400, 7, 3, 10, "*", "CpG", "s1"),
            ("chr1", 600, 4, 6, 10, "*", "CpG", "s1"),   # gap > 1 from 400
        ])
        out = _merge_cpg_pairs_by_position(df).sort("pos")
        assert out.height == 2
        assert out["pos"].to_list() == [400, 600]

    def test_multi_chrom_no_cross_pairing(self):
        """Positions on different chromosomes must not be paired together."""
        df = self._make_df([
            # chr1 pos 100 has no +1 neighbour on chr1
            ("chr1", 100, 5, 5, 10, "*", "CpG", "s1"),
            # chr2 pos 101 must NOT pair with chr1 pos 100
            ("chr2", 101, 3, 7, 10, "*", "CpG", "s1"),
        ])
        out = _merge_cpg_pairs_by_position(df).sort("chrom", "pos")
        assert out.height == 2

    def test_canonical_output_schema(self):
        """Output must contain the canonical set of columns."""
        CANONICAL = {"chrom", "pos", "strand", "context",
                     "N_meth", "N_unmeth", "coverage", "sample"}
        df = self._make_df([
            ("chr1", 10, 8, 2, 10, "*", "CpG", "s1"),
            ("chr1", 11, 6, 4, 10, "*", "CpG", "s1"),
        ])
        out = _merge_cpg_pairs_by_position(df)
        assert set(out.columns) == CANONICAL

    def test_empty_frame_roundtrip(self):
        """Empty input returns empty output without error."""
        df = self._make_df([])
        out = _merge_cpg_pairs_by_position(df)
        assert out.height == 0

    def test_warning_fires_on_real_merge(self, caplog):
        """A WARNING mentioning reference_fasta is emitted exactly once when dyads are merged."""
        df = self._make_df([
            ("chr1", 100, 8, 2, 10, "*", "CpG", "s1"),  # dyad pair-start
            ("chr1", 101, 6, 4, 10, "*", "CpG", "s1"),  # dyad partner
        ])
        with caplog.at_level(logging.WARNING, logger="epykit.convert"):
            _merge_cpg_pairs_by_position(df)
        warning_records = [
            r for r in caplog.records
            if r.levelno == logging.WARNING and "reference_fasta" in r.message
        ]
        assert len(warning_records) == 1, (
            f"Expected exactly 1 reference_fasta WARNING, got {len(warning_records)}: "
            f"{[r.message for r in warning_records]}"
        )

    def test_no_warning_when_nothing_merged(self, caplog):
        """No WARNING is emitted when all sites are singletons (no dyads merged)."""
        df = self._make_df([
            ("chr1", 100, 7, 3, 10, "*", "CpG", "s1"),  # isolated — no ±1 neighbour
            ("chr1", 500, 4, 6, 10, "*", "CpG", "s1"),  # isolated — no ±1 neighbour
        ])
        with caplog.at_level(logging.WARNING, logger="epykit.convert"):
            _merge_cpg_pairs_by_position(df)
        warning_records = [
            r for r in caplog.records
            if r.levelno == logging.WARNING and "reference_fasta" in r.message
        ]
        assert len(warning_records) == 0, (
            f"Expected no WARNING when nothing is merged, got: "
            f"{[r.message for r in warning_records]}"
        )


# ---------------------------------------------------------------------------
# End-to-end tests via read_bismark / convert_sample
# ---------------------------------------------------------------------------

class TestMergeStrandsNoReferenceBismark:
    """End-to-end tests exercising read_bismark with merge_strands=True
    and no reference_fasta."""

    def test_basic_dyad_halves_site_count(self, tmp_path):
        """Two rows per dyad in the .cov become one merged site in the store."""
        cov = tmp_path / "s1.cov"
        # Two CpG dyads at 1-based (101,102) and (201,202).
        # 0-based store pos: (100,101) -> merged at 100; (200,201) -> merged at 200.
        _write_cov(cov, [
            ("chr1", 101, 101, 80.0, 8, 2),   # + strand, dyad 1
            ("chr1", 102, 102, 60.0, 6, 4),   # - strand, dyad 1
            ("chr1", 201, 201, 70.0, 7, 3),   # + strand, dyad 2
            ("chr1", 202, 202, 50.0, 5, 5),   # - strand, dyad 2
        ])
        sheet = tmp_path / "sheet.csv"
        _write_sheet(sheet, {"s1": ("treatment", cov), "s2": ("control", cov)})

        md = ep.read_bismark(
            str(sheet),
            treatment_group="treatment",
            control_group="control",
            store_dir=str(tmp_path / "store"),
        )
        df = _read_store(md, "s1", "chr1")

        # Should be exactly 2 sites (one per dyad), not 4
        assert df.height == 2, (
            f"Expected 2 merged dyad sites, got {df.height}. "
            "merge_strands=True is likely still a no-op without reference_fasta."
        )

    def test_merged_coverage_is_summed(self, tmp_path):
        """Merged site coverage equals the sum of both strand rows."""
        cov = tmp_path / "s1.cov"
        # Dyad at 1-based (301, 302): + strand 8 meth / 2 unmeth = 10 cov
        #                              - strand 6 meth / 4 unmeth = 10 cov
        # Expected merged: N_meth=14, N_unmeth=6, coverage=20
        _write_cov(cov, [
            ("chr1", 301, 301, 80.0, 8, 2),
            ("chr1", 302, 302, 60.0, 6, 4),
        ])
        sheet = tmp_path / "sheet.csv"
        _write_sheet(sheet, {"s1": ("treatment", cov), "s2": ("control", cov)})

        md = ep.read_bismark(
            str(sheet),
            treatment_group="treatment",
            control_group="control",
            store_dir=str(tmp_path / "store"),
        )
        df = _read_store(md, "s1", "chr1")

        assert df.height == 1
        row = df.row(0, named=True)
        assert row["pos"] == 300, f"Expected merged pos 300, got {row['pos']}"
        assert row["N_meth"] == 14
        assert row["N_unmeth"] == 6
        assert row["coverage"] == 20

    def test_cgcg_adjacency_end_to_end(self, tmp_path):
        """CGCG sequence (4 consecutive rows) must produce exactly 2 merged sites."""
        P = 500  # 0-based start of first dyad
        cov = tmp_path / "s1.cov"
        # 1-based: P+1=501, P+2=502, P+3=503, P+4=504
        _write_cov(cov, [
            ("chr1", P + 1, P + 1, 80.0, 8, 2),   # dyad1 + strand
            ("chr1", P + 2, P + 2, 80.0, 8, 2),   # dyad1 - strand
            ("chr1", P + 3, P + 3, 60.0, 6, 4),   # dyad2 + strand
            ("chr1", P + 4, P + 4, 60.0, 6, 4),   # dyad2 - strand
        ])
        sheet = tmp_path / "sheet.csv"
        _write_sheet(sheet, {"s1": ("treatment", cov), "s2": ("control", cov)})

        md = ep.read_bismark(
            str(sheet),
            treatment_group="treatment",
            control_group="control",
            store_dir=str(tmp_path / "store"),
        )
        df = _read_store(md, "s1", "chr1")

        assert df.height == 2, (
            f"CGCG adjacency: expected 2 merged sites, got {df.height}"
        )
        assert df["pos"].to_list() == [P, P + 2]

    def test_singleton_survives(self, tmp_path):
        """A site with no ±1 neighbour passes through unmodified."""
        cov = tmp_path / "s1.cov"
        # Two sites far apart (0-based 100 and 300); neither has a ±1 neighbour
        _write_cov(cov, [
            ("chr1", 101, 101, 80.0, 8, 2),   # 0-based 100
            ("chr1", 301, 301, 70.0, 7, 3),   # 0-based 300
        ])
        sheet = tmp_path / "sheet.csv"
        _write_sheet(sheet, {"s1": ("treatment", cov), "s2": ("control", cov)})

        md = ep.read_bismark(
            str(sheet),
            treatment_group="treatment",
            control_group="control",
            store_dir=str(tmp_path / "store"),
        )
        df = _read_store(md, "s1", "chr1")
        assert df.height == 2
        assert df["pos"].to_list() == [100, 300]

    def test_merged_strand_is_star(self, tmp_path):
        """Merged sites must carry strand='*' (strand genuinely unknown)."""
        cov = tmp_path / "s1.cov"
        _write_cov(cov, [
            ("chr1", 101, 101, 80.0, 8, 2),
            ("chr1", 102, 102, 60.0, 6, 4),
        ])
        sheet = tmp_path / "sheet.csv"
        _write_sheet(sheet, {"s1": ("treatment", cov), "s2": ("control", cov)})

        md = ep.read_bismark(
            str(sheet),
            treatment_group="treatment",
            control_group="control",
            store_dir=str(tmp_path / "store"),
        )
        df = _read_store(md, "s1", "chr1")
        assert df.height == 1
        assert df["strand"][0] == "*", (
            f"Merged strand should be '*' (unknown), got {df['strand'][0]!r}"
        )


class TestMergeStrandsDirectConvertSample:
    """Direct convert_sample tests so we can also check merge_strands=False."""

    def test_merge_strands_false_leaves_two_rows(self, tmp_path):
        """With merge_strands=False the two dyad rows are preserved as-is."""
        cov = tmp_path / "s1.cov"
        _write_cov(cov, [
            ("chr1", 101, 101, 80.0, 8, 2),
            ("chr1", 102, 102, 60.0, 6, 4),
        ])
        out = tmp_path / "store"
        convert_sample(str(cov), "s1", str(out), merge_strands=False)
        part = out / "sample=s1" / "chrom=chr1" / "part-0.parquet"
        df = pl.read_parquet(part).sort("pos")
        assert df.height == 2

    def test_merge_strands_true_no_ref_merges(self, tmp_path):
        """convert_sample with merge_strands=True and no reference merges."""
        cov = tmp_path / "s1.cov"
        _write_cov(cov, [
            ("chr1", 101, 101, 80.0, 8, 2),
            ("chr1", 102, 102, 60.0, 6, 4),
        ])
        out = tmp_path / "store"
        convert_sample(str(cov), "s1", str(out), merge_strands=True)
        part = out / "sample=s1" / "chrom=chr1" / "part-0.parquet"
        df = pl.read_parquet(part).sort("pos")
        assert df.height == 1
        assert df["N_meth"][0] == 14
        assert df["coverage"][0] == 20
