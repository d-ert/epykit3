"""Smoke tests for run_epykit_gse.py (Phase 4 Task 6).

The real DMC run on whole-genome WGBS (~30M CpGs across 6 samples) takes
~1-2 hours, so the test suite does NOT exercise the full pipeline end to
end. Instead it covers the pieces that benefit from CI coverage:

* ``_resolve_samplesheet`` correctly rewrites the relative ``path`` column
  to absolute, raising on missing files.
* The concordance helpers (``_interval_jaccard_dmr``,
  ``_direction_agreement``, ``_overlap_join``) behave correctly on hand-
  built fixtures with known overlap geometry.
* ``build_concordance`` end-to-end on a tiny synthetic DMR set,
  producing the expected output parquets + SUMMARY.md.
* ``main(--skip-dmc --skip-dmr --skip-concordance)`` is wired up and
  exits cleanly (no engine work required).

A real-data ingestion smoke test (``test_real_ingest_skipif``) runs
``ingest()`` only when the raw ``../epykit2/GSE263850_RAW/`` directory is
present, exercises the 12-col BED reader, and asserts the 6-sample,
3-vs-3 obs frame. It is skipped on machines that don't have the raw
data on disk.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import polars as pl
import pytest


_SCRIPTS = Path(__file__).resolve().parents[1]
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))


@pytest.fixture(scope="module")
def runner():
    mod = importlib.import_module("run_epykit_gse")
    return importlib.reload(mod)


# ---------------------------------------------------------------------------
# _resolve_samplesheet
# ---------------------------------------------------------------------------


def test_resolve_samplesheet_rewrites_to_absolute(runner, tmp_path):
    """Relative paths in the input sheet become absolute in the output."""
    f1 = tmp_path / "a.cov.gz"
    f2 = tmp_path / "b.cov.gz"
    f1.write_bytes(b"")
    f2.write_bytes(b"")

    raw = tmp_path / "sheet.csv"
    raw.write_text(
        "sample_id,group,path\n"
        f"s1,treat,{f1.name}\n"      # bare filename -> resolves vs REPO root
        f"s2,ctrl,{f2}\n",            # already absolute -> pass through
        encoding="utf-8",
    )

    # The bare filename resolves against REPO root; copy the files there
    # for the existence check to pass. Use a sibling tmp dir under REPO.
    repo = runner.REPO
    placed_f1 = repo / f1.name
    placed_f1.write_bytes(b"")
    try:
        out_sheet = tmp_path / "resolved.csv"
        result = runner._resolve_samplesheet(raw, out_sheet)
        assert result == out_sheet
        df = pl.read_csv(out_sheet)
        assert df.height == 2
        # Both paths now absolute and existing.
        for p in df["path"].to_list():
            assert Path(p).is_absolute(), p
            assert Path(p).exists(), p
    finally:
        placed_f1.unlink(missing_ok=True)


def test_resolve_samplesheet_raises_on_missing(runner, tmp_path):
    raw = tmp_path / "sheet.csv"
    raw.write_text(
        "sample_id,group,path\n"
        "s1,treat,/definitely/not/a/real/file.cov.gz\n",
        encoding="utf-8",
    )
    with pytest.raises(FileNotFoundError):
        runner._resolve_samplesheet(raw, tmp_path / "resolved.csv")


# ---------------------------------------------------------------------------
# Concordance helpers
# ---------------------------------------------------------------------------


def _mk_dmr_df(rows: list[tuple[str, int, int, float, float]]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={
            "chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
            "qvalue": pl.Float64, "meth_diff": pl.Float64,
        },
        orient="row",
    )


def test_interval_jaccard_dmr_perfect_overlap(runner):
    a = _mk_dmr_df([("chr1", 100, 200, 0.01, 0.2)])
    b = _mk_dmr_df([("chr1", 100, 200, 0.02, 0.3)])
    res = runner._interval_jaccard_dmr(a, b)
    assert res["jaccard"] == 1.0
    assert res["bp_intersection"] == 100
    assert res["bp_union"] == 100


def test_interval_jaccard_dmr_partial_overlap(runner):
    a = _mk_dmr_df([("chr1", 100, 300, 0.01, 0.2)])
    b = _mk_dmr_df([("chr1", 200, 400, 0.02, 0.3)])
    res = runner._interval_jaccard_dmr(a, b)
    # intersection [200, 300) = 100bp, union [100, 400) = 300bp
    assert res["bp_intersection"] == 100
    assert res["bp_union"] == 300
    assert res["jaccard"] == pytest.approx(100 / 300)


def test_interval_jaccard_dmr_disjoint(runner):
    a = _mk_dmr_df([("chr1", 100, 200, 0.01, 0.2)])
    b = _mk_dmr_df([("chr2", 100, 200, 0.02, 0.3)])
    res = runner._interval_jaccard_dmr(a, b)
    assert res["jaccard"] == 0.0
    assert res["bp_intersection"] == 0


def test_direction_agreement(runner):
    """Two overlapping pairs: one sign-matched, one mismatched -> 50%."""
    a = _mk_dmr_df([
        ("chr1", 100, 200, 0.01, 0.4),    # hyper
        ("chr1", 500, 700, 0.01, -0.3),   # hypo
    ])
    b = _mk_dmr_df([
        ("chr1", 150, 250, 0.02, 0.5),    # hyper -- matches a[0]
        ("chr1", 600, 800, 0.02, 0.2),    # hyper -- mismatches a[1]
    ])
    res = runner._direction_agreement(a, b)
    assert res["n_overlapping_pairs"] == 2
    assert res["n_direction_agree"] == 1
    assert res["direction_agree_frac"] == 0.5


def test_overlap_join_schema(runner):
    a = _mk_dmr_df([("chr1", 100, 300, 0.01, 0.4)])
    b = _mk_dmr_df([("chr1", 200, 400, 0.02, 0.3)])
    out = runner._overlap_join(a, b)
    assert out.height == 1
    r = out.row(0, named=True)
    assert r["chrom"] == "chr1"
    assert r["bp_inter"] == 100
    assert r["bp_union"] == 300
    assert r["direction_match"] is True
    # Required columns for downstream consumers.
    for c in ("q_a", "q_b", "meth_diff_a", "meth_diff_b",
              "tool_a"):
        # tool_a is added by build_concordance, not the bare helper
        if c == "tool_a":
            continue
        assert c in out.columns


def test_overlap_join_empty_returns_typed_empty(runner):
    """Empty overlap still gives a valid (zero-row) DataFrame with the
    expected schema -- callers concat across pairs without dtype drift.
    """
    a = _mk_dmr_df([("chr1", 100, 200, 0.01, 0.4)])
    b = _mk_dmr_df([("chr2", 100, 200, 0.02, 0.3)])
    out = runner._overlap_join(a, b)
    assert out.height == 0
    for c in ("a_idx", "b_idx", "chrom", "q_a", "q_b",
              "meth_diff_a", "meth_diff_b", "direction_match"):
        assert c in out.columns


# ---------------------------------------------------------------------------
# build_concordance (synthetic DMR sets)
# ---------------------------------------------------------------------------


def test_build_concordance_smoke(runner, tmp_path, monkeypatch):
    """End-to-end concordance on tiny synthetic DMR sets.

    Writes a fake epykit chain_merge parquet, points the DSS + methylKit
    paths at small synthetic CSVs, and asserts the output parquets +
    SUMMARY.md are produced with sensible row counts.
    """
    # 1) Fake epykit chain_merge.parquet
    ek_dmr = pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [100, 1000, 100],
        "end": [300, 1200, 200],
        # last row filtered out by the q<0.05 threshold
        "qvalue": [0.01, 0.02, 0.5],
        "mean_meth_diff": [0.3, -0.2, 0.4],
    })
    tmp_dmr_dir = tmp_path / "epykit_dmr"
    tmp_dmr_dir.mkdir()
    (tmp_dmr_dir / "chain_merge.parquet").write_bytes(b"")  # placeholder
    ek_dmr.write_parquet(tmp_dmr_dir / "chain_merge.parquet")

    # 2) Fake DSS csv -- diff_Methy_DSSfit column + the same coord shape
    dss_csv = tmp_path / "dmr_dss.csv"
    pl.DataFrame({
        "chrom": ["chr1", "chr2"],
        "start": [150, 100],
        "end": [350, 200],
        "diff_Methy_DSSfit": [0.25, 0.35],
    }).write_csv(dss_csv)

    # 3) Fake methylKit csv -- meth_diff in percent, qvalue column
    mk_csv = tmp_path / "dmr_all_tiles.csv"
    pl.DataFrame({
        "chrom": ["chr1", "chr3"],
        "start": [200, 100],
        "end": [400, 200],
        "qvalue": [0.001, 0.5],   # second row filtered out (q > 0.05)
        "meth_diff": [25.0, 10.0],
    }).write_csv(mk_csv)

    # Redirect runner paths
    out_dir = tmp_path / "comparisons_post_phase3"
    monkeypatch.setattr(runner, "OUT_DMR", tmp_dmr_dir)
    monkeypatch.setattr(runner, "DSS_DMR_CSV", dss_csv)
    monkeypatch.setattr(runner, "MK_DMR_CSV", mk_csv)

    headline = runner.build_concordance(out_dir, engines_done=("lr",))

    # Outputs exist
    assert (out_dir / "dmr_iou.parquet").exists()
    assert (out_dir / "per_dmr_stat_concordance.parquet").exists()
    assert (out_dir / "SUMMARY.md").exists()

    iou = pl.read_parquet(out_dir / "dmr_iou.parquet")
    # 3 tools => C(3,2)=3 pairs
    assert iou.height == 3
    pairs = set(zip(iou["tool_a"].to_list(), iou["tool_b"].to_list()))
    assert pairs == {
        ("dss", "epykit_chain_merge"),
        ("dss", "methylkit"),
        ("epykit_chain_merge", "methylkit"),
    }
    # All Jaccards in [0, 1]
    for j in iou["jaccard"].to_list():
        assert 0.0 <= j <= 1.0

    # headline carries set_sizes
    assert "set_sizes" in headline
    assert headline["set_sizes"]["epykit_chain_merge"] == 2  # both q<0.05
    assert headline["set_sizes"]["methylkit"] == 1           # one filtered out

    # SUMMARY.md is non-empty + mentions the call sizes
    text = (out_dir / "SUMMARY.md").read_text(encoding="utf-8")
    assert "Phase-3 cross-tool DMR concordance" in text
    assert "bp-Jaccard" in text


# ---------------------------------------------------------------------------
# main() driver: --skip-* paths
# ---------------------------------------------------------------------------


def test_main_all_skipped(runner, tmp_path, monkeypatch):
    """``--skip-dmc --skip-dmr --skip-concordance`` writes only the manifest."""
    monkeypatch.setattr(runner, "OUT_BASE", tmp_path / "epykit_post_phase3")
    monkeypatch.setattr(runner, "OUT_DMC", tmp_path / "epykit_post_phase3" / "dmc")
    monkeypatch.setattr(runner, "OUT_DMR", tmp_path / "epykit_post_phase3" / "dmr")
    monkeypatch.setattr(
        runner, "MANIFEST_PATH",
        tmp_path / "epykit_post_phase3" / "MANIFEST.txt",
    )
    monkeypatch.setattr(runner, "CMP_OUT", tmp_path / "cmp")

    rc = runner.main([
        "--skip-dmc", "--skip-dmr", "--skip-concordance",
    ])
    assert rc == 0
    assert (tmp_path / "epykit_post_phase3" / "MANIFEST.txt").exists()


# ---------------------------------------------------------------------------
# Real-data ingestion smoke (skipped when raw data absent)
# ---------------------------------------------------------------------------


def _raw_data_available() -> bool:
    expected = Path(
        "../epykit2/GSE263850_RAW/"
        "GSM8200106_Clone16_untreated.readset_sorted.dedup.filtered.bed.gz"
    ).resolve()
    return expected.exists()


@pytest.mark.slow
@pytest.mark.skipif(
    not _raw_data_available(),
    reason="GSE263850 raw .bed.gz not present at ../epykit2/GSE263850_RAW/",
)
def test_real_ingest_obs_shape(runner, tmp_path):
    """When raw data is present, ingest() returns 6 samples in 3v3 design.

    This is a real (slow) smoke test: it does the full Parquet
    partitioning of all six samples, but does NOT run any DMC engines.
    Wallclock is dominated by the parquet conversion (a few minutes).
    """
    md = runner.ingest(store_dir=tmp_path / "store")
    assert md.obs.height == 6
    groups = sorted(md.obs["group"].to_list())
    assert groups.count("sbp009") == 3
    assert groups.count("clone") == 3
    assert sorted(md.obs["treatment"].to_list()) == [0, 0, 0, 1, 1, 1]
