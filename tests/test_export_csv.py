"""Tests for CSV/TSV export of epykit result tables."""
from __future__ import annotations

import csv
from pathlib import Path

import polars as pl
import pytest
import epykit as ep
from epykit import export


def _make_dmr_frame() -> pl.DataFrame:
    return pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [200, 100, 50],
        "end":   [300, 200, 150],
        "meth_diff": [0.3, -0.4, 0.2],
        "qvalue":  [0.01, 0.02, 0.5],
        "dmr_type": ["hyper", "hypo", "hyper"],
    })


def _stub_md_with_dmr(tmp_path: Path):
    """Build a minimal MethylData carrying a DMR table in md.uns['dmr']."""
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.uns["dmr"] = _make_dmr_frame()
    return md


def test_dmr_to_tsv_writes_full_table_chrom_start_sorted(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.tsv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    assert len(rows) == 3
    assert [(r["chrom"], int(r["start"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr2", 50),
    ]


def test_dmr_to_csv_uses_comma_for_csv_suffix(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.csv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    # Header line must contain commas, no tabs.
    header = text.splitlines()[0]
    assert "," in header
    assert "\t" not in header


def _make_dmc_frame(with_combined: bool = False) -> pl.DataFrame:
    rows = {
        "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "pos":   [100, 200, 300, 50, 150],
        "meth_diff": [0.3, -0.4, 0.05, 0.2, -0.1],
        "pvalue": [1e-6, 1e-5, 0.5, 1e-3, 0.4],
        "qvalue": [1e-5, 1e-4, 0.6, 1e-2, 0.5],
    }
    if with_combined:
        # Flip combined values so significance differs from raw qvalue.
        rows["pvalue_combined"] = [0.5, 1e-7, 0.5, 0.5, 1e-9]
        rows["qvalue_combined"] = [0.6, 1e-6, 0.6, 0.6, 1e-8]
    return pl.DataFrame(rows)


def _stub_md_with_dmc(tmp_path: Path, *, key: str = "dmc_lr",
                       with_combined: bool = False):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"], "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm[key] = _make_dmc_frame(with_combined=with_combined)
    md.uns["dmc"] = {"last_key": key}
    return md


def test_dmc_to_tsv_significant_only_qvalue_asc(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.significant.tsv"
    export.dmc_to_tsv(md, str(out))  # default alpha=0.05, full=False

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # 3 rows have qvalue < 0.05 (1e-5, 1e-4, 1e-2); the other 2 (q=0.5, 0.6) are dropped.
    assert len(rows) == 3
    # qvalue ascending
    qvalues = [float(r["qvalue"]) for r in rows]
    assert qvalues == sorted(qvalues)


def test_dmc_to_tsv_full_writes_all_rows_genomic_order(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.tsv"
    export.dmc_to_tsv(md, str(out), full=True)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 5
    assert [(r["chrom"], int(r["pos"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr1", 300),
        ("chr2", 50),  ("chr2", 150),
    ]


def test_dmc_to_tsv_alpha_override(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.strict.tsv"
    export.dmc_to_tsv(md, str(out), alpha=1e-3)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # Only rows with qvalue < 1e-3 -> qvalue 1e-5 and 1e-4
    assert len(rows) == 2
    assert all(float(r["qvalue"]) < 1e-3 for r in rows)


def test_dmc_to_tsv_uses_qvalue_combined_when_present(tmp_path):
    md = _stub_md_with_dmc(tmp_path, with_combined=True)
    out = tmp_path / "dmc.combined.tsv"
    export.dmc_to_tsv(md, str(out))

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # qvalue_combined < 0.05 for rows with combined values 1e-6 and 1e-8.
    assert len(rows) == 2
    qc = sorted(float(r["qvalue_combined"]) for r in rows)
    assert qc == [1e-8, 1e-6]


def test_dmc_to_tsv_flattens_nested_list_columns(tmp_path):
    """Annotated DMC frames carry List(String) columns (overlapping genes /
    features). polars' CSV writer rejects nested data, so the writer must
    flatten them to '; '-joined strings instead of raising ComputeError."""
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({"sample_id": ["s1", "s2"], "group": ["case", "ctrl"]})
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm["dmc_lr"] = pl.DataFrame({
        "chrom": ["chr1", "chr2"],
        "pos": [100, 200],
        "meth_diff": [0.3, -0.4],
        "qvalue": [1e-5, 1e-4],
        "all_overlapping_genes": [["GENEA", "GENEB"], ["GENEC"]],
        "all_overlapping_features": [["intron"], ["exon", "utr"]],
    })
    md.uns["dmc"] = {"last_key": "dmc_lr"}

    out = tmp_path / "dmc.annotated.tsv"
    export.dmc_to_tsv(md, str(out), full=True)  # must not raise

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    genes = {r["all_overlapping_genes"] for r in rows}
    assert genes == {"GENEA; GENEB", "GENEC"}
    feats = {r["all_overlapping_features"] for r in rows}
    assert feats == {"intron", "exon; utr"}


def _stub_md_with_dvc(tmp_path: Path):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"], "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm["dvc"] = pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "pos":   [100, 200, 50],
        "var_log_ratio": [1.2, 0.1, 0.9],
        "p_variance": [1e-5, 0.6, 1e-3],
        "q_variance": [1e-4, 0.7, 1e-2],
        "is_dvc": [True, False, True],
    })
    return md


def test_dvc_to_tsv_significant_only(tmp_path):
    md = _stub_md_with_dvc(tmp_path)
    out = tmp_path / "dvc.tsv"
    export.dvc_to_tsv(md, str(out))

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 2
    qs = [float(r["q_variance"]) for r in rows]
    assert qs == sorted(qs)


def test_dvc_to_tsv_full(tmp_path):
    md = _stub_md_with_dvc(tmp_path)
    out = tmp_path / "dvc.full.tsv"
    export.dvc_to_tsv(md, str(out), full=True)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 3


def test_qc_to_tsv_writes_md_obs(tmp_path):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
        "mean_coverage": [12.3, 8.1],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))

    out = tmp_path / "qc.tsv"
    export.qc_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert set(header) == {"sample_id", "group", "mean_coverage"}
    assert "s1\t" in text and "s2\t" in text


def _stub_md_with_all_tables(tmp_path: Path):
    """A MethylData carrying DMC + DMR + DVC tables and a populated obs."""
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
        "mean_coverage": [12.3, 8.1],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm["dmc_lr"] = _make_dmc_frame()
    md.uns["dmc"] = {"last_key": "dmc_lr"}
    md.uns["dmr"] = _make_dmr_frame()
    md.varm["dvc"] = pl.DataFrame({
        "chrom": ["chr1", "chr1"],
        "pos": [100, 200],
        "p_variance": [1e-5, 0.6],
        "q_variance": [1e-4, 0.7],
        "is_dvc": [True, False],
    })
    return md


def test_export_tables_writes_all_present_tables(tmp_path):
    md = _stub_md_with_all_tables(tmp_path)
    out = tmp_path / "tables"
    written = md.export_tables(str(out))

    assert set(written) == {"dmc_significant", "dmr", "dvc_significant", "qc"}
    assert (out / "dmc_lr.significant.tsv").exists()
    assert (out / "dmr.tsv").exists()
    assert (out / "dvc.significant.tsv").exists()
    assert (out / "qc_summary.tsv").exists()
    # Returned paths are the real files on disk.
    assert all(Path(p).exists() for p in written.values())


def test_export_tables_full_adds_full_dmc_and_dvc(tmp_path):
    md = _stub_md_with_all_tables(tmp_path)
    out = tmp_path / "tables"
    written = md.export_tables(str(out), full=True)

    assert "dmc_full" in written and "dvc_full" in written
    assert (out / "dmc_lr.tsv").exists()
    assert (out / "dvc.tsv").exists()
    # Full DMC keeps every row; significant-only drops the q>=0.05 ones.
    full_rows = len(Path(out / "dmc_lr.tsv").read_text(encoding="utf-8").splitlines()) - 1
    sig_rows = len(Path(out / "dmc_lr.significant.tsv").read_text(encoding="utf-8").splitlines()) - 1
    assert full_rows == 5 and sig_rows == 3


def test_export_tables_skips_missing_tables(tmp_path):
    """Only a DMC table present -> only DMC (+ qc from obs) is written, no raise."""
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "tables"
    written = md.export_tables(str(out))

    assert "dmc_significant" in written
    assert "dmr" not in written and "dvc_significant" not in written
    assert not (out / "dmr.tsv").exists()


def test_export_tables_csv_fmt_uses_comma_and_csv_suffix(tmp_path):
    md = _stub_md_with_all_tables(tmp_path)
    out = tmp_path / "tables"
    written = md.export_tables(str(out), fmt="csv")

    assert (out / "dmr.csv").exists()
    header = Path(written["dmr"]).read_text(encoding="utf-8").splitlines()[0]
    assert "," in header and "\t" not in header


def test_tl_dmc_tsv_kwarg_writes_file(tmp_path, synth_md_filtered):
    out = tmp_path / "dmc.significant.tsv"
    ep.tl.dmc(synth_md_filtered, test="lr", tsv=str(out))

    assert out.exists()
    text = out.read_text(encoding="utf-8")
    # Header from md.varm['dmc_lr']: chrom, pos, ... qvalue ...
    header = text.splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header and "qvalue" in header


def test_tl_dmc_tsv_full_writes_every_row(tmp_path, synth_md_filtered):
    full_out = tmp_path / "dmc.tsv"
    ep.tl.dmc(synth_md_filtered, test="lr", tsv=str(full_out), tsv_full=True)

    n_full = len(full_out.read_text(encoding="utf-8").splitlines()) - 1
    n_varm = len(synth_md_filtered.varm["dmc_lr"])
    assert n_full == n_varm


def test_tl_dvc_tsv_kwarg_writes_file(tmp_path, synth_md_filtered):
    ep.tl.dvc(synth_md_filtered, test="bartlett")
    out = tmp_path / "dvc.significant.tsv"
    ep.tl.dvc(synth_md_filtered, test="bartlett", tsv=str(out))
    # File should exist; may be empty (no significant DVCs in the fixture)
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header


def test_tl_dmr_tsv_kwarg_writes_file(tmp_path, synth_md_filtered):
    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "dmr.tsv"
    ep.tl.dmr(synth_md_filtered, method="chain_merge", tsv=str(out))

    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "start" in header and "end" in header


def test_tl_qc_tsv_kwarg_writes_md_obs(tmp_path, synth_md_filtered):
    out = tmp_path / "qc.tsv"
    ep.tl.qc(synth_md_filtered, tsv=str(out))

    assert out.exists()
    text = out.read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert "sample_id" in header


def test_tl_dmc_csv_kwarg_is_deprecated_alias(tmp_path, synth_md_filtered):
    """The legacy `csv=` kwarg still works but emits a DeprecationWarning."""
    out = tmp_path / "dmc.significant.tsv"
    with pytest.warns(DeprecationWarning, match="deprecated"):
        ep.tl.dmc(synth_md_filtered, test="lr", csv=str(out))
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "qvalue" in header


def test_tl_dmr_csv_kwarg_is_deprecated_alias(tmp_path, synth_md_filtered):
    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "dmr.tsv"
    with pytest.warns(DeprecationWarning, match="deprecated"):
        ep.tl.dmr(synth_md_filtered, method="chain_merge", csv=str(out))
    assert out.exists()


def test_tl_dmc_tsv_wins_when_both_given(tmp_path, synth_md_filtered):
    """If both tsv= and csv= are passed, tsv= takes precedence (still warns)."""
    tsv_out = tmp_path / "from_tsv.tsv"
    csv_out = tmp_path / "from_csv.tsv"
    with pytest.warns(DeprecationWarning):
        ep.tl.dmc(synth_md_filtered, test="lr", tsv=str(tsv_out), csv=str(csv_out))
    assert tsv_out.exists()
    assert not csv_out.exists()


def test_tl_dmc_auto_emits_to_analysis_root_by_default(synth_md_filtered):
    """Default-on: bare ep.tl.dmc(md) writes <analysis_root>/results/dmc.significant.tsv."""
    auto = Path(synth_md_filtered.analysis_root) / "results" / "dmc.significant.tsv"
    if auto.exists():
        auto.unlink()
    ep.tl.dmc(synth_md_filtered, test="lr")  # no tsv= -> auto-emit
    assert auto.exists(), "default-on auto-emit should write the significant DMC TSV"


def test_tl_dmc_tsv_false_disables_auto_emit(synth_md_filtered):
    auto = Path(synth_md_filtered.analysis_root) / "results" / "dmc.significant.tsv"
    if auto.exists():
        auto.unlink()
    ep.tl.dmc(synth_md_filtered, test="lr", tsv=False)
    assert not auto.exists()


def test_tl_dmr_auto_emits_to_analysis_root_by_default(synth_md_filtered):
    ep.tl.dmc(synth_md_filtered, test="lr", tsv=False)
    auto = Path(synth_md_filtered.analysis_root) / "results" / "dmr.tsv"
    if auto.exists():
        auto.unlink()
    ep.tl.dmr(synth_md_filtered, method="chain_merge")  # no tsv= -> auto-emit
    assert auto.exists()


def test_cli_dmc_auto_emits_sibling_significant_tsv(tmp_path, synth_bundle, monkeypatch):
    """`epykit dmc --output X.parquet` writes X.significant.tsv next to it."""
    import sys
    import epykit as ep
    from epykit.cli import main

    # Populate the methylstore so the CLI dmc command has data to process.
    # read_bismark returns a MethylData whose .store points at the raw parquet
    # partition tree (sample=.../chrom=.../part-0.parquet) — that is the path
    # process_chromosomes_dmc expects as --methylstore.
    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store  # e.g. <store_dir>/.cache/raw/

    out_parquet = tmp_path / "dmc.parquet"
    sibling_tsv = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
    ])
    main()

    assert out_parquet.exists()
    assert sibling_tsv.exists()
    # Should contain a header line and some data rows
    lines = sibling_tsv.read_text(encoding="utf-8").splitlines()
    assert len(lines) >= 1
    assert "chrom" in lines[0]


def test_cli_dmc_no_tsv_suppresses_sibling(tmp_path, synth_bundle, monkeypatch):
    import sys
    import epykit as ep
    from epykit.cli import main

    # Populate the store the same way as the auto-emit test.
    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    out_parquet = tmp_path / "dmc.parquet"
    sibling_tsv = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
        "--no-tsv",
    ])
    main()

    assert out_parquet.exists()
    assert not sibling_tsv.exists()


def test_cli_dmc_no_csv_still_suppresses_deprecated(tmp_path, synth_bundle, monkeypatch):
    """The deprecated --no-csv flag is still honoured (alias for --no-tsv)."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    out_parquet = tmp_path / "dmc.parquet"
    sibling_tsv = tmp_path / "dmc.significant.tsv"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", md_setup.store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
        "--no-csv",
    ])
    main()

    assert out_parquet.exists()
    assert not sibling_tsv.exists()


def test_cli_dmr_auto_emits_sibling_tsv(tmp_path, synth_bundle, monkeypatch):
    """`epykit dmr --output X.parquet` writes X.tsv next to it."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    # First make a DMC parquet, since chain_merge consumes one.
    dmc_parquet = tmp_path / "dmc.parquet"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(dmc_parquet),
        "--test", "lr",
        "--no-tsv",  # don't pollute tmp_path with the dmc sibling
    ])
    main()

    dmr_parquet = tmp_path / "dmr.parquet"
    sibling = tmp_path / "dmr.tsv"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmr",
        "--method", "chain_merge",
        "--dmc-results", str(dmc_parquet),
        "--output", str(dmr_parquet),
        "--preset", "permissive",
    ])
    main()

    assert dmr_parquet.exists()
    assert sibling.exists()
    header = sibling.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "start" in header


def test_cli_annotate_auto_emits_sibling_tsv(tmp_path, synth_bundle, monkeypatch):
    """`epykit annotate --output X.parquet` writes X.tsv next to it."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    dmc_parquet = tmp_path / "dmc.parquet"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(dmc_parquet),
        "--test", "lr",
        "--no-csv",
    ])
    main()

    annotated = tmp_path / "annotated.parquet"
    sibling = tmp_path / "annotated.tsv"
    # No --gtf / --cpg-islands -> annotate is a pass-through, but the sibling
    # TSV must still be written.
    monkeypatch.setattr(sys, "argv", [
        "epykit", "annotate",
        "--input", str(dmc_parquet),
        "--output", str(annotated),
    ])
    main()

    assert annotated.exists()
    assert sibling.exists()
    header = sibling.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header


def test_cli_qc_report_auto_emits_sibling_tsvs(tmp_path, synth_bundle, monkeypatch):
    """`epykit qc-report --output-dir DIR` writes sibling .tsv files in DIR."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    out_dir = tmp_path / "qc_out"
    samples = ",".join(synth_bundle.treatment_ids + synth_bundle.control_ids)

    monkeypatch.setattr(sys, "argv", [
        "epykit", "qc-report",
        "--methylstore", populated_store,
        "--samples", samples,
        "--output-dir", str(out_dir),
    ])
    main()

    assert (out_dir / "global_methylation.parquet").exists()
    assert (out_dir / "global_methylation.tsv").exists()
    assert (out_dir / "coverage_uniformity.parquet").exists()
    assert (out_dir / "coverage_uniformity.tsv").exists()


def test_env_var_suppresses_cli_auto_emit(tmp_path, synth_bundle, monkeypatch):
    """EPYKIT_NO_AUTO_CSV=1 suppresses the sibling write across the CLI."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    out_parquet = tmp_path / "dmc.parquet"
    sibling = tmp_path / "dmc.significant.tsv"

    monkeypatch.setenv("EPYKIT_NO_AUTO_CSV", "1")
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
    ])
    main()

    assert out_parquet.exists()
    assert not sibling.exists(), (
        "EPYKIT_NO_AUTO_CSV=1 must suppress the sibling write"
    )


def test_explicit_csv_path_wins_over_auto_emit_name(tmp_path, synth_bundle, monkeypatch):
    """`--csv` flag overrides the derived `<stem>.significant.tsv` name and
    picks the delimiter from the explicit path suffix."""
    import sys
    import epykit as ep
    from epykit.cli import main

    md_setup = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )
    populated_store = md_setup.store

    out_parquet = tmp_path / "dmc.parquet"
    explicit = tmp_path / "my_hits.csv"          # .csv suffix -> comma delim
    default_sibling = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", populated_store,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
        "--csv", str(explicit),
    ])
    main()

    assert explicit.exists()
    assert not default_sibling.exists()
    header = explicit.read_text(encoding="utf-8").splitlines()[0]
    assert "," in header and "\t" not in header
