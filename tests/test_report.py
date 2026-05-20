"""Tests for the HTML report generator."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest


pytest.importorskip("jinja2")
pytest.importorskip("plotly")


def test_report_partial_pipeline(synth_md_filtered, tmp_path):
    """A MethylData with only DMC (no DMR/annotation) still renders an HTML
    report and reports the missing sections gracefully.
    """
    import epykit as ep

    ep.tl.dmc(synth_md_filtered, test="lr")

    out = tmp_path / "report_partial.html"
    synth_md_filtered.report(str(out), title="partial pipeline")
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert "<title>partial pipeline</title>" in html
    # DMC section populated, DMR section says "not yet called"
    assert "DMC volcano" in html or "volcano" in html.lower()
    assert "ep.tl.dmr" in html  # the "not run yet" notice for DMR
    # Plotly was bundled via CDN
    assert "plotly" in html.lower()


def test_report_full_pipeline(synth_md_filtered, tmp_path):
    """Full pipeline (dmc + dmr + injected annotation) yields all sections."""
    import epykit as ep
    import numpy as np

    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=3,
              min_mean_qvalue=1.0)
    # Inject fake annotations so the annotated section renders
    dmc = md.varm["dmc_lr"]
    rng = np.random.default_rng(0)
    features = rng.choice(["promoter", "exon", "intron", "intergenic"], len(dmc))
    contexts = rng.choice(["island", "shore", "shelf", "open_sea"], len(dmc))
    md.varm["dmc_lr_annotated"] = dmc.with_columns([
        pl.Series("feature_type", features.tolist()),
        pl.Series("cpg_context", contexts.tolist()),
    ])

    out = tmp_path / "report_full.html"
    md.report(str(out), title="full pipeline", alpha=0.05, min_abs_diff=0.1)
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert "full pipeline" in html
    # Big enough to be a real report
    assert out.stat().st_size > 30_000, f"report too small: {out.stat().st_size} bytes"
    # Sample IDs should appear
    for sid in md.obs.get_column("sample_id").to_list()[:2]:
        assert sid in html
    # Provenance present
    assert "Provenance" in html


def test_report_skips_unrun_sections(synth_md_filtered, tmp_path):
    """No DMC / DMR yet -- every analysis section should render a `not run`
    notice instead of crashing."""
    import epykit as ep

    out = tmp_path / "report_empty.html"
    synth_md_filtered.report(str(out))
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert "DMC not yet called" in html or "ep.tl.dmc" in html
    assert "DMR not yet called" in html or "ep.tl.dmr" in html
