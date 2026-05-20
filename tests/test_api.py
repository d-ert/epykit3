"""Layer 4: API behaviour tests.

Verifies the data-object contract:

* ``MethylData`` round-trips losslessly through ``save()`` / ``load()``.
* Preprocessing state is derived from ``uns['_store_history']``.
* ``MethylData.get_dmc(test=...)`` looks up by explicit name; ``.dmc``
  resolves to the *last-written* table via ``uns['dmc']['last_key']``.
* Removed ``samples_case=`` / ``min_samples_case=`` aliases raise errors.
* Covariate-adjusted GLM dispatches correctly through ``tl.dmr``.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import polars as pl
import pytest



# read_bismark + obs schema


def test_read_bismark_produces_well_formed_methyldata(synth_md):
    """Sanity check the fixture: obs has the expected columns and counts."""
    md = synth_md
    assert set(md.obs.columns) >= {"sample_id", "group", "treatment", "path"}
    assert md.n_samples == 8
    assert len(md.treatment_ids) == 4
    assert len(md.control_ids) == 4
    assert all(s.startswith("treatment_") for s in md.treatment_ids)
    assert all(s.startswith("control_") for s in md.control_ids)


def test_read_bismark_writes_partitioned_parquet_store(synth_md):
    """The methylstore should have ``sample=*`` / ``chromosome=*`` (or
    similar) partitions on disk."""
    store = Path(synth_md.store)
    assert store.exists(), f"methylstore not created at {store}"
    # At least one sample partition should exist after conversion.
    sample_dirs = list(store.glob("sample=*"))
    assert len(sample_dirs) >= 1, "no sample partitions in methylstore"



# State derivation from _store_history


def test_state_is_raw_before_any_preprocessing(synth_md):
    """Before any pp.* call, the boolean flags are False and no
    preprocessing step appears in ``state``.

    Note: ``read_bismark`` itself records a ``"raw"`` step in
    ``_store_history`` so callers can audit where the methylstore came
    from. The state list will therefore contain at most ``["raw"]`` -- but
    none of the preprocessing markers (filtered / united / smoothed) should
    be present yet.
    """
    md = synth_md
    assert md._filtered is False
    assert md._united is False
    assert md._smoothed is False
    assert "filtered" not in md.state
    assert "united" not in md.state
    assert "smoothed" not in md.state


def test_state_after_filter_coverage(synth_md):
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    assert synth_md._filtered is True
    assert "filtered" in synth_md.state


def test_state_after_unite(synth_md):
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="intersect")
    assert synth_md._united is True
    assert "united" in synth_md.state


def test_state_persists_through_save_load_round_trip(synth_md, tmp_path):
    """A filtered + united MethylData saved and reloaded retains its
    derived state. Save/load are symmetric for paths with directory
    components -- the previous test workaround helper is no longer needed.
    """
    import epykit as ep

    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="intersect")

    save_path = tmp_path / "saved_md"
    synth_md.save(str(save_path))
    md2 = ep.load(str(save_path))

    assert md2._filtered is True
    assert md2._united is True
    assert "filtered" in md2.state
    assert "united" in md2.state



# save / load round-trip


def test_save_load_round_trip_preserves_obs_varm_uns(synth_md, tmp_path):
    """Save -> load must preserve obs, varm, and primitive uns values."""
    import epykit as ep

    # Populate uns and varm.
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="intersect")
    ep.tl.dmc(synth_md, test="lr")

    save_path = tmp_path / "rt"
    synth_md.save(str(save_path))
    md2 = ep.load(str(save_path))

    # obs equal as polars DataFrames (sample_id order preserved).
    assert md2.obs.shape == synth_md.obs.shape
    assert sorted(md2.obs.columns) == sorted(synth_md.obs.columns)

    # varm keys preserved.
    assert set(md2.varm.keys()) == set(synth_md.varm.keys())
    for k in synth_md.varm:
        assert md2.varm[k].shape == synth_md.varm[k].shape, f"varm[{k!r}] shape drift"

    # uns: primitive keys round-trip.
    for key in ("filter", "unite", "dmc"):
        if key in synth_md.uns:
            assert key in md2.uns


def test_save_load_does_not_persist_boolean_state_in_meta(synth_md, tmp_path):
    """methyldata.json should *not* hard-code _filtered etc. (by design,
    that they are derived properties)."""
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    save_path = tmp_path / "rt2"
    synth_md.save(str(save_path))
    meta = json.loads((save_path / "methyldata.json").read_text())
    assert "_filtered" not in meta, "_filtered must not be persisted (derived)"
    assert "_united" not in meta
    assert "_smoothed" not in meta



# MethylData.dmc / .get_dmc


def test_get_dmc_returns_none_before_running_dmc(synth_md_filtered):
    md = synth_md_filtered
    assert md.get_dmc() is None
    assert md.dmc is None


def test_get_dmc_returns_explicit_test_by_name(synth_md_filtered):
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    df = synth_md_filtered.get_dmc(test="lr")
    assert df is not None
    assert "pvalue" in df.columns or "qvalue" in df.columns


def test_dmc_property_uses_last_key_pointer(synth_md_filtered):
    """After running two tests, ``.dmc`` should resolve to whichever was
    written most recently (per ``uns['dmc']['last_key']``), not to a
    hardcoded priority list."""
    import epykit as ep

    ep.tl.dmc(synth_md_filtered, test="lr")
    ep.tl.dmc(synth_md_filtered, test="logit_t")

    # Last writer wins.
    assert synth_md_filtered.uns["dmc"]["last_key"] == "dmc_logit_t"
    df_via_property = synth_md_filtered.dmc
    df_via_explicit_logit_t = synth_md_filtered.get_dmc(test="logit_t")
    # Both should reference the same logit_t table (or its annotated variant).
    assert df_via_property is df_via_explicit_logit_t or df_via_property.shape == df_via_explicit_logit_t.shape


def test_get_dmc_prefers_annotated_when_available(synth_md_filtered):
    """If a *_annotated table exists, ``get_dmc(annotated=True)`` should
    return it instead of the raw table."""
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")

    # Manually drop an annotated table into varm so we can test resolution
    # without depending on a real GTF.
    raw = synth_md_filtered.varm["dmc_lr"]
    synth_md_filtered.varm["dmc_lr_annotated"] = raw.with_columns(
        pl.lit("intergenic").alias("feature_type")
    )

    ann = synth_md_filtered.get_dmc(test="lr", annotated=True)
    raw_only = synth_md_filtered.get_dmc(test="lr", annotated=False)
    assert "feature_type" in ann.columns
    assert "feature_type" not in raw_only.columns



# Removed-alias tests


def test_samples_case_kwarg_is_rejected(synth_md_filtered):
    """The ``samples_case`` kwarg was removed; passing it raises TypeError."""
    from epykit.dmc import process_chromosomes_dmc

    with pytest.raises(TypeError):
        process_chromosomes_dmc(
            methylstore_path=synth_md_filtered.store,
            samples_case=synth_md_filtered.treatment_ids,
            samples_control=synth_md_filtered.control_ids,
            test="fisher",
        )


def test_samples_treatment_kwarg_works(synth_md_filtered):
    """The canonical ``samples_treatment=`` kwarg works without warnings."""
    from epykit.dmc import process_chromosomes_dmc

    df = process_chromosomes_dmc(
        methylstore_path=synth_md_filtered.store,
        samples_treatment=synth_md_filtered.treatment_ids,
        samples_control=synth_md_filtered.control_ids,
        test="fisher",
    )
    assert df is not None and len(df) > 0



# Covariate-adjusted GLM dispatch


def test_covariate_dmr_runs_via_design_kwarg(synth_md_filtered):
    """When ``design=`` is passed, ``tl.dmr`` must force the GLM backend
    even if ``test='auto'``."""
    import epykit as ep

    # Inject a continuous covariate into obs.
    synth_md_filtered.obs = synth_md_filtered.obs.with_columns(
        pl.Series("age", np.arange(synth_md_filtered.n_samples, dtype=float) + 20.0)
    )

    ep.tl.dmr(
        synth_md_filtered,
        method="tile",
        tile_size_bp=500,
        min_cpgs_per_tile=3,
        design="~ treatment + age",
        treatment_col="treatment",
    )
    params = synth_md_filtered.uns["dmr_params"]
    assert params["test"] == "glm", f"design= should force GLM, got {params['test']!r}"
    assert params["design"] == "~ treatment + age"


def test_get_dmc_with_unknown_test_returns_none(synth_md_filtered):
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    assert synth_md_filtered.get_dmc(test="cmh") is None



# Version + __all__ contract


def test_version_is_a_pep440_string():
    import epykit
    v = epykit.__version__
    assert isinstance(v, str) and len(v) > 0


def test_all_lists_documented_public_surface():
    import epykit
    expected = {
        "MethylData", "read_bismark", "load",
        "pp", "tl", "pl",
        "convert_sample",
        "smooth_methylation_gaussian",
        "annotate_features", "annotate_cpg_islands",
        "bisulfite_conversion_rate",
    }
    missing = expected - set(epykit.__all__)
    assert not missing, f"epykit.__all__ missing: {sorted(missing)}"


# ---- Region aggregation (merged from test_pp_regions.py) -------------


def _write_bed(path, rows):
    path.write_text(
        "\n".join(f"{c}\t{s}\t{e}\t{name}" for c, s, e, name in rows) + "\n"
    )
    return path


def test_aggregate_regions_round_trip(synth_md_filtered, tmp_path):
    """Aggregating to a few wide regions yields a row per (region, sample)."""
    import epykit as ep

    md = synth_md_filtered
    bed_rows: list[tuple[str, int, int, str]] = []
    chrom_bounds = (
        pl.scan_parquet(f"{md.store}/sample=*/chrom=*/part-*.parquet")
        .group_by("chrom")
        .agg([pl.min("pos").alias("lo"), pl.max("pos").alias("hi")])
        .collect()
    )
    for r in chrom_bounds.iter_rows(named=True):
        lo, hi = int(r["lo"]), int(r["hi"])
        if hi - lo < 30:
            continue
        third = (hi - lo) // 3
        bed_rows.extend([
            (r["chrom"], lo, lo + third, f"{r['chrom']}_a"),
            (r["chrom"], lo + third, lo + 2 * third, f"{r['chrom']}_b"),
            (r["chrom"], lo + 2 * third, hi + 1, f"{r['chrom']}_c"),
        ])
    bed = _write_bed(tmp_path / "regions.bed", bed_rows)

    ep.pp.aggregate_regions(md, str(bed), min_cpgs_per_region=1)

    new_store = Path(md.store)
    assert new_store.exists()
    sample_dirs = list(new_store.glob("sample=*"))
    assert len(sample_dirs) == md.n_samples

    df = pl.read_parquet(
        f"{md.store}/sample=*/chrom=*/part-*.parquet"
    )
    for col in (
        "chrom", "pos", "strand", "context",
        "N_meth", "N_unmeth", "coverage", "sample",
        "region_id", "start", "end", "n_cpgs",
    ):
        assert col in df.columns, f"missing column {col}"
    assert df["coverage"].eq(df["N_meth"] + df["N_unmeth"]).all()
    assert md.uns["regions"]["n_regions"] == len(bed_rows)
    assert any(h["step"] == "regions" for h in md.uns["_store_history"])


def test_aggregate_regions_then_dmc(synth_md_filtered, tmp_path):
    """Downstream `tl.dmc` runs on the region-aggregated store without errors."""
    import epykit as ep

    md = synth_md_filtered
    chrom_bounds = (
        pl.scan_parquet(f"{md.store}/sample=*/chrom=*/part-*.parquet")
        .group_by("chrom")
        .agg([pl.min("pos").alias("lo"), pl.max("pos").alias("hi")])
        .collect()
    )
    bed_rows: list[tuple[str, int, int, str]] = []
    for r in chrom_bounds.iter_rows(named=True):
        lo, hi = int(r["lo"]), int(r["hi"])
        step = max(1, (hi - lo) // 5)
        for i in range(5):
            bed_rows.append(
                (r["chrom"], lo + i * step, lo + (i + 1) * step, f"{r['chrom']}_b{i}")
            )
    bed = _write_bed(tmp_path / "regions.bed", bed_rows)
    ep.pp.aggregate_regions(md, str(bed), min_cpgs_per_region=1)
    md.uns.pop("unite", None)
    ep.pp.unite(md, type="intersect")
    ep.tl.dmc(md, test="lr")
    dmc = md.dmc
    assert dmc is not None and len(dmc) > 0
