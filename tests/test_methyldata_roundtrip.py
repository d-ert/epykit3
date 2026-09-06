"""``MethylData.save()`` / ``load()`` round-trip contract.

Uses a small synthetic cohort (3 vs 3, one chromosome, 200 CpGs) so the
whole module stays in the fast tier.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit.methyldata import MethylData
from tests.fixtures.synth import SimConfig, generate


@pytest.fixture(scope="module")
def samplesheet(tmp_path_factory) -> str:
    cfg = SimConfig(
        n_per_group=3,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    return generate(cfg, tmp_path_factory.mktemp("roundtrip_synth"))["samplesheet"]


def _filtered_md(samplesheet: str, store_dir: Path) -> MethylData:
    md = ep.read_bismark(
        samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(store_dir),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


@dataclass
class RoundTrip:
    original: MethylData
    loaded: MethylData
    save_dir: Path
    meta: dict


@pytest.fixture(scope="module")
def roundtrip(samplesheet, tmp_path_factory) -> RoundTrip:
    """One lr DMC run, plus a DataFrame and a transient cache on ``uns``,
    saved and loaded once for the module."""
    root = tmp_path_factory.mktemp("roundtrip")
    md = _filtered_md(samplesheet, root / "store")
    ep.tl.dmc(md, test="lr", tsv=False)
    md.uns["dmr"] = pl.DataFrame(
        {"chrom": ["chr1", "chr1"], "start": [1000, 5000], "end": [1500, 5200]}
    )
    md.uns["_report_cache"] = {"pca": object()}

    save_dir = root / "saved"
    md.save(str(save_dir))
    meta = json.loads((save_dir / "methyldata.json").read_text())
    return RoundTrip(original=md, loaded=ep.load(str(save_dir)), save_dir=save_dir, meta=meta)


def test_dmcstore_backed_varm_round_trips_through_linked_chrom_files(roundtrip):
    """The DMC table that ``uns["dmc"]["last_key"]`` points at is saved by
    linking the DMCStore's per-chromosome parquets (no re-encode), flagged
    ``varm_format == "dmcstore"``, and loads back equal."""
    assert roundtrip.meta["varm_keys"] == ["dmc_lr"]
    assert roundtrip.meta["varm_format"] == {"dmc_lr": "dmcstore"}
    linked = roundtrip.save_dir / "varm_dmc_lr"
    assert sorted(p.name for p in linked.iterdir()) == [
        ".epykit_dmc_manifest.json",
        "chrom=chr1.parquet",
    ]
    assert not (roundtrip.save_dir / "varm_dmc_lr.parquet").exists()

    original = roundtrip.original.varm["dmc_lr"]
    loaded = roundtrip.loaded.varm["dmc_lr"]
    assert loaded.columns == original.columns
    assert loaded.equals(original)


def test_pvalue_combined_table_falls_back_to_a_single_parquet(samplesheet, tmp_path):
    """``neighbour_combine=True`` adds ``pvalue_combined`` (and friends) to the
    in-memory table only, after the DMCStore files were written. ``save()``
    must not link those files; it writes one parquet so the extra columns
    survive, and flags ``varm_format == "parquet"``."""
    md = _filtered_md(samplesheet, tmp_path / "store")
    ep.tl.dmc(md, test="lr", neighbour_combine=True, tsv=False)
    assert "pvalue_combined" in md.varm["dmc_lr"].columns

    save_dir = tmp_path / "saved"
    md.save(str(save_dir))
    meta = json.loads((save_dir / "methyldata.json").read_text())
    assert meta["varm_format"] == {"dmc_lr": "parquet"}
    assert (save_dir / "varm_dmc_lr.parquet").exists()
    assert not (save_dir / "varm_dmc_lr").exists()

    loaded = ep.load(str(save_dir)).varm["dmc_lr"]
    assert loaded.columns == md.varm["dmc_lr"].columns
    assert loaded.equals(md.varm["dmc_lr"])


def test_uns_dataframes_become_parquet_sidecars_and_json_values_round_trip(roundtrip):
    original, loaded, meta = roundtrip.original, roundtrip.loaded, roundtrip.meta

    assert meta["uns"]["dmr"] == {"__parquet__": "uns_dmr.parquet"}
    assert (roundtrip.save_dir / "uns_dmr.parquet").exists()
    assert loaded.uns["dmr"].equals(original.uns["dmr"])

    for key in ("filter", "unite", "dmc", "_store_history", "n_sites_filtered"):
        assert loaded.uns[key] == original.uns[key], key


def test_report_cache_is_not_persisted(roundtrip):
    assert "_report_cache" in roundtrip.original.uns  # still there in memory
    assert "_report_cache" not in roundtrip.meta["uns"]
    assert "_report_cache" not in roundtrip.loaded.uns


def test_state_dmc_pointers_and_sample_ids_survive_reload(roundtrip):
    original, loaded = roundtrip.original, roundtrip.loaded

    assert loaded.state == original.state == ["raw", "filtered", "united"]
    assert loaded.treatment_ids == original.treatment_ids
    assert loaded.control_ids == original.control_ids
    assert loaded.assembly == "synth"

    # .dmc resolves through uns["dmc"]["last_key"] to the reloaded table.
    assert loaded.dmc is loaded.varm["dmc_lr"]
    # .dmc_store re-opens the on-disk DMCStore recorded in uns["dmc"]["store_path"].
    store = loaded.dmc_store
    assert store is not None
    assert store.path == Path(original.uns["dmc"]["store_path"])
    assert store.total_sites == original.uns["dmc"]["n_sites"] == len(loaded.dmc)
