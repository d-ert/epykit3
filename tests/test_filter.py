"""Characterization tests for ``epykit.filter``.

Pins three behaviours a refactor of the store-level filter step could
break: blacklist interval semantics, per-sample coverage normalisation,
and the per-sample filter manifest that lets a rerun skip finished work.

The fixture is a tiny raw methylstore (2 vs 2 samples, two chromosomes,
50 CpGs each) built once per module through ``read_bismark``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit import _cache
from epykit.filter import FILTER_MANIFEST_NAME, filter_sites, normalize_coverage_store
from tests.fixtures.synth import SimConfig, generate


@pytest.fixture(scope="module")
def raw_store(tmp_path_factory) -> Path:
    cfg = SimConfig(
        n_per_group=2,
        chromosomes=("chr1", "chr2"),
        cpgs_per_chrom=50,
        n_dmrs=1,
        n_scattered_dmcs=5,
    )
    out_dir = tmp_path_factory.mktemp("filter_synth")
    res = generate(cfg, out_dir)
    md = ep.read_bismark(
        res["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(out_dir / "store"),
    )
    return Path(md.store)


def _samples(store: Path) -> list[str]:
    return sorted(d.name.removeprefix("sample=") for d in store.glob("sample=*"))


def _read_sample(store: Path, sample: str) -> pl.DataFrame:
    parts = sorted((store / f"sample={sample}").glob("chrom=*/part-*.parquet"))
    return pl.concat([pl.read_parquet(str(p)) for p in parts]).sort(["chrom", "pos"])


def _sites(df: pl.DataFrame) -> set[tuple[str, int]]:
    return set(zip(df["chrom"].to_list(), df["pos"].to_list()))


# Identity coverage bounds: every synthetic site has coverage >= 1, and
# quantile 1.0 is the per-sample maximum, so only the blacklist can drop rows.
_KEEP_ALL = {"min_coverage": 1, "max_coverage_quantile": 1.0}


def test_blacklist_drops_exactly_the_sites_inside_bed_intervals(raw_store, tmp_path):
    """A BED interval is 0-based half-open: site ``pos`` is dropped iff
    ``start <= pos < end``. The interval end itself is kept, an interval
    between two sites hits nothing, and chromosomes are matched by name."""
    samples = _samples(raw_store)
    chr1 = _read_sample(raw_store, samples[0]).filter(pl.col("chrom") == "chr1")["pos"].to_list()
    chr2 = _read_sample(raw_store, samples[0]).filter(pl.col("chrom") == "chr2")["pos"].to_list()
    p0, p1, p2, p3, p4, p5 = chr1[:6]
    q0, q1, q2, q3 = chr2[:4]

    bed = tmp_path / "blacklist.bed"
    bed.write_text(
        f"chr1\t{p0}\t{p0 + 1}\n"  # single-site interval -> drops p0 only
        f"chr1\t{p2}\t{p4}\n"  # covers p2, p3; ends AT p4 -> p4 kept
        f"chr1\t{p4 + 1}\t{p5}\n"  # falls between two sites -> no hit
        f"chr2\t0\t{q3}\n"  # chromosome scoping: drops q0..q2 on chr2 only
    )
    expected_dropped = {
        ("chr1", p0),
        ("chr1", p2),
        ("chr1", p3),
        ("chr2", q0),
        ("chr2", q1),
        ("chr2", q2),
    }
    assert ("chr1", p1) not in expected_dropped  # the neighbour between two hits is kept

    out = tmp_path / "filtered"
    filter_sites(str(raw_store), str(out), blacklist_bed=str(bed), **_KEEP_ALL)

    for sample in samples:
        before = _read_sample(raw_store, sample)
        after = _read_sample(out, sample)
        assert _sites(before) - _sites(after) == expected_dropped, sample
        # Kept rows are carried over untouched (same columns, same values).
        keep = pl.Series([(c, p) not in expected_dropped for c, p in _sites_ordered(before)])
        assert after.equals(before.filter(keep)), sample


def _sites_ordered(df: pl.DataFrame) -> list[tuple[str, int]]:
    return list(zip(df["chrom"].to_list(), df["pos"].to_list()))


def test_normalize_coverage_store_aligns_medians_and_rescales_counts(raw_store, tmp_path):
    """``method="median"``: factor_i = median(per-sample medians) / median_i,
    and every row's ``N_meth`` and ``coverage`` are the original counts
    scaled by that factor, rounded to int, with ``coverage`` rebuilt as the
    sum of the two scaled parts (so it stays consistent after rounding)."""
    out = tmp_path / "normalized"
    factors = normalize_coverage_store(str(raw_store), str(out), method="median")

    samples = _samples(raw_store)
    assert set(factors) == set(samples)
    medians = {s: float(_read_sample(raw_store, s)["coverage"].median()) for s in samples}
    target = float(np.median(list(medians.values())))
    assert factors == pytest.approx({s: target / medians[s] for s in samples})

    for sample in samples:
        f = factors[sample]
        before = _read_sample(raw_store, sample)
        after = _read_sample(out, sample)
        n_meth = (before["N_meth"].cast(pl.Float64) * f).round().cast(pl.Int32)
        n_unmeth = ((before["coverage"] - before["N_meth"]).cast(pl.Float64) * f).round()
        expected_cov = n_meth + n_unmeth.cast(pl.Int32)
        assert after["N_meth"].to_list() == n_meth.to_list(), sample
        assert after["coverage"].to_list() == expected_cov.to_list(), sample
        assert after.select("chrom", "pos").equals(before.select("chrom", "pos")), sample


@pytest.mark.xfail(
    strict=True,
    reason=(
        "normalize_coverage_store rebuilds coverage from the scaled parts but "
        "leaves the store's N_unmeth column at its original value, so the "
        "documented invariant coverage == N_meth + N_unmeth does not hold on "
        "the normalised store. Fixing that is a behaviour change (out of scope "
        "for the tl.dmc refactor series); this xfail turns into a failure when "
        "it is fixed, so the marker gets removed with the fix."
    ),
)
def test_normalize_keeps_n_unmeth_consistent_with_coverage(raw_store, tmp_path):
    out = tmp_path / "normalized"
    normalize_coverage_store(str(raw_store), str(out), method="median")
    for sample in _samples(raw_store):
        after = _read_sample(out, sample)
        assert after.filter(pl.col("coverage") != pl.col("N_meth") + pl.col("N_unmeth")).is_empty()


_FILTER_PARAMS = {"min_coverage": 3, "max_coverage_quantile": 0.99}


def test_filter_manifest_skips_rerun_with_identical_params(raw_store, tmp_path, caplog):
    """The per-sample ``.epykit_filter_manifest.json`` fingerprints the
    upstream sample and the params; a rerun with the same params leaves every
    output file untouched and logs each sample as cached."""
    out = tmp_path / "filtered"
    filter_sites(str(raw_store), str(out), **_FILTER_PARAMS)

    samples = _samples(raw_store)
    manifests = {s: _cache.load_json(out / f"sample={s}" / FILTER_MANIFEST_NAME) for s in samples}
    for sample, manifest in manifests.items():
        assert manifest is not None, sample
        assert manifest["params"] == {
            "min_coverage": 3,
            "max_coverage_quantile": 0.99,
            "blacklist_bed_sig": None,
        }
        assert manifest["chroms"] == ["chrom=chr1", "chrom=chr2"]
    parts = sorted(out.rglob("part-*.parquet"))
    assert len(parts) == len(samples) * 2
    mtimes = {p: p.stat().st_mtime_ns for p in parts}

    with caplog.at_level(logging.INFO, logger="epykit.filter"):
        filter_sites(str(raw_store), str(out), **_FILTER_PARAMS)

    assert {p: p.stat().st_mtime_ns for p in parts} == mtimes
    assert {
        s: _cache.load_json(out / f"sample={s}" / FILTER_MANIFEST_NAME) for s in samples
    } == manifests
    cached = [r.getMessage() for r in caplog.records if r.getMessage().endswith(": cached")]
    assert len(cached) == len(samples), caplog.text


def test_filter_manifest_invalidates_when_params_change(raw_store, tmp_path, caplog):
    """A different ``min_coverage`` misses the manifest and rewrites the sample."""
    out = tmp_path / "filtered"
    filter_sites(str(raw_store), str(out), **_FILTER_PARAMS)
    sample = _samples(raw_store)[0]
    manifest_path = out / f"sample={sample}" / FILTER_MANIFEST_NAME
    first = _cache.load_json(manifest_path)

    with caplog.at_level(logging.INFO, logger="epykit.filter"):
        filter_sites(str(raw_store), str(out), min_coverage=4, max_coverage_quantile=0.99)

    second = _cache.load_json(manifest_path)
    assert second["params"]["min_coverage"] == 4
    assert second["source"] == first["source"]  # same upstream, only params moved
    assert not any(r.getMessage().endswith(": cached") for r in caplog.records), caplog.text
