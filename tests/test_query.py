"""Tabix-on-Parquet random-access query tests (0.4.0).

The contract:

  1. ``query_region`` returns the same rows as a naive full-scan + filter
     on the same store. (Correctness.)
  2. ``query_regions`` is row-wise equivalent to looping ``query_region``
     over each row of the input frame.
  3. ``query_sites`` returns only rows at the exact requested
     (chrom, pos) -- never bystanders from the partition.
  4. Querying a chromosome / sample that doesn't exist returns an empty
     frame, not an error.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit.query import query_region, query_regions, query_sites


def _naive_full_scan(store: str, chrom: str, start: int, end: int) -> pl.DataFrame:
    """Reference implementation: read every partition for the chrom, filter."""
    store_p = Path(store)
    rows = []
    for sample_dir in sorted(store_p.glob("sample=*")):
        sample = sample_dir.name.removeprefix("sample=")
        part = sample_dir / f"chrom={chrom}" / "part-0.parquet"
        if not part.exists():
            continue
        df = pl.read_parquet(str(part), columns=["pos", "strand", "N_meth", "coverage"])
        df = df.filter((pl.col("pos") >= start) & (pl.col("pos") < end))
        if df.height == 0:
            continue
        rows.append(df.with_columns(pl.lit(sample).alias("sample_id"),
                                     pl.lit(chrom).alias("chrom")))
    if not rows:
        return pl.DataFrame(schema={"sample_id": pl.Utf8, "chrom": pl.Utf8, "pos": pl.Int32})
    return pl.concat(rows).sort(["sample_id", "pos"])


def test_query_region_matches_naive_scan(synth_md_filtered):
    """query_region must produce the same rows as a full-scan-then-filter."""
    md = synth_md_filtered
    # Pick a chrom + range that has data.
    chrom = "chr1"
    # Build a window over the middle of the chrom's pos range.
    sample = (md.obs.get_column("sample_id").to_list())[0]
    part = Path(md.store) / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
    if not part.exists():
        pytest.skip(f"synth store has no {chrom} for sample {sample}")
    positions = pl.read_parquet(str(part), columns=["pos"])["pos"].to_numpy()
    if len(positions) < 10:
        pytest.skip("synth chrom too small for a meaningful range query")
    start = int(positions[len(positions) // 4])
    end = int(positions[3 * len(positions) // 4])

    fast = query_region(md.store, chrom, start, end)
    naive = _naive_full_scan(md.store, chrom, start, end)

    # Compare on the common columns. The fast path adds beta + reorders;
    # the naive scan omits beta.
    fast_keys = fast.select(["sample_id", "chrom", "pos"]).sort(["sample_id", "pos"])
    naive_keys = naive.select(["sample_id", "chrom", "pos"]).sort(["sample_id", "pos"])
    assert fast_keys.equals(naive_keys), (
        f"query_region returned different rows than naive scan: "
        f"{fast_keys.height} vs {naive_keys.height}"
    )


def test_query_region_empty_for_unknown_chrom(synth_md_filtered):
    md = synth_md_filtered
    result = query_region(md.store, "chrUnknown", 0, 1_000_000)
    assert result.height == 0


def test_query_region_rejects_bad_interval(synth_md_filtered):
    md = synth_md_filtered
    with pytest.raises(ValueError, match="must be > start"):
        query_region(md.store, "chr1", 100, 100)


def test_query_region_filters_samples(synth_md_filtered):
    """Passing samples= restricts the result to those samples only."""
    md = synth_md_filtered
    all_samples = sorted(md.obs.get_column("sample_id").to_list())
    if len(all_samples) < 2:
        pytest.skip("need >= 2 samples")
    target = all_samples[0]
    result = query_region(md.store, "chr1", 0, 10_000_000_000, samples=[target])
    if result.height == 0:
        pytest.skip("synth store empty on chr1")
    returned_samples = set(result["sample_id"].to_list())
    assert returned_samples == {target}


def test_query_regions_matches_loop(synth_md_filtered):
    """Batched query_regions == concat of per-region query_region calls."""
    md = synth_md_filtered
    regions = pl.DataFrame({
        "chrom": ["chr1", "chr2"],
        "start": [0, 0],
        "end":   [1_000_000_000, 1_000_000_000],
    })
    batched = query_regions(md.store, regions).sort(["region_id", "sample_id", "pos"])
    looped_parts = []
    for region_id, (chrom, start, end) in enumerate(zip(
        regions["chrom"], regions["start"], regions["end"]
    )):
        df = query_region(md.store, chrom, int(start), int(end))
        if df.height:
            looped_parts.append(df.with_columns(
                pl.lit(region_id).cast(pl.Int32).alias("region_id")
            ))
    if not looped_parts:
        assert batched.height == 0
        return
    looped = pl.concat(looped_parts).sort(["region_id", "sample_id", "pos"])
    assert batched.equals(looped)


def test_query_sites_returns_only_requested_positions(synth_md_filtered):
    """query_sites must return rows at exactly the requested (chrom, pos)."""
    md = synth_md_filtered
    # Take 5 known positions from the store.
    sample = md.obs.get_column("sample_id").to_list()[0]
    part = Path(md.store) / f"sample={sample}" / "chrom=chr1" / "part-0.parquet"
    if not part.exists():
        pytest.skip("no chr1 partition for first sample")
    pos_some = pl.read_parquet(str(part), columns=["pos"])["pos"].to_list()[:5]
    if len(pos_some) < 5:
        pytest.skip("need >= 5 positions")
    sites = pl.DataFrame({
        "chrom": ["chr1"] * len(pos_some),
        "pos":   pos_some,
    })
    result = query_sites(md.store, sites)
    assert result.height > 0
    # Every returned pos must be in the requested set.
    returned_pos = set(result["pos"].to_list())
    assert returned_pos.issubset(set(pos_some))


def test_query_module_is_exposed_on_epykit():
    """ep.query.query_region works without an explicit import."""
    assert hasattr(ep, "query")
    assert callable(ep.query.query_region)
    assert callable(ep.query.query_regions)
    assert callable(ep.query.query_sites)
