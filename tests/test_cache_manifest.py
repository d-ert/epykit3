"""Characterization tests for the pipeline-manifest primitives in
``epykit._cache`` that ``tl.dmc(resumable=True)`` relies on."""

from __future__ import annotations

import os

import polars as pl

from epykit._cache import count_store_rows, input_signature, manifest_append, manifest_read


def test_manifest_append_replaces_same_stage_entry(tmp_path):
    """Recording a stage twice keeps one entry, carrying the latest params,
    moved to the end of ``stages`` (most recent completion order). Other
    stages are untouched."""
    manifest_append(tmp_path, "filtered", params={"lo": 5}, input_sig="f1", output_path="f")
    manifest_append(tmp_path, "dmc_lr", params={"v": 1}, input_sig="A", output_path="A.parquet")
    manifest_append(tmp_path, "united", params={}, input_sig="u1", output_path="u")
    manifest_append(
        tmp_path,
        "dmc_lr",
        params={"v": 2},
        input_sig="B",
        output_path="B.parquet",
        extra={"n_sites": 7},
    )

    stages = manifest_read(tmp_path)["stages"]
    assert [s["name"] for s in stages] == ["filtered", "united", "dmc_lr"]
    (entry,) = [s for s in stages if s["name"] == "dmc_lr"]
    assert entry["params"] == {"v": 2}
    assert entry["input_sig"] == "B"
    assert entry["output_path"] == "B.parquet"
    assert entry["extra"] == {"n_sites": 7}


def test_input_signature_ignores_dict_key_order(tmp_path):
    assert input_signature("stage", {"a": 1, "b": [1, 2]}) == input_signature(
        "stage", {"b": [1, 2], "a": 1}
    )
    assert input_signature("stage", {"a": 1}) != input_signature("stage", {"a": 2})


def test_input_signature_tracks_referenced_file_size_and_mtime(tmp_path):
    """A path argument that exists is fingerprinted by (path, size, mtime_ns):
    appending bytes or touching the file changes the signature, and the
    signature is otherwise stable across calls."""
    data = tmp_path / "input.parquet"
    data.write_bytes(b"abc")
    sig = input_signature(str(data), {"k": 1})
    assert sig == input_signature(str(data), {"k": 1})

    data.write_bytes(b"abcd")
    sig_size = input_signature(str(data), {"k": 1})
    assert sig_size != sig

    stat = data.stat()
    os.utime(data, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    assert input_signature(str(data), {"k": 1}) != sig_size


def test_count_store_rows_sums_all_part_files(tmp_path):
    store = tmp_path / "store"
    for sample, chrom, n in (("a", "chr1", 3), ("a", "chr2", 5), ("b", "chr1", 4)):
        part_dir = store / f"sample={sample}" / f"chrom={chrom}"
        part_dir.mkdir(parents=True)
        pl.DataFrame({"pos": list(range(n))}).write_parquet(part_dir / "part-0.parquet")
    (store / "sample=a" / "chrom=chr1" / "notes.parquet").write_bytes(b"")  # not a part file

    assert count_store_rows(str(store)) == 12
