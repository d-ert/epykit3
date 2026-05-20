"""Bismark .cov -> partitioned Parquet converter.

Coordinate system: ``start`` is treated as **0-based** (BED-format), which
is what ``bismark2bedGraph`` and nf-core/methylseq emit. The 1-based
output of ``bismark_methylation_extractor --comprehensive`` /
``coverage2cytosine`` must be pre-shifted by -1 before being passed here;
otherwise every CpG ends up offset by 1 bp relative to GTF / CpG-island
annotations. (Quick check: a 0-based file has ``start = end - 1``; a
1-based file has ``start == end``.)

Strand: not present in Bismark merged .cov files. When ``reference_fasta``
is provided, strand is inferred from the reference base at ``pos`` (``+``
for C, ``-`` otherwise); without a reference it defaults to ``*``.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from . import _cache


RAW_MANIFEST_NAME = ".epykit_raw_manifest.json"

# Bismark .cov column order. MethylDackel's extract --mergeContext output
# uses the same six columns (chrom, start, end, percent, M, U) but prepends
# a single ``track type="bedGraph" ...`` header line. _FORMAT_SKIP_ROWS
# encodes that difference so a polars scan_csv can ingest both without
# format-specific parsing logic.
_COV_COLUMNS = ["chrom", "start", "end", "methyl_percent", "N_meth", "N_unmeth"]

_COV_SCHEMA: dict[str, type[pl.DataType]] = {
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "methyl_percent": pl.Float32,
    "N_meth": pl.Int32,
    "N_unmeth": pl.Int32,
}

_FORMAT_SKIP_ROWS: dict[str, int] = {
    "bismark": 0,
    "methyldackel": 1,
    # 12-col methylation BED (strand-collapsed CpG dyad summary)
    # Layout: chrom, start, end, fwd_M, fwd_T, fwd_%, rev_M, rev_T, rev_%, M, T, %
    "combined_strand_bed": 0,
}


# Wider schema for the 12-col combined-strand BED. We read the full 12
# columns then project to the canonical 6-col Bismark-equivalent layout.
_COMBINED_BED_COLUMNS = [
    "chrom", "start", "end",
    "fwd_M", "fwd_T", "fwd_pct",
    "rev_M", "rev_T", "rev_pct",
    "M", "T", "methyl_percent",
]
_COMBINED_BED_SCHEMA: dict[str, type[pl.DataType]] = {
    "chrom":         pl.Utf8,
    "start":         pl.Int32,
    "end":           pl.Int32,
    "fwd_M":         pl.Int32,
    "fwd_T":         pl.Int32,
    "fwd_pct":       pl.Float32,
    "rev_M":         pl.Int32,
    "rev_T":         pl.Int32,
    "rev_pct":       pl.Float32,
    "M":             pl.Int32,
    "T":             pl.Int32,
    "methyl_percent": pl.Float32,
}


# Manifest helpers

@dataclass(frozen=True)
class _SampleManifest:
    sample_name: str
    source: dict[str, object]
    chroms: list[str]
    row_group_size: int
    format: str = "bismark"


_file_signature = _cache.file_signature
_load_json = _cache.load_json
_write_json = _cache.write_json
_sample_dir = _cache.sample_dir
_expected_chrom_dirs = _cache.expected_chrom_dirs
_sample_is_complete = _cache.sample_is_complete


def _manifest_path(sample_dir: Path) -> Path:
    return sample_dir / RAW_MANIFEST_NAME


def _manifest_payload(manifest: _SampleManifest) -> dict[str, object]:
    return {
        "sample_name": manifest.sample_name,
        "source": manifest.source,
        "chroms": manifest.chroms,
        "row_group_size": manifest.row_group_size,
        "format": manifest.format,
    }


def _can_reuse_sample(
    input_path: Path, sample_dir: Path, row_group_size: int,
    format: str = "bismark",
) -> bool:
    manifest = _load_json(_manifest_path(sample_dir))
    if not manifest:
        return False
    if manifest.get("source") != _file_signature(input_path):
        return False
    if manifest.get("row_group_size") != row_group_size:
        return False
    # Reject cached conversions made under a different source format. The
    # default "bismark" preserves cache compatibility for stores converted
    # before the format key existed.
    if manifest.get("format", "bismark") != format:
        return False
    chroms = manifest.get("chroms")
    if not isinstance(chroms, list) or not all(
        isinstance(c, str) for c in chroms
    ):
        return False
    return _sample_is_complete(sample_dir, chroms)


def _promote_sample_dir(temp_sample_dir: Path, final_sample_dir: Path) -> None:
    backup_dir = final_sample_dir.with_name(f"{final_sample_dir.name}.bak")
    if backup_dir.exists():
        shutil.rmtree(backup_dir)
    if final_sample_dir.exists():
        final_sample_dir.rename(backup_dir)
    try:
        temp_sample_dir.rename(final_sample_dir)
    except Exception:
        if backup_dir.exists() and not final_sample_dir.exists():
            backup_dir.rename(final_sample_dir)
        raise
    else:
        if backup_dir.exists():
            shutil.rmtree(backup_dir)


# CpG strand merging

def _merge_cpg_pairs(df: pl.DataFrame) -> pl.DataFrame:
    """Merge + and - strand CpG pairs into single sites at the + strand position.
    
    When Bismark .cov files contain both strands, a CpG dinucleotide appears as:
      - + strand at position N (C position)
      - - strand at position N+1 (G position on reverse strand)
    
    This function merges them by:
      1. Shifting - strand positions back by 1 (N+1 -> N)
      2. Grouping by (chrom, pos) and summing counts
      3. Setting all merged sites to + strand
    
    Validates pairing and warns if unpaired sites are found (may indicate
    incomplete bisulfite conversion or quality issues).
    
    If the input already has strand-merged data (e.g., from bismark2bedGraph),
    this function is a no-op.
    """
    if "strand" not in df.columns:
        # No strand information, return as-is
        return df
    
    # Separate + and - strands
    plus = df.filter(pl.col("strand") == "+")
    minus = df.filter(pl.col("strand") == "-")
    
    if len(minus) == 0:
        # No - strand data, already merged or only + strand present
        return df
    
    # Shift - strand positions to + strand coordinate (N+1 -> N)
    minus = minus.with_columns(
        (pl.col("pos") - 1).alias("pos")
    )
    
    # VALIDATION: Check for proper pairing. Polars semi-join on
    # (chrom, pos) is cheaper than materialising two Python int sets
    # and is also correct in the (rare) multi-chrom case where the
    # same `pos` value exists on different chromosomes -- the old
    # set-diff conflated them.
    total_plus_orig = len(plus)
    total_minus_orig = len(minus)
    key = ["chrom", "pos"]
    n_paired = plus.join(minus, on=key, how="semi").height
    n_unpaired_plus = total_plus_orig - n_paired
    n_unpaired_minus = total_minus_orig - n_paired

    if n_unpaired_plus or n_unpaired_minus:
        import logging
        logger = logging.getLogger(__name__)
        logger.warning(
            f"CpG strand pairing incomplete:\n"
            f"  Paired sites: {n_paired:,}\n"
            f"  Unpaired + strand: {n_unpaired_plus:,}/{total_plus_orig:,}\n"
            f"  Unpaired - strand: {n_unpaired_minus:,}/{total_minus_orig:,}\n"
            f"  Merging will sum paired sites and keep unpaired sites as-is.\n"
            f"  This may indicate incomplete bisulfite conversion or data quality issues."
        )
    
    # Combine and merge by position
    combined = pl.concat([plus, minus])
    
    # Group by chrom and pos, summing methylation counts
    merged = combined.group_by(["chrom", "pos"], maintain_order=True).agg([
        pl.sum("N_meth").alias("N_meth"),
        pl.sum("N_unmeth").alias("N_unmeth"),
        pl.sum("coverage").alias("coverage"),
        pl.first("sample").alias("sample"),
        pl.first("context").alias("context"),
    ]).with_columns(
        pl.lit("+").alias("strand")
    )
    
    return merged.sort("pos")


# Optional strand inference

def _infer_strand(df: pl.DataFrame, reference_fasta: str) -> pl.Series:
    """Infer strand from reference sequence for each CpG position.

    A cytosine on the + strand sits at position `start` in the reference.
    Its complement on the - strand is at `start + 1`. Bismark merged .cov
    coordinates are 0-based start, 1-based end (BED-like).

    Requires pyfaidx:  pip install pyfaidx

    Parameters
    ----------
    df : pl.DataFrame
        Must contain columns: chrom (str), start (Int32)
    reference_fasta : str
        Path to indexed reference FASTA (.fai index must exist)

    Returns
    -------
    pl.Series (Utf8)
        "+" where the reference base at `start` is C (or c),
        "-" where it is G (complement C on the - strand),
        "*" for anything else (non-CpG context or N base).
    """
    try:
        from pyfaidx import Fasta  # optional dependency
    except ImportError as exc:
        raise ImportError(
            "pyfaidx is required for strand inference. "
            "Install it with: pip install pyfaidx"
        ) from exc

    import numpy as np

    fasta = Fasta(reference_fasta, as_raw=True)

    # Per-chromosome vectorised lookup: load each chromosome sequence
    # once, index by all of its `start` positions in bulk via numpy.
    # Preserves the per-row "+"/"-"/"*" mapping and the KeyError/
    # IndexError -> "*" fallback of the original Python loop.
    chrom_arr = df["chrom"].to_numpy()
    start_arr = df["start"].to_numpy()
    out = np.full(df.height, ord("*"), dtype=np.uint8)

    C, c_lower = ord("C"), ord("c")
    G, g_lower = ord("G"), ord("g")
    PLUS, MINUS = ord("+"), ord("-")

    for chrom in np.unique(chrom_arr):
        try:
            seq_obj = fasta[chrom][:]
        except (KeyError, IndexError):
            continue
        # ``as_raw=True`` makes pyfaidx return either ``str`` or ``bytes``
        # depending on version; normalise to bytes for numpy.
        seq_bytes = seq_obj.encode("ascii") if isinstance(seq_obj, str) \
            else bytes(seq_obj)
        seq = np.frombuffer(seq_bytes, dtype=np.uint8)

        mask = chrom_arr == chrom
        positions = start_arr[mask]
        in_bounds = positions < seq.size
        safe_pos = np.where(in_bounds, positions, 0)
        bases = seq[safe_pos]
        result = np.full(positions.size, ord("*"), dtype=np.uint8)
        result[(bases == C) | (bases == c_lower)] = PLUS
        result[(bases == G) | (bases == g_lower)] = MINUS
        result[~in_bounds] = ord("*")
        out[mask] = result

    return pl.Series(
        "strand",
        [chr(b) for b in out.tolist()],
        dtype=pl.Utf8,
    )


# Public API

def convert_sample(
    input_path: str,
    sample_name: str,
    output_dir: str,
    row_group_size: int = 1_000_000,
    context: str = "CpG",
    reference_fasta: str | None = None,
    merge_strands: bool = True,
    format: str = "bismark",
) -> None:
    """Convert a Bismark .cov or MethylDackel .bedGraph file into a
    partitioned Parquet store.

    Parameters
    ----------
    input_path : str
        Path to the .cov / .cov.gz / .bedGraph / .bedGraph.gz file.
    sample_name : str
        Sample identifier written into the `sample` column
    output_dir : str
        Directory where Parquet partitions will be written
    row_group_size : int
        Approximate Parquet row-group size (default 1 000 000)
    context : str
        Methylation context label stored in the `context` column
        ("CpG", "CHG", "CHH"). Default "CpG".
    reference_fasta : str, optional
        Path to an indexed reference FASTA. When provided, strand is inferred
        from the reference base at each position . Without this
        argument, strand defaults to "*".
    merge_strands : bool
        If True (default), merge + and - strand CpG pairs into single sites
        at the + strand position . This is appropriate for .cov files
        from bismark_methylation_extractor with both strands. Files from
        bismark2bedGraph are typically already strand-merged.
    format : {"bismark", "methyldackel"}
        Source format. Both use the same 6-column layout
        (chrom, start, end, methylation_percent, count_methylated,
        count_unmethylated); MethylDackel prepends a one-line ``track``
        header that is skipped automatically. Default "bismark".

    Output schema
    -------------
    chrom   Utf8
    pos     Int32   (0-based, == Bismark start)
    strand  Utf8    ("+" | "-" | "*")
    context Utf8    ("CpG" | "CHG" | "CHH")
    N_meth  Int32
    N_unmeth Int32
    coverage Int32
    sample  Utf8
    """
    if format not in _FORMAT_SKIP_ROWS:
        raise ValueError(
            f"Unknown format {format!r}; expected one of "
            f"{sorted(_FORMAT_SKIP_ROWS)}"
        )

    p = Path(input_path)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if format == "combined_strand_bed":
        # 12-col methylation BED: use the combined-strand triplet (cols 10-12).
        # Project to the canonical Bismark layout so downstream is unchanged.
        lf = (
            pl.scan_csv(
                str(p),
                separator="\t",
                has_header=False,
                skip_rows=_FORMAT_SKIP_ROWS[format],
                new_columns=_COMBINED_BED_COLUMNS,
                schema_overrides=_COMBINED_BED_SCHEMA,
            )
            .with_columns([
                pl.col("M").alias("N_meth"),
                (pl.col("T") - pl.col("M")).cast(pl.Int32).alias("N_unmeth"),
                pl.col("T").alias("coverage"),
                pl.lit(sample_name).alias("sample"),
                pl.col("start").alias("pos"),
                pl.lit(context).alias("context"),
            ])
            .select(["chrom", "pos", "context", "N_meth", "N_unmeth",
                     "coverage", "sample", "start"])
        )
    else:
        lf = pl.scan_csv(
            str(p),
            separator="\t",
            has_header=False,
            skip_rows=_FORMAT_SKIP_ROWS[format],
            new_columns=_COV_COLUMNS,
            schema_overrides=_COV_SCHEMA,
        ).with_columns(
            [
                (pl.col("N_meth") + pl.col("N_unmeth")).alias("coverage"),
                pl.lit(sample_name).alias("sample"),
                pl.col("start").alias("pos"),
                pl.lit(context).alias("context"),
            ]
        ).select(
            ["chrom", "pos", "context", "N_meth", "N_unmeth", "coverage", "sample",
             "start"]   # keep start temporarily for strand inference
        )

    df = lf.collect()

    # Strand inference : requires reference FASTA via pyfaidx
    if reference_fasta is not None:
        strand_series = _infer_strand(df, reference_fasta)
    else:
        strand_series = pl.Series("strand", ["*"] * len(df), dtype=pl.Utf8)

    df = df.with_columns(strand_series).drop("start").select(
        ["chrom", "pos", "strand", "context", "N_meth", "N_unmeth", "coverage",
         "sample"]
    )

    # CpG strand merging : merge + and - strand pairs if requested
    if merge_strands and reference_fasta is not None:
        df = _merge_cpg_pairs(df)

    # Write one Parquet file per chromosome. partition_by is a single
    # hash-partition pass; the prior unique()+filter() loop scanned the
    # frame once per chromosome.
    for key, sub in df.partition_by(
        "chrom", as_dict=True, maintain_order=False
    ).items():
        chrom = key[0] if isinstance(key, tuple) else key
        part_dir = out / f"sample={sample_name}" / f"chrom={chrom}"
        part_dir.mkdir(parents=True, exist_ok=True)
        sub.write_parquet(
            str(part_dir / "part-0.parquet"),
            compression="zstd",
            row_group_size=row_group_size,
        )


def ensure_converted_sample(
    input_path: str,
    sample_name: str,
    output_dir: str,
    row_group_size: int = 1_000_000,
    context: str = "CpG",
    reference_fasta: str | None = None,
    format: str = "bismark",
) -> bool:
    """Convert a sample unless a valid on-disk conversion already exists.

    Returns True when a fresh conversion was performed, False when the
    existing partitioned store was reused without changes.
    """
    source_path = Path(input_path)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    final_sample_dir = _sample_dir(output_root, sample_name)
    if _can_reuse_sample(
        source_path, final_sample_dir, row_group_size, format=format,
    ):
        return False

    temp_root = output_root.parent / f".{output_root.name}.{sample_name}.tmp"
    if temp_root.exists():
        shutil.rmtree(temp_root)

    try:
        convert_sample(
            input_path,
            sample_name,
            str(temp_root),
            row_group_size=row_group_size,
            context=context,
            reference_fasta=reference_fasta,
            format=format,
        )
        temp_sample_dir = _sample_dir(temp_root, sample_name)
        chroms = _expected_chrom_dirs(temp_sample_dir)
        manifest = _SampleManifest(
            sample_name=sample_name,
            source=_file_signature(source_path),
            chroms=chroms,
            row_group_size=row_group_size,
            format=format,
        )
        _write_json(_manifest_path(temp_sample_dir), _manifest_payload(manifest))
        _promote_sample_dir(temp_sample_dir, final_sample_dir)
    finally:
        if temp_root.exists():
            shutil.rmtree(temp_root)

    return True


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Convert Bismark .cov to partitioned Parquet"
    )
    ap.add_argument("--input", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument(
        "--context",
        default="CpG",
        choices=["CpG", "CHG", "CHH"],
        help="Methylation context (default: CpG)",
    )
    ap.add_argument(
        "--reference-fasta",
        default=None,
        help="Optional reference FASTA for strand inference (requires pyfaidx)",
    )
    args = ap.parse_args()
    convert_sample(
        args.input,
        args.sample_id,
        args.output_dir,
        context=args.context,
        reference_fasta=args.reference_fasta,
    )