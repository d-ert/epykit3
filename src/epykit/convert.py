"""Bismark .cov / MethylDackel .bedGraph -> partitioned Parquet converter.

Coordinate system: the internal methylstore is **0-based** (``pos`` is the
0-based coordinate of the cytosine). Inputs differ in convention:

* Standard Bismark coverage (``.bismark.cov[.gz]`` from
  ``bismark_methylation_extractor`` / ``coverage2cytosine``) is **1-based**:
  each row is a single cytosine with ``start == end``.
* MethylDackel ``.bedGraph`` and the 12-column combined-strand BED are
  **0-based** half-open (``end == start + 1``).

``convert_sample(..., coordinate_base="auto")`` (the default) inspects
``start``/``end`` and shifts 1-based Bismark input by ``-1`` so every source
lands on the same 0-based ``pos``. Pass ``coordinate_base="one_based"`` or
``"zero_based"`` to force the convention. (Quick check: a 0-based file has
``start == end - 1``; a 1-based file has ``start == end``.)

Strand: not present in Bismark merged .cov files. When ``reference_fasta``
is provided, strand is inferred from the reference base; without a
reference it defaults to ``*``.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from . import _cache


RAW_MANIFEST_NAME = ".epykit_raw_manifest.json"

logger = logging.getLogger(__name__)

# Bumped to 2 when the Bismark .cov 1-based coordinate fix (C1) landed.
# Stores written under manifest_version < 2 may carry +1-shifted positions
# for real (1-based) Bismark .cov input, so they fail _can_reuse_sample and
# are rebuilt under the corrected, coordinate-aware converter.
_MANIFEST_VERSION = 2

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
    coordinate_base: str = "auto"            # requested convention
    resolved_coordinate_base: str = "zero_based"  # convention actually applied


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
        "manifest_version": _MANIFEST_VERSION,
        "sample_name": manifest.sample_name,
        "source": manifest.source,
        "chroms": manifest.chroms,
        "row_group_size": manifest.row_group_size,
        "format": manifest.format,
        "coordinate_base": manifest.coordinate_base,
        "resolved_coordinate_base": manifest.resolved_coordinate_base,
    }


def _can_reuse_sample(
    input_path: Path, sample_dir: Path, row_group_size: int,
    format: str = "bismark",
    coordinate_base: str = "auto",
) -> bool:
    manifest = _load_json(_manifest_path(sample_dir))
    if not manifest:
        return False
    # Reject stores written before the coordinate fix (C1): they may carry
    # +1-shifted positions for real 1-based Bismark .cov, so rebuild them.
    if manifest.get("manifest_version") != _MANIFEST_VERSION:
        return False
    if manifest.get("coordinate_base", "auto") != coordinate_base:
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


def _merge_cpg_pairs_by_position(df: pl.DataFrame) -> pl.DataFrame:
    """Merge CpG dyad rows by position when strand labels are unavailable.

    A two-strand Bismark ``.cov`` emits two rows per CpG dinucleotide — the
    + strand cytosine at position P and the − strand cytosine at position P+1.
    When no reference FASTA was supplied, all rows carry ``strand='*'`` and
    ``_merge_cpg_pairs`` cannot be used (it needs explicit +/- labels).

    This function recovers the dyad structure from position alone, operating
    **per chromosome** so positions are only compared within a chromosome:

    1. Sort positions.
    2. Greedy left-to-right scan: take the lowest unconsumed position ``p``.
       If ``p+1`` is present in the chromosome's position set, pair them and
       consume both; otherwise emit ``p`` alone.
    3. For a merged pair: sum ``N_meth``, ``N_unmeth``, ``coverage``; keep the
       lower position as the site coordinate; set ``strand='*'`` (genuinely
       unknown without a reference).

    This is **heuristic** — a leading unpaired − strand site (where the +
    strand row is absent or zero-coverage) can cause a mis-pair.  A one-time
    warning advises users to pass ``reference_fasta=`` for guaranteed accuracy.

    Parameters
    ----------
    df : pl.DataFrame
        Full-genome frame with columns
        ``["chrom","pos","strand","context","N_meth","N_unmeth","coverage","sample"]``.
        All rows should carry ``strand='*'`` (the no-reference path).

    Returns
    -------
    pl.DataFrame
        Same column set and sort order (sorted by chrom then pos within each
        chromosome partition; final frame sorted by the original chrom order
        then pos).
    """
    if df.is_empty():
        return df

    import numpy as np

    chrom_col   = df["chrom"].to_numpy()
    pos_col     = df["pos"].to_numpy()
    nm_col      = df["N_meth"].to_numpy()
    nu_col      = df["N_unmeth"].to_numpy()
    cov_col     = df["coverage"].to_numpy()
    strand_col  = df["strand"].to_numpy()
    context_col = df["context"].to_numpy()
    sample_col  = df["sample"].to_numpy()

    out_chrom:   list = []
    out_pos:     list = []
    out_strand:  list = []
    out_context: list = []
    out_nm:      list = []
    out_nu:      list = []
    out_cov:     list = []
    out_sample:  list = []

    n_paired_total   = 0
    n_unpaired_total = 0

    # Process per-chromosome so positions are only compared within a chrom.
    unique_chroms, chrom_first_idx = np.unique(chrom_col, return_index=True)
    # Preserve original chromosome order (np.unique sorts lexicographically).
    chrom_order = unique_chroms[np.argsort(chrom_first_idx)]

    for chrom in chrom_order:
        mask = chrom_col == chrom
        idx  = np.where(mask)[0]
        # Sort by position within chromosome.
        order = np.argsort(pos_col[idx], kind="stable")
        idx   = idx[order]

        positions = pos_col[idx]
        pos_set   = set(positions.tolist())

        consumed = np.zeros(len(idx), dtype=bool)
        # Map position -> local index in `idx` for fast partner lookup.
        pos_to_local: dict[int, int] = {int(p): i for i, p in enumerate(positions)}

        for i, p in enumerate(positions.tolist()):
            if consumed[i]:
                continue
            p1 = p + 1
            if p1 in pos_set:
                j = pos_to_local[p1]
                if not consumed[j]:
                    # Pair (i, j)
                    gi = idx[i]
                    gj = idx[j]
                    out_chrom.append(chrom)
                    out_pos.append(p)
                    out_strand.append("*")
                    out_context.append(context_col[gi])
                    out_nm.append(int(nm_col[gi]) + int(nm_col[gj]))
                    out_nu.append(int(nu_col[gi]) + int(nu_col[gj]))
                    out_cov.append(int(cov_col[gi]) + int(cov_col[gj]))
                    out_sample.append(sample_col[gi])
                    consumed[i] = True
                    consumed[j] = True
                    n_paired_total += 1
                    continue
            # Singleton — no ±1 neighbour (or neighbour already consumed).
            gi = idx[i]
            out_chrom.append(chrom)
            out_pos.append(p)
            out_strand.append(str(strand_col[gi]))
            out_context.append(context_col[gi])
            out_nm.append(int(nm_col[gi]))
            out_nu.append(int(nu_col[gi]))
            out_cov.append(int(cov_col[gi]))
            out_sample.append(sample_col[gi])
            n_unpaired_total += 1

    logger.warning(
        "convert: merge_strands=True but no reference_fasta given — "
        "CpG dyads are merged by position (strand-free heuristic). "
        "Paired %d dyads, left %d positions unpaired. "
        "For guaranteed strand-aware merging pass reference_fasta=.",
        n_paired_total,
        n_unpaired_total,
    )

    return pl.DataFrame({
        "chrom":    out_chrom,
        "pos":      pl.Series(out_pos, dtype=pl.Int32),
        "strand":   out_strand,
        "context":  out_context,
        "N_meth":   pl.Series(out_nm,  dtype=pl.Int32),
        "N_unmeth": pl.Series(out_nu,  dtype=pl.Int32),
        "coverage": pl.Series(out_cov, dtype=pl.Int32),
        "sample":   out_sample,
    })


# Optional strand inference

def _infer_strand(df: pl.DataFrame, reference_fasta: str) -> pl.Series:
    """Infer strand from reference sequence for each CpG position.

    The lookup uses the internal **0-based ``pos``** coordinate (the
    cytosine's own position), NOT the raw input ``start``. Bismark ``.cov``
    is 1-based, so indexing the 0-based reference array by ``start`` would
    read the base one position 3' of the cytosine and systematically
    mislabel + strand CpGs as - (it would land on the forward G of the
    dinucleotide). Reading ``pos`` is correct for every source format
    because ``pos`` is normalised to 0-based at ingestion (see
    ``_resolve_coordinate_offset``).

    A cytosine measured on the + strand has reference base C at its own
    ``pos``; the cytosine measured on the - strand pairs with a forward G,
    so its forward reference base at ``pos`` is G.

    Requires pyfaidx:  pip install pyfaidx

    Parameters
    ----------
    df : pl.DataFrame
        Must contain columns: chrom (str), pos (Int32, 0-based)
    reference_fasta : str
        Path to indexed reference FASTA (.fai index must exist)

    Returns
    -------
    pl.Series (Utf8)
        "+" where the reference base at ``pos`` is C (or c),
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
    # once, index by all of its 0-based `pos` coordinates in bulk via numpy.
    # Preserves the per-row "+"/"-"/"*" mapping and the KeyError/
    # IndexError -> "*" fallback of the original Python loop.
    chrom_arr = df["chrom"].to_numpy()
    pos_arr = df["pos"].to_numpy()
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
        positions = pos_arr[mask]
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

def _resolve_coordinate_offset(
    raw_lf: pl.LazyFrame, format: str, coordinate_base: str
) -> tuple[int, str]:
    """Resolve the ``start -> pos`` offset and the convention actually applied.

    ``pos = start + offset``. Standard Bismark .cov is 1-based (``start ==
    end``) and is shifted by ``-1``; MethylDackel bedGraph is 0-based
    half-open (no shift). With ``coordinate_base="auto"`` the convention is
    detected from a sample of ``start``/``end`` rows; ``"one_based"`` /
    ``"zero_based"`` force it. (C1)
    """
    if coordinate_base == "one_based":
        return -1, "one_based"
    if coordinate_base == "zero_based":
        return 0, "zero_based"
    if coordinate_base != "auto":
        raise ValueError(
            "coordinate_base must be 'auto', 'one_based' or 'zero_based', "
            f"got {coordinate_base!r}"
        )
    # auto-detect
    if format != "bismark":
        # MethylDackel bedGraph is 0-based half-open.
        return 0, "zero_based"
    sample = raw_lf.select(["start", "end"]).head(2000).collect()
    if sample.height == 0:
        return 0, "zero_based"
    n = sample.height
    frac_eq = int((sample["start"] == sample["end"]).sum()) / n
    frac_half = int((sample["end"] == sample["start"] + 1).sum()) / n
    if frac_eq >= 0.9:
        logger.info(
            "convert: detected 1-based Bismark .cov (start == end in %.0f%% "
            "of sampled rows); shifting pos by -1.", 100 * frac_eq,
        )
        return -1, "one_based"
    if frac_half >= 0.9:
        logger.info(
            "convert: detected 0-based input (end == start + 1 in %.0f%% of "
            "sampled rows); no shift.", 100 * frac_half,
        )
        return 0, "zero_based"
    logger.warning(
        "convert: ambiguous coordinate convention (%.0f%% start==end, "
        "%.0f%% end==start+1); assuming 0-based (no shift). Pass "
        "coordinate_base='one_based' or 'zero_based' to override.",
        100 * frac_eq, 100 * frac_half,
    )
    return 0, "zero_based"


def convert_sample(
    input_path: str,
    sample_name: str,
    output_dir: str,
    row_group_size: int = 1_000_000,
    context: str = "CpG",
    reference_fasta: str | None = None,
    merge_strands: bool = True,
    format: str = "bismark",
    coordinate_base: str = "auto",
) -> str:
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
    coordinate_base : {"auto", "one_based", "zero_based"}
        Input coordinate convention. ``"auto"`` (default) detects 1-based
        Bismark .cov (``start == end``) vs 0-based bedGraph and shifts so
        ``pos`` is always 0-based. Override to force the convention. (C1)

    Returns
    -------
    str
        The coordinate convention actually applied
        ("one_based" / "zero_based").

    Output schema
    -------------
    chrom   Utf8
    pos     Int32   (0-based; 1-based Bismark .cov start is shifted -1)
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
        # This BED is genuinely 0-based half-open, so pos = start (no shift).
        resolved_base = "zero_based"
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
        raw_lf = pl.scan_csv(
            str(p),
            separator="\t",
            has_header=False,
            skip_rows=_FORMAT_SKIP_ROWS[format],
            new_columns=_COV_COLUMNS,
            schema_overrides=_COV_SCHEMA,
        )
        # Bismark .cov is 1-based (start == end); MethylDackel bedGraph is
        # 0-based. Resolve the offset so pos is always 0-based. (C1)
        offset, resolved_base = _resolve_coordinate_offset(
            raw_lf, format, coordinate_base
        )
        lf = raw_lf.with_columns(
            [
                (pl.col("N_meth") + pl.col("N_unmeth")).alias("coverage"),
                pl.lit(sample_name).alias("sample"),
                (pl.col("start") + offset).alias("pos"),
                pl.lit(context).alias("context"),
            ]
        ).select(
            ["chrom", "pos", "context", "N_meth", "N_unmeth", "coverage", "sample",
             "start"]   # start retained only to be dropped below; strand
                        # inference uses the 0-based `pos` column
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

    # CpG strand merging : merge + and - strand pairs if requested.
    # When reference_fasta is given, strand labels (+/-) are available and
    # _merge_cpg_pairs uses them directly.  Without a reference, strand is
    # unknown ("*") for every row; _merge_cpg_pairs_by_position recovers
    # dyad structure from position alone via a greedy per-chrom pairwise pass
    # (heuristic -- a one-time warning advises users to pass reference_fasta=
    # for guaranteed accuracy).
    if merge_strands and reference_fasta is not None:
        df = _merge_cpg_pairs(df)
    elif merge_strands and reference_fasta is None:
        df = _merge_cpg_pairs_by_position(df)

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

    return resolved_base


def ensure_converted_sample(
    input_path: str,
    sample_name: str,
    output_dir: str,
    row_group_size: int = 1_000_000,
    context: str = "CpG",
    reference_fasta: str | None = None,
    format: str = "bismark",
    coordinate_base: str = "auto",
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
        coordinate_base=coordinate_base,
    ):
        return False

    temp_root = output_root.parent / f".{output_root.name}.{sample_name}.tmp"
    if temp_root.exists():
        shutil.rmtree(temp_root)

    try:
        resolved_base = convert_sample(
            input_path,
            sample_name,
            str(temp_root),
            row_group_size=row_group_size,
            context=context,
            reference_fasta=reference_fasta,
            format=format,
            coordinate_base=coordinate_base,
        )
        temp_sample_dir = _sample_dir(temp_root, sample_name)
        chroms = _expected_chrom_dirs(temp_sample_dir)
        manifest = _SampleManifest(
            sample_name=sample_name,
            source=_file_signature(source_path),
            chroms=chroms,
            row_group_size=row_group_size,
            format=format,
            coordinate_base=coordinate_base,
            resolved_coordinate_base=resolved_base,
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