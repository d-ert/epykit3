"""IGV / UCSC-friendly exports for a MethylData object.

Functions
---------
* ``to_bedgraph``  : single-sample beta (or coverage) -> BedGraph (chrom start end value).
* ``to_bigwig``    : same, BigWig via pyBigWig (optional dep).
* ``dmcs_to_bed``  : DMC table -> BED, optionally filtered by alpha / |Deltabeta|.
* ``dmrs_to_bed``  : DMR table -> BED.
* ``dmr_to_tsv``   : DMR table -> TSV/CSV (full table, sorted by chrom/start).
* ``dmc_to_tsv``   : DMC table -> TSV/CSV (significant or full, with lr+ combined-p support).
* ``dvc_to_tsv``   : DVC table -> TSV/CSV (significant or full, filtered on q_variance).
* ``qc_to_tsv``    : per-sample QC summary (md.obs) -> TSV/CSV verbatim.

All BedGraph / BED writers work with no extra dependencies. ``to_bigwig``
imports ``pyBigWig`` lazily and raises a friendly error if missing
(``pip install 'epykit[export]'``). There is no Windows wheel on PyPI for
pyBigWig at the time of writing; bedgraph + BED cover that path.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

import polars as pl

from .methyldata import MethylData

ValueKind = Literal["beta", "coverage", "N_meth"]

logger = logging.getLogger(__name__)


def _chromosome_sort_key(chrom: str) -> tuple[int, str | int]:
    """Sort chromosomes the way IGV/UCSC expect: 1-22, X, Y, M, then the rest."""
    stripped = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if stripped.isdigit():
        return (0, int(stripped))
    if stripped.upper() in ("X", "Y", "M", "MT"):
        return (1, stripped.upper())
    return (2, stripped)


def _value_lf(lf: pl.LazyFrame, value: ValueKind) -> pl.LazyFrame:
    """Attach the requested ``value`` column to a lazy frame of CpG rows."""
    if value == "beta":
        return lf.filter(pl.col("coverage") > 0).with_columns(
            (pl.col("N_meth").cast(pl.Float64) / pl.col("coverage")).alias("value")
        )
    if value == "coverage":
        return lf.with_columns(pl.col("coverage").cast(pl.Float64).alias("value"))
    if value == "N_meth":
        return lf.with_columns(pl.col("N_meth").cast(pl.Float64).alias("value"))
    raise ValueError(f"value must be one of 'beta', 'coverage', 'N_meth'; got {value!r}")


def _validate_sample(md: MethylData, sample: str) -> None:
    if sample not in md.obs.get_column("sample_id").to_list():
        raise ValueError(
            f"Sample {sample!r} not in md.obs.sample_id. Known samples: "
            f"{md.obs.get_column('sample_id').to_list()}"
        )


def _iter_sample_chrom_value(md: MethylData, sample: str, value: ValueKind):
    """Yield ``(chrom, DataFrame[pos, value])`` one chromosome at a time.

    Streams the methylstore partition tree so peak memory stays
    O(largest chromosome) rather than materialising the whole sample
    (~22-28M CpGs) and then full-genome Python lists (M12). Chromosomes
    are yielded in IGV/UCSC order.
    """
    _validate_sample(md, sample)
    # ``value`` is validated up front so an unknown kind raises before any I/O.
    _value_lf(pl.LazyFrame({"N_meth": [], "coverage": []}), value)
    sample_dir = Path(md.store) / f"sample={sample}"
    chrom_dirs = sorted(
        sample_dir.glob("chrom=*"),
        key=lambda d: _chromosome_sort_key(d.name.removeprefix("chrom=")),
    )
    for chrom_dir in chrom_dirs:
        chrom = chrom_dir.name.removeprefix("chrom=")
        lf = pl.scan_parquet(str(chrom_dir / "part-*.parquet")).select(
            ["pos", "N_meth", "coverage"]
        )
        df = _value_lf(lf, value).select(["pos", "value"]).collect().sort("pos")
        if len(df) > 0:
            yield chrom, df


def to_bedgraph(
    md: MethylData,
    sample: str,
    output: str,
    *,
    value: ValueKind = "beta",
) -> str:
    """Write a single sample's beta (or coverage) as a 4-column BedGraph.

    Output rows: ``chrom\\tstart\\tend\\tvalue`` (1 bp interval per CpG;
    end = start + 1). 0-based half-open, matching the methylstore.

    Parameters
    ----------
    md : MethylData
    sample : str
        Sample identifier (must appear in ``md.obs.sample_id``).
    output : str
        Output file path. Suffix ``.bedgraph`` / ``.bg`` is conventional
        but not enforced.
    value : {"beta", "coverage", "N_meth"}
        Column to emit. beta is in [0, 1].

    Returns
    -------
    str
        The output path (absolute), for convenience.
    """
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_rows = 0
    with out.open("w", encoding="utf-8") as f:
        f.write(f'track type=bedGraph name="{sample}_{value}" description="{value} from epykit"\n')
        # Stream one chromosome at a time so we never hold the whole genome
        # as Python lists (M12). Peak memory is O(largest chromosome).
        for chrom, df in _iter_sample_chrom_value(md, sample, value):
            starts = df["pos"].to_list()
            vals = df["value"].to_list()
            for s, v in zip(starts, vals, strict=True):
                if v is None:
                    continue
                f.write(f"{chrom}\t{s}\t{s + 1}\t{v:.6g}\n")
                n_rows += 1
    logger.info("Wrote bedgraph: %s (%d rows)", out, n_rows)
    return str(out.resolve())


def _infer_chrom_sizes(md: MethylData) -> dict[str, int]:
    """Compute the maximum CpG position + 1 per chromosome from the store."""
    pattern = f"{md.store}/sample=*/chrom=*/part-*.parquet"
    sizes_df = (
        pl.scan_parquet(pattern).group_by("chrom").agg(pl.max("pos").alias("max_pos")).collect()
    )
    return {row["chrom"]: int(row["max_pos"]) + 1 for row in sizes_df.iter_rows(named=True)}


def to_bigwig(
    md: MethylData,
    sample: str,
    output: str,
    *,
    value: ValueKind = "beta",
    chrom_sizes: dict[str, int] | None = None,
) -> str:
    """Write a single sample's beta (or coverage) as a BigWig.

    Requires ``pyBigWig`` -- install with ``pip install 'epykit[export]'``.
    There is no Windows wheel on PyPI at the time of writing; use
    :func:`to_bedgraph` instead on Windows.

    Parameters
    ----------
    md, sample, output, value : see :func:`to_bedgraph`.
    chrom_sizes : dict[str, int], optional
        Mapping ``{chrom: length}``. If None, sizes are inferred from the
        maximum CpG position observed per chromosome in ``md.store``.
        Inference is approximate (max-pos + 1), which is sufficient for
        BigWig validity but does not extend past the last CpG.
    """
    try:
        import pyBigWig
    except ImportError as exc:
        raise ImportError(
            "pyBigWig is required for BigWig export. "
            "Install it with: pip install 'epykit[export]'  "
            "(no Windows wheel -- use to_bedgraph instead on Windows)."
        ) from exc

    _validate_sample(md, sample)
    if chrom_sizes is None:
        chrom_sizes = _infer_chrom_sizes(md)

    # Header must list every chromosome up front, but we read CpG data one
    # chromosome at a time so peak memory stays O(largest chromosome) (M12).
    # Chrom names for the header come from the partition directory listing
    # (no data scan); entries are streamed per chromosome below.
    sample_dir = Path(md.store) / f"sample={sample}"
    observed_chroms = sorted(
        (d.name.removeprefix("chrom=") for d in sample_dir.glob("chrom=*")),
        key=_chromosome_sort_key,
    )
    if not observed_chroms:
        raise ValueError(f"No data for sample {sample!r}")

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    bw = pyBigWig.open(str(out), "w")
    try:
        header = [(c, chrom_sizes[c]) for c in observed_chroms if c in chrom_sizes]
        if not header:
            raise ValueError("chrom_sizes does not cover any chromosome present in the sample")
        bw.addHeader(header)
        for c, sub in _iter_sample_chrom_value(md, sample, value):
            if c not in chrom_sizes:
                continue
            starts = sub["pos"].cast(pl.Int64).to_list()
            ends = [s + 1 for s in starts]
            vals = [float(v) if v is not None else 0.0 for v in sub["value"].to_list()]
            bw.addEntries([c] * len(starts), starts, ends=ends, values=vals)
    finally:
        bw.close()

    logger.info("Wrote BigWig: %s", out)
    return str(out.resolve())


def _resolve_dmc_table(md: MethylData, test: str | None) -> pl.DataFrame:
    df = md.get_dmc(test=test, annotated=True)
    if df is None:
        raise ValueError("No DMC results on this MethylData. Run ep.tl.dmc(md) first.")
    return df


def dmcs_to_bed(
    md: MethylData,
    output: str,
    *,
    alpha: float = 0.05,
    min_abs_diff: float = 0.0,
    test: str | None = None,
) -> str:
    """Write significant DMCs as a 6-column BED (chrom start end name score strand).

    * ``name``   = ``"dmc_<i>"`` (1-indexed by genomic order)
    * ``score``  = ``round(1000 * (1 - qvalue))`` clipped to [0, 1000]
    * ``strand`` = ``+`` for hyper, ``-`` for hypo (so IGV can colour them)

    Parameters
    ----------
    md : MethylData
    output : str
    alpha : float
        q-value (or p-value if no q-value) threshold.
    min_abs_diff : float
        Minimum |meth_diff| to include.
    test : str, optional
        Specific DMC test backend. Defaults to the most recent.
    """
    df = _resolve_dmc_table(md, test)
    p_col = "qvalue" if "qvalue" in df.columns else "pvalue"
    filt = df.filter(
        (pl.col(p_col) < alpha)
        & (pl.col("meth_diff").abs() >= min_abs_diff)
        & pl.col(p_col).is_not_nan()
    ).sort(["chrom", "pos"])

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", encoding="utf-8") as f:
        f.write('track name=epykit_dmcs description="epykit DMC calls"\n')
        for i, row in enumerate(filt.iter_rows(named=True), start=1):
            chrom = row["chrom"]
            pos = int(row["pos"])
            q = row.get(p_col)
            md_diff = row.get("meth_diff") or 0.0
            score = 0 if q is None or q != q else int(max(0.0, min(1000.0, 1000.0 * (1.0 - q))))
            strand = "+" if md_diff > 0 else "-"
            name = f"dmc_{i}"
            f.write(f"{chrom}\t{pos}\t{pos + 1}\t{name}\t{score}\t{strand}\n")
    logger.info("Wrote DMC BED: %s (%d rows)", out, len(filt))
    return str(out.resolve())


def _resolve_dmr_table(md: MethylData) -> pl.DataFrame:
    df = md.uns.get("dmr")
    if df is None or not isinstance(df, pl.DataFrame) or len(df) == 0:
        raise ValueError("No DMR results on this MethylData. Run ep.tl.dmr(md) first.")
    return df


def dmrs_to_bed(
    md: MethylData,
    output: str,
) -> str:
    """Write the DMR table as a 6-column BED.

    Works with both the tile-based DMR schema (``meth_diff``, ``qvalue``,
    ``dmr_type``) and the sliding-window schema (``mean_meth_diff``,
    ``combined_qvalue``). ``name`` is set to ``dmr_<i>``, score to the
    ``-log10`` of the q-value clipped to [0, 1000], strand to ``+`` for
    ``hyper`` and ``-`` for ``hypo`` (``.`` for ``mixed``).
    """
    df = _resolve_dmr_table(md)
    diff_col = "meth_diff" if "meth_diff" in df.columns else "mean_meth_diff"
    q_col = (
        "qvalue"
        if "qvalue" in df.columns
        else "combined_qvalue"
        if "combined_qvalue" in df.columns
        else "combined_pvalue"
        if "combined_pvalue" in df.columns
        else None
    )

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    df_sorted = df.sort(["chrom", "start"])
    with out.open("w", encoding="utf-8") as f:
        f.write('track name=epykit_dmrs description="epykit DMR calls"\n')
        for i, row in enumerate(df_sorted.iter_rows(named=True), start=1):
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            diff = row.get(diff_col) or 0.0
            if "dmr_type" in row and row.get("dmr_type"):
                t = row["dmr_type"]
                strand = "+" if t == "hyper" else "-" if t == "hypo" else "."
            else:
                strand = "+" if diff > 0 else "-"
            score = 0
            if q_col is not None:
                q = row.get(q_col)
                if q is not None and q == q and q > 0:
                    import math

                    score = int(max(0.0, min(1000.0, -100.0 * math.log10(q))))
            name = f"dmr_{i}"
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    logger.info("Wrote DMR BED: %s (%d rows)", out, len(df_sorted))
    return str(out.resolve())


def _separator_for(path: str) -> str:
    """Tab unless the path ends in .csv (case-insensitive)."""
    return "," if str(path).lower().endswith(".csv") else "\t"


def _flatten_nested_for_csv(df: pl.DataFrame) -> pl.DataFrame:
    """Stringify List/Struct columns so the CSV writer can serialise them.

    polars' ``write_csv`` raises ``CSV format does not support nested data`` on
    List or Struct columns. Annotated DMC tables carry ``List(String)`` columns
    (``all_overlapping_genes`` / ``all_overlapping_features``), so without this
    every ``*_to_tsv`` writer would fail on a post-annotation frame. List
    columns are joined with ``; ``; Struct columns are JSON-encoded.
    """
    exprs = []
    for name, dtype in df.schema.items():
        if isinstance(dtype, pl.List):
            exprs.append(pl.col(name).cast(pl.List(pl.String)).list.join("; ").alias(name))
        elif isinstance(dtype, pl.Struct):
            exprs.append(pl.col(name).struct.json_encode().alias(name))
    return df.with_columns(exprs) if exprs else df


def _write_table(df: pl.DataFrame, path: str) -> str:
    """Write `df` to `path` using the suffix-derived delimiter.

    Returns the resolved absolute path. Creates parent directories. Nested
    (List/Struct) columns are flattened to strings first so annotated frames
    survive the CSV writer.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    _flatten_nested_for_csv(df).write_csv(str(out), separator=_separator_for(path))
    logger.info("Wrote table: %s (%d rows)", out, len(df))
    return str(out.resolve())


def dmr_to_tsv(md: MethylData, path: str) -> str:
    """Write the DMR table (md.uns['dmr']) as TSV/CSV.

    Full table, sorted by (chrom, start). Delimiter is derived from the
    path suffix (.csv -> comma, otherwise tab).
    """
    df = _resolve_dmr_table(md)
    sort_cols = ["chrom", "start"] if "start" in df.columns else ["chrom"]
    return _write_table(df.sort(sort_cols), path)


def dmc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
    test: str | None = None,
) -> str:
    """Write the DMC table as TSV/CSV.

    Default: significant-only (qvalue < alpha) sorted by qvalue ascending.
    full=True: every row, sorted by (chrom, pos). When the frame carries
    `qvalue_combined` (from the lr+ neighbour-combine knob), that column
    drives the significance filter and the sort; otherwise `qvalue`.

    Delimiter is derived from the path suffix (.csv -> comma, else tab).
    """
    df = _resolve_dmc_table(md, test)

    if full:
        out_df = df.sort(["chrom", "pos"])
    else:
        # Significance gate: prefer qvalue_combined (lr+ neighbour-combine),
        # then raw qvalue, then pvalue_combined, then raw pvalue.
        if "qvalue_combined" in df.columns:
            gate_col = "qvalue_combined"
        elif "qvalue" in df.columns:
            gate_col = "qvalue"
        elif "pvalue_combined" in df.columns:
            gate_col = "pvalue_combined"
        else:
            gate_col = "pvalue"

        out_df = df.filter(
            pl.col(gate_col).is_not_null()
            & pl.col(gate_col).is_not_nan()
            & (pl.col(gate_col) < alpha)
        ).sort(gate_col)
    return _write_table(out_df, path)


def dvc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
) -> str:
    """Write the DVC table (md.varm['dvc']) as TSV/CSV.

    Default: significant-only (q_variance < alpha) sorted by q_variance
    ascending. full=True keeps every row in (chrom, pos) order.
    Delimiter is derived from the path suffix.
    """
    df = md.varm.get("dvc")
    if df is None or len(df) == 0:
        raise ValueError("No DVC results on this MethylData. Run ep.tl.dvc(md) first.")
    if full:
        sort_cols = ["chrom", "pos"] if "pos" in df.columns else ["chrom"]
        out_df = df.sort(sort_cols)
    else:
        gate = "q_variance" if "q_variance" in df.columns else "p_variance"
        out_df = df.filter(
            pl.col(gate).is_not_null() & pl.col(gate).is_not_nan() & (pl.col(gate) < alpha)
        ).sort(gate)
    return _write_table(out_df, path)


def qc_to_tsv(md: MethylData, path: str) -> str:
    """Write the per-sample QC summary (md.obs) as TSV/CSV.

    After ep.tl.qc(md), md.obs carries the per-sample metrics joined onto
    the existing samplesheet columns. This writer dumps it verbatim.
    Delimiter is derived from the path suffix.
    """
    return _write_table(md.obs, path)


def export_tables(
    md: MethylData,
    out_dir: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
    fmt: Literal["tsv", "csv"] = "tsv",
    dmc: bool = True,
    dmr: bool = True,
    dvc: bool = True,
    qc: bool = True,
) -> dict[str, str]:
    """Dump every available result table to ``out_dir`` in one call.

    A convenience over the individual ``*_to_tsv`` writers: it writes whichever
    of the DMC / DMR / DVC / QC tables are present on ``md`` and skips the rest
    silently (nothing is raised for a missing table). Files written:

    * ``<dmc_key>.significant.<ext>`` -- significant DMCs (``qvalue < alpha``).
      With ``full=True``, also ``<dmc_key>.<ext>`` (every row of the resolved
      DMC frame; after annotation that frame is the annotated subset).
    * ``dmr.<ext>``                   -- the DMR table (``md.uns['dmr']``).
    * ``dvc.significant.<ext>``       -- significant DVCs (``q_variance < alpha``).
      With ``full=True``, also ``dvc.<ext>``.
    * ``qc_summary.<ext>``            -- per-sample QC summary (``md.obs``).

    Parameters
    ----------
    md : MethylData
    out_dir : str
        Directory to write into (created if absent).
    alpha : float
        Significance threshold for the significant-only DMC / DVC tables.
    full : bool
        Also emit the full (all-rows) DMC and DVC tables alongside the
        significant ones.
    fmt : {"tsv", "csv"}
        ``"tsv"`` -> tab-delimited ``.tsv`` files; ``"csv"`` -> comma ``.csv``.
    dmc, dmr, dvc, qc : bool
        Per-table switches; set any to False to skip it.

    Returns
    -------
    dict[str, str]
        Logical table name (``"dmc_significant"``, ``"dmc_full"``, ``"dmr"``,
        ``"dvc_significant"``, ``"dvc_full"``, ``"qc"``) -> absolute path written.
        Keys are only present for tables that existed and were exported.
    """
    ext = "csv" if fmt == "csv" else "tsv"
    out = Path(out_dir)
    written: dict[str, str] = {}

    if dmc:
        dmc_key = md.uns.get("dmc", {}).get("last_key")
        if dmc_key is not None and md.get_dmc() is not None:
            written["dmc_significant"] = dmc_to_tsv(
                md, str(out / f"{dmc_key}.significant.{ext}"), alpha=alpha
            )
            if full:
                written["dmc_full"] = dmc_to_tsv(md, str(out / f"{dmc_key}.{ext}"), full=True)

    if dmr:
        dmr_df = md.uns.get("dmr")
        if isinstance(dmr_df, pl.DataFrame) and len(dmr_df):
            written["dmr"] = dmr_to_tsv(md, str(out / f"dmr.{ext}"))

    if dvc:
        dvc_df = md.varm.get("dvc")
        if isinstance(dvc_df, pl.DataFrame) and len(dvc_df):
            written["dvc_significant"] = dvc_to_tsv(
                md, str(out / f"dvc.significant.{ext}"), alpha=alpha
            )
            if full:
                written["dvc_full"] = dvc_to_tsv(md, str(out / f"dvc.{ext}"), full=True)

    if qc and md.obs is not None and len(md.obs):
        written["qc"] = qc_to_tsv(md, str(out / f"qc_summary.{ext}"))

    logger.info(
        "export_tables wrote %d table(s) to %s: %s",
        len(written),
        out,
        ", ".join(sorted(written)) or "(none)",
    )
    return written


__all__ = [
    "dmc_to_tsv",
    "dmcs_to_bed",
    "dmr_to_tsv",
    "dmrs_to_bed",
    "dvc_to_tsv",
    "export_tables",
    "qc_to_tsv",
    "to_bedgraph",
    "to_bigwig",
]
