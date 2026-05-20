"""IGV / UCSC-friendly exports for a MethylData object.

Functions
---------
* ``to_bedgraph``  : single-sample beta (or coverage) -> BedGraph (chrom start end value).
* ``to_bigwig``    : same, BigWig via pyBigWig (optional dep).
* ``dmcs_to_bed``  : DMC table -> BED, optionally filtered by alpha / |Deltabeta|.
* ``dmrs_to_bed``  : DMR table -> BED.

All BedGraph / BED writers work with no extra dependencies. ``to_bigwig``
imports ``pyBigWig`` lazily and raises a friendly error if missing
(``pip install 'epykit[export]'``). There is no Windows wheel on PyPI for
pyBigWig at the time of writing; bedgraph + BED cover that path.
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl

from .methyldata import MethylData

logger = logging.getLogger(__name__)


def _chromosome_sort_key(chrom: str) -> tuple[int, str | int]:
    """Sort chromosomes the way IGV/UCSC expect: 1-22, X, Y, M, then the rest."""
    stripped = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if stripped.isdigit():
        return (0, int(stripped))
    if stripped.upper() in ("X", "Y", "M", "MT"):
        return (1, stripped.upper())
    return (2, stripped)


def _load_sample_beta(md: MethylData, sample: str, value: str) -> pl.DataFrame:
    """Read a single sample's per-CpG rows. Returns chrom, pos, value."""
    if sample not in md.obs.get_column("sample_id").to_list():
        raise ValueError(
            f"Sample {sample!r} not in md.obs.sample_id. Known samples: "
            f"{md.obs.get_column('sample_id').to_list()}"
        )

    pattern = f"{md.store}/sample={sample}/chrom=*/part-*.parquet"
    lf = pl.scan_parquet(pattern).select(
        ["chrom", "pos", "N_meth", "coverage"]
    )
    if value == "beta":
        lf = lf.filter(pl.col("coverage") > 0).with_columns(
            (pl.col("N_meth").cast(pl.Float64) / pl.col("coverage")).alias("value")
        )
    elif value == "coverage":
        lf = lf.with_columns(pl.col("coverage").cast(pl.Float64).alias("value"))
    elif value == "N_meth":
        lf = lf.with_columns(pl.col("N_meth").cast(pl.Float64).alias("value"))
    else:
        raise ValueError(
            f"value must be one of 'beta', 'coverage', 'N_meth'; got {value!r}"
        )
    return lf.select(["chrom", "pos", "value"]).collect().sort(["chrom", "pos"])


def to_bedgraph(
    md: MethylData,
    sample: str,
    output: str,
    *,
    value: str = "beta",
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
    df = _load_sample_beta(md, sample, value)
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    chroms = df["chrom"].to_list()
    starts = df["pos"].to_list()
    vals = df["value"].to_list()
    with out.open("w", encoding="utf-8") as f:
        f.write(
            f'track type=bedGraph name="{sample}_{value}" '
            f'description="{value} from epykit"\n'
        )
        for c, s, v in zip(chroms, starts, vals):
            if v is None:
                continue
            f.write(f"{c}\t{s}\t{s + 1}\t{v:.6g}\n")
    logger.info("Wrote bedgraph: %s (%d rows)", out, len(chroms))
    return str(out.resolve())


def _infer_chrom_sizes(md: MethylData) -> dict[str, int]:
    """Compute the maximum CpG position + 1 per chromosome from the store."""
    pattern = f"{md.store}/sample=*/chrom=*/part-*.parquet"
    sizes_df = (
        pl.scan_parquet(pattern)
        .group_by("chrom")
        .agg(pl.max("pos").alias("max_pos"))
        .collect()
    )
    return {
        row["chrom"]: int(row["max_pos"]) + 1
        for row in sizes_df.iter_rows(named=True)
    }


def to_bigwig(
    md: MethylData,
    sample: str,
    output: str,
    *,
    value: str = "beta",
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

    df = _load_sample_beta(md, sample, value)
    if len(df) == 0:
        raise ValueError(f"No data for sample {sample!r}")

    if chrom_sizes is None:
        chrom_sizes = _infer_chrom_sizes(md)

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    bw = pyBigWig.open(str(out), "w")
    try:
        observed_chroms = sorted(df["chrom"].unique().to_list(), key=_chromosome_sort_key)
        header = [(c, chrom_sizes[c]) for c in observed_chroms if c in chrom_sizes]
        if not header:
            raise ValueError(
                "chrom_sizes does not cover any chromosome present in the sample"
            )
        bw.addHeader(header)
        for c in observed_chroms:
            if c not in chrom_sizes:
                continue
            sub = df.filter(pl.col("chrom") == c)
            starts = sub["pos"].cast(pl.Int64).to_list()
            ends = [s + 1 for s in starts]
            vals = [
                float(v) if v is not None else 0.0
                for v in sub["value"].to_list()
            ]
            bw.addEntries([c] * len(starts), starts, ends=ends, values=vals)
    finally:
        bw.close()

    logger.info("Wrote BigWig: %s", out)
    return str(out.resolve())


def _resolve_dmc_table(
    md: MethylData, test: str | None
) -> pl.DataFrame:
    df = md.get_dmc(test=test, annotated=True)
    if df is None:
        raise ValueError(
            "No DMC results on this MethylData. Run ep.tl.dmc(md) first."
        )
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
        f.write("track name=epykit_dmcs description=\"epykit DMC calls\"\n")
        for i, row in enumerate(filt.iter_rows(named=True), start=1):
            chrom = row["chrom"]
            pos = int(row["pos"])
            q = row.get(p_col)
            md_diff = row.get("meth_diff") or 0.0
            score = 0 if q is None or q != q else int(max(0.0, min(1000.0, 1000.0 * (1.0 - q))))
            strand = "+" if md_diff > 0 else "-"
            name = f"dmc_{i}"
            f.write(
                f"{chrom}\t{pos}\t{pos + 1}\t{name}\t{score}\t{strand}\n"
            )
    logger.info("Wrote DMC BED: %s (%d rows)", out, len(filt))
    return str(out.resolve())


def _resolve_dmr_table(md: MethylData) -> pl.DataFrame:
    df = md.uns.get("dmr")
    if df is None or not isinstance(df, pl.DataFrame) or len(df) == 0:
        raise ValueError(
            "No DMR results on this MethylData. Run ep.tl.dmr(md) first."
        )
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
        "qvalue" if "qvalue" in df.columns
        else "combined_qvalue" if "combined_qvalue" in df.columns
        else "combined_pvalue" if "combined_pvalue" in df.columns
        else None
    )

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    df_sorted = df.sort(["chrom", "start"])
    with out.open("w", encoding="utf-8") as f:
        f.write("track name=epykit_dmrs description=\"epykit DMR calls\"\n")
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


__all__ = [
    "to_bedgraph",
    "to_bigwig",
    "dmcs_to_bed",
    "dmrs_to_bed",
]
