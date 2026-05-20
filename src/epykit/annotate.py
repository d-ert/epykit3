"""Genomic annotation for DMC / DMR results.

Public API: ``annotate_features`` (gene-feature overlap from a GTF) and
``annotate_cpg_islands`` (island / shore / shelf / open-sea context from a
UCSC CpG-island BED).

By default ``annotate_features`` returns a single best gene per site
(intronic-host-first priority). Pass ``multi_annotation=True`` to also
populate:

  - ``nearest_tss_gene`` / ``nearest_tss_distance`` -- HOMER-style nearest
    TSS assignment (signed distance, ``-`` strand flipped so positive is
    downstream of the TSS in transcription direction). Independent of the
    feature-overlap pick -- answers a different biological question (likely
    regulated promoter) than the intronic-host gene.
  - ``all_overlapping_genes`` / ``all_overlapping_features`` -- annotatr-style
    one-to-many: every gene whose feature interval overlaps the site, and
    every feature class that overlaps. Useful when a site sits inside one
    gene's intron while also being in another gene's promoter window --
    something the single-best pick necessarily hides.

The per-chromosome join loop bounds peak memory by the largest single
chromosome rather than the whole genome. GTFs are parsed once per process
and cached in a bounded LRU (``_GTF_CACHE``, default 2 slots; override via
``EPYKIT_GTF_CACHE_SIZE`` or :func:`set_gtf_cache_size`) keyed on the
canonical file path. Interval overlaps go through bioframe, which is pure
Python (pandas + numpy) and avoids the C-extension install pain that comes
with pyranges/ncls/sorted_nearest.
"""

from __future__ import annotations

import gc
import gzip
import logging
import os
import re
import time
import traceback
from collections import OrderedDict
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

_FEATURE_PRIORITY: dict[str, int] = {
    # HOMER-style priority order. Sites are tagged with the highest-priority
    # (lowest number) feature that overlaps them. Coarse features
    # (promoter / exon / intron / intergenic) keep their historic priorities
    # so existing tests and downstream consumers don't shift; the
    # fine-grained HOMER additions (5UTR / 3UTR / TTS / noncoding) slot in
    # between them.
    "promoter":   0,
    "5UTR":       1,
    "3UTR":       2,
    "TTS":        3,
    "exon":       4,
    "intron":     5,
    "noncoding":  6,
    "intergenic": 7,
}

# HOMER-equivalent feature catalog: pass this tuple to ``annotate_features``
# to reproduce a HOMER-style 8-category breakdown
# (matches the typical methylation-paper pie chart).
HOMER_FEATURES: tuple[str, ...] = (
    "promoter", "5UTR", "exon", "intron", "3UTR", "TTS", "noncoding",
)

_FEAT_COLS = ["Chromosome", "Start", "End", "Strand", "Feature", "gene_id", "gene_name"]

# GTF cache: bounded LRU of {path -> (genes_pd, exons_pd)}.
#
# A parsed human GTF can be ~1.5 GB resident in pandas form. The default of 2
# slots is enough to keep two genomes hot (e.g. annotating DMC + DMR against
# the same GTF, or comparing mouse + human in one session) without unbounded
# growth in long-running notebooks. Override via the ``EPYKIT_GTF_CACHE_SIZE``
# env var or :func:`set_gtf_cache_size`.
_GTF_CACHE_MAX_SIZE: int = max(1, int(os.environ.get("EPYKIT_GTF_CACHE_SIZE", "2")))
_GTF_CACHE: "OrderedDict[str, tuple[Any, Any]]" = OrderedDict()

# Built-feature-index cache: bounded LRU of {key -> (features_by_chrom, tss_series, strand_lut)}.
#
# ``annotate_features`` calls Steps 2-5 (dedup exons, build feature intervals,
# group by chromosome, build TSS / strand lookups) on every invocation even
# when the GTF was already parsed via _GTF_CACHE. On a human GENCODE GTF the
# rebuild is ~5-6 s, which is meaningful when annotate is called twice in a
# row (once on the DMC table, once on the DMR table -- see ``ep.tl.annotate``).
# Caching the built index lets the second call skip straight to the per-site
# overlap loop. Key includes promoter window + feature tuple so changing
# either invalidates the cache automatically.
_BUILT_FEATURES_CACHE_MAX_SIZE: int = max(
    1, int(os.environ.get("EPYKIT_BUILT_FEATURES_CACHE_SIZE", "2"))
)
# Bundle: (features_by_chrom, tss_series, strand_lut, tss_by_chrom)
_BUILT_FEATURES_CACHE: "OrderedDict[tuple, tuple[Any, Any, Any, Any]]" = OrderedDict()


def set_gtf_cache_size(max_size: int) -> None:
    """Set the maximum number of parsed GTFs held in memory.

    The cache is keyed by canonical file path; one slot per distinct GTF.
    Decreasing the size evicts the least-recently-used entries immediately.
    """
    global _GTF_CACHE_MAX_SIZE
    if max_size < 1:
        raise ValueError(f"max_size must be >= 1, got {max_size}")
    _GTF_CACHE_MAX_SIZE = int(max_size)
    while len(_GTF_CACHE) > _GTF_CACHE_MAX_SIZE:
        evicted, _ = _GTF_CACHE.popitem(last=False)
        logger.debug("[annotate] GTF cache evicted (resize): %s", evicted)


def _gtf_cache_get(key: str) -> tuple[Any, Any] | None:
    val = _GTF_CACHE.get(key)
    if val is not None:
        _GTF_CACHE.move_to_end(key)
    return val


def _gtf_cache_put(key: str, value: tuple[Any, Any]) -> None:
    if key in _GTF_CACHE:
        _GTF_CACHE.move_to_end(key)
        _GTF_CACHE[key] = value
        return
    _GTF_CACHE[key] = value
    while len(_GTF_CACHE) > _GTF_CACHE_MAX_SIZE:
        evicted, _ = _GTF_CACHE.popitem(last=False)
        logger.debug("[annotate] GTF cache evicted (LRU): %s", evicted)


def _built_features_cache_get(key: tuple) -> tuple[Any, Any, Any, Any] | None:
    val = _BUILT_FEATURES_CACHE.get(key)
    if val is not None:
        _BUILT_FEATURES_CACHE.move_to_end(key)
    return val


def _built_features_cache_put(key: tuple, value: tuple[Any, Any, Any, Any]) -> None:
    if key in _BUILT_FEATURES_CACHE:
        _BUILT_FEATURES_CACHE.move_to_end(key)
        _BUILT_FEATURES_CACHE[key] = value
        return
    _BUILT_FEATURES_CACHE[key] = value
    while len(_BUILT_FEATURES_CACHE) > _BUILT_FEATURES_CACHE_MAX_SIZE:
        evicted, _ = _BUILT_FEATURES_CACHE.popitem(last=False)
        logger.debug("[annotate] built-features cache evicted (LRU): %s", evicted[0])


def _log(msg: str) -> None:
    """Debug-level annotation log (silent at default INFO)."""
    logger.debug("[annotate] %s", msg)


def _df_info(name: str, df) -> str:
    try:
        return f"{name}: {len(df):,} rows"
    except Exception:
        return f"{name}: (unknown shape)"


def _sites_to_df(sites: pl.DataFrame) -> "pd.DataFrame":
    """Build a pandas DataFrame of single-base (or explicit-range) site intervals
    using the Capitalized column convention (``Chromosome``/``Start``/``End``)
    shared with the GTF-derived feature DataFrames downstream.
    """
    import pandas as pd
    if "start" in sites.columns and "end" in sites.columns:
        return pd.DataFrame({
            "Chromosome": sites["chrom"].to_list(),
            "Start":      sites["start"].to_list(),
            "End":        sites["end"].to_list(),
        })
    pos = sites["pos"].to_list()
    return pd.DataFrame({
        "Chromosome": sites["chrom"].to_list(),
        "Start":      pos,
        "End":        [p + 1 for p in pos],
    })


def _build_promoter_df(genes_pd, upstream_bp: int, downstream_bp: int) -> "pd.DataFrame":
    import pandas as pd
    plus  = genes_pd[genes_pd["Strand"] == "+"].copy()
    minus = genes_pd[genes_pd["Strand"] == "-"].copy()
    plus["End"]    = plus["Start"] + downstream_bp
    plus["Start"]  = (plus["Start"] - upstream_bp).clip(lower=0)
    tss_minus      = minus["End"].copy()
    minus["Start"] = (tss_minus - downstream_bp).clip(lower=0)
    minus["End"]   = tss_minus + upstream_bp
    combined = pd.concat([plus, minus], ignore_index=True)
    combined["Feature"] = "promoter"
    return combined


def _build_intron_df(exons_pd, genes_pd) -> "pd.DataFrame":
    import pandas as pd
    if len(exons_pd) == 0 or len(genes_pd) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    gene_meta = (
        genes_pd
        .drop_duplicates("gene_id")
        .set_index("gene_id")[["Chromosome", "Start", "End", "Strand", "gene_name"]]
        .rename(columns={"Chromosome": "_g_chrom", "Start": "_g_start", "End": "_g_end"})
    )
    ex = exons_pd[["gene_id", "Start", "End"]].join(gene_meta, on="gene_id", how="inner").copy()
    ex["Start"] = ex[["Start", "_g_start"]].max(axis=1).astype(np.int64)
    ex["End"]   = ex[["End",   "_g_end"]  ].min(axis=1).astype(np.int64)
    ex = ex[ex["Start"] < ex["End"]].copy()
    if len(ex) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    ex = ex.sort_values(["gene_id", "Start"]).reset_index(drop=True)
    ex["_prev_end"] = ex.groupby("gene_id", sort=False)["End"].shift(1)
    first_mask = ex["_prev_end"].isna()
    ex.loc[first_mask, "_prev_end"] = ex.loc[first_mask, "_g_start"]
    ex["_prev_end"] = ex["_prev_end"].astype(np.int64)
    introns = ex[ex["_prev_end"] < ex["Start"]].copy()
    if len(introns) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    introns["_intron_end"] = introns["Start"]
    introns["Start"]       = introns["_prev_end"]
    introns["End"]         = introns["_intron_end"]
    introns["Feature"]     = "intron"
    introns                = introns.rename(columns={"_g_chrom": "Chromosome"})
    return introns[_FEAT_COLS].reset_index(drop=True)


_UTR_PARSE_CACHE: "OrderedDict[str, Any]" = OrderedDict()


def _utr_parse_cache_get(key: str):
    if key in _UTR_PARSE_CACHE:
        _UTR_PARSE_CACHE.move_to_end(key)
        return _UTR_PARSE_CACHE[key]
    return None


def _utr_parse_cache_put(key: str, value) -> None:
    _UTR_PARSE_CACHE[key] = value
    while len(_UTR_PARSE_CACHE) > _GTF_CACHE_MAX_SIZE:
        _UTR_PARSE_CACHE.popitem(last=False)


def _parse_gtf_utrs(gtf_path: str) -> "pd.DataFrame":
    """Stream-parse the same GTF for ``five_prime_utr`` / ``three_prime_utr``
    rows, returning a single DataFrame tagged by ``Feature``.

    Kept in a separate cache so the primary :func:`_parse_gtf_streaming`
    (genes + exons only) stays cheap and 2-tuple-shaped for the callers
    that don't care about UTRs. Returns an empty frame when the GTF
    doesn't emit UTR rows (older GENCODE versions, custom assemblies);
    callers should treat absence as 'fall through to other features'.
    """
    import pandas as pd

    cache_key = "utrs::" + str(Path(gtf_path).resolve())
    cached = _utr_parse_cache_get(cache_key)
    if cached is not None:
        return cached

    rows: list[dict] = []
    attr_re = re.compile(r'(\w+)\s+"([^"]+)"')
    is_gzip = gtf_path.endswith(".gz")
    open_fn = gzip.open if is_gzip else open

    try:
        with open_fn(gtf_path, "rt") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                feature = parts[2]
                if feature not in ("five_prime_utr", "three_prime_utr"):
                    continue
                chrom = parts[0]
                start = int(parts[3]) - 1
                end = int(parts[4])
                strand = parts[6]
                attrs = {}
                for m in attr_re.finditer(parts[8]):
                    attrs[m.group(1)] = m.group(2)
                gene_id = attrs.get("gene_id", "")
                gene_name = attrs.get("gene_name", gene_id)
                rows.append({
                    "Chromosome": chrom, "Start": start, "End": end,
                    "Strand": strand, "gene_id": gene_id,
                    "gene_name": gene_name,
                    "Feature": "5UTR" if feature == "five_prime_utr" else "3UTR",
                })
    except Exception:
        _log(f"  warning: UTR parse failed for {gtf_path}; "
             "falling back to no-UTR mode")
        rows = []

    cols = ["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name", "Feature"]
    out = pd.DataFrame(rows) if rows else pd.DataFrame(columns=cols)
    _utr_parse_cache_put(cache_key, out)
    return out


def _build_tts_df(
    genes_pd,
    *,
    upstream_bp: int = 100,
    downstream_bp: int = 1000,
):
    """Build a transcription-termination-site (TES) window DataFrame.

    Analogous to :func:`_build_promoter_df` but anchored at the gene's
    3' end. HOMER's defaults are ``-100/+1000`` around the TES.
    """
    import pandas as pd
    if len(genes_pd) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    plus  = genes_pd[genes_pd["Strand"] == "+"].copy()
    minus = genes_pd[genes_pd["Strand"] == "-"].copy()
    # + strand: TES = End -> window [End-upstream, End+downstream)
    plus["Start"] = (plus["End"] - upstream_bp).clip(lower=0)
    plus["End"]   = plus["End"] + downstream_bp
    # - strand: TES = Start -> window [Start-downstream, Start+upstream)
    tts_minus = minus["Start"].copy()
    minus["Start"] = (tts_minus - downstream_bp).clip(lower=0)
    minus["End"]   = tts_minus + upstream_bp
    combined = pd.concat([plus, minus], ignore_index=True)
    combined["Feature"] = "TTS"
    if "gene_name" not in combined.columns:
        combined["gene_name"] = combined.get("gene_id", "")
    return combined[_FEAT_COLS]


def _build_noncoding_df(genes_pd):
    """Whole gene-body intervals for genes flagged as non-protein-coding.

    HOMER's "noncoding" bucket gets a hit when a site overlaps a
    transcript whose biotype is not ``protein_coding`` (lincRNA,
    antisense, miRNA, snoRNA, ...). When the source has no ``gene_type``
    column, returns an empty frame so the priority chain falls through.
    """
    import pandas as pd
    if len(genes_pd) == 0 or "gene_type" not in genes_pd.columns:
        return pd.DataFrame(columns=_FEAT_COLS)
    nc = genes_pd[genes_pd["gene_type"].fillna("") != "protein_coding"].copy()
    if len(nc) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    nc["Feature"] = "noncoding"
    if "gene_name" not in nc.columns:
        nc["gene_name"] = nc.get("gene_id", "")
    return nc[_FEAT_COLS].reset_index(drop=True)


def _build_utr_df_from_gtf_utrs(utr_pd, side: str):
    """Filter the cached GTF-derived UTR frame to one side (5UTR / 3UTR)."""
    import pandas as pd
    if len(utr_pd) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    sub = utr_pd[utr_pd["Feature"] == side].copy()
    if len(sub) == 0:
        return pd.DataFrame(columns=_FEAT_COLS)
    return sub[_FEAT_COLS].reset_index(drop=True)


def _build_utr_df_from_refgene(genes_pd, exons_pd, side: str):
    """Derive UTR intervals from refGene-style CDS coordinates.

    refGene stores ``cdsStart`` / ``cdsEnd`` on each transcript but we
    only kept tx-level Start/End in ``genes_pd``. For refGene-sourced
    UTRs the cleanest path is to fall through to "no UTRs" so the
    priority chain demotes those sites to ``exon`` / ``intron``. Users
    who need UTR resolution should annotate with a GTF.
    """
    import pandas as pd
    return pd.DataFrame(columns=_FEAT_COLS)


def _parse_gtf_streaming(gtf_path: str) -> tuple["pd.DataFrame", "pd.DataFrame"]:
    """Stream-parse a GTF file, extracting only gene and exon rows.

    Results are cached in _GTF_CACHE keyed by canonical path so that
    repeated calls (e.g. annotating DMC then DMR) pay the I/O cost once.
    """
    import pandas as pd

    cache_key = str(Path(gtf_path).resolve())
    cached = _gtf_cache_get(cache_key)
    if cached is not None:
        _log(f"  GTF cache hit for {cache_key}")
        return cached

    gene_rows: list[dict] = []
    exon_rows: list[dict] = []
    attr_re = re.compile(r'(\w+)\s+"([^"]+)"')

    is_gzip = gtf_path.endswith('.gz')
    open_fn = gzip.open if is_gzip else open

    lines_read = 0
    try:
        with open_fn(gtf_path, 'rt') as f:
            for line in f:
                lines_read += 1
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                chrom   = parts[0]
                feature = parts[2]
                if feature not in ('gene', 'exon'):
                    continue
                # GTF 1-based closed -> 0-based half-open: subtract
                # 1 from start; end already correct.
                start  = int(parts[3]) - 1
                end    = int(parts[4])
                strand = parts[6]
                attrs  = {}
                for m in attr_re.finditer(parts[8]):
                    attrs[m.group(1)] = m.group(2)
                gene_id   = attrs.get('gene_id', '')
                gene_name = attrs.get('gene_name', attrs.get('gene_id', ''))
                # GENCODE uses ``gene_type``; Ensembl uses ``gene_biotype``.
                # Accept either so the same parser handles both vendor GTFs.
                # Defaults to "" for files that omit it entirely (the
                # gene_type_filter path treats "" as "unknown" -> excluded
                # when a filter is in effect, included when no filter).
                gene_type = attrs.get('gene_type', attrs.get('gene_biotype', ''))
                row = {
                    'Chromosome': chrom, 'Start': start, 'End': end,
                    'Strand': strand, 'gene_id': gene_id, 'gene_name': gene_name,
                    'gene_type': gene_type,
                }
                if feature == 'gene':
                    gene_rows.append(row)
                else:
                    exon_rows.append(row)
    except Exception as e:
        _log(f"  ERROR parsing GTF (read {lines_read:,} lines): {e}")
        raise

    _log(f"  GTF streaming complete: {lines_read:,} lines read")
    _log(f"  Extracted {len(gene_rows):,} gene rows, {len(exon_rows):,} exon rows")

    _empty_cols = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_name', 'gene_type']
    genes_pd = pd.DataFrame(gene_rows) if gene_rows else pd.DataFrame(columns=_empty_cols)
    exons_pd = pd.DataFrame(exon_rows) if exon_rows else pd.DataFrame(columns=_empty_cols)

    result = (genes_pd, exons_pd)
    _gtf_cache_put(cache_key, result)
    return result


def _parse_refgene_streaming(refgene_path: str) -> tuple["pd.DataFrame", "pd.DataFrame"]:
    """Stream-parse a UCSC refGene.txt(.gz) file into (genes_pd, exons_pd).

    Produces the same DataFrame schema as :func:`_parse_gtf_streaming` so
    every downstream consumer (feature interval builders, TSS map, overlap
    join, nearest-TSS lookup) works unchanged.

    refGene schema (relevant cols): ``bin, name, chrom, strand, txStart,
    txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score,
    name2``. Coords are already 0-based half-open. ``name`` is the RefSeq
    accession (NM_*/NR_*); ``name2`` is the curated gene symbol. Each
    transcript becomes one ``genes_pd`` row (``gene_id`` = accession,
    ``gene_name`` = symbol) and one ``exons_pd`` row per exon. ``gene_type``
    is derived from the accession prefix: ``NM_`` -> ``protein_coding``,
    ``NR_`` -> ``non-coding``.

    Cached in the same ``_GTF_CACHE`` keyed with a ``refgene::`` prefix
    so RefSeq and GTF sources never collide.
    """
    import pandas as pd

    cache_key = "refgene::" + str(Path(refgene_path).resolve())
    cached = _gtf_cache_get(cache_key)
    if cached is not None:
        _log(f"  refGene cache hit for {cache_key}")
        return cached

    gene_rows: list[dict] = []
    exon_rows: list[dict] = []
    is_gzip = refgene_path.endswith('.gz')
    open_fn = gzip.open if is_gzip else open

    lines_read = 0
    try:
        with open_fn(refgene_path, 'rt') as f:
            for line in f:
                lines_read += 1
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 13:
                    continue
                acc      = parts[1]
                chrom    = parts[2]
                strand   = parts[3]
                tx_start = int(parts[4])   # refGene is already 0-based half-open
                tx_end   = int(parts[5])
                ex_starts = [int(x) for x in parts[9].rstrip(',').split(',') if x]
                ex_ends   = [int(x) for x in parts[10].rstrip(',').split(',') if x]
                symbol   = parts[12]
                gene_type = "protein_coding" if acc.startswith("NM_") else "non-coding"

                # One genes_pd row per transcript: TSS = txStart on + strand,
                # txEnd on - strand. Multiple transcripts of the same symbol
                # become multiple entries, which is correct: a gene's
                # alternative TSSs are real biology and the nearest-TSS
                # rule should consider each one.
                gene_rows.append({
                    'Chromosome': chrom, 'Start': tx_start, 'End': tx_end,
                    'Strand': strand, 'gene_id': acc, 'gene_name': symbol,
                    'gene_type': gene_type,
                })
                for es, ee in zip(ex_starts, ex_ends):
                    exon_rows.append({
                        'Chromosome': chrom, 'Start': es, 'End': ee,
                        'Strand': strand, 'gene_id': acc, 'gene_name': symbol,
                        'gene_type': gene_type,
                    })
    except Exception as e:
        _log(f"  ERROR parsing refGene (read {lines_read:,} lines): {e}")
        raise

    _log(f"  refGene streaming complete: {lines_read:,} transcripts read")
    _log(f"  Extracted {len(gene_rows):,} gene rows, {len(exon_rows):,} exon rows")

    _empty_cols = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'gene_name', 'gene_type']
    genes_pd = pd.DataFrame(gene_rows) if gene_rows else pd.DataFrame(columns=_empty_cols)
    exons_pd = pd.DataFrame(exon_rows) if exon_rows else pd.DataFrame(columns=_empty_cols)

    result = (genes_pd, exons_pd)
    _gtf_cache_put(cache_key, result)
    return result


def _pick_best_overlap(joined_df) -> "pd.DataFrame":
    df = joined_df.copy()
    feat_col = "Feature_b" if "Feature_b" in df.columns else "Feature"
    df["_priority"] = df[feat_col].map(_FEATURE_PRIORITY).fillna(99)
    return (
        df.sort_values("_priority")
          .groupby("_row_idx", as_index=False)
          .first()
          .drop(columns=["_priority"])
    )


def _annotate_chromosome_chunk(
    chrom: str,
    chrom_sites: pl.DataFrame,
    chrom_features_df: "pd.DataFrame",
    multi_annotation: bool = False,
) -> "pd.DataFrame":
    """Run overlap + best-pick for one chromosome. Returns pandas DataFrame.

    When ``multi_annotation`` is True, also adds two object columns:
    ``all_genes`` (list[str] per site) and ``all_features`` (list[str] per
    site). They are aggregated from every overlap row before the best-pick
    reduction, so a site that lies in one gene's intron AND another gene's
    promoter window is faithfully represented.
    """
    import bioframe
    import pandas as pd

    COLS = ("Chromosome", "Start", "End")

    chunk_n   = len(chrom_sites)
    orig_idxs = chrom_sites["_orig_idx"].to_numpy()

    result = pd.DataFrame({
        "_orig_idx":    orig_idxs,
        "gene_id":      np.full(chunk_n, "", dtype=object),
        "gene_name":    np.full(chunk_n, "", dtype=object),
        "feature_type": np.full(chunk_n, "intergenic", dtype=object),
    })
    if multi_annotation:
        # Default empty lists for sites with no feature overlaps.
        result["all_genes"]    = [[] for _ in range(chunk_n)]
        result["all_features"] = [[] for _ in range(chunk_n)]

    if chrom_features_df.empty:
        _log(f"  {chrom}: no features -> all intergenic")
        return result

    _log(f"  {chrom}: building sites DataFrame ({chunk_n:,} sites)")
    t0 = time.time()
    try:
        sites_pd = _sites_to_df(chrom_sites)
        sites_pd["_row_idx"] = np.arange(chunk_n, dtype=np.int32)
        _log(f"  {chrom}: sites DataFrame built in {time.time()-t0:.1f}s")
    except Exception:
        _log(f"  {chrom}: ERROR building sites DataFrame:\n{traceback.format_exc()}")
        return result

    _log(f"  {chrom}: features DataFrame ({len(chrom_features_df):,} features)")
    feat_df = chrom_features_df

    _log(f"  {chrom}: running overlap ...")
    t0 = time.time()
    try:
        joined = bioframe.overlap(
            sites_pd, feat_df,
            how="left",
            cols1=COLS, cols2=COLS,
            suffixes=("", "_b"),
        )
        join_rows = len(joined)
        _log(f"  {chrom}: overlap done in {time.time()-t0:.1f}s  -> {join_rows:,} rows")
    except Exception:
        _log(f"  {chrom}: ERROR during overlap:\n{traceback.format_exc()}")
        return result

    if join_rows == 0:
        del joined
        gc.collect()
        _log(f"  {chrom}: overlap returned 0 rows -> all intergenic")
        return result

    try:
        joined_df = joined
        del joined
        gc.collect()

        # ----- annotatr-style multi-annotation aggregation -----
        # Done BEFORE _pick_best_overlap reduces joined_df, so every gene /
        # feature that overlaps each site is captured (intronic-host AND
        # promoter-of-neighbour, etc.). Empty/sentinel values from the
        # left-join's no-match rows are filtered out so a no-overlap site
        # stays with the empty-list defaults seeded above.
        if multi_annotation:
            feat_col_j = "Feature_b"   if "Feature_b"   in joined_df.columns else "Feature"
            gnm_col_j  = "gene_name_b" if "gene_name_b" in joined_df.columns else "gene_name"
            valid = (
                joined_df[gnm_col_j].notna()
                & (joined_df[gnm_col_j].astype(str) != "")
                & (joined_df[gnm_col_j].astype(str) != "-1")
            )
            if valid.any():
                multi = joined_df.loc[valid, ["_row_idx", gnm_col_j, feat_col_j]].copy()
                multi[gnm_col_j]  = multi[gnm_col_j].astype(str)
                multi[feat_col_j] = multi[feat_col_j].astype(str)
                # Sorted unique per row -- stable, dedup'd, deterministic output
                gene_lists = (
                    multi.groupby("_row_idx", sort=False)[gnm_col_j]
                         .apply(lambda s: sorted(set(s)))
                )
                feat_lists = (
                    multi.groupby("_row_idx", sort=False)[feat_col_j]
                         .apply(lambda s: sorted(set(s)))
                )
                for row_idx, glist in gene_lists.items():
                    result.at[int(row_idx), "all_genes"] = glist
                for row_idx, flist in feat_lists.items():
                    result.at[int(row_idx), "all_features"] = flist

        _log(f"  {chrom}: picking best overlaps ...")
        best = _pick_best_overlap(joined_df)
        _log(f"  {chrom}: {_df_info('best', best)}")
        del joined_df
        gc.collect()

        feat_col = "Feature_b"   if "Feature_b"   in best.columns else "Feature"
        gid_col  = "gene_id_b"   if "gene_id_b"   in best.columns else "gene_id"
        gnm_col  = "gene_name_b" if "gene_name_b" in best.columns else "gene_name"
        best_slim = (
            best[["_row_idx", gid_col, gnm_col, feat_col]]
            .rename(columns={
                gid_col:  "gene_id",
                gnm_col:  "gene_name",
                feat_col: "feature_type",
            })
        )
        del best

        local_df = (
            pd.DataFrame({"_row_idx": np.arange(chunk_n, dtype=np.int32)})
            .merge(best_slim, on="_row_idx", how="left")
        )
        result["gene_id"]      = local_df["gene_id"].fillna("").astype(str).replace("-1", "").to_numpy()
        result["gene_name"]    = local_df["gene_name"].fillna("").astype(str).replace("-1", "").to_numpy()
        result["feature_type"] = local_df["feature_type"].fillna("intergenic").astype(str).replace("-1", "intergenic").to_numpy()

        n_annotated = int((result["gene_id"] != "").sum())
        _log(f"  {chrom}: {n_annotated:,}/{chunk_n:,} sites annotated")

    except Exception:
        _log(f"  {chrom}: ERROR during post-join assembly:\n{traceback.format_exc()}")

    return result


# Public API

def _build_features_index(
    annotation_path: str,
    features: tuple[str, ...],
    promoter_upstream_bp: int,
    promoter_downstream_bp: int,
    source: str = "gtf",
    gene_type_filter: tuple[str, ...] | None = None,
) -> tuple[dict[str, Any], Any, Any, dict[str, tuple[Any, Any, Any]]]:
    """Build (or fetch from cache) the per-chromosome feature index + lookups.

    Parameters
    ----------
    annotation_path : str
        Path to the gene model file. Format determined by ``source``.
    source : {"gtf", "refgene"}, default "gtf"
        Annotation file format. ``"gtf"`` parses GENCODE/Ensembl GTFs;
        ``"refgene"`` parses UCSC ``refGene.txt`` (HOMER's default catalog).
    gene_type_filter : tuple of str or None, default None
        If set, only genes whose ``gene_type`` matches one of these strings
        are used to build feature intervals and the nearest-TSS index.
        Typical: ``("protein_coding",)`` to drop lincRNAs / pseudogenes /
        novel ENSG predictions and match HOMER+RefSeq's effective behavior.

    Returns
    -------
    (features_by_chrom, tss_series, strand_lut, tss_by_chrom)
        ``tss_by_chrom[chrom] = (sorted_positions, gene_names, strands)`` --
        three parallel numpy arrays sorted by TSS position, used for
        annotatr-style nearest-TSS lookup via bisect.

    Caches the bundle in ``_BUILT_FEATURES_CACHE`` keyed on the resolved
    path, source, feature tuple, promoter window, and gene-type filter so
    different combinations don't collide.
    """
    import pandas as pd

    feature_key = tuple(sorted(set(features)))
    gtf_key = (
        tuple(sorted(set(gene_type_filter)))
        if gene_type_filter is not None else None
    )
    cache_key = (
        str(Path(annotation_path).resolve()),
        source,
        feature_key,
        int(promoter_upstream_bp),
        int(promoter_downstream_bp),
        gtf_key,
    )
    cached = _built_features_cache_get(cache_key)
    if cached is not None:
        _log(
            f"  built-features cache hit for {cache_key[0]} "
            f"(features={feature_key}, prom={promoter_upstream_bp}/{promoter_downstream_bp})"
        )
        return cached

    # ------------------------------------------------------------------
    # Step 1: Parse annotation source (uses _GTF_CACHE after first call)
    # ------------------------------------------------------------------
    if source == "gtf":
        _log("Step 1/8: stream-parsing GTF (gene and exon rows only) ...")
    elif source == "refgene":
        _log("Step 1/8: stream-parsing UCSC refGene ...")
    else:
        raise ValueError(f"Unknown source: {source!r} (expected 'gtf' or 'refgene')")

    t0 = time.time()
    try:
        if source == "gtf":
            genes_pd, exons_pd = _parse_gtf_streaming(annotation_path)
        else:
            genes_pd, exons_pd = _parse_refgene_streaming(annotation_path)
        _log(f"  parsed in {time.time()-t0:.1f}s")
        _log(f"  {_df_info('genes_pd', genes_pd)}")
        _log(f"  {_df_info('exons_pd (raw)', exons_pd)}")
        gc.collect()
        _log("  Intermediate data freed")
    except Exception:
        _log(f"FATAL: error parsing annotation source:\n{traceback.format_exc()}")
        raise

    if "gene_id" not in genes_pd.columns:
        raise ValueError("Annotation source missing 'gene_id' attribute column")
    if "gene_name" not in genes_pd.columns:
        genes_pd["gene_name"] = genes_pd["gene_id"]
    if "gene_name" not in exons_pd.columns:
        exons_pd = exons_pd.merge(
            genes_pd[["gene_id", "gene_name"]].drop_duplicates(),
            on="gene_id", how="left",
        )

    # Apply gene_type filter (if requested). Drops both gene rows and any
    # exons belonging to those genes, so downstream feature intervals and
    # the nearest-TSS index only see the kept genes.
    if gene_type_filter is not None:
        allow = set(gene_type_filter)
        n_before = len(genes_pd)
        if "gene_type" not in genes_pd.columns:
            _log(f"  WARNING: gene_type_filter={allow} requested but source "
                 f"didn't expose gene_type; falling through (no filter applied)")
        else:
            kept_gene_ids = set(genes_pd.loc[genes_pd["gene_type"].isin(allow), "gene_id"])
            genes_pd = genes_pd[genes_pd["gene_id"].isin(kept_gene_ids)].reset_index(drop=True)
            if "gene_id" in exons_pd.columns:
                exons_pd = exons_pd[exons_pd["gene_id"].isin(kept_gene_ids)].reset_index(drop=True)
            _log(f"  gene_type filter {allow}: {n_before:,} -> {len(genes_pd):,} genes "
                 f"({len(genes_pd)/max(n_before,1):.1%} retained)")

    # ------------------------------------------------------------------
    # Step 2: Deduplicate exons
    # ------------------------------------------------------------------
    _log("Step 2/8: deduplicating exons ...")
    _exon_key = ["Chromosome", "Start", "End", "Strand", "gene_id"]
    if all(c in exons_pd.columns for c in _exon_key):
        n_before = len(exons_pd)
        extra = [c for c in ["gene_name"] if c in exons_pd.columns]
        exons_pd = (
            exons_pd[_exon_key + extra]
            .drop_duplicates(subset=["Chromosome", "Start", "End", "gene_id"])
            .reset_index(drop=True)
        )
        gc.collect()
        _log(f"  exons: {n_before:,} -> {len(exons_pd):,} (removed {n_before - len(exons_pd):,} duplicates)")
    else:
        _log("  WARNING: expected exon columns not all present; skipping dedup.")

    # ------------------------------------------------------------------
    # Step 3: Build combined feature DataFrame
    # ------------------------------------------------------------------
    _log("Step 3/8: building feature intervals ...")
    feature_dfs: list[pd.DataFrame] = []

    if "promoter" in features:
        t0 = time.time()
        prom_df = _build_promoter_df(genes_pd, promoter_upstream_bp, promoter_downstream_bp)
        _log(f"  {_df_info('promoters', prom_df)}  ({time.time()-t0:.1f}s)")
        feature_dfs.append(prom_df[_FEAT_COLS])

    if "exon" in features and len(exons_pd) > 0:
        ex = exons_pd[["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]].copy()
        ex["Feature"] = "exon"
        _log(f"  {_df_info('exons (feature)', ex)}")
        feature_dfs.append(ex[_FEAT_COLS])

    if "intron" in features and len(exons_pd) > 0 and len(genes_pd) > 0:
        t0 = time.time()
        _log("  building introns (vectorised) ...")
        try:
            intron_df = _build_intron_df(exons_pd, genes_pd)
            _log(f"  {_df_info('introns', intron_df)}  ({time.time()-t0:.1f}s)")
            if len(intron_df) > 0:
                feature_dfs.append(intron_df[_FEAT_COLS])
        except Exception:
            _log(f"  ERROR building introns:\n{traceback.format_exc()}")

    # HOMER-style additions. UTRs come from a second GTF pass (cheap with
    # the dedicated UTR cache); TTS is derived from gene-body coordinates;
    # noncoding filters genes by biotype.
    if "5UTR" in features or "3UTR" in features:
        if source == "gtf":
            utr_pd = _parse_gtf_utrs(annotation_path)
            if "5UTR" in features:
                u5 = _build_utr_df_from_gtf_utrs(utr_pd, "5UTR")
                _log(f"  {_df_info('5UTRs', u5)}")
                if len(u5) > 0:
                    feature_dfs.append(u5[_FEAT_COLS])
            if "3UTR" in features:
                u3 = _build_utr_df_from_gtf_utrs(utr_pd, "3UTR")
                _log(f"  {_df_info('3UTRs', u3)}")
                if len(u3) > 0:
                    feature_dfs.append(u3[_FEAT_COLS])
        else:
            _log("  WARNING: 5UTR/3UTR requested but source is refGene; "
                 "refGene's CDS coordinates aren't carried through the parser. "
                 "Use a GTF source for UTR-level resolution.")

    if "TTS" in features and len(genes_pd) > 0:
        tts_df = _build_tts_df(genes_pd)
        _log(f"  {_df_info('TTS', tts_df)}")
        if len(tts_df) > 0:
            feature_dfs.append(tts_df[_FEAT_COLS])

    if "noncoding" in features and len(genes_pd) > 0:
        nc_df = _build_noncoding_df(genes_pd)
        _log(f"  {_df_info('noncoding', nc_df)}")
        if len(nc_df) > 0:
            feature_dfs.append(nc_df[_FEAT_COLS])

    if feature_dfs:
        t0 = time.time()
        all_features_df = (
            pd.concat(feature_dfs, ignore_index=True)
            [_FEAT_COLS]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        _log(f"  {_df_info('all_features_df', all_features_df)}  ({time.time()-t0:.1f}s)")
    else:
        _log("  WARNING: no feature DataFrames built; all sites will be intergenic")
        all_features_df = pd.DataFrame(columns=_FEAT_COLS)

    del feature_dfs
    gc.collect()

    # ------------------------------------------------------------------
    # Step 4: Group features by chromosome
    # ------------------------------------------------------------------
    _log("Step 4/8: grouping features by chromosome ...")
    features_by_chrom: dict[str, pd.DataFrame] = {}
    for chrom_name, grp in all_features_df.groupby("Chromosome", sort=False):
        features_by_chrom[str(chrom_name)] = grp.reset_index(drop=True)
    n_feat_chroms = len(features_by_chrom)
    _log(f"  features grouped across {n_feat_chroms} chromosomes")
    for c, df in sorted(features_by_chrom.items()):
        _log(f"    {c}: {len(df):,} feature intervals")
    del all_features_df
    gc.collect()

    # ------------------------------------------------------------------
    # Step 5: Build TSS map and strand lookup (both keyed by gene_id;
    # used in Step 8 for TSS-distance + sign).
    # ------------------------------------------------------------------
    _log("Step 5/8: building TSS map (per-gene_id) and per-chrom TSS arrays ...")
    _g = (
        genes_pd[["gene_id", "Chromosome", "Start", "End", "Strand", "gene_name"]]
        .drop_duplicates("gene_id")
    )
    tss_values = np.where(
        _g["Strand"].to_numpy() != "-",
        _g["Start"].to_numpy(),
        _g["End"].to_numpy(),
    ).astype(np.int64)
    tss_series = pd.Series(tss_values, index=_g["gene_id"].to_numpy(), dtype="Int64")
    strand_lut = pd.Series(
        _g["Strand"].to_numpy(),
        index=_g["gene_id"].to_numpy(),
        dtype=object,
    )

    # Per-chromosome sorted TSS arrays for bisect-based nearest-TSS lookup
    # (annotatr/HOMER-style). Each entry holds three parallel arrays already
    # sorted ascending by TSS position so a chunk can do
    # ``np.searchsorted(positions, center)`` for O(log N) lookup. Built
    # unconditionally -- the cost is negligible (a sort per chromosome) and
    # keeping it in the same cache means ``multi_annotation`` toggling
    # doesn't invalidate the bundle.
    _g_tss = _g.assign(_tss=tss_values).sort_values(["Chromosome", "_tss"])
    tss_by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom_name, grp in _g_tss.groupby("Chromosome", sort=False):
        tss_by_chrom[str(chrom_name)] = (
            grp["_tss"].to_numpy(dtype=np.int64),
            grp["gene_name"].to_numpy(dtype=object),
            grp["Strand"].to_numpy(dtype=object),
        )
    _log(f"  TSS map built: {len(tss_series):,} genes; "
         f"per-chrom TSS arrays across {len(tss_by_chrom)} chromosomes")
    del _g, _g_tss, genes_pd, exons_pd
    gc.collect()

    bundle = (features_by_chrom, tss_series, strand_lut, tss_by_chrom)
    _built_features_cache_put(cache_key, bundle)
    return bundle


def _detect_annotation_source(path: str) -> str:
    """Infer ``"gtf"`` or ``"refgene"`` from a file path.

    Accepts the canonical extensions and a fallback heuristic on the
    filename. Raises ``ValueError`` with a useful message if neither
    matches so the caller knows to pass ``source=`` explicitly.
    """
    p = path.lower()
    if p.endswith(".gtf") or p.endswith(".gtf.gz"):
        return "gtf"
    if p.endswith(".refgene.txt") or p.endswith(".refgene.txt.gz"):
        return "refgene"
    # Loose heuristic: "refgene" anywhere in the basename catches UCSC's
    # canonical "refGene.txt[.gz]" plus most user-renamed variants.
    base = Path(path).name.lower()
    if "refgene" in base:
        return "refgene"
    raise ValueError(
        f"Cannot auto-detect annotation source from path {path!r}. "
        f"Expected a .gtf[.gz] or refGene.txt[.gz] file, or pass "
        f"source='gtf'|'refgene' explicitly."
    )


def annotate_features(
    sites: pl.DataFrame,
    annotation: str,
    *,
    source: str = "auto",
    features: list[str] | tuple[str, ...] = ("promoter", "exon", "intron"),
    promoter_upstream_bp: int = 2000,
    promoter_downstream_bp: int = 200,
    multi_annotation: bool = True,
    gene_type_filter: str | list[str] | None = None,
) -> pl.DataFrame:
    """Annotate DMC / DMR sites with gene-level genomic features.

    Parameters
    ----------
    sites : polars DataFrame
        Per-site or per-region table with at least ``chrom`` and either
        ``pos`` (single-base sites) or ``start`` + ``end`` (regions).
    annotation : str
        Path to a gene-model file. Format is auto-detected from the
        extension: ``.gtf`` / ``.gtf.gz`` -> GENCODE/Ensembl GTF;
        ``refGene.txt`` / ``refGene.txt.gz`` -> UCSC RefSeq (HOMER's
        default catalog). Pass ``source=`` to override if your filename
        doesn't follow the convention.
    source : {"auto", "gtf", "refgene"}, keyword-only, default "auto"
        Format override for the annotation file. ``"auto"`` infers from
        the extension; the explicit values are an escape hatch for
        unusually-named files.
    features : sequence of str, keyword-only
        Feature classes to build overlap intervals for. Sites that don't
        hit any are reported as ``feature_type="intergenic"`` automatically
        -- "intergenic" is the fallback, not something you opt into.
        Default builds promoter / exon / intron.
    promoter_upstream_bp, promoter_downstream_bp : int, keyword-only
        Promoter window around each gene's TSS. Default ``(-2000, +200)``
        matches the conventional "core promoter" definition.
    multi_annotation : bool, keyword-only, default True
        If True (default), also adds four annotatr-style columns alongside
        the single-best pick. Set to False to get only the legacy four-
        column output (``gene_id``, ``gene_name``, ``feature_type``,
        ``distance_to_tss``) -- useful for narrow storage or when downstream
        code asserts on exact column sets. New columns when True:

          - ``nearest_tss_gene`` (Utf8) -- gene whose TSS is closest to the
            site center (HOMER's rule).
          - ``nearest_tss_distance`` (Int32) -- signed bp distance from
            site center to that TSS, flipped for ``-`` strand so positive
            is downstream in transcription direction.
          - ``all_overlapping_genes`` (List[Utf8]) -- every gene whose
            feature interval overlaps the site (one-to-many).
          - ``all_overlapping_features`` (List[Utf8]) -- every feature type
            that overlaps (e.g. ``["intron", "promoter"]``).
    gene_type_filter : str or sequence of str or None, keyword-only
        If set, only genes whose ``gene_type`` matches are used to build
        overlap intervals and the nearest-TSS index. Typical:
        ``"protein_coding"`` to drop lincRNAs / pseudogenes / novel
        predictions. Works on both GTF (via ``gene_type`` /
        ``gene_biotype`` attribute) and refGene (NM_* -> protein_coding,
        NR_* -> non-coding). A bare string is treated as a 1-element list.

    Examples
    --------
    >>> # GTF, auto-detected
    >>> annotate_features(sites, "genes.gtf.gz", gene_type_filter="protein_coding")
    >>> # UCSC RefSeq (HOMER's default catalog), auto-detected
    >>> annotate_features(sites, "refGene.txt.gz")
    >>> # Unusually-named file; explicit override
    >>> annotate_features(sites, "my_genes.tsv", source="refgene")
    """
    # Resolve source (auto-detect if needed)
    if source == "auto":
        source = _detect_annotation_source(annotation)
    elif source not in ("gtf", "refgene"):
        raise ValueError(
            f"source must be 'auto', 'gtf', or 'refgene'; got {source!r}"
        )

    # Normalize gene_type_filter to tuple or None
    if isinstance(gene_type_filter, str):
        gene_type_filter_norm: tuple[str, ...] | None = (gene_type_filter,)
    elif gene_type_filter is None:
        gene_type_filter_norm = None
    else:
        gene_type_filter_norm = tuple(gene_type_filter)
    try:
        import bioframe  # noqa: F401  (presence-check; used inside _annotate_chromosome_chunk)
    except ImportError as exc:
        raise ImportError("bioframe is required. pip install bioframe") from exc

    import pandas as pd

    _log("=" * 60)
    _log("annotate_features START")
    _log(f"  sites input: {_df_info('sites', sites)}")
    _log(f"  source: {source} -> {annotation}")
    _log(f"  features requested: {list(features)}")
    _log(f"  promoter window: -{promoter_upstream_bp} / +{promoter_downstream_bp}")
    _log(f"  multi_annotation: {multi_annotation}")
    _log(f"  gene_type_filter: {gene_type_filter_norm}")

    n = len(sites)
    t_total = time.time()

    # ------------------------------------------------------------------
    # Steps 1-5: source parse + dedup + feature intervals + groupby +
    # TSS/strand maps. All five outputs are pure functions of (source path,
    # source format, feature tuple, promoter window, gene_type filter) --
    # independent of the input ``sites`` -- so we cache the bundle in
    # _BUILT_FEATURES_CACHE keyed on all of those. A second call inside the
    # same script (e.g. annotate(md) which annotates DMC then DMR) hits the
    # cache and skips the 5-6 s of dedup / interval construction.
    # ------------------------------------------------------------------
    features_by_chrom, tss_series, strand_lut, tss_by_chrom = _build_features_index(
        annotation,
        tuple(features),
        promoter_upstream_bp,
        promoter_downstream_bp,
        source=source,
        gene_type_filter=gene_type_filter_norm,
    )

    # ------------------------------------------------------------------
    # Step 6: Tag sites with original row index
    # ------------------------------------------------------------------
    _log("Step 6/8: tagging sites with original row index ...")
    sites_with_idx = sites.with_columns(
        pl.Series("_orig_idx", np.arange(n, dtype=np.int32))
    )
    chromosomes = sorted(sites["chrom"].unique().to_list())
    _log(f"  {n:,} sites across {len(chromosomes)} chromosomes: {chromosomes}")

    # ------------------------------------------------------------------
    # Step 7: Per-chromosome annotation loop
    # ------------------------------------------------------------------
    _log("Step 7/8: per-chromosome annotation loop ...")
    annot_parts: list[pd.DataFrame] = []

    for i, chrom in enumerate(chromosomes, 1):
        chrom_sites    = sites_with_idx.filter(pl.col("chrom") == chrom)
        chunk_n        = len(chrom_sites)
        chrom_features = features_by_chrom.get(chrom, pd.DataFrame(columns=_FEAT_COLS))

        _log(f"[{i}/{len(chromosomes)}] {chrom}: {chunk_n:,} sites, "
             f"{len(chrom_features):,} features")

        if chunk_n == 0:
            _log(f"  {chrom}: 0 sites, skipping")
            continue

        t0 = time.time()
        try:
            part = _annotate_chromosome_chunk(
                chrom, chrom_sites, chrom_features,
                multi_annotation=multi_annotation,
            )
            annot_parts.append(part)
            _log(f"  {chrom}: done in {time.time()-t0:.1f}s")
        except Exception:
            _log(f"  {chrom}: UNHANDLED ERROR:\n{traceback.format_exc()}")
            part = pd.DataFrame({
                "_orig_idx":    chrom_sites["_orig_idx"].to_numpy(),
                "gene_id":      np.full(chunk_n, "", dtype=object),
                "gene_name":    np.full(chunk_n, "", dtype=object),
                "feature_type": np.full(chunk_n, "intergenic", dtype=object),
            })
            if multi_annotation:
                part["all_genes"]    = [[] for _ in range(chunk_n)]
                part["all_features"] = [[] for _ in range(chunk_n)]
            annot_parts.append(part)

        gc.collect()

    # ------------------------------------------------------------------
    # Step 8: Reassemble + TSS distance
    # ------------------------------------------------------------------
    _log("Step 8/8: reassembling results ...")
    if annot_parts:
        annot_all = (
            pd.concat(annot_parts, ignore_index=True)
            .sort_values("_orig_idx")
            .reset_index(drop=True)
        )
    else:
        _log("  WARNING: annot_parts is empty -- returning all-intergenic")
        annot_all = pd.DataFrame({
            "_orig_idx":    np.arange(n, dtype=np.int32),
            "gene_id":      np.full(n, "", dtype=object),
            "gene_name":    np.full(n, "", dtype=object),
            "feature_type": np.full(n, "intergenic", dtype=object),
        })
        if multi_annotation:
            annot_all["all_genes"]    = [[] for _ in range(n)]
            annot_all["all_features"] = [[] for _ in range(n)]

    _log(f"  {_df_info('annot_all (reassembled)', annot_all)}")

    gene_ids      = annot_all["gene_id"].to_numpy(dtype=object)
    gene_names    = annot_all["gene_name"].to_numpy(dtype=object)
    feature_types = annot_all["feature_type"].to_numpy(dtype=object)

    n_annotated = int((gene_ids != "").sum())
    ft_counts   = {k: int((feature_types == k).sum()) for k in _FEATURE_PRIORITY}
    ft_counts["intergenic"] = int((feature_types == "intergenic").sum())
    _log(f"  annotation summary: {n_annotated:,}/{n:,} sites have a gene  | {ft_counts}")

    if "pos" in sites.columns:
        site_mids = sites["pos"].to_numpy().astype(np.float64)
    else:
        site_mids = (
            (sites["start"].to_numpy() + sites["end"].to_numpy()) / 2.0
        ).astype(np.float64)

    tss_positions = (
        pd.Series(gene_ids.tolist())
        .map(tss_series)
        .to_numpy(dtype=np.float64, na_value=np.nan)
    )

    # TSS distance: positive = downstream. On - strand, TSS sits at End and a
    # higher genomic coordinate is upstream, so flip the sign. ``strand_lut``
    # came from the cached feature index, so this is O(1) on a rerun.
    strand_arr  = pd.Series(gene_ids.tolist()).map(strand_lut).to_numpy(dtype=object)
    strand_sign = np.where(strand_arr == "-", -1.0, 1.0).astype(np.float64)
    dist_to_tss              = (strand_sign * (site_mids - tss_positions)).astype(np.float32)
    dist_to_tss[gene_ids == ""] = np.nan

    # ------------------------------------------------------------------
    # Multi-annotation extras (annotatr-style): nearest-TSS lookup +
    # one-to-many gene/feature lists. Only built when requested.
    # ------------------------------------------------------------------
    multi_columns: list[pl.Series] = []
    if multi_annotation:
        _log("Step 8b/8: computing nearest-TSS (annotatr/HOMER-style) ...")
        site_chroms = sites["chrom"].to_list()
        # Reuse site_mids computed just above
        nearest_genes    = np.full(n, "", dtype=object)
        nearest_dist     = np.full(n, np.iinfo(np.int32).min, dtype=np.int64)
        for i in range(n):
            chrom = site_chroms[i]
            entry = tss_by_chrom.get(chrom)
            if entry is None:
                continue
            sorted_pos, sorted_names, sorted_strands = entry
            center = int(site_mids[i])
            # Two candidates around bisect_left, pick whichever has smaller
            # |center - tss|. Handles both edges via list indexing guards.
            idx = int(np.searchsorted(sorted_pos, center, side="left"))
            cand_idxs = []
            if idx > 0: cand_idxs.append(idx - 1)
            if idx < len(sorted_pos): cand_idxs.append(idx)
            if not cand_idxs:
                continue
            best_i = min(cand_idxs, key=lambda k: abs(int(sorted_pos[k]) - center))
            sign = -1 if sorted_strands[best_i] == "-" else 1
            nearest_genes[i] = str(sorted_names[best_i])
            nearest_dist[i]  = sign * (center - int(sorted_pos[best_i]))

        # Clip to int32 range to keep the polars dtype small; sentinel for
        # chroms missing from the GTF is "" gene + NaN distance.
        sentinel_min = np.iinfo(np.int32).min
        missing = nearest_dist == sentinel_min
        nearest_dist_clip = np.clip(nearest_dist, np.iinfo(np.int32).min + 1,
                                    np.iinfo(np.int32).max).astype(np.int32)
        # Polars Int32 doesn't carry NaN -- encode "no TSS found" as int32.min
        nearest_dist_clip[missing] = np.iinfo(np.int32).min

        multi_columns = [
            pl.Series("nearest_tss_gene",     nearest_genes.tolist(), dtype=pl.Utf8),
            pl.Series("nearest_tss_distance", nearest_dist_clip,       dtype=pl.Int32),
            pl.Series("all_overlapping_genes",
                      annot_all["all_genes"].tolist(),    dtype=pl.List(pl.Utf8)),
            pl.Series("all_overlapping_features",
                      annot_all["all_features"].tolist(), dtype=pl.List(pl.Utf8)),
        ]

    _log(f"annotate_features DONE  total elapsed {time.time()-t_total:.1f}s")
    _log("=" * 60)

    return sites.with_columns([
        pl.Series("gene_id",         gene_ids.tolist(),      dtype=pl.Utf8),
        pl.Series("gene_name",       gene_names.tolist(),    dtype=pl.Utf8),
        pl.Series("feature_type",    feature_types.tolist(), dtype=pl.Utf8),
        pl.Series("distance_to_tss", dist_to_tss,            dtype=pl.Float32),
        *multi_columns,
    ])


def annotate_cpg_islands(
    sites: pl.DataFrame,
    cpg_island_bed: str,
) -> pl.DataFrame:
    """Classify each CpG site by CpG-island context."""
    _log("=" * 60)
    _log("annotate_cpg_islands START")
    _log(f"  sites: {_df_info('sites', sites)}")
    _log(f"  BED: {cpg_island_bed}")

    if len(sites) == 0:
        _log("  sites is empty -- returning early with no cpg_context column")
        return sites

    try:
        import bioframe
    except ImportError as exc:
        raise ImportError("bioframe is required. pip install bioframe") from exc

    import pandas as pd

    COLS = ("Chromosome", "Start", "End")

    t_total = time.time()
    n = len(sites)

    _log("Step 1/3: loading BED ...")
    try:
        t0 = time.time()
        islands_df = bioframe.read_table(
            cpg_island_bed, schema="bed3", usecols=[0, 1, 2]
        ).rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
        _log(f"  BED loaded in {time.time()-t0:.1f}s: {len(islands_df):,} islands")
    except Exception:
        _log(f"FATAL: error loading BED:\n{traceback.format_exc()}")
        raise

    if len(islands_df) == 0:
        _log("  WARNING: BED is empty -> all sites open_sea")
        return sites.with_columns(pl.lit("open_sea").alias("cpg_context"))

    SHORE_DIST = 2_000
    SHELF_DIST = 4_000

    def _flanks(df: pd.DataFrame, inner: int, outer: int, label: str) -> pd.DataFrame:
        up = df[["Chromosome", "Start", "End"]].copy()
        up["End"]   = (up["Start"] - inner).clip(lower=0)
        up["Start"] = (up["Start"] - outer).clip(lower=0)
        up["_ctx"]  = label
        dn = df[["Chromosome", "Start", "End"]].copy()
        dn["Start"] = dn["End"] + inner
        dn["End"]   = dn["End"] + outer
        dn["_ctx"]  = label
        return pd.concat([up, dn], ignore_index=True)

    shore_df           = _flanks(islands_df, 0,          SHORE_DIST, "shore")
    shelf_df           = _flanks(islands_df, SHORE_DIST, SHELF_DIST, "shelf")
    islands_df["_ctx"] = "island"
    _log(f"  flanks built: {len(shore_df):,} shore intervals, {len(shelf_df):,} shelf intervals")

    _log("Step 2/3: building sites DataFrame ...")
    try:
        sites_pr_df = _sites_to_df(sites)
        sites_pr_df["_row_idx"] = np.arange(n, dtype=np.int32)
        _log(f"  sites DataFrame: {n:,} rows")
    except Exception:
        _log(f"FATAL: error building sites DataFrame:\n{traceback.format_exc()}")
        raise

    cpg_context = np.full(n, "open_sea", dtype=object)

    _log("Step 3/3: overlapping shelf / shore / island ...")
    for ctx_df, ctx_label in [
        (shelf_df,                                           "shelf"),
        (shore_df,                                           "shore"),
        (islands_df[["Chromosome", "Start", "End", "_ctx"]], "island"),
    ]:
        t0 = time.time()
        try:
            overlap = bioframe.overlap(
                sites_pr_df, ctx_df,
                how="inner",
                cols1=COLS, cols2=COLS,
                suffixes=("", "_b"),
            )
            n_hits  = len(overlap)
            _log(f"  {ctx_label}: {n_hits:,} hits in {time.time()-t0:.1f}s")
            if n_hits == 0:
                continue
            hit_idxs = overlap["_row_idx"].drop_duplicates().to_numpy(dtype=np.int32)
            cpg_context[hit_idxs] = ctx_label
        except Exception:
            _log(f"  ERROR during {ctx_label} overlap:\n{traceback.format_exc()}")

    counts = {lbl: int((cpg_context == lbl).sum())
              for lbl in ["island", "shore", "shelf", "open_sea"]}
    _log(f"  context summary: {counts}")
    _log(f"annotate_cpg_islands DONE  total elapsed {time.time()-t_total:.1f}s")
    _log("=" * 60)

    return sites.with_columns(
        pl.Series("cpg_context", list(cpg_context), dtype=pl.Utf8)
    )