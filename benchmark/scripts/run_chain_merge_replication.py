"""Clean chain-merge replication of GSE263850 with paper-matched parameters.

Reference paper (Farhangdoost et al., 2024) — Methods, DMR section:

    Differential methylation analyses were performed using the R package DSS
    (DMLfit.multiFactor function with parameter smoothing = TRUE). DMRs were
    called using the function callDMR at default parameters with a minimum
    length of 50 base pairs and 3 CpG sites (delta = 0, p.threshold = 1e-5,
    minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5). Annotation
    determined using Homer (hg38). Methylation values and P-values are
    averaged within each DMR. For DMR-gene expression correlations, we only
    included genes associated within a 100 kb distance from the DMR
    (distance of 100 kb from the TSS to the middle of the DMR).

Parameter mapping (DSS -> epykit):

    DSS::DMLfit.multiFactor(smoothing=TRUE)
        -> ep.tl.dmc(test='lr', dispersion='site',
                     smoothing=True, smoothing_span_bp=500)
           (quasi-binomial LR, per-site Pearson dispersion, Gaussian
           smoothing of counts; the cheapest DSS-style configuration.)

    DSS::callDMR(p.threshold=1e-5, delta=0, minlen=50,
                 minCG=3, dis.merge=100, pct.sig=0.5)
        -> ep.tl.dmr(method='chain_merge',
                     alpha=1e-5,           # p.threshold
                     min_abs_meth_diff=0,  # delta
                     minlen_bp=50,         # minlen
                     min_cpgs=3,           # minCG
                     dis_merge_bp=100,     # dis.merge
                     pct_sig=0.5,          # pct.sig
                     min_mean_qvalue=0.05) # safety BH-q filter on DMR table

    Annotation: HOMER nearest-TSS via UCSC refGene catalog
        -> ep.tl.annotate(refgene=<refGene.txt.gz>)

    100 kb DMR-gene table (TSS to DMR midpoint <= 100 kb)
        -> built directly here from the refGene transcripts.

Outputs (in FINAL_REPORT/data/study3/chain_merge/):
    dmr_chain_merge.parquet / .csv      - 813-row-style annotated DMR table
    dmr_gene_links_100kb.csv            - every (DMR, gene) pair with
                                          |TSS - midpoint| <= 100 kb
    parameters.json                     - exact paper-mapped parameters
    summary.md                          - paper-comparable headline numbers
    run_log.txt                         - INFO-level log of the run
    README.md                           - file index

Single-shot script. Re-running with the same inputs yields the same outputs
because every step is deterministic and per-CpG counts are cached in a
fresh per-run epykit store under data/study3/chain_merge/_store/.
"""

from __future__ import annotations

import bisect
import gzip
import json
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import polars as pl


# ---- Paths ------------------------------------------------------------------

REPO_ROOT  = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
EPYKIT_SRC = Path(r"D:/Coding/Projeler/methyl_lib/epykit3/src")
RAW_DIR    = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

OUT_DIR  = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STORE_DIR    = OUT_DIR / "_store"
TMP_DIR      = OUT_DIR / "_tmp"
LOG_PATH     = OUT_DIR / "run_log.txt"
SAMPLESHEET  = RAW_DIR / "samplesheet.csv"
REFGENE      = RAW_DIR / "refseq" / "refGene.txt.gz"
CPG_ISLANDS  = RAW_DIR / "hg38_cpg_islands.bed"
PAPER_DMR    = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"


# ---- Paper-mapped parameters ------------------------------------------------

DMC_KWARGS = dict(
    test="lr",
    dispersion="site",
    smoothing=True,
    smoothing_span_bp=500,
)

DMR_KWARGS = dict(
    method="chain_merge",
    alpha=1e-5,            # DSS p.threshold
    min_abs_meth_diff=0.0, # DSS delta
    minlen_bp=50,          # DSS minlen
    min_cpgs=3,            # DSS minCG
    dis_merge_bp=100,      # DSS dis.merge
    pct_sig=0.5,           # DSS pct.sig
    min_mean_qvalue=0.05,  # post-hoc BH-q filter on the DMR table
)

GENE_LINK_MAX_BP = 100_000  # paper: 100 kb from TSS to DMR midpoint
COVERAGE_MIN = 5            # paper: keep CpGs with >= 5x coverage


# ---- Logging ----------------------------------------------------------------

LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
log = logging.getLogger("chain_merge_replication")


def _configure_logging(mode: str = "w") -> None:
    """Configure root logging. Called by main(); not at import time, so
    importing this module from a helper does not truncate the log file."""
    handlers = [
        logging.FileHandler(LOG_PATH, mode=mode, encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=handlers,
        force=True,
    )


# ---- Step 1: ingest + filter + unite ---------------------------------------

def run_methylation_pipeline():
    """Returns a MethylData object with DMC table populated."""
    sys.path.insert(0, str(EPYKIT_SRC))
    import epykit as ep

    ep.set_tmp_dir(str(TMP_DIR))

    log.info("[1/5] Reading combined-strand BEDs from %s", SAMPLESHEET)
    md = ep.read_combined_strand_bed(
        str(SAMPLESHEET),
        treatment_group="Het_AKAP11_KO",
        control_group="WT",
        assembly="hg38",
        store_dir=str(STORE_DIR),
    )
    log.info("    %d samples loaded; store at %s", md.n_samples, md.store)
    log.info("    Treatment: %s", md.treatment_ids)
    log.info("    Control:   %s", md.control_ids)

    log.info("[2/5] Coverage filter (>= %dx)", COVERAGE_MIN)
    ep.pp.filter_coverage(md, lo_count=COVERAGE_MIN, hi_perc=99.9)

    log.info("[3/5] Intersect CpGs across samples (pp.unite type=intersect)")
    ep.pp.unite(md, type="intersect")

    log.info("[4/5] DMC with %s", DMC_KWARGS)
    t0 = time.time()
    ep.tl.dmc(md, **DMC_KWARGS)
    dt_dmc = time.time() - t0
    dmc_key = md.uns["dmc"]["last_key"]
    dmc = md.varm[dmc_key]
    n_sig = dmc.filter(pl.col("qvalue") < 0.05).height
    log.info("    DMC: %s CpGs tested; %s sig at q<0.05  (%.1fs)",
             f"{dmc.height:,}", f"{n_sig:,}", dt_dmc)

    log.info("[5/5] Chain-merge DMR with %s", DMR_KWARGS)
    t0 = time.time()
    ep.tl.dmr(md, **DMR_KWARGS)
    dt_dmr = time.time() - t0
    dmr_df = md.uns["dmr"]
    log.info("    DMR: %s candidates after BH q<0.05  (%.1fs)",
             f"{dmr_df.height:,}", dt_dmr)

    log.info("[anno] HOMER (refGene) annotation via ep.tl.annotate")
    ep.tl.annotate(
        md,
        refgene=str(REFGENE),
        cpg_islands=str(CPG_ISLANDS) if CPG_ISLANDS.exists() else None,
    )
    annotated = md.uns["dmr"].clone()
    log.info("    Annotated DMR columns: %s", annotated.columns)

    return md, annotated, dict(
        n_cpgs_tested=int(dmc.height),
        n_dmc_sig=int(n_sig),
        n_dmr=int(annotated.height),
        dmc_seconds=round(dt_dmc, 1),
        dmr_seconds=round(dt_dmr, 1),
    )


# ---- Step 2: 100 kb DMR-gene linkage (TSS -> midpoint) ----------------------

def parse_refgene_tss(refgene_path: Path) -> pd.DataFrame:
    """Return one row per (transcript, gene) with chrom / strand / tss.

    HOMER's 100 kb DMR-gene rule is *gene*-level — we collapse transcripts
    to one canonical TSS per gene by taking the most upstream TSS on each
    strand (matches HOMER's behaviour when there are multiple isoforms).
    """
    rows = []
    with gzip.open(refgene_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            chrom = p[2]
            if not chrom.startswith("chr") or "_" in chrom:
                continue
            strand   = p[3]
            tx_start = int(p[4])
            tx_end   = int(p[5])
            gene     = p[12]
            tss = tx_start if strand == "+" else tx_end
            rows.append((chrom, strand, tss, gene))
    df = pd.DataFrame(rows, columns=["chrom", "strand", "tss", "gene"])
    # Collapse to one TSS per gene per chromosome (most upstream on strand)
    # to match HOMER's nearest-TSS-of-gene semantics.
    plus = (df[df.strand == "+"]
            .groupby(["chrom", "gene"], as_index=False)["tss"].min())
    minus = (df[df.strand == "-"]
             .groupby(["chrom", "gene"], as_index=False)["tss"].max())
    plus["strand"] = "+"
    minus["strand"] = "-"
    canonical = pd.concat([plus, minus], ignore_index=True)
    return canonical


def build_100kb_gene_links(
    dmr_df: pl.DataFrame,
    gene_tss: pd.DataFrame,
    max_bp: int = GENE_LINK_MAX_BP,
) -> pd.DataFrame:
    """For every DMR, list every gene whose TSS is within `max_bp` of the
    DMR midpoint. Matches the paper's "100 kb from TSS to middle of DMR".

    Returns long-form table: one row per (DMR, gene) pair.
    """
    by_chrom: dict[str, tuple[list[int], list[str]]] = {}
    for ch, grp in gene_tss.groupby("chrom"):
        s = grp.sort_values("tss")
        by_chrom[str(ch)] = (s["tss"].astype(int).tolist(),
                             s["gene"].astype(str).tolist())

    pdf = dmr_df.to_pandas()
    pdf["chrom"] = pdf["chrom"].astype(str)
    pdf["midpoint"] = ((pdf["start"].astype(int) +
                        pdf["end"].astype(int)) // 2).astype(int)

    out_rows = []
    for r in pdf.itertuples(index=True):
        chrom = r.chrom
        if chrom not in by_chrom:
            continue
        positions, genes = by_chrom[chrom]
        lo = bisect.bisect_left(positions, r.midpoint - max_bp)
        hi = bisect.bisect_right(positions, r.midpoint + max_bp)
        seen = set()
        for j in range(lo, hi):
            gene = genes[j]
            if gene in seen:
                continue
            seen.add(gene)
            dist = positions[j] - r.midpoint  # signed: + means TSS downstream
            out_rows.append(dict(
                dmr_index=r.Index,
                chrom=chrom,
                dmr_start=int(r.start),
                dmr_end=int(r.end),
                dmr_midpoint=int(r.midpoint),
                dmr_type=getattr(r, "dmr_type", ""),
                mean_meth_diff=float(getattr(r, "mean_meth_diff", float("nan"))),
                gene=gene,
                tss=int(positions[j]),
                distance_tss_to_midpoint=int(dist),
                abs_distance=int(abs(dist)),
            ))
    out = pd.DataFrame(out_rows).sort_values(
        ["dmr_index", "abs_distance"]
    ).reset_index(drop=True)
    return out


# ---- Step 3: comparison to paper (headline numbers) ------------------------

def _interval_index(df: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    by_ch: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for _, r in df.iterrows():
        by_ch[str(r["chrom"])].append((int(r["start"]), int(r["end"])))
    for ch in by_ch:
        by_ch[ch].sort()
    return by_ch


def _any_overlap(intervals, s, e) -> bool:
    for ps, pe in intervals:
        if pe < s:
            continue
        if ps > e:
            break
        return True
    return False


def compare_to_paper(dmr_df: pl.DataFrame) -> dict:
    if not PAPER_DMR.exists():
        log.warning("Paper DMR table missing: %s — skipping comparison", PAPER_DMR)
        return {"paper_table_available": False}

    paper = pd.read_excel(PAPER_DMR, sheet_name=0)
    paper_norm = paper[["chr", "start", "end"]].rename(columns={"chr": "chrom"}).copy()
    paper_norm["chrom"] = paper_norm["chrom"].astype(str)
    paper_norm["start"] = paper_norm["start"].astype(int)
    paper_norm["end"]   = paper_norm["end"].astype(int)

    ours = dmr_df.select(["chrom", "start", "end"]).to_pandas()
    ours["chrom"] = ours["chrom"].astype(str)
    paper_by_ch = _interval_index(paper_norm)
    ours_by_ch  = _interval_index(ours)

    n_ours = len(ours)
    n_paper = len(paper_norm)

    ours_overlap = sum(
        _any_overlap(paper_by_ch.get(str(r["chrom"]), []),
                     int(r["start"]), int(r["end"]))
        for _, r in ours.iterrows()
    )
    paper_overlap = sum(
        _any_overlap(ours_by_ch.get(str(r["chrom"]), []),
                     int(r["start"]), int(r["end"]))
        for _, r in paper_norm.iterrows()
    )

    n_hyper = int((paper["diff.meth_mean"] > 0).sum())
    n_hypo  = int((paper["diff.meth_mean"] < 0).sum())
    median_len = int((paper["end"] - paper["start"]).median())

    return dict(
        paper_table_available=True,
        paper_n_dmr=n_paper,
        paper_n_hyper=n_hyper,
        paper_n_hypo=n_hypo,
        paper_median_length_bp=median_len,
        our_n_dmr=n_ours,
        coord_overlap_count=int(ours_overlap),
        coord_recall_of_paper=round(paper_overlap / max(n_paper, 1), 4),
        coord_precision=round(ours_overlap / max(n_ours, 1), 4),
    )


# ---- Step 4: write structured outputs --------------------------------------

def write_outputs(annotated: pl.DataFrame, links: pd.DataFrame,
                  pipeline_stats: dict, paper_stats: dict) -> None:

    # DMR table (parquet + csv)
    dmr_parquet = OUT_DIR / "dmr_chain_merge.parquet"
    dmr_csv     = OUT_DIR / "dmr_chain_merge.csv"
    annotated.write_parquet(str(dmr_parquet))
    # Flatten any list-typed annotation columns (multi_annotation=True
    # gives all_overlapping_genes / all_overlapping_features as List[str])
    # so CSV write works.
    flat = annotated
    for col, dtype in zip(annotated.columns, annotated.dtypes):
        if dtype == pl.List(pl.Utf8):
            flat = flat.with_columns(
                pl.col(col).list.join(";").alias(col)
            )
    flat.write_csv(str(dmr_csv))
    log.info("Wrote %s (%d rows, %d cols)",
             dmr_parquet.name, annotated.height, len(annotated.columns))

    # 100 kb DMR-gene links (csv only; parquet would be redundant for a flat table)
    links_csv = OUT_DIR / "dmr_gene_links_100kb.csv"
    links.to_csv(links_csv, index=False)
    log.info("Wrote %s (%d rows; %d unique genes)",
             links_csv.name, len(links), links["gene"].nunique())

    # Parameters
    params = {
        "paper_reference": (
            "Farhangdoost et al., 2024. AKAP11 heterozygous LOF in human "
            "iPSC-derived cortical neurons (GSE263850), hg38."
        ),
        "dataset": "GSE263850 (n=6: 3 Het-AKAP11-KO + 3 WT)",
        "ingest": {
            "samplesheet": str(SAMPLESHEET),
            "assembly": "hg38",
            "coverage_min": COVERAGE_MIN,
            "unite_mode": "intersect",
        },
        "dmc": {
            "epykit_call": "ep.tl.dmc(...)",
            "kwargs": DMC_KWARGS,
            "dss_equivalent": "DSS::DMLfit.multiFactor(smoothing=TRUE)",
        },
        "dmr": {
            "epykit_call": "ep.tl.dmr(...)",
            "kwargs": DMR_KWARGS,
            "dss_equivalent": (
                "DSS::callDMR(p.threshold=1e-5, delta=0, minlen=50, "
                "minCG=3, dis.merge=100, pct.sig=0.5)"
            ),
        },
        "annotation": {
            "tool": "ep.tl.annotate(refgene=...)",
            "catalog": str(REFGENE),
            "paper_tool": "HOMER annotatePeaks.pl (hg38, RefSeq)",
        },
        "gene_linkage_100kb": {
            "rule": (
                "Include every gene whose canonical TSS lies within 100 kb "
                "of the DMR midpoint (paper's exact rule). Multiple genes "
                "per DMR allowed; one row per (DMR, gene) pair."
            ),
            "max_bp": GENE_LINK_MAX_BP,
            "tss_source": "UCSC refGene.txt.gz, one canonical TSS per "
                          "gene per chromosome (most-upstream on strand).",
        },
        "outputs": [
            "dmr_chain_merge.parquet",
            "dmr_chain_merge.csv",
            "dmr_gene_links_100kb.csv",
            "summary.md",
            "parameters.json",
            "run_log.txt",
            "README.md",
        ],
    }
    (OUT_DIR / "parameters.json").write_text(
        json.dumps(params, indent=2), encoding="utf-8"
    )
    log.info("Wrote parameters.json")

    # Summary markdown
    lines = []
    lines.append("# Chain-merge replication of GSE263850 (paper-matched parameters)\n")
    lines.append("Per-CpG: `ep.tl.dmc(test='lr', dispersion='site', "
                 "smoothing=True, smoothing_span_bp=500)`\n\n"
                 "DMR aggregation: `ep.tl.dmr(method='chain_merge', "
                 "alpha=1e-5, min_abs_meth_diff=0, minlen_bp=50, "
                 "min_cpgs=3, dis_merge_bp=100, pct_sig=0.5)`\n\n"
                 "Annotation: HOMER-equivalent (UCSC refGene) via "
                 "`ep.tl.annotate(refgene=...)`\n")
    lines.append("\n## Pipeline\n")
    lines.append(f"- CpGs tested: **{pipeline_stats['n_cpgs_tested']:,}**")
    lines.append(f"- Significant DMCs (q < 0.05): **{pipeline_stats['n_dmc_sig']:,}**")
    lines.append(f"- DMRs called (chain_merge, alpha=1e-5, BH q < 0.05): "
                 f"**{pipeline_stats['n_dmr']:,}**")
    lines.append(f"- DMC runtime: {pipeline_stats['dmc_seconds']} s")
    lines.append(f"- DMR runtime: {pipeline_stats['dmr_seconds']} s")
    lines.append("")

    # DMR morphology
    dmr_pdf = annotated.to_pandas()
    lengths = (dmr_pdf["end"].astype(int) - dmr_pdf["start"].astype(int))
    n_hyper_ours = int((dmr_pdf.get("dmr_type", "") == "hyper").sum())
    n_hypo_ours  = int((dmr_pdf.get("dmr_type", "") == "hypo").sum())
    lines.append("\n## DMR morphology\n")
    lines.append(f"- Hyper / hypo: **{n_hyper_ours} / {n_hypo_ours}** "
                 f"({100*n_hyper_ours/max(n_hyper_ours+n_hypo_ours,1):.1f}% hyper)")
    lines.append(f"- Median length: **{int(lengths.median())} bp** "
                 f"(IQR {int(lengths.quantile(0.25))}-{int(lengths.quantile(0.75))} bp)")
    lines.append(f"- Mean length: {lengths.mean():.0f} bp")
    lines.append(f"- Range: {int(lengths.min())}-{int(lengths.max())} bp")
    lines.append("")

    # Paper comparison
    if paper_stats.get("paper_table_available"):
        lines.append("\n## Coordinate-level comparison to paper (Supp Table 5)\n")
        lines.append(f"- Paper DMRs: **{paper_stats['paper_n_dmr']}** "
                     f"(hyper {paper_stats['paper_n_hyper']} / "
                     f"hypo {paper_stats['paper_n_hypo']}; "
                     f"median {paper_stats['paper_median_length_bp']} bp)")
        lines.append(f"- Our DMRs:   **{paper_stats['our_n_dmr']}**")
        lines.append(f"- Coord recall (>=1 bp overlap of paper DMRs by ours): "
                     f"**{paper_stats['coord_recall_of_paper']*100:.1f}%** "
                     f"({int(paper_stats['coord_recall_of_paper']*paper_stats['paper_n_dmr'])} "
                     f"/ {paper_stats['paper_n_dmr']})")
        lines.append(f"- Coord precision (our DMRs hitting any paper DMR): "
                     f"**{paper_stats['coord_precision']*100:.1f}%** "
                     f"({paper_stats['coord_overlap_count']} / "
                     f"{paper_stats['our_n_dmr']})")
        lines.append("")

    # Gene linkage
    n_links = len(links)
    n_genes = links["gene"].nunique()
    n_dmrs_with_link = links["dmr_index"].nunique()
    lines.append("\n## 100 kb DMR-gene linkage\n")
    lines.append(f"- (DMR, gene) pairs: **{n_links:,}**")
    lines.append(f"- Unique genes within 100 kb of any DMR: **{n_genes:,}**")
    lines.append(f"- DMRs with >= 1 linked gene: "
                 f"**{n_dmrs_with_link:,} / {pipeline_stats['n_dmr']:,}**")
    if n_links:
        per_dmr = links.groupby("dmr_index").size()
        lines.append(f"- Genes per DMR: median {int(per_dmr.median())}, "
                     f"mean {per_dmr.mean():.1f}, max {int(per_dmr.max())}")
    lines.append("")

    (OUT_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    log.info("Wrote summary.md")

    # README
    readme = f"""# Chain-merge replication of GSE263850 (paper-matched parameters)

This folder contains a clean, single-run reproduction of the differential
methylation analysis in Farhangdoost et al. (2024) on GSE263850, using
epykit's `dmr_chain_merge` engine with parameters mapped one-for-one to
the paper's DSS::callDMR call. See `parameters.json` for the full mapping.

## Files

| File | Description |
|---|---|
| `dmr_chain_merge.parquet` | All DMRs (chain_merge, alpha=1e-5, BH q < 0.05), HOMER-style annotated. Polars/Pandas-friendly. |
| `dmr_chain_merge.csv` | Same as above, CSV. |
| `dmr_gene_links_100kb.csv` | Long-form (DMR, gene) pairs where the gene's canonical TSS is within 100 kb of the DMR midpoint. Matches the paper's exact gene-linkage rule. |
| `parameters.json` | Exact parameters used, with DSS <-> epykit mapping. |
| `summary.md` | Headline numbers: DMR counts, morphology, paper-coord overlap, gene linkage. |
| `run_log.txt` | INFO log of the run. |
| `_store/` | Per-run epykit MethylStore (CpG counts, DMC parquet). Safe to delete; will be regenerated. |
| `_tmp/` | Scratch space; safe to delete. |

## Reproducing

```
py {Path(__file__).relative_to(REPO_ROOT)}
```

Single command. Deterministic (modulo floating-point summation order).
Total wall time: ~8-10 minutes on the local box.

## Headline numbers

See `summary.md` for the full table. Highlights:

- DMRs (chain_merge, alpha=1e-5): **{pipeline_stats['n_dmr']:,}**
- Paper DMRs (DSS, p.threshold=1e-5): **{paper_stats.get('paper_n_dmr', 'NA')}**
- Coordinate recall of paper DMRs: **{paper_stats.get('coord_recall_of_paper', 'NA') if isinstance(paper_stats.get('coord_recall_of_paper'), str) else f"{paper_stats.get('coord_recall_of_paper', 0)*100:.1f}%"}**
- Genes within 100 kb of any DMR: **{links['gene'].nunique():,}**
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")
    log.info("Wrote README.md")


# ---- Main -------------------------------------------------------------------

def main() -> None:
    _configure_logging(mode="w")
    log.info("=" * 72)
    log.info("GSE263850 chain_merge replication (paper-matched parameters)")
    log.info("Output dir: %s", OUT_DIR)
    log.info("=" * 72)

    md, annotated, pipeline_stats = run_methylation_pipeline()

    log.info("[gene-link] Building 100 kb DMR-gene table from %s", REFGENE)
    gene_tss = parse_refgene_tss(REFGENE)
    log.info("    %d (chrom, gene) TSS entries", len(gene_tss))
    links = build_100kb_gene_links(annotated, gene_tss, GENE_LINK_MAX_BP)
    log.info("    %d (DMR, gene) links built; %d unique genes",
             len(links), links["gene"].nunique())

    paper_stats = compare_to_paper(annotated)
    if paper_stats.get("paper_table_available"):
        log.info(
            "[compare] paper=%d, ours=%d, recall=%.1f%%, precision=%.1f%%",
            paper_stats["paper_n_dmr"],
            paper_stats["our_n_dmr"],
            paper_stats["coord_recall_of_paper"] * 100,
            paper_stats["coord_precision"] * 100,
        )

    write_outputs(annotated, links, pipeline_stats, paper_stats)

    log.info("Done. Outputs in %s", OUT_DIR)


if __name__ == "__main__":
    main()
