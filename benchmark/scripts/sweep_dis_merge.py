"""dis.merge sweep on the chain_merge replication.

Re-uses the cached DMC table at FINAL_REPORT/data/study3/chain_merge/_store/
to avoid re-running the ~220 s DMC step. Sweeps dis_merge_bp across
{100, 150, 200, 250, 500} (paper used 100), holding every other parameter
identical to run_chain_merge_replication.py.

For each dis_merge value, writes:

    sweep/dis_merge_<n>/dmr.parquet            - annotated DMRs
    sweep/dis_merge_<n>/dmr_gene_links_100kb.csv
    sweep/dis_merge_<n>/comparison.json        - paper-coord overlap metrics

Plus a top-level summary:

    sweep_summary.csv / .md  - one row per dis_merge value

Output root: FINAL_REPORT/data/study3/chain_merge_dis_merge_sweep/
"""

from __future__ import annotations

import bisect
import gzip
import json
import logging
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

# ---- Paths -----------------------------------------------------------------

REPO_ROOT  = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
EPYKIT_SRC = Path(r"D:/Coding/Projeler/methyl_lib/epykit3/src")
RAW_DIR    = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

CM_DIR     = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
OUT_DIR    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" \
             / "chain_merge_dis_merge_sweep"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STORE_DIR    = CM_DIR / "_store"       # reuse chain_merge cache
TMP_DIR      = CM_DIR / "_tmp"
SAMPLESHEET  = RAW_DIR / "samplesheet.csv"
REFGENE      = RAW_DIR / "refseq" / "refGene.txt.gz"
CPG_ISLANDS  = RAW_DIR / "hg38_cpg_islands.bed"
PAPER_DMR    = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"
PANEL_E      = REPO_ROOT / "FINAL_REPORT" / "shinygo_lists" / "outputs" \
               / "reactome" / "table8.xlsx"

DIS_MERGE_GRID = [100, 150, 200, 250, 500]

# Same as run_chain_merge_replication, just dis_merge_bp varies:
DMR_KWARGS_BASE = dict(
    method="chain_merge",
    alpha=1e-5,
    min_abs_meth_diff=0.0,
    minlen_bp=50,
    min_cpgs=3,
    pct_sig=0.5,
    min_mean_qvalue=0.05,
)
DMC_KWARGS = dict(
    test="lr",
    dispersion="site",
    smoothing=True,
    smoothing_span_bp=500,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[
        logging.FileHandler(OUT_DIR / "sweep_log.txt", mode="w",
                            encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
    force=True,
)
log = logging.getLogger("sweep_dis_merge")


# ---- Pipeline helpers -------------------------------------------------------

def setup_methyldata():
    sys.path.insert(0, str(EPYKIT_SRC))
    import epykit as ep
    ep.set_tmp_dir(str(TMP_DIR))

    log.info("Re-using existing store: %s", STORE_DIR)
    md = ep.read_combined_strand_bed(
        str(SAMPLESHEET),
        treatment_group="Het_AKAP11_KO",
        control_group="WT",
        assembly="hg38",
        store_dir=str(STORE_DIR),
    )
    log.info("  %d samples loaded", md.n_samples)
    log.info("Coverage filter (>= 5x), unite intersect")
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")

    log.info("DMC (should hit cache): %s", DMC_KWARGS)
    t0 = time.time()
    ep.tl.dmc(md, **DMC_KWARGS)
    log.info("  DMC done in %.1fs (cache hit if fast)", time.time() - t0)
    return md, ep


# ---- HOMER-style nearest-TSS for 100 kb gene linkage -----------------------

def parse_refgene_tss(refgene_path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(refgene_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            chrom = p[2]
            if not chrom.startswith("chr") or "_" in chrom:
                continue
            strand = p[3]
            tx_start = int(p[4]); tx_end = int(p[5]); gene = p[12]
            tss = tx_start if strand == "+" else tx_end
            rows.append((chrom, strand, tss, gene))
    df = pd.DataFrame(rows, columns=["chrom", "strand", "tss", "gene"])
    plus  = df[df.strand == "+"].groupby(["chrom", "gene"], as_index=False)["tss"].min()
    minus = df[df.strand == "-"].groupby(["chrom", "gene"], as_index=False)["tss"].max()
    plus["strand"] = "+"; minus["strand"] = "-"
    return pd.concat([plus, minus], ignore_index=True)


def build_100kb_links(dmrs: pl.DataFrame, gene_tss: pd.DataFrame,
                      max_bp: int = 100_000) -> pd.DataFrame:
    by_chrom: dict[str, tuple[list[int], list[str]]] = {}
    for ch, grp in gene_tss.groupby("chrom"):
        s = grp.sort_values("tss")
        by_chrom[str(ch)] = (s["tss"].astype(int).tolist(),
                             s["gene"].astype(str).tolist())
    pdf = dmrs.to_pandas()
    pdf["chrom"] = pdf["chrom"].astype(str)
    pdf["midpoint"] = ((pdf["start"].astype(int) + pdf["end"].astype(int)) // 2).astype(int)

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
            g = genes[j]
            if g in seen:
                continue
            seen.add(g)
            out_rows.append(dict(
                dmr_index=r.Index, chrom=chrom,
                dmr_start=int(r.start), dmr_end=int(r.end),
                dmr_midpoint=int(r.midpoint),
                dmr_type=getattr(r, "dmr_type", ""),
                mean_meth_diff=float(getattr(r, "mean_meth_diff", float("nan"))),
                gene=g, tss=int(positions[j]),
                distance_tss_to_midpoint=int(positions[j] - r.midpoint),
                abs_distance=int(abs(positions[j] - r.midpoint)),
            ))
    return pd.DataFrame(out_rows)


# ---- Paper comparison helpers ----------------------------------------------

def _interval_index(df: pd.DataFrame, start_col="start", end_col="end",
                    idx_col=None) -> dict[str, list]:
    out: dict[str, list] = defaultdict(list)
    for i, r in df.iterrows():
        idx = i if idx_col is None else r[idx_col]
        out[str(r["chrom"])].append(
            (int(r[start_col]), int(r[end_col]), int(idx))
        )
    for ch in out:
        out[ch].sort(key=lambda t: t[0])
    return out


def _jaccard(a_s, a_e, b_s, b_e):
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def _compare(dmrs: pl.DataFrame, paper: pd.DataFrame,
             panel_genes: list[str], links_df: pd.DataFrame) -> dict:
    ours = dmrs.to_pandas().copy()
    ours["chrom"] = ours["chrom"].astype(str)
    ours["start"] = ours["start"].astype(int)
    ours["end"]   = ours["end"].astype(int)
    ours["our_index"] = np.arange(len(ours))
    ours["length"]    = ours["end"] - ours["start"]
    ours["direction"] = np.where(ours["mean_meth_diff"] > 0, "hyper",
                         np.where(ours["mean_meth_diff"] < 0, "hypo", "none"))

    paper_idx = _interval_index(paper, "start", "end", "paper_index")
    ours_idx  = _interval_index(ours,  "start", "end", "our_index")

    # Recall: each paper DMR's best Jaccard against ours
    paper_jacs = []
    paper_dir_match = []
    paper_dir_lookup = dict(zip(paper["paper_index"], paper["direction"]))
    ours_dir_lookup  = dict(zip(ours["our_index"],   ours["direction"]))
    for _, r in paper.iterrows():
        ch = str(r["chrom"]); s = int(r["start"]); e = int(r["end"])
        cands = ours_idx.get(ch, [])
        best_j = 0.0; best_i = -1
        for ts, te, ti in cands:
            if te < s: continue
            if ts > e: break
            j = _jaccard(s, e, ts, te)
            if j > best_j:
                best_j = j; best_i = ti
        paper_jacs.append(best_j)
        if best_i == -1:
            paper_dir_match.append(None)
        else:
            paper_dir_match.append(paper_dir_lookup[r["paper_index"]]
                                   == ours_dir_lookup[best_i])
    paper_jacs = np.array(paper_jacs)
    dir_total = sum(x is not None for x in paper_dir_match)
    dir_agree = sum(x is True for x in paper_dir_match)

    # Precision: each our DMR overlapping ANY paper DMR
    our_hits = 0
    for _, r in ours.iterrows():
        ch = str(r["chrom"]); s = int(r["start"]); e = int(r["end"])
        cands = paper_idx.get(ch, [])
        for ps, pe, pi in cands:
            if pe < s: continue
            if ps > e: break
            our_hits += 1; break

    # Gene-level
    paper_genes = paper["gene_u"].dropna().unique().tolist()
    paper_genes = [g for g in paper_genes if g and g != "NAN"]
    our_nearest = set(
        ours.get("nearest_tss_gene", pd.Series(dtype=str))
            .dropna().astype(str).str.upper().unique()
    ) - {""}
    our_100kb = set()
    if len(links_df):
        our_100kb = set(links_df["gene"].astype(str).str.upper().unique()) - {""}

    t5_near = len(set(paper_genes) & our_nearest)
    t5_100  = len(set(paper_genes) & our_100kb)
    pe_near = len(set(panel_genes) & our_nearest)
    pe_100  = len(set(panel_genes) & our_100kb)

    return dict(
        n_dmr=int(len(ours)),
        n_hyper=int((ours["direction"] == "hyper").sum()),
        n_hypo =int((ours["direction"] == "hypo").sum()),
        pct_hyper=round(100 * (ours["direction"] == "hyper").sum() / max(len(ours),1), 1),
        median_length_bp=int(ours["length"].median()) if len(ours) else 0,
        mean_length_bp=int(ours["length"].mean()) if len(ours) else 0,
        max_length_bp=int(ours["length"].max()) if len(ours) else 0,
        recall_anybp=round(float((paper_jacs > 0).mean()), 4),
        recall_J_0_25=round(float((paper_jacs >= 0.25).mean()), 4),
        recall_J_0_5 =round(float((paper_jacs >= 0.5).mean()), 4),
        recall_J_0_75=round(float((paper_jacs >= 0.75).mean()), 4),
        precision_anybp=round(our_hits / max(len(ours), 1), 4),
        direction_agree_n=int(dir_agree),
        direction_agree_total=int(dir_total),
        direction_agree_frac=round(dir_agree / max(dir_total, 1), 4),
        gene_recall_table5_nearest_tss=round(t5_near / max(len(paper_genes), 1), 4),
        gene_recall_table5_100kb     =round(t5_100  / max(len(paper_genes), 1), 4),
        panel_e_recall_nearest_tss   =round(pe_near / max(len(panel_genes), 1), 4),
        panel_e_recall_100kb         =round(pe_100  / max(len(panel_genes), 1), 4),
    )


# ---- Main ------------------------------------------------------------------

def main() -> None:
    log.info("=" * 72)
    log.info("dis.merge sweep on chain_merge replication")
    log.info("Grid: %s", DIS_MERGE_GRID)
    log.info("=" * 72)

    md, ep = setup_methyldata()
    gene_tss = parse_refgene_tss(REFGENE)
    log.info("refGene TSS rows: %d", len(gene_tss))

    paper = pd.read_excel(PAPER_DMR, sheet_name=0).rename(columns={"chr":"chrom"})
    paper["chrom"] = paper["chrom"].astype(str)
    paper["start"] = paper["start"].astype(int)
    paper["end"]   = paper["end"].astype(int)
    paper["paper_index"] = np.arange(len(paper))
    paper["direction"]   = np.where(paper["diff.meth_mean"] > 0, "hyper",
                            np.where(paper["diff.meth_mean"] < 0, "hypo", "none"))
    paper["gene_u"] = paper["Gene.Name"].fillna("").astype(str).str.strip().str.upper()
    log.info("paper DMRs: %d", len(paper))

    panel_e = pd.read_excel(PANEL_E, sheet_name=0)
    panel_genes = (panel_e["Gene"].astype(str).str.strip().str.upper()
                                   .dropna().unique().tolist())
    panel_genes = [g for g in panel_genes if g and g != "NAN"]
    log.info("panel-E genes: %d", len(panel_genes))

    rows = []
    for dis_merge in DIS_MERGE_GRID:
        log.info("\n--- dis_merge = %d ---", dis_merge)
        sub_dir = OUT_DIR / f"dis_merge_{dis_merge}"
        sub_dir.mkdir(parents=True, exist_ok=True)

        kwargs = {**DMR_KWARGS_BASE, "dis_merge_bp": dis_merge}
        t0 = time.time()
        ep.tl.dmr(md, **kwargs)
        dt_dmr = time.time() - t0
        log.info("  DMR (%d): %d candidates  (%.1fs)",
                 dis_merge, md.uns["dmr"].height, dt_dmr)

        ep.tl.annotate(md, refgene=str(REFGENE),
                       cpg_islands=str(CPG_ISLANDS) if CPG_ISLANDS.exists() else None)
        annotated = md.uns["dmr"].clone()
        annotated.write_parquet(str(sub_dir / "dmr.parquet"))

        links = build_100kb_links(annotated, gene_tss, 100_000)
        links.to_csv(sub_dir / "dmr_gene_links_100kb.csv", index=False)

        cmp = _compare(annotated, paper, panel_genes, links)
        cmp.update(dict(dis_merge_bp=dis_merge,
                        dmr_seconds=round(dt_dmr, 1),
                        n_100kb_links=int(len(links)),
                        n_100kb_unique_genes=int(links["gene"].nunique()
                                                  if len(links) else 0)))
        (sub_dir / "comparison.json").write_text(json.dumps(cmp, indent=2),
                                                 encoding="utf-8")
        rows.append(cmp)
        log.info("  recall=%s, prec=%s, dir_agree=%s, panel_E_nearest=%s",
                 cmp["recall_anybp"], cmp["precision_anybp"],
                 cmp["direction_agree_frac"],
                 cmp["panel_e_recall_nearest_tss"])

    # Summary table
    df = pd.DataFrame(rows)
    cols = ["dis_merge_bp", "n_dmr", "median_length_bp", "pct_hyper",
            "recall_anybp", "recall_J_0_25", "recall_J_0_5", "recall_J_0_75",
            "precision_anybp", "direction_agree_frac",
            "gene_recall_table5_nearest_tss", "gene_recall_table5_100kb",
            "panel_e_recall_nearest_tss", "panel_e_recall_100kb",
            "n_100kb_unique_genes", "dmr_seconds"]
    df = df[cols]
    df.to_csv(OUT_DIR / "sweep_summary.csv", index=False)

    # Markdown
    L = ["# dis.merge sweep on chain_merge replication\n",
         "All other parameters held at the paper-matched values "
         "(`alpha=1e-5, delta=0, minlen=50, minCG=3, pct.sig=0.5`).\n",
         "Paper value: `dis.merge = 100`. Larger values = more aggressive "
         "merging of adjacent chains.\n",
         "## Summary\n",
         "| dis.merge | n DMRs | median bp | %hyper | recall any-bp | "
         "recall J≥0.25 | recall J≥0.5 | recall J≥0.75 | precision | "
         "dir agree | Table 5 gene (nearest) | Table 5 gene (100kb) | "
         "Panel E (nearest) | Panel E (100kb) | DMR s |",
         "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for _, r in df.iterrows():
        L.append(
            f"| **{int(r['dis_merge_bp'])}** | {int(r['n_dmr']):,} | "
            f"{int(r['median_length_bp'])} | {r['pct_hyper']:.1f}% | "
            f"{r['recall_anybp']*100:.1f}% | {r['recall_J_0_25']*100:.1f}% | "
            f"{r['recall_J_0_5']*100:.1f}% | {r['recall_J_0_75']*100:.1f}% | "
            f"{r['precision_anybp']*100:.1f}% | "
            f"{r['direction_agree_frac']*100:.1f}% | "
            f"{r['gene_recall_table5_nearest_tss']*100:.1f}% | "
            f"{r['gene_recall_table5_100kb']*100:.1f}% | "
            f"{r['panel_e_recall_nearest_tss']*100:.1f}% | "
            f"{r['panel_e_recall_100kb']*100:.1f}% | "
            f"{r['dmr_seconds']:.1f} |"
        )
    L.append("")
    L.append("Paper reference morphology: 813 DMRs, median 240 bp, 78.5% hyper.")
    L.append("")
    (OUT_DIR / "sweep_summary.md").write_text("\n".join(L), encoding="utf-8")

    log.info("\nDone. Outputs in %s", OUT_DIR)
    log.info("Summary:\n%s", df.to_string(index=False))


if __name__ == "__main__":
    main()
