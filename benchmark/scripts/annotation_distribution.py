"""Build a long-form annotation feature distribution table across callers.

For each caller (paper, methylKit-tile, ek-tile, ek-chain_merge-100,
ek-chain_merge-250, DSS) produce (feature_type, count, fraction) and
compute chi-square vs the paper distribution.

Annotation source: UCSC refGene (HOMER-equivalent), same catalog HOMER
ships, so all callers share the same feature definitions.

Outputs (FINAL_REPORT/data/study3/comparisons/):
  annotation_distribution.csv   long form: caller × feature_type → count, fraction
  annotation_distribution_summary.json   chi-square vs paper + bias commentary
  methylkit_dmrs_annotated.csv  side-effect: methylKit DMRs with feature_type

Feeds figure F5.
"""

from __future__ import annotations

import bisect
import gzip
import json
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy.stats import chi2_contingency

warnings.filterwarnings("ignore")

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
RAW_DIR   = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

CM_PQ      = (REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
              / "dmr_chain_merge.parquet")
DM250_PQ   = (REPO_ROOT / "FINAL_REPORT" / "data" / "study3"
              / "chain_merge_dis_merge_sweep" / "dis_merge_250" / "dmr.parquet")
DSS_CSV    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss" / "dmr_dss.csv"
EK_TILE    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dmr_significant_lenient.csv"
MK_TILE    = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/"
                  r"scripts_and_results/methylkit_results/"
                  r"dmr_significant_lenient.csv")
PAPER_T5   = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"
REFGENE    = RAW_DIR / "refseq" / "refGene.txt.gz"

OUT_DIR    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "comparisons"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# HOMER feature priority (lower = higher priority)
FEAT_PRIO = {
    "promoter-TSS": 0, "TTS": 1, "5UTR": 2, "3UTR": 3,
    "exon": 4, "intron": 5, "non-coding": 6, "intergenic": 7,
}
FEAT_ORDER = ["promoter-TSS", "5UTR", "exon", "intron", "3UTR",
              "TTS", "non-coding", "intergenic"]

PROMOTER_UP, PROMOTER_DOWN = 1000, 100
TTS_UP, TTS_DOWN = 100, 1000


# ---- refGene loading (same algorithm as compare_homer_refseq.py) -----------

def load_refgene(path: Path):
    print(f"Loading refGene {path.name} …", flush=True)
    genes = []
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom = p[2]
            if not chrom.startswith("chr") or "_" in chrom:
                continue
            rec = dict(
                acc=p[1], chrom=chrom, strand=p[3],
                txStart=int(p[4]), txEnd=int(p[5]),
                cdsStart=int(p[6]), cdsEnd=int(p[7]),
                exonStarts=[int(x) for x in p[9].rstrip(",").split(",")],
                exonEnds=[int(x) for x in p[10].rstrip(",").split(",")],
                gene=p[12],
                coding=p[1].startswith("NM_"),
            )
            rec["tss"] = rec["txStart"] if rec["strand"] == "+" else rec["txEnd"]
            rec["tts"] = rec["txEnd"]   if rec["strand"] == "+" else rec["txStart"]
            genes.append(rec)
    by_chr: dict[str, list] = defaultdict(list)
    for g in genes:
        by_chr[g["chrom"]].append(g)
    print(f"  {len(genes):,} transcripts across {len(by_chr)} chromosomes",
          flush=True)
    return by_chr


def classify_center(center: int, chrom: str, by_chr: dict) -> str:
    cands = by_chr.get(chrom, [])
    if not cands:
        return "intergenic"
    best_prio = 8; best_feat = "intergenic"
    for g in cands:
        if g["strand"] == "+":
            prom_lo, prom_hi = g["tss"] - PROMOTER_UP, g["tss"] + PROMOTER_DOWN
            tts_lo,  tts_hi  = g["tts"] - TTS_UP,     g["tts"] + TTS_DOWN
        else:
            prom_lo, prom_hi = g["tss"] - PROMOTER_DOWN, g["tss"] + PROMOTER_UP
            tts_lo,  tts_hi  = g["tts"] - TTS_DOWN,     g["tts"] + TTS_UP
        feat = None
        if prom_lo <= center <= prom_hi:
            feat = "promoter-TSS"
        elif tts_lo <= center <= tts_hi:
            feat = "TTS"
        elif g["txStart"] <= center <= g["txEnd"]:
            if not g["coding"]:
                feat = "non-coding"
            else:
                in_exon = any(es <= center <= ee
                              for es, ee in zip(g["exonStarts"], g["exonEnds"]))
                if not in_exon:
                    feat = "intron"
                elif g["cdsStart"] == g["cdsEnd"]:
                    feat = "non-coding"
                elif g["strand"] == "+":
                    if center < g["cdsStart"]:
                        feat = "5UTR"
                    elif center >= g["cdsEnd"]:
                        feat = "3UTR"
                    else:
                        feat = "exon"
                else:
                    if center >= g["cdsEnd"]:
                        feat = "5UTR"
                    elif center < g["cdsStart"]:
                        feat = "3UTR"
                    else:
                        feat = "exon"
        if feat is not None:
            p = FEAT_PRIO[feat]
            if p < best_prio:
                best_prio = p; best_feat = feat
    return best_feat


# ---- Per-source loaders ---------------------------------------------------

def annotate_via_refgene(df: pd.DataFrame, by_chr: dict) -> pd.DataFrame:
    """Annotate a (chrom, start, end) table with HOMER feature_type."""
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    feats = []
    for _, r in df.iterrows():
        center = (int(r["start"]) + int(r["end"])) // 2
        feats.append(classify_center(center, str(r["chrom"]), by_chr))
    df["feature_type"] = feats
    return df


def normalize_paper_annotation(s: str) -> str:
    """Map paper's Annotation strings to our 8-feature schema."""
    if pd.isna(s):
        return "intergenic"
    s = str(s).strip().lower()
    # Strip HOMER multi-anno suffixes like ".2", ".7"
    if "." in s:
        s = s.split(".")[0].strip()
    if "promoter" in s: return "promoter-TSS"
    if "tts" in s:      return "TTS"
    if "5" in s and "utr" in s: return "5UTR"
    if "3" in s and "utr" in s: return "3UTR"
    if "exon" in s:     return "exon"
    if "intron" in s:   return "intron"
    if "non" in s and "coding" in s: return "non-coding"
    if "intergenic" in s: return "intergenic"
    return "intergenic"


def normalize_epykit_feature(s: str) -> str:
    """epykit's feature_type is already in this schema (case-sensitive)."""
    if pd.isna(s) or s == "":
        return "intergenic"
    s = str(s).strip()
    # Handle some variants
    if s.lower().startswith("promoter"): return "promoter-TSS"
    if s.lower().startswith("tts"):       return "TTS"
    if s.lower() in {"5utr", "5'utr"}:    return "5UTR"
    if s.lower() in {"3utr", "3'utr"}:    return "3UTR"
    return s if s in FEAT_PRIO else "intergenic"


# ---- Main pipeline --------------------------------------------------------

def main() -> None:
    by_chr = load_refgene(REFGENE)

    # ---- Paper ----------------------------------------------------------
    paper = pd.read_excel(PAPER_T5, sheet_name=0)
    paper["feature_type"] = paper["Annotation"].apply(normalize_paper_annotation)

    # ---- methylKit tile (needs annotation) ------------------------------
    mk = pd.read_csv(MK_TILE)[["chrom", "start", "end"]]
    print(f"methylKit-tile {len(mk)} DMRs → annotating …", flush=True)
    mk = annotate_via_refgene(mk, by_chr)
    mk.to_csv(OUT_DIR / "methylkit_dmrs_annotated.csv", index=False)

    # ---- epykit-tile (already in feature_type if any) — re-annotate via refGene
    ek_tile_df = pd.read_csv(EK_TILE)[["chrom", "start", "end"]]
    print(f"epykit-tile {len(ek_tile_df)} DMRs → annotating …", flush=True)
    ek_tile_df = annotate_via_refgene(ek_tile_df, by_chr)

    # ---- epykit chain_merge (use existing feature_type, normalize) -----
    ek100 = pl.read_parquet(CM_PQ).to_pandas()
    ek100["feature_type"] = ek100["feature_type"].apply(normalize_epykit_feature)

    ek250 = pl.read_parquet(DM250_PQ).to_pandas()
    ek250["feature_type"] = ek250["feature_type"].apply(normalize_epykit_feature)

    # ---- DSS (use existing feature_type, normalize) --------------------
    dss = pd.read_csv(DSS_CSV)
    dss["feature_type"] = dss["feature_type"].apply(normalize_epykit_feature)

    # ---- Build long-form table -----------------------------------------
    sources = [
        ("paper-DSS (Supp Table 5)", paper),
        ("methylKit-tile",            mk),
        ("epykit-tile",               ek_tile_df),
        ("epykit-chain_merge-100",    ek100),
        ("epykit-chain_merge-250",    ek250),
        ("DSS-from-scratch",          dss),
    ]

    rows = []
    for caller, df in sources:
        n = len(df)
        counts = df["feature_type"].value_counts()
        for feat in FEAT_ORDER:
            c = int(counts.get(feat, 0))
            rows.append(dict(
                caller=caller,
                feature_type=feat,
                count=c,
                fraction=round(c / max(n, 1), 4),
                n_total=n,
            ))
    long = pd.DataFrame(rows)
    long.to_csv(OUT_DIR / "annotation_distribution.csv", index=False)

    # ---- Chi-square vs paper ------------------------------------------
    paper_counts = np.array([
        int(paper["feature_type"].value_counts().get(f, 0))
        for f in FEAT_ORDER
    ])
    summary = {"paper_total": int(paper_counts.sum())}
    for caller, df in sources:
        if caller.startswith("paper"):
            continue
        our_counts = np.array([
            int(df["feature_type"].value_counts().get(f, 0))
            for f in FEAT_ORDER
        ])
        ct = np.array([paper_counts, our_counts])
        # Add 0.5 to avoid zero-cell warnings? Use raw — chi2_contingency
        # tolerates zeros.
        chi2, pval, dof, _ = chi2_contingency(ct + 0)
        # Pearson residuals per feature for diagnostics
        expected = ct.sum(axis=0, keepdims=True) * \
                   (ct.sum(axis=1, keepdims=True) / ct.sum())
        resid = (ct - expected) / np.sqrt(expected + 1e-9)
        summary[caller] = {
            "n": int(our_counts.sum()),
            "chi2": round(float(chi2), 2),
            "pvalue": float(pval),
            "dof": int(dof),
            "per_feature_residual_vs_paper": {
                FEAT_ORDER[i]: round(float(resid[1, i]), 2)
                for i in range(len(FEAT_ORDER))
            },
        }
    (OUT_DIR / "annotation_distribution_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    # ---- Print preview ------------------------------------------------
    print("\n--- ANNOTATION DISTRIBUTION (% of total per caller) ---")
    pivot = long.pivot(index="feature_type", columns="caller",
                       values="fraction")
    pivot = pivot.reindex(FEAT_ORDER)
    cols = ["paper-DSS (Supp Table 5)", "methylKit-tile", "epykit-tile",
            "epykit-chain_merge-100", "epykit-chain_merge-250",
            "DSS-from-scratch"]
    pivot = pivot[[c for c in cols if c in pivot.columns]]
    print((pivot * 100).round(1).to_string())
    print("\nChi-square vs paper (smaller chi2 ≈ closer to paper):")
    for caller, info in summary.items():
        if isinstance(info, dict):
            print(f"  {caller}: chi2={info['chi2']:>8.1f}  "
                  f"p={info['pvalue']:.2e}  n={info['n']}")


if __name__ == "__main__":
    main()
