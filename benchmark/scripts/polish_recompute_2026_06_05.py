"""Recompute paper §3.3.5 / §3.3.8 / §3.3.9 numbers against post-rerun call sets.

Validates and (if needed) refreshes:
- §3.3.5 annotation distribution for ek-cm-100 (852) and ek-cm-250 (1139).
- §3.3.8 Panel-E gene capture (46 genes from panel_e_capture_dss.csv).
- §3.3.9 per-DMR effect-size concordance (Pearson r, Spearman rho, direction)
  at J >= 0.5 between epykit chain_merge and DSS-from-scratch.

Writes a single JSON report to:
  benchmark/data/study3/comparisons/epykit_vs_dss/polish_recompute_2026_06_05.json
"""

from __future__ import annotations
import json
import sys
import io
from pathlib import Path
import polars as pl
import numpy as np
from scipy.stats import pearsonr, spearmanr

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]  # benchmark/
EK100_CSV = ROOT / "data" / "study3" / "chain_merge" / "dmr_chain_merge.csv"
EK250_CSV = ROOT / "data" / "multi_thread_and_chain_sweep" / "chain_merge_dis_merge_sweep" / "dis_merge_250" / "dmr.csv"
DSS_CSV = ROOT / "data" / "study3" / "dss" / "dmr_dss.csv"
MK_TILE_CSV = ROOT / "data" / "study3" / "comparisons" / "epykit_vs_dss" / "coord_overlap_per_our_dmr.csv"  # only ek100
PANEL_E_CSV = ROOT / "data" / "study3" / "comparisons" / "epykit_vs_dss" / "panel_e_capture_dss.csv"
EK100_GENE_LINKS = ROOT / "data" / "study3" / "chain_merge" / "dmr_gene_links_100kb.csv"
DSS_GENE_LINKS = ROOT / "data" / "study3" / "dss" / "dmr_gene_links_100kb.csv"
EK100_OVERLAP_CSV = ROOT / "data" / "study3" / "comparisons" / "epykit_vs_dss" / "coord_overlap_per_our_dmr.csv"
OUT_JSON = ROOT / "data" / "study3" / "comparisons" / "epykit_vs_dss" / "polish_recompute_2026_06_05.json"


def annotation_pie(csv_path, label):
    df = pl.read_csv(csv_path)
    total = df.height
    counts = df.group_by("feature_type").len().sort("len", descending=True)
    counts_dict = {row["feature_type"]: row["len"] for row in counts.iter_rows(named=True)}
    pct = {k: round(100.0 * v / total, 2) for k, v in counts_dict.items()}
    return {"label": label, "n": total, "counts": counts_dict, "pct": pct}


def panel_e_capture(dmr_csv, gene_links_csv, panel_genes, label):
    """Two definitions: nearest-TSS (single gene per DMR) and 100kb (any gene linked within 100kb).
    """
    dmr_df = pl.read_csv(dmr_csv)
    nearest = set(g for g in dmr_df["nearest_tss_gene"].drop_nulls().to_list() if g)
    n_panel = len(panel_genes)
    captured_nearest = sum(1 for g in panel_genes if g in nearest)
    out = {
        "label": label,
        "n_panel_e": n_panel,
        "captured_nearest_tss": captured_nearest,
        "recall_nearest_tss": round(captured_nearest / n_panel, 4),
    }
    if gene_links_csv is not None and gene_links_csv.exists():
        gl = pl.read_csv(gene_links_csv)
        # column name is typically "gene_name" or "gene"
        gcol = "gene_name" if "gene_name" in gl.columns else ("gene" if "gene" in gl.columns else None)
        if gcol is not None:
            linked = set(g for g in gl[gcol].drop_nulls().to_list() if g)
            captured_100kb = sum(1 for g in panel_genes if g in linked)
            out["captured_100kb"] = captured_100kb
            out["recall_100kb"] = round(captured_100kb / n_panel, 4)
            out["gene_links_n_pairs"] = gl.height
            out["gene_links_n_unique_genes"] = len(linked)
    return out


def interval_overlap(ek_df, dss_df):
    """For each ek DMR, find DSS DMRs with any-bp overlap, return list of (ek_idx, dss_idx, jaccard, recip_frac)."""
    # Group DSS by chrom for speed
    dss_by_chrom = {}
    for i, row in enumerate(dss_df.iter_rows(named=True)):
        dss_by_chrom.setdefault(row["chrom"], []).append((i, row["start"], row["end"], row["diff_Methy_fromCounts"]))

    matches = []  # (ek_idx, ek_meth_diff, dss_idx, dss_meth_diff, jaccard, recip_frac, direction_match)
    for ek_i, ek in enumerate(ek_df.iter_rows(named=True)):
        if ek["chrom"] not in dss_by_chrom:
            continue
        es, ee = ek["start"], ek["end"]
        ek_md = ek["mean_meth_diff"]
        best_j = 0.0
        best_recip = 0.0
        best_dss_idx = -1
        best_dss_md = None
        for dss_i, ds, de, dss_md in dss_by_chrom[ek["chrom"]]:
            if min(ee, de) <= max(es, ds):
                continue  # no overlap
            inter = min(ee, de) - max(es, ds)
            union = max(ee, de) - min(es, ds)
            j = inter / union if union > 0 else 0.0
            recip = inter / min(ee - es, de - ds) if min(ee - es, de - ds) > 0 else 0.0
            if j > best_j:
                best_j = j
                best_recip = recip
                best_dss_idx = dss_i
                best_dss_md = dss_md
        if best_dss_idx >= 0:
            dir_match = (np.sign(ek_md) == np.sign(best_dss_md))
            matches.append((ek_i, ek_md, best_dss_idx, best_dss_md, best_j, best_recip, dir_match))
    return matches


def concordance_at_threshold(matches, j_threshold):
    sel = [m for m in matches if m[4] >= j_threshold]
    if len(sel) < 3:
        return {"n_pairs": len(sel), "pearson_r": None, "spearman_rho": None, "direction_agree_n": None, "direction_agree_total": None}
    ek_md = np.array([m[1] for m in sel])
    dss_md = np.array([m[3] for m in sel])
    dir_match = sum(1 for m in sel if m[6])
    pr = pearsonr(ek_md, dss_md)
    sr = spearmanr(ek_md, dss_md)
    return {
        "n_pairs": len(sel),
        "pearson_r": round(float(pr.statistic), 4),
        "spearman_rho": round(float(sr.correlation), 4),
        "direction_agree_n": dir_match,
        "direction_agree_total": len(sel),
        "direction_agree_frac": round(dir_match / len(sel), 4),
    }


def main():
    panel_genes = pl.read_csv(PANEL_E_CSV)["panel_e_gene"].to_list()
    print(f"Panel-E genes: {len(panel_genes)}")
    print(f"  e.g. {panel_genes[:5]}")

    out = {
        "generated_at": "2026-06-05",
        "context": "Phase B recompute for paper §3.3.5 / §3.3.8 / §3.3.9 against post-rerun call sets (ek-cm-100=852, ek-cm-250=1139, DSS=922).",
    }

    # § 3.3.5 annotation distribution
    out["annotation"] = {
        "ek_cm_100": annotation_pie(EK100_CSV, "ek-cm-100 (852 DMRs)"),
        "ek_cm_250": annotation_pie(EK250_CSV, "ek-cm-250 (1139 DMRs)"),
        "dss": annotation_pie(DSS_CSV, "DSS-from-scratch (922 DMRs)"),
    }

    # § 3.3.8 panel-E capture
    out["panel_e"] = {
        "ek_cm_100": panel_e_capture(EK100_CSV, EK100_GENE_LINKS, panel_genes, "ek-cm-100"),
        "ek_cm_250": panel_e_capture(EK250_CSV, None, panel_genes, "ek-cm-250"),  # gene_links table not generated for dis=250
        "dss": panel_e_capture(DSS_CSV, DSS_GENE_LINKS, panel_genes, "DSS-from-scratch"),
    }

    # § 3.3.9 per-DMR effect-size concordance
    ek100 = pl.read_csv(EK100_CSV, columns=["chrom", "start", "end", "mean_meth_diff"])
    ek250 = pl.read_csv(EK250_CSV, columns=["chrom", "start", "end", "mean_meth_diff"])
    dss = pl.read_csv(DSS_CSV, columns=["chrom", "start", "end", "diff_Methy_fromCounts"])

    print("Computing ek100 vs DSS overlaps ...")
    matches_100 = interval_overlap(ek100, dss)
    print(f"  ek100 any-bp matches: {len(matches_100)}")

    print("Computing ek250 vs DSS overlaps ...")
    matches_250 = interval_overlap(ek250, dss)
    print(f"  ek250 any-bp matches: {len(matches_250)}")

    out["per_dmr_concordance"] = {
        "ek_cm_100_vs_dss_J025": concordance_at_threshold(matches_100, 0.25),
        "ek_cm_100_vs_dss_J050": concordance_at_threshold(matches_100, 0.50),
        "ek_cm_100_vs_dss_anybp": concordance_at_threshold(matches_100, 0.0),
        "ek_cm_250_vs_dss_J025": concordance_at_threshold(matches_250, 0.25),
        "ek_cm_250_vs_dss_J050": concordance_at_threshold(matches_250, 0.50),
        "ek_cm_250_vs_dss_anybp": concordance_at_threshold(matches_250, 0.0),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_JSON}")
    print()
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
