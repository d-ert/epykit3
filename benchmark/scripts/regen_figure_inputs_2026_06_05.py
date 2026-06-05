"""Regenerate the intermediate comparison files the three-way figure
scripts (figures_v2/f*.py) consume, from the CURRENT post-rerun data tree.

The figure scripts were authored against an external pre-rerun tree
(benchmarkin_merges/FINAL_REPORT) whose chain_merge call sets were
702 / 940. The committed epykit3/benchmark tree is post-rerun
(852 / 1,139 vs DSS-922). This script rebuilds the three intermediate
files that drifted, so the figures can be regenerated against current
data:

  1. comparisons/annotation_distribution.csv      (F5, F10)
  2. comparisons/per_dmr_stat_concordance.csv      (F9)
  3. chain_merge_dis_merge_sweep/sweep_summary.csv (F2, F10)

Architecturally-unchanged caller rows (paper-DSS, methylKit-tile,
epykit-tile) are carried forward verbatim from the external file;
only the chain_merge-100 / chain_merge-250 / DSS rows are recomputed
from current data.
"""

from __future__ import annotations
import sys
import io
from pathlib import Path
import numpy as np
import polars as pl

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

BENCH = Path(__file__).resolve().parents[1]
STUDY3 = BENCH / "data" / "study3"
SWEEP = BENCH / "data" / "multi_thread_and_chain_sweep" / "chain_merge_dis_merge_sweep"
COMP = STUDY3 / "comparisons"
EXT_ANNOT = Path(
    r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges/FINAL_REPORT/data/study3/comparisons/annotation_distribution.csv"
)

EK100 = STUDY3 / "chain_merge" / "dmr_chain_merge.csv"
EK250 = SWEEP / "dis_merge_250" / "dmr.csv"
DSS = STUDY3 / "dss" / "dmr_dss.csv"
PANEL_E = COMP / "epykit_vs_dss" / "panel_e_capture_dss.csv"
SENS = SWEEP / "dis_merge_vs_dss_sensitivity.csv"

# Feature-type label normalisation: current DMR CSVs use short labels,
# the figure's FEAT_ORDER uses HOMER-style labels.
FEAT_MAP = {
    "promoter": "promoter-TSS",
    "promoter-TSS": "promoter-TSS",
    "5UTR": "5UTR",
    "exon": "exon",
    "intron": "intron",
    "3UTR": "3UTR",
    "TTS": "TTS",
    "noncoding": "non-coding",
    "non-coding": "non-coding",
    "intergenic": "intergenic",
}
FEAT_ORDER = ["promoter-TSS", "5UTR", "exon", "intron", "3UTR",
              "TTS", "non-coding", "intergenic"]


def annotation_rows(dmr_csv, caller, n_total):
    df = pl.read_csv(dmr_csv)
    counts = {}
    for row in df.group_by("feature_type").len().iter_rows(named=True):
        key = FEAT_MAP.get(row["feature_type"], row["feature_type"])
        counts[key] = counts.get(key, 0) + row["len"]
    rows = []
    for feat in FEAT_ORDER:
        c = counts.get(feat, 0)
        rows.append({
            "caller": caller,
            "feature_type": feat,
            "count": c,
            "fraction": round(c / n_total, 4),
            "n_total": n_total,
        })
    return rows


def regen_annotation():
    ext = pl.read_csv(EXT_ANNOT)
    # Carry forward the rows whose DMR sets did not change in the rerun.
    carry = ext.filter(
        pl.col("caller").is_in([
            "paper-DSS (Supp Table 5)", "methylKit-tile", "epykit-tile",
        ])
    )
    new_rows = []
    new_rows += annotation_rows(EK100, "epykit-chain_merge-100", 852)
    new_rows += annotation_rows(EK250, "epykit-chain_merge-250", 1139)
    new_rows += annotation_rows(DSS, "DSS-from-scratch", 922)
    new_df = pl.DataFrame(new_rows, schema=carry.schema)
    out = pl.concat([carry, new_df])
    COMP.mkdir(parents=True, exist_ok=True)
    out.write_csv(COMP / "annotation_distribution.csv")
    print(f"wrote {COMP / 'annotation_distribution.csv'} ({out.height} rows)")


def overlap_pairs(ek_csv, dss_df, caller):
    """Per-(ek DMR, best DSS DMR) any-bp overlap rows with Δβ + jaccard."""
    ek_df = pl.read_csv(ek_csv, columns=["chrom", "start", "end", "mean_meth_diff"])
    dss_by_chrom = {}
    for i, row in enumerate(dss_df.iter_rows(named=True)):
        dss_by_chrom.setdefault(row["chrom"], []).append(
            (row["start"], row["end"], row["diff_Methy_fromCounts"]))
    rows = []
    for ek in ek_df.iter_rows(named=True):
        if ek["chrom"] not in dss_by_chrom:
            continue
        es, ee, ek_md = ek["start"], ek["end"], ek["mean_meth_diff"]
        best_j, best_md = 0.0, None
        for ds, de, dss_md in dss_by_chrom[ek["chrom"]]:
            if min(ee, de) <= max(es, ds):
                continue
            inter = min(ee, de) - max(es, ds)
            union = max(ee, de) - min(es, ds)
            j = inter / union if union > 0 else 0.0
            if j > best_j:
                best_j, best_md = j, dss_md
        if best_md is not None:
            rows.append({
                "caller": caller,
                "jaccard": round(best_j, 6),
                "ek_mean_meth_diff": ek_md,
                "dss_diff_Methy_fromCounts": best_md,
            })
    return rows


def regen_concordance():
    dss_df = pl.read_csv(DSS, columns=["chrom", "start", "end", "diff_Methy_fromCounts"])
    rows = []
    rows += overlap_pairs(EK100, dss_df, "ek100")
    rows += overlap_pairs(EK250, dss_df, "ek250")
    out = pl.DataFrame(rows)
    out.write_csv(COMP / "per_dmr_stat_concordance.csv")
    n100 = sum(1 for r in rows if r["caller"] == "ek100")
    n250 = sum(1 for r in rows if r["caller"] == "ek250")
    print(f"wrote {COMP / 'per_dmr_stat_concordance.csv'} (ek100={n100}, ek250={n250} pairs)")


def panel_e_recall(dmr_csv, panel_genes):
    df = pl.read_csv(dmr_csv)
    nearest = set(g for g in df["nearest_tss_gene"].drop_nulls().to_list() if g)
    return sum(1 for g in panel_genes if g in nearest) / len(panel_genes)


def regen_sweep_summary():
    sens = pl.read_csv(SENS)
    panel_genes = pl.read_csv(PANEL_E)["panel_e_gene"].to_list()
    rows = []
    for r in sens.iter_rows(named=True):
        dm = r["dis_merge_bp"]
        dmr_csv = (EK100 if dm == 100 else SWEEP / f"dis_merge_{dm}" / "dmr.csv")
        df = pl.read_csv(dmr_csv)
        lengths = (df["end"] - df["start"]).to_numpy()
        pct_hyper = (df["dmr_type"] == "hyper").mean() if "dmr_type" in df.columns else None
        rows.append({
            "dis_merge_bp": dm,
            "n_dmr": r["n_epykit"],
            "median_length_bp": float(np.median(lengths)),
            "pct_hyper": round(float(pct_hyper) * 100, 1) if pct_hyper is not None else None,
            "recall_anybp": r["recall_anybp"],
            "recall_J_0_25": r["recall_J_0_25"],
            "recall_J_0_5": r["recall_J_0_5"],
            "recall_J_0_75": r["recall_J_0_75"],
            "precision_anybp": r["precision_anybp"],
            "direction_agree_frac": r["direction_agree_frac"],
            "panel_e_recall_nearest_tss": round(panel_e_recall(dmr_csv, panel_genes), 4),
            "n_dss": r["n_dss"],
        })
    out = pl.DataFrame(rows)
    out.write_csv(SWEEP / "sweep_summary.csv")
    print(f"wrote {SWEEP / 'sweep_summary.csv'} ({out.height} rows)")
    print(out)


def main():
    regen_annotation()
    regen_concordance()
    regen_sweep_summary()


if __name__ == "__main__":
    main()
