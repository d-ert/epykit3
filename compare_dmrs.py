#!/usr/bin/env python3
"""
Compare DMR results between the published paper (DMR_total_list.xlsx)
and epykit pipeline output (results/tables/dmr.tsv).

Calculates agreement and recovery against the published reference.
The reference is not exhaustive biological truth.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd


def load_data(paper_file: str, epykit_file: str):
    """Load paper DMRs and epykit DMRs."""
    if not os.path.exists(paper_file):
        raise FileNotFoundError(f"Paper DMR file not found: {paper_file}")
    if not os.path.exists(epykit_file):
        raise FileNotFoundError(f"Epykit DMR file not found: {epykit_file}")

    paper_df = pd.read_excel(paper_file)
    epykit_df = pd.read_csv(epykit_file, sep="\t")
    return paper_df, epykit_df


def standardize_columns(paper_df: pd.DataFrame, epykit_df: pd.DataFrame):
    """Standardize column names and types for easy comparison."""
    p_df = paper_df.copy()
    e_df = epykit_df.copy()

    p_df["chr"] = p_df["chr"].astype(str)
    p_df["start"] = p_df["start"].astype(int)
    p_df["end"] = p_df["end"].astype(int)
    # Paper/DSS: one-based inclusive. Internal comparison: zero-based half-open.
    p_df["start"] -= 1
    p_df["paper_length"] = p_df["end"] - p_df["start"]
    p_df["paper_direction"] = np.where(p_df["diff.meth_mean"] > 0, "hyper", "hypo")

    e_df["chrom"] = e_df["chrom"].astype(str)
    e_df["start"] = e_df["start"].astype(int)
    e_df["end"] = e_df["end"].astype(int)
    e_df["epykit_length"] = e_df["end"] - e_df["start"]
    if "dmr_type" not in e_df.columns:
        e_df["dmr_type"] = np.where(e_df["mean_meth_diff"] > 0, "hyper", "hypo")

    return p_df, e_df


def find_overlaps(
    p_df: pd.DataFrame,
    e_df: pd.DataFrame,
    min_overlap_bp: int = 50,
    min_overlap_frac: float = 0.0,
):
    """Find overlapping regions between epykit and paper DMRs."""
    overlap_records = []
    p_grouped = {c: g for c, g in p_df.groupby("chr")}

    for e_idx, e_row in e_df.iterrows():
        chrom = e_row["chrom"]
        if chrom not in p_grouped:
            continue

        p_chrom_df = p_grouped[chrom]
        ov_mask = (p_chrom_df["start"] < e_row["end"]) & (p_chrom_df["end"] > e_row["start"])
        ov_matches = p_chrom_df[ov_mask]

        for p_idx, p_row in ov_matches.iterrows():
            ov_start = max(e_row["start"], p_row["start"])
            ov_end = min(e_row["end"], p_row["end"])
            ov_len = max(0, ov_end - ov_start)

            e_len = max(1, e_row["end"] - e_row["start"])
            p_len = max(1, p_row["end"] - p_row["start"])
            min_frac = min(ov_len / e_len, ov_len / p_len)

            if ov_len >= min_overlap_bp and min_frac >= min_overlap_frac:
                same_direction = (e_row["mean_meth_diff"] * p_row["diff.meth_mean"]) > 0

                overlap_records.append({
                    "chrom": chrom,
                    "epykit_start": e_row["start"],
                    "epykit_end": e_row["end"],
                    "epykit_len": e_len,
                    "epykit_cpgs": e_row.get("n_cpgs", np.nan),
                    "epykit_meth_diff": e_row["mean_meth_diff"],
                    "epykit_type": e_row["dmr_type"],
                    "epykit_qvalue": e_row.get("combined_qvalue", np.nan),
                    "paper_start": p_row["start"],
                    "paper_end": p_row["end"],
                    "paper_len": p_len,
                    "paper_cpgs": p_row.get("nCG", np.nan),
                    "paper_meth_diff": p_row["diff.meth_mean"],
                    "paper_direction": p_row["paper_direction"],
                    "paper_gene": p_row.get("Gene.Name", ""),
                    "paper_annotation": p_row.get("Annotation", ""),
                    "overlap_bp": ov_len,
                    "concordant_direction": same_direction,
                    "epykit_index": e_idx,
                    "paper_index": p_idx,
                })

    return pd.DataFrame(overlap_records)


def merge_genomic_intervals(df: pd.DataFrame, chrom_col: str, start_col: str, end_col: str):
    """Merge overlapping intervals by chromosome."""
    merged = {}
    for chrom, group in df.groupby(chrom_col):
        sorted_ints = group[[start_col, end_col]].sort_values(by=start_col).values.tolist()
        m = []
        for s, e in sorted_ints:
            if not m or m[-1][1] < s:
                m.append([s, e])
            else:
                m[-1][1] = max(m[-1][1], e)
        merged[chrom] = m
    return merged


def calculate_ground_truth_metrics(p_df: pd.DataFrame, e_df: pd.DataFrame, ov_df: pd.DataFrame):
    """Calculate reference agreement/recovery; legacy column names are retained."""
    n_truth_total = len(p_df)
    n_pred_total = len(e_df)

    tp_pred = ov_df["epykit_index"].nunique() if not ov_df.empty else 0
    tp_truth = ov_df["paper_index"].nunique() if not ov_df.empty else 0

    prec = tp_pred / n_pred_total if n_pred_total > 0 else 0.0
    rec = tp_truth / n_truth_total if n_truth_total > 0 else 0.0
    f1 = (2 * prec * rec) / (prec + rec) if (prec + rec) > 0 else 0.0

    # Hyper
    p_hyper = p_df[p_df["paper_direction"] == "hyper"]
    e_hyper = e_df[e_df["dmr_type"] == "hyper"]
    ov_hyper = ov_df[(ov_df["paper_direction"] == "hyper") & (ov_df["epykit_type"] == "hyper")] if not ov_df.empty else pd.DataFrame()
    tp_pred_hyper = ov_hyper["epykit_index"].nunique() if not ov_hyper.empty else 0
    tp_truth_hyper = ov_hyper["paper_index"].nunique() if not ov_hyper.empty else 0
    prec_hyper = tp_pred_hyper / len(e_hyper) if len(e_hyper) > 0 else 0.0
    rec_hyper = tp_truth_hyper / len(p_hyper) if len(p_hyper) > 0 else 0.0
    f1_hyper = (2 * prec_hyper * rec_hyper) / (prec_hyper + rec_hyper) if (prec_hyper + rec_hyper) > 0 else 0.0

    # Hypo
    p_hypo = p_df[p_df["paper_direction"] == "hypo"]
    e_hypo = e_df[e_df["dmr_type"] == "hypo"]
    ov_hypo = ov_df[(ov_df["paper_direction"] == "hypo") & (ov_df["epykit_type"] == "hypo")] if not ov_df.empty else pd.DataFrame()
    tp_pred_hypo = ov_hypo["epykit_index"].nunique() if not ov_hypo.empty else 0
    tp_truth_hypo = ov_hypo["paper_index"].nunique() if not ov_hypo.empty else 0
    prec_hypo = tp_pred_hypo / len(e_hypo) if len(e_hypo) > 0 else 0.0
    rec_hypo = tp_truth_hypo / len(p_hypo) if len(p_hypo) > 0 else 0.0
    f1_hypo = (2 * prec_hypo * rec_hypo) / (prec_hypo + rec_hypo) if (prec_hypo + rec_hypo) > 0 else 0.0

    # Basepair-level
    p_merged = merge_genomic_intervals(p_df, "chr", "start", "end")
    e_merged = merge_genomic_intervals(e_df, "chrom", "start", "end")
    total_truth_bp = sum(sum(e - s for s, e in ints) for ints in p_merged.values())
    total_pred_bp = sum(sum(e - s for s, e in ints) for ints in e_merged.values())

    shared_bp = 0
    for chrom, e_ints in e_merged.items():
        if chrom in p_merged:
            for es, ee in e_ints:
                for ps, pe in p_merged[chrom]:
                    ov_s = max(es, ps)
                    ov_e = min(ee, pe)
                    if ov_s < ov_e:
                        shared_bp += (ov_e - ov_s)

    prec_bp = shared_bp / total_pred_bp if total_pred_bp > 0 else 0.0
    rec_bp = shared_bp / total_truth_bp if total_truth_bp > 0 else 0.0
    f1_bp = (2 * prec_bp * rec_bp) / (prec_bp + rec_bp) if (prec_bp + rec_bp) > 0 else 0.0

    metrics_rows = [
        {"Scope": "Region-Level", "Truth": n_truth_total, "Pred": n_pred_total, "Precision": prec, "Recall": rec, "F1": f1},
        {"Scope": "Stratified: Hyper", "Truth": len(p_hyper), "Pred": len(e_hyper), "Precision": prec_hyper, "Recall": rec_hyper, "F1": f1_hyper},
        {"Scope": "Stratified: Hypo", "Truth": len(p_hypo), "Pred": len(e_hypo), "Precision": prec_hypo, "Recall": rec_hypo, "F1": f1_hypo},
        {"Scope": "Basepair-Level (bp)", "Truth": total_truth_bp, "Pred": total_pred_bp, "Precision": prec_bp, "Recall": rec_bp, "F1": f1_bp},
    ]
    return pd.DataFrame(metrics_rows)


def sweep_thresholds(p_df: pd.DataFrame, e_df: pd.DataFrame):
    """Compute sensitivity across standard thresholds."""
    configs = [
        ("≥ 1 bp (any)", 1, 0.0),
        ("≥ 30 bp", 30, 0.0),
        ("≥ 50 bp (default)", 50, 0.0),
        ("≥ 50% reciprocal", 1, 0.50),
    ]
    sweep_results = []
    n_p, n_e = len(p_df), len(e_df)
    p_grouped = {c: g for c, g in p_df.groupby("chr")}

    for label, min_bp, min_frac in configs:
        tp_e = 0
        tp_p = set()
        conc = 0

        for _, e_row in e_df.iterrows():
            c = e_row["chrom"]
            if c not in p_grouped:
                continue
            matches = p_grouped[c][(p_grouped[c]["start"] < e_row["end"]) & (p_grouped[c]["end"] > e_row["start"])]

            matched_e = False
            matched_conc = False
            for _, p_row in matches.iterrows():
                ov_len = max(0, min(e_row["end"], p_row["end"]) - max(e_row["start"], p_row["start"]))
                frac = min(ov_len / max(1, e_row["end"] - e_row["start"]), ov_len / max(1, p_row["end"] - p_row["start"]))
                if ov_len >= min_bp and frac >= min_frac:
                    matched_e = True
                    tp_p.add(p_row.name)
                    if (e_row["mean_meth_diff"] * p_row["diff.meth_mean"]) > 0:
                        matched_conc = True

            if matched_e:
                tp_e += 1
            if matched_conc:
                conc += 1

        prec = tp_e / n_e if n_e > 0 else 0
        rec = len(tp_p) / n_p if n_p > 0 else 0
        f1 = (2 * prec * rec) / (prec + rec) if (prec + rec) > 0 else 0
        conc_pct = (conc / tp_e * 100) if tp_e > 0 else 0.0

        sweep_results.append({
            "Threshold": label,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,
            "Direction_Agreement": conc_pct,
        })
    return pd.DataFrame(sweep_results)


def print_clean_report(
    p_df: pd.DataFrame,
    e_df: pd.DataFrame,
    ov_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    sweep_df: pd.DataFrame,
    threshold_str: str,
    verbose: bool = False,
):
    """Print a concise, punchy, on-point benchmark report."""
    n_p = len(p_df)
    n_e = len(e_df)
    p_hyper = (p_df["paper_direction"] == "hyper").sum()
    p_hypo = (p_df["paper_direction"] == "hypo").sum()
    e_hyper = (e_df["dmr_type"] == "hyper").sum()
    e_hypo = (e_df["dmr_type"] == "hypo").sum()

    n_ov = len(ov_df)
    n_conc = ov_df["concordant_direction"].sum() if not ov_df.empty else 0
    corr = ov_df[["epykit_meth_diff", "paper_meth_diff"]].corr().iloc[0, 1] if not ov_df.empty else 0.0

    print("\n" + "=" * 68)
    print(f"            EPYKIT vs PAPER DMR BENCHMARK ({threshold_str})")
    print("=" * 68)

    print("\n[COUNTS & DIRECTION]")
    print(f"  • Paper (Reference): {n_p:>4} DMRs ({p_hyper} hyper, {p_hypo} hypo)")
    print(f"  • Epykit (Pred): {n_e:>4} DMRs ({e_hyper} hyper, {e_hypo} hypo)")

    print(f"\n[CONCORDANCE: {threshold_str}]")
    print(f"  • Direction agreement: {n_conc}/{n_ov} pairs ({(n_conc/n_ov*100 if n_ov else 0):.1f}%)")
    print(f"  • Methylation diff correlation: Pearson r = {corr:.4f}")

    print(f"\n[REFERENCE AGREEMENT (does not estimate biological FDR)]")
    hdr_m = f"  {'Scope':<23} {'Truth':<10} {'Pred':<8} {'Prec (%)':<10} {'Rec (%)':<10} {'F1':<6}"
    print(hdr_m)
    print("  " + "-" * (len(hdr_m) - 2))
    for _, r in metrics_df.iterrows():
        t_str = f"{int(r['Truth']):,}"
        p_str = f"{int(r['Pred']):,}"
        prec_str = f"{r['Precision']*100:.2f}%"
        rec_str = f"{r['Recall']*100:.2f}%"
        f1_str = f"{r['F1']:.4f}"
        print(f"  {r['Scope']:<23} {t_str:<10} {p_str:<8} {prec_str:<10} {rec_str:<10} {f1_str:<6}")

    print("\n[THRESHOLD SENSITIVITY]")
    hdr_s = f"  {'Threshold':<20} {'Precision':<12} {'Recall':<12} {'F1 Score':<10} {'Direction Agreement':<18}"
    print(hdr_s)
    print("  " + "-" * (len(hdr_s) - 2))
    for _, r in sweep_df.iterrows():
        print(f"  {r['Threshold']:<20} {r['Precision']*100:>6.2f}%      {r['Recall']*100:>6.2f}%      {r['F1']:>6.4f}     {r['Direction_Agreement']:>5.1f}%")

    if verbose and not ov_df.empty:
        print(f"\n[TOP OVERLAPPING DMRS (Sample of 5)]")
        top_ov = ov_df.sort_values(by="overlap_bp", ascending=False).head(5)
        for _, row in top_ov.iterrows():
            print(f"  • {row['chrom']}:{row['epykit_start']}-{row['epykit_end']} ({row['paper_gene']}): Epy={row['epykit_meth_diff']:+.2f}, Pap={row['paper_meth_diff']:+.2f}")

    print("\n" + "=" * 68 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare paper DMRs with epykit DMRs (concise benchmark report)"
    )
    default_paper = "/Users/dert/Documents/Coding/Epykit/epykit3/paper_dmrs/DMR_total_list.xlsx"
    default_epykit_dir = "/Users/dert/Documents/Coding/Epykit/epykit3/results/tables"
    default_epykit = os.path.join(default_epykit_dir, "dmr.tsv")
    default_output_overlaps = os.path.join(default_epykit_dir, "dmr_comparison_with_paper.tsv")
    default_output_metrics = os.path.join(default_epykit_dir, "dmr_benchmark_metrics.tsv")

    parser.add_argument("--paper", default=default_paper, help="Path to paper DMR excel file")
    parser.add_argument("--epykit", default=default_epykit, help="Path to epykit dmr.tsv file")
    parser.add_argument("--output-overlaps", default=default_output_overlaps, help="Path to save overlapping DMRs TSV")
    parser.add_argument("--output-metrics", default=default_output_metrics, help="Path to save benchmark metrics TSV")
    parser.add_argument("--min-overlap-bp", type=int, default=50, help="Minimum bp overlap (default: 50 bp)")
    parser.add_argument("--min-overlap-frac", type=float, default=0.0, help="Minimum reciprocal fraction (default: 0.0)")
    parser.add_argument("--verbose", action="store_true", help="Print detailed individual DMR list")

    args = parser.parse_args()

    paper_df, epykit_df = load_data(args.paper, args.epykit)
    p_df, e_df = standardize_columns(paper_df, epykit_df)

    ov_df = find_overlaps(p_df, e_df, min_overlap_bp=args.min_overlap_bp, min_overlap_frac=args.min_overlap_frac)

    threshold_str = f"≥ {args.min_overlap_bp} bp"
    if args.min_overlap_frac > 0:
        threshold_str += f" & ≥ {int(args.min_overlap_frac*100)}% reciprocal"

    metrics_df = calculate_ground_truth_metrics(p_df, e_df, ov_df)
    sweep_df = sweep_thresholds(p_df, e_df)

    print_clean_report(p_df, e_df, ov_df, metrics_df, sweep_df, threshold_str, verbose=args.verbose)

    if not ov_df.empty and args.output_overlaps:
        save_df = ov_df.drop(columns=["epykit_index", "paper_index"])
        save_df.to_csv(args.output_overlaps, sep="\t", index=False)
        print(f"Saved: {args.output_overlaps}")

    if not metrics_df.empty and args.output_metrics:
        metrics_df.to_csv(args.output_metrics, sep="\t", index=False)
        print(f"Saved: {args.output_metrics}")


if __name__ == "__main__":
    main()
