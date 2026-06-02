"""Per-DMR statistic concordance: epykit chain_merge vs DSS.

For each DMR that overlaps with Jaccard >= 0.5 between epykit (any of
ek-100 or ek-250) and DSS, join the statistics from both sides:

    epykit:  combined_qvalue, mean_meth_diff, n_cpgs, length
    DSS:     areaStat, diff_Methy_DSSfit, diff_Methy_fromCounts, length,
             pvals_mean_in_DMR (NOT YET — derive from cached DMLtest if
             desired in a follow-up; for now use areaStat as the
             significance summary).

Outputs (FINAL_REPORT/data/study3/comparisons/):
  per_dmr_stat_concordance.csv   one row per matched (ek, DSS) DMR pair
                                  with both sides' statistics
  per_dmr_stat_concordance_summary.json  Pearson r on Δβ, Spearman ρ on
                                          (epykit_q, DSS_areaStat),
                                          counts at each Jaccard cutoff.

Feeds figure F9.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats as sps

warnings.filterwarnings("ignore")

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
CM_DIR    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
DM250_PQ  = (REPO_ROOT / "FINAL_REPORT" / "data" / "study3"
             / "chain_merge_dis_merge_sweep" / "dis_merge_250" / "dmr.parquet")
DSS_CSV   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss" / "dmr_dss.csv"
OUT_DIR   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "comparisons"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def jaccard(a_s, a_e, b_s, b_e) -> float:
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def best_match_pairs(ek: pd.DataFrame, dss: pd.DataFrame,
                     min_j: float = 0.0) -> pd.DataFrame:
    """For each ek DMR, find the DSS DMR with the highest Jaccard overlap
    (if any). Returns long-form table joining both sides."""
    by_chrom: dict[str, list] = {}
    for r in dss.itertuples():
        by_chrom.setdefault(str(r.chrom), []).append(
            (int(r.start), int(r.end), int(r.idx))
        )
    for ch in by_chrom:
        by_chrom[ch].sort(key=lambda t: t[0])

    rows = []
    for r in ek.itertuples():
        ch = str(r.chrom); s = int(r.start); e = int(r.end)
        cands = by_chrom.get(ch, [])
        best_j = 0.0; best_i = -1
        for ts, te, ti in cands:
            if te < s:
                continue
            if ts > e:
                break
            j = jaccard(s, e, ts, te)
            if j > best_j:
                best_j = j; best_i = ti
        if best_i == -1 or best_j < min_j:
            continue
        d = dss[dss["idx"] == best_i].iloc[0]
        rows.append(dict(
            ek_idx=int(r.idx),
            ek_chrom=ch, ek_start=s, ek_end=e,
            ek_length=int(e - s),
            ek_n_cpgs=int(getattr(r, "n_cpgs", 0)),
            ek_n_significant=int(getattr(r, "n_significant", 0)),
            ek_mean_meth_diff=float(r.mean_meth_diff),
            ek_combined_qvalue=float(getattr(r, "combined_qvalue", float("nan"))),
            ek_combined_pvalue=float(getattr(r, "combined_pvalue", float("nan"))),
            ek_dmr_type=str(getattr(r, "dmr_type", "")),
            jaccard=round(best_j, 4),
            dss_idx=int(best_i),
            dss_chrom=str(d["chrom"]),
            dss_start=int(d["start"]),
            dss_end=int(d["end"]),
            dss_length=int(d["length"]),
            dss_n_cpgs=int(d.get("n_cpgs", d.get("nCG", 0))),
            dss_areaStat=float(d["areaStat"]),
            dss_diff_Methy_DSSfit=float(d.get("diff_Methy_DSSfit", float("nan"))),
            dss_diff_Methy_fromCounts=float(d.get("diff_Methy_fromCounts",
                                                   float("nan"))),
            dss_dmr_type=str(d.get("dmr_type", "")),
        ))
    return pd.DataFrame(rows)


def load_ek(parquet_path: Path, label: str) -> pd.DataFrame:
    df = pl.read_parquet(parquet_path).to_pandas()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["idx"]   = np.arange(len(df))
    df["caller"] = label
    return df


def load_dss(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["idx"]   = np.arange(len(df))
    return df


def summary_stats(df: pd.DataFrame, label: str) -> dict:
    if len(df) == 0:
        return {f"{label}_n_matched": 0}
    # Effect-size agreement
    r_ds_fit  = sps.pearsonr(df["ek_mean_meth_diff"],
                              df["dss_diff_Methy_DSSfit"]).statistic
    r_ds_from = sps.pearsonr(df["ek_mean_meth_diff"],
                              df["dss_diff_Methy_fromCounts"]).statistic
    rho_ds_fit  = sps.spearmanr(df["ek_mean_meth_diff"],
                                 df["dss_diff_Methy_DSSfit"]).statistic
    # Significance rank concordance: epykit -log10(q) vs DSS |areaStat|
    ek_sig = -np.log10(df["ek_combined_qvalue"].clip(lower=1e-300))
    dss_sig = df["dss_areaStat"].abs()
    rho_sig = sps.spearmanr(ek_sig, dss_sig).statistic
    # Direction agreement
    dir_agree = int(((df["ek_mean_meth_diff"].fillna(0)
                       * df["dss_diff_Methy_fromCounts"].fillna(0)) > 0).sum())
    return {
        f"{label}_n_matched": int(len(df)),
        f"{label}_pearson_dB_ek_vs_DSSfit": round(float(r_ds_fit), 4),
        f"{label}_pearson_dB_ek_vs_DSSfromCounts": round(float(r_ds_from), 4),
        f"{label}_spearman_dB_ek_vs_DSSfit": round(float(rho_ds_fit), 4),
        f"{label}_spearman_eksig_vs_DSSareaStat": round(float(rho_sig), 4),
        f"{label}_direction_agree_n": dir_agree,
        f"{label}_direction_agree_frac":
            round(dir_agree / max(len(df), 1), 4),
        f"{label}_jaccard_p25_p50_p75":
            [round(float(df["jaccard"].quantile(q)), 4)
             for q in (0.25, 0.5, 0.75)],
    }


def main() -> None:
    print("Loading inputs …")
    ek100 = load_ek(CM_DIR / "dmr_chain_merge.parquet", "ek100")
    ek250 = load_ek(DM250_PQ, "ek250")
    dss = load_dss(DSS_CSV)
    print(f"  ek100={len(ek100)}, ek250={len(ek250)}, dss={len(dss)}")

    print("\nMatching ek100 vs DSS …")
    df100 = best_match_pairs(ek100, dss, min_j=0.0)
    df100["caller"] = "ek100"
    print(f"  pairs: {len(df100)}")

    print("Matching ek250 vs DSS …")
    df250 = best_match_pairs(ek250, dss, min_j=0.0)
    df250["caller"] = "ek250"
    print(f"  pairs: {len(df250)}")

    combined = pd.concat([df100, df250], ignore_index=True)
    out_csv = OUT_DIR / "per_dmr_stat_concordance.csv"
    combined.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv.name}  ({len(combined)} rows)")

    # Summary stats at multiple Jaccard cutoffs
    summary = {}
    for jcut, jlab in [(0.0, "anybp"), (0.25, "J025"), (0.5, "J05"),
                       (0.75, "J075")]:
        sub100 = df100[df100["jaccard"] >= jcut]
        sub250 = df250[df250["jaccard"] >= jcut]
        summary[f"ek100_{jlab}"] = summary_stats(sub100, f"ek100_{jlab}")
        summary[f"ek250_{jlab}"] = summary_stats(sub250, f"ek250_{jlab}")
    out_json = OUT_DIR / "per_dmr_stat_concordance_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {out_json.name}")

    print("\n--- SUMMARY (J>=0.5) ---")
    for k in (f"ek100_J05", f"ek250_J05"):
        print(f"  {k}:")
        print("    " + json.dumps(summary[k], indent=2).replace("\n", "\n    "))


if __name__ == "__main__":
    main()
