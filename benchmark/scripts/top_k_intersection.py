"""Top-K intersection of epykit vs methylKit calls on real GSE263850 data.

Computes empirical overlap at K = 5, 10, 25, 50, 100 for both DMCs and DMRs,
annotates each hit with the nearest gene (from epykit's tss_distance_* files),
and writes a markdown report to FINAL_REPORT/top_k_report.md.

Coordinate convention: methylKit .cov is 1-based; epykit BED is 0-based, so
epykit positions are shifted by +1 before comparison.
"""

from __future__ import annotations

import os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"
MK_SOURCE = ROOT.parent.parent / "methylkıt_realResults" / "scripts_and_results" / "methylkit_results"

OUT = ROOT / "top_k_report.md"


def load_calls(path: Path, kind: str) -> pd.DataFrame:
    """kind in {'dmc', 'dmr'}. Returns sorted by q-value ascending."""
    df = pd.read_csv(path)
    df = df.sort_values("qvalue", ascending=True).reset_index(drop=True)
    return df


def normalise_dmc_key(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Build a comparable (chrom, pos_1based) key.

    methylKit positions are 1-based; epykit are 0-based. Shift epykit by +1.
    """
    df = df.copy()
    if source == "epykit":
        df["pos_1based"] = df["pos"].astype(int) + 1
    else:
        df["pos_1based"] = df["pos"].astype(int)
    df["key"] = df["chrom"].astype(str) + ":" + df["pos_1based"].astype(str)
    return df


def normalise_dmr_key(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Build a comparable tile key. methylKit tiles are 1-based, epykit 0-based."""
    df = df.copy()
    if source == "epykit":
        df["start_1based"] = df["start"].astype(int) + 1
    else:
        df["start_1based"] = df["start"].astype(int)
    df["key"] = (
        df["chrom"].astype(str)
        + ":"
        + df["start_1based"].astype(str)
        + "-"
        + df["end"].astype(int).astype(str)
    )
    return df


def top_k_intersection(mk: pd.DataFrame, ep: pd.DataFrame, K: int):
    mk_top = set(mk.head(K)["key"])
    ep_top = set(ep.head(K)["key"])
    inter = mk_top & ep_top
    return len(inter), len(mk_top), len(ep_top), inter


def load_gene_annotation() -> dict[str, str]:
    """Return chrom:pos_1based -> gene_name dictionary from epykit tss_distance files.

    tss_distance_hyper.csv is indexed by row position in the HYPER subset of
    the significant DMC table; tss_distance_hypo.csv into the HYPO subset.
    Both subsets are derived from dmc_significant_qval05.csv in its native
    file order, then split by sign of meth_diff.
    """
    annot = {}
    sig_dmc = pd.read_csv(SOURCE / "epykit_results" / "dmc_significant_qval05.csv")

    # Recreate the hyper/hypo split that epykit's pipeline uses (sign of meth_diff)
    hyper_subset = sig_dmc[sig_dmc["meth_diff"] > 0].reset_index(drop=True)
    hypo_subset = sig_dmc[sig_dmc["meth_diff"] < 0].reset_index(drop=True)

    for fn, subset in [
        ("tss_distance_hyper.csv", hyper_subset),
        ("tss_distance_hypo.csv", hypo_subset),
    ]:
        p = SOURCE / "epykit_results" / fn
        if not p.exists():
            continue
        t = pd.read_csv(p)
        if len(t) != len(subset):
            print(
                f"  WARNING: {fn} has {len(t):,} rows but subset has "
                f"{len(subset):,}; alignment may be off"
            )
        for i, row in t.iterrows():
            ridx = int(row["target.row"]) - 1
            if 0 <= ridx < len(subset):
                site = subset.iloc[ridx]
                key = f"{site['chrom']}:{int(site['pos']) + 1}"
                gene = row["feature.name"]
                if pd.isna(gene) or str(gene) == "nan":
                    continue
                annot[key] = f"{gene} ({row['dist.to.feature']:+d} bp)"
    return annot


def annotate_dmrs(ep_top: pd.DataFrame, annot: dict[str, str]) -> list[str]:
    """For each DMR row, find the nearest annotated DMC inside that tile."""
    out = []
    for _, row in ep_top.iterrows():
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        # search annot for a key in this tile
        hits = []
        for key, gene in annot.items():
            ch, pos = key.split(":")
            if ch == chrom and start + 1 <= int(pos) <= end:
                hits.append(gene)
        if hits:
            # show the first 2 unique genes
            uniq = []
            for h in hits:
                gn = h.split(" ")[0]
                if gn not in [u.split(" ")[0] for u in uniq]:
                    uniq.append(h)
                if len(uniq) >= 2:
                    break
            out.append(", ".join(uniq))
        else:
            out.append("(no DMC-level annotation in tile)")
    return out


def fmt_p(p: float) -> str:
    if p < 1e-10:
        return f"{p:.2e}"
    elif p < 0.001:
        return f"{p:.3e}"
    else:
        return f"{p:.4f}"


def main():
    lines = []
    lines.append("# Top-K intersection — epykit vs methylKit on GSE263850\n")
    lines.append(
        "Computed empirically from the real call files. Coordinates aligned "
        "(epykit 0-based + 1 → methylKit 1-based).\n"
    )

    # === DMCs ============================================================
    print("Loading DMC files...")
    ep_dmc = normalise_dmc_key(
        load_calls(ROOT / "data" / "study3" / "dmc_significant_qval05.csv", "dmc"),
        "epykit",
    )
    mk_dmc = normalise_dmc_key(
        load_calls(MK_SOURCE / "dmc_significant_qval05.csv", "dmc"),
        "methylkit",
    )
    print(f"  epykit DMCs: {len(ep_dmc):,}")
    print(f"  methylKit DMCs: {len(mk_dmc):,}")

    lines.append("## DMC intersection by top-K\n")
    lines.append(
        f"- epykit significant DMCs: **{len(ep_dmc):,}**\n"
        f"- methylKit significant DMCs: **{len(mk_dmc):,}**\n"
    )
    lines.append("\n| K  | epykit ∩ methylKit | recall (∩/K) |")
    lines.append("|---:|---:|---:|")
    for K in (5, 10, 25, 50, 100, 250, 500):
        n_inter, _, _, _ = top_k_intersection(mk_dmc, ep_dmc, K)
        lines.append(f"| {K} | {n_inter} | {n_inter / K * 100:.1f}% |")
    lines.append("")

    # Gene annotation
    print("Loading gene annotation...")
    annot = load_gene_annotation()
    print(f"  annotations loaded: {len(annot):,}")

    # === Top 10 DMC table ================================================
    K = 10
    mk_top10 = mk_dmc.head(K).copy()
    ep_top10 = ep_dmc.head(K).copy()
    mk_keys = set(mk_top10["key"])
    ep_keys = set(ep_top10["key"])

    lines.append("## Top-10 DMCs by q-value — methylKit\n")
    lines.append(
        "| # | chrom:pos | gene (Δbp from TSS) | meth_diff | q-value | also in epykit top-10? |"
    )
    lines.append("|---:|---|---|---:|---:|:---:|")
    for i, (_, row) in enumerate(mk_top10.iterrows(), 1):
        gene = annot.get(row["key"], "(not in epykit sig set)")
        also = "✓" if row["key"] in ep_keys else " "
        lines.append(
            f"| {i} | {row['chrom']}:{row['pos_1based']:,} | "
            f"{gene} | {row['meth_diff']:+.1f}% | "
            f"{fmt_p(row['qvalue'])} | {also} |"
        )
    lines.append("")

    lines.append("## Top-10 DMCs by q-value — epykit\n")
    lines.append(
        "| # | chrom:pos | gene (Δbp from TSS) | meth_diff | q-value | also in methylKit top-10? |"
    )
    lines.append("|---:|---|---|---:|---:|:---:|")
    for i, (_, row) in enumerate(ep_top10.iterrows(), 1):
        gene = annot.get(row["key"], "(no annotation)")
        also = "✓" if row["key"] in mk_keys else " "
        lines.append(
            f"| {i} | {row['chrom']}:{row['pos_1based']:,} | "
            f"{gene} | {row['meth_diff']:+.1f}% | "
            f"{fmt_p(row['qvalue'])} | {also} |"
        )
    lines.append("")

    # === DMRs ============================================================
    print("Loading DMR files...")
    ep_dmr = normalise_dmr_key(
        load_calls(ROOT / "data" / "study3" / "dmr_significant_lenient.csv", "dmr"),
        "epykit",
    )
    mk_dmr = normalise_dmr_key(
        load_calls(MK_SOURCE / "dmr_significant_lenient.csv", "dmr"),
        "methylkit",
    )
    print(f"  epykit DMRs (lenient): {len(ep_dmr):,}")
    print(f"  methylKit DMRs (lenient): {len(mk_dmr):,}")

    lines.append("## DMR intersection by top-K (lenient threshold)\n")
    lines.append(
        f"- epykit significant DMRs: **{len(ep_dmr):,}**\n"
        f"- methylKit significant DMRs: **{len(mk_dmr):,}**\n"
    )
    lines.append("\n| K  | epykit ∩ methylKit | recall (∩/K) |")
    lines.append("|---:|---:|---:|")
    for K in (5, 10, 25, 50, 100, 250, 500):
        n_inter, _, _, _ = top_k_intersection(mk_dmr, ep_dmr, K)
        lines.append(f"| {K} | {n_inter} | {n_inter / K * 100:.1f}% |")
    lines.append("")

    # === Top 10 DMR tables ===============================================
    K = 10
    mk_top10 = mk_dmr.head(K).copy()
    ep_top10 = ep_dmr.head(K).copy()
    mk_keys = set(mk_top10["key"])
    ep_keys = set(ep_top10["key"])

    # Annotate each DMR with nearby gene names
    mk_genes = annotate_dmrs(mk_top10, annot)
    ep_genes = annotate_dmrs(ep_top10, annot)

    lines.append("## Top-10 DMRs by q-value — methylKit\n")
    lines.append(
        "| # | chrom:start-end | gene(s) in tile | meth_diff | q-value | also in epykit top-10? |"
    )
    lines.append("|---:|---|---|---:|---:|:---:|")
    for i, ((_, row), gene) in enumerate(zip(mk_top10.iterrows(), mk_genes), 1):
        also = "✓" if row["key"] in ep_keys else " "
        lines.append(
            f"| {i} | {row['chrom']}:{int(row['start_1based']):,}–{int(row['end']):,} | "
            f"{gene} | {row['meth_diff']:+.1f}% | "
            f"{fmt_p(row['qvalue'])} | {also} |"
        )
    lines.append("")

    lines.append("## Top-10 DMRs by q-value — epykit\n")
    lines.append(
        "| # | chrom:start-end | gene(s) in tile | meth_diff | q-value | also in methylKit top-10? |"
    )
    lines.append("|---:|---|---|---:|---:|:---:|")
    for i, ((_, row), gene) in enumerate(zip(ep_top10.iterrows(), ep_genes), 1):
        also = "✓" if row["key"] in mk_keys else " "
        lines.append(
            f"| {i} | {row['chrom']}:{int(row['start_1based']):,}–{int(row['end']):,} | "
            f"{gene} | {row['meth_diff']:+.1f}% | "
            f"{fmt_p(row['qvalue'])} | {also} |"
        )
    lines.append("")

    # === Summary =========================================================
    n5_dmc, _, _, _ = top_k_intersection(mk_dmc, ep_dmc, 5)
    n10_dmc, _, _, _ = top_k_intersection(mk_dmc, ep_dmc, 10)
    n5_dmr, _, _, _ = top_k_intersection(mk_dmr, ep_dmr, 5)
    n10_dmr, _, _, _ = top_k_intersection(mk_dmr, ep_dmr, 10)

    lines.append("## Summary\n")
    lines.append(
        f"- **Top 5 DMC overlap:** {n5_dmc} / 5 = {n5_dmc / 5 * 100:.0f}%\n"
        f"- **Top 10 DMC overlap:** {n10_dmc} / 10 = {n10_dmc / 10 * 100:.0f}%\n"
        f"- **Top 5 DMR overlap:** {n5_dmr} / 5 = {n5_dmr / 5 * 100:.0f}%\n"
        f"- **Top 10 DMR overlap:** {n10_dmr} / 10 = {n10_dmr / 10 * 100:.0f}%\n"
    )

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {OUT}")
    print(f"  Top-5 DMC overlap:  {n5_dmc}/5")
    print(f"  Top-10 DMC overlap: {n10_dmc}/10")
    print(f"  Top-5 DMR overlap:  {n5_dmr}/5")
    print(f"  Top-10 DMR overlap: {n10_dmr}/10")


if __name__ == "__main__":
    main()
