"""F1 — 3-way DMR overlap, hand-rolled UpSet-style figure.

Sets: paper-Table5 (813), epykit-chain_merge-100 (702),
DSS-from-scratch (922). For each DMR we record which sets it belongs
to (via any-bp adjacency OR J>=0.5). Two panels: any-bp and J>=0.5.

Custom render avoids the upsetplot/matplotlib 3.10 incompatibility.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, setup, save_dual)

PAPER_T5 = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/"
                r"Paper resources/DMR_total_list.xlsx")

SOURCE_ORDER = ["paper", "epykit", "DSS"]
SOURCE_COLORS = {"paper": "#2c3e50", "epykit": "#3498db", "DSS": "#27ae60"}
SOURCE_TOTALS = {}  # filled later: total count per source


def load_intervals(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["source"] = label
    df["idx"] = np.arange(len(df))
    return df[["source", "chrom", "start", "end", "idx"]]


def jaccard(a_s, a_e, b_s, b_e) -> float:
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def memberships(p: pd.DataFrame, e: pd.DataFrame, d: pd.DataFrame,
                min_j: float = 0.0) -> dict[tuple[str, ...], int]:
    """For each DMR (from any source), determine which sources it shares
    overlap with. Returns Counter keyed by sorted membership tuple.

    Each membership tuple is counted once per *source-instance* DMR — i.e.,
    a paper DMR matched to one epykit DMR counts once under (paper,epykit)
    and that same epykit DMR's pair also contributes once under
    (paper,epykit). Net effect: paired regions are over-counted relative
    to "unique regions" but UpSet's standard interpretation is per-source-
    row, which is what's plotted.
    """
    sources = {"paper": p, "epykit": e, "DSS": d}
    idx = {}
    for s, df in sources.items():
        by_ch = defaultdict(list)
        for _, r in df.iterrows():
            by_ch[r["chrom"]].append((r["start"], r["end"]))
        for ch in by_ch:
            by_ch[ch].sort(key=lambda t: t[0])
        idx[s] = by_ch

    counts: Counter = Counter()
    for s, df in sources.items():
        for _, r in df.iterrows():
            ch = r["chrom"]; rs = r["start"]; re_ = r["end"]
            sets_in = {s}
            for other in SOURCE_ORDER:
                if other == s:
                    continue
                for ts, te in idx[other].get(ch, []):
                    if te < rs:
                        continue
                    if ts > re_:
                        break
                    if min_j == 0.0:
                        if min(re_, te) > max(rs, ts):
                            sets_in.add(other); break
                    else:
                        if jaccard(rs, re_, ts, te) >= min_j:
                            sets_in.add(other); break
            counts[tuple(sorted(sets_in))] += 1
    return counts


def draw_upset(counts: Counter, ax_bar, ax_matrix, ax_left, title: str) -> None:
    """Custom UpSet-style figure with three axes:
      ax_bar:    top intersection bar chart
      ax_matrix: bottom dot matrix
      ax_left:   left set-size bar chart
    """
    # Sort intersections by descending count
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    keys = [k for k, _ in items]
    vals = [v for _, v in items]

    # Bar chart of intersection sizes
    xs = np.arange(len(keys))
    bars = ax_bar.bar(xs, vals, color="#34495e", edgecolor="black",
                       linewidth=0.6, width=0.7)
    for b, v in zip(bars, vals):
        ax_bar.text(b.get_x() + b.get_width() / 2, v * 1.02,
                     f"{v:,}", ha="center", va="bottom", fontsize=8.5)
    ax_bar.set_xticks(xs); ax_bar.set_xticklabels([])
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_ylim(0, max(vals) * 1.15)
    ax_bar.set_title(title, fontsize=11, pad=10)
    ax_bar.grid(axis="y", alpha=0.2)
    for spine in ("right", "top"):
        ax_bar.spines[spine].set_visible(False)

    # Dot matrix
    set_y = {s: i for i, s in enumerate(SOURCE_ORDER)}
    ax_matrix.set_xlim(-0.5, len(keys) - 0.5)
    ax_matrix.set_ylim(-0.5, len(SOURCE_ORDER) - 0.5)
    ax_matrix.invert_yaxis()
    for j, k in enumerate(keys):
        present = set(k)
        for s, y in set_y.items():
            color = SOURCE_COLORS[s] if s in present else "#d0d3d4"
            ax_matrix.scatter([j], [y], s=80, color=color, zorder=3,
                               edgecolor="black", linewidth=0.5)
        # Connect dots in this intersection with a line
        present_y = sorted(set_y[s] for s in present)
        if len(present_y) > 1:
            ax_matrix.plot([j, j], [present_y[0], present_y[-1]],
                            color="#2c3e50", lw=1.5, zorder=2)
    ax_matrix.set_yticks(list(set_y.values()))
    ax_matrix.set_yticklabels(
        [f"{s} ({SOURCE_TOTALS[s]:,})" for s in SOURCE_ORDER],
        fontsize=9)
    ax_matrix.set_xticks([])
    ax_matrix.set_xticklabels([])
    for spine in ("right", "top", "bottom"):
        ax_matrix.spines[spine].set_visible(False)
    ax_matrix.tick_params(left=False)

    # Left set-size bar
    set_totals = [SOURCE_TOTALS[s] for s in SOURCE_ORDER]
    ax_left.barh(list(set_y.values()), set_totals,
                  color=[SOURCE_COLORS[s] for s in SOURCE_ORDER],
                  edgecolor="black", linewidth=0.6, height=0.6)
    ax_left.invert_xaxis()
    ax_left.invert_yaxis()
    ax_left.set_yticks([]); ax_left.set_xlabel("Set size")
    for spine in ("top", "right", "left"):
        ax_left.spines[spine].set_visible(False)
    ax_left.tick_params(left=False)
    for i, v in enumerate(set_totals):
        ax_left.text(v * 0.5, i, f"{v:,}", ha="center", va="center",
                       fontsize=8.5, color="white", fontweight="bold")


def main() -> None:
    setup()

    paper = pd.read_excel(PAPER_T5).rename(columns={"chr": "chrom"})
    p = load_intervals(paper, "paper")
    e = load_intervals(
        pl.read_parquet(DATA_DIR / "chain_merge" / "dmr_chain_merge.parquet")
          .to_pandas(),
        "epykit")
    d = load_intervals(
        pd.read_csv(DATA_DIR / "dss" / "dmr_dss.csv"),
        "DSS")
    SOURCE_TOTALS.update({"paper": len(p), "epykit": len(e), "DSS": len(d)})

    cnt_any = memberships(p, e, d, min_j=0.0)
    cnt_j5  = memberships(p, e, d, min_j=0.5)

    # Save the underlying counts
    def to_df(c: Counter, lbl: str) -> pd.DataFrame:
        return pd.DataFrame([
            {"intersection": "+".join(k), "count": v, "criterion": lbl}
            for k, v in sorted(c.items(), key=lambda kv: -kv[1])
        ])
    pd.concat([to_df(cnt_any, "any_bp"), to_df(cnt_j5, "J0.5")]) \
       .to_csv(THREE_WAY / "F1_upset_data.csv", index=False)

    # ---- Two-panel figure ----------------------------------------------
    for cnt, label, fname in [
        (cnt_any, "Any-bp overlap criterion",
         "F1a_upset_any_bp"),
        (cnt_j5,  "Strict overlap criterion (Jaccard ≥ 0.5)",
         "F1b_upset_J05"),
    ]:
        fig = plt.figure(figsize=(8.5, 5.2))
        # Layout: top bar, middle matrix, left set-size
        gs = fig.add_gridspec(
            2, 2,
            width_ratios=[1.2, 4],
            height_ratios=[2.5, 1.3],
            hspace=0.05, wspace=0.05,
        )
        ax_top    = fig.add_subplot(gs[0, 1])
        ax_matrix = fig.add_subplot(gs[1, 1], sharex=ax_top)
        ax_left   = fig.add_subplot(gs[1, 0], sharey=ax_matrix)
        ax_empty  = fig.add_subplot(gs[0, 0]); ax_empty.axis("off")
        draw_upset(cnt, ax_top, ax_matrix, ax_left,
                    f"F1 · 3-way DMR overlap — {label}")
        fig.subplots_adjust(top=0.91)
        save_dual(fig, THREE_WAY / fname)
        plt.close(fig)


if __name__ == "__main__":
    main()
