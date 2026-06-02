"""F9 — Per-DMR Δβ scatter: epykit chain_merge vs DSS, for matched pairs.

Two panels:
  A. epykit (ek-100) vs DSS Δβ, J ≥ 0.25 matched pairs only,
     colored by Jaccard quantile, identity line.
  B. ek-250 vs DSS, same.
  Inset / title bar: Pearson r, Spearman ρ, 100% direction agreement.

Source: comparisons/per_dmr_stat_concordance.csv from P1.2.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, setup, save_dual)


def main() -> None:
    setup()
    df = pd.read_csv(DATA_DIR / "comparisons" / "per_dmr_stat_concordance.csv")
    df = df[df["jaccard"] >= 0.0]   # already filtered min_j=0.0

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), constrained_layout=True)

    for ax, caller in zip(axes, ("ek100", "ek250")):
        sub = df[(df["caller"] == caller) & (df["jaccard"] >= 0.25)]
        ek = sub["ek_mean_meth_diff"].values
        ds = sub["dss_diff_Methy_fromCounts"].values
        jq = sub["jaccard"].values
        scatter = ax.scatter(ek, ds, c=jq, cmap="viridis",
                              vmin=0.25, vmax=1.0, s=22, alpha=0.85,
                              edgecolor="white", linewidth=0.4)
        ax.plot([-1, 1], [-1, 1], "k--", lw=0.8, alpha=0.5)
        ax.axhline(0, color="#bdc3c7", lw=0.6, alpha=0.7)
        ax.axvline(0, color="#bdc3c7", lw=0.6, alpha=0.7)
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1)
        ax.set_xlabel("epykit chain_merge  meth_diff (KO − WT)")
        ax.set_ylabel("DSS  diff.Methy_fromCounts (KO − WT)")
        if len(sub) > 1:
            r = sps.pearsonr(ek, ds).statistic
            rho = sps.spearmanr(ek, ds).statistic
            dir_agree = ((ek * ds) > 0).mean() * 100
            title = (f"{caller.replace('ek', 'ek-chain_merge-')} vs DSS\n"
                      f"n = {len(sub)}  ·  Pearson r = {r:.4f}  ·  "
                      f"Spearman ρ = {rho:.3f}  ·  "
                      f"direction agreement = {dir_agree:.0f}%")
        else:
            title = caller
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)

    cax = fig.add_axes([0.30, -0.04, 0.4, 0.02])
    cb = fig.colorbar(scatter, cax=cax, orientation="horizontal")
    cb.set_label("Jaccard of matched-pair overlap", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.suptitle("F9 · Per-DMR effect-size concordance (epykit chain_merge vs DSS)",
                 fontsize=12, y=1.04)
    save_dual(fig, THREE_WAY / "F9_per_dmr_concordance")
    plt.close(fig)


if __name__ == "__main__":
    main()
