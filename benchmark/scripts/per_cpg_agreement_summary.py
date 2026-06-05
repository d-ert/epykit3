#!/usr/bin/env python
"""Per-CpG agreement summary: epykit vs methylKit on the Study-3 (GSE263850)
genome-wide per-CpG DMC output.

Backs paper.md §3.3.1 / REPORT.md §3.5 (Pearson r, Spearman rho, direction
agreement on the ~15.6 M shared CpGs). The raw per-CpG frames are large
(~1 GB CSV each) and live on the analysis host, not in the repo; this script
distils them to a compact, committable artifact:

  * a JSON with the headline scalars (n_shared, pearson_r, spearman_rho,
    direction-agreement count + fraction), and
  * a 2-D histogram parquet (binned epykit vs methylKit meth_diff) so the
    concordance figure can be regenerated without the 15.6 M-row inputs.

Both inputs are `dmc_all_sites.csv` with columns
`chrom,pos,end,strand,pvalue,qvalue,meth_diff,mean_beta_control,mean_beta_case,meth_diff_prop`.
`meth_diff` is on the percent scale (-100..100); `meth_diff_prop` is fractional
(-1..1). epykit `pos` is 0-based; methylKit `pos` is 1-based, so the join is
`epykit.pos + 1 == methylKit.pos` on the same chromosome (Methods §2.2/§4.3).

Example (paths are on the analysis host, outside the repo):

    uv run python benchmark/scripts/per_cpg_agreement_summary.py \
        --epykit-csv ".../epykit_results/dmc_all_sites.csv" \
        --methylkit-csv ".../methylkit_results/dmc_all_sites.csv" \
        --out-json benchmark/data/study3/comparisons/per_cpg_agreement.json \
        --out-hist benchmark/data/study3/comparisons/per_cpg_agreement_hist2d.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import polars as pl


def _load(path: str) -> pl.DataFrame:
    """Lazy-scan a dmc_all_sites.csv, keep only the join keys + effect size."""
    return (
        pl.scan_csv(path)
        .select(
            pl.col("chrom").cast(pl.Utf8),
            pl.col("pos").cast(pl.Int64),
            pl.col("meth_diff_prop").cast(pl.Float64),
        )
        .collect()
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epykit-csv", required=True)
    ap.add_argument("--methylkit-csv", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-hist", required=True)
    ap.add_argument("--bins", type=int, default=100)
    args = ap.parse_args()

    ek = _load(args.epykit_csv).rename({"meth_diff_prop": "ek"})
    mk = _load(args.methylkit_csv).rename({"meth_diff_prop": "mk"})
    n_ek, n_mk = ek.height, mk.height

    # epykit 0-based -> methylKit 1-based; inner join on (chrom, pos+1 == pos).
    ek = ek.with_columns((pl.col("pos") + 1).alias("pos"))
    j = ek.join(mk, on=["chrom", "pos"], how="inner")
    n_shared = j.height

    ek_v = j["ek"].to_numpy()
    mk_v = j["mk"].to_numpy()

    pearson = float(np.corrcoef(ek_v, mk_v)[0, 1])
    # Spearman = Pearson on ranks (average ranks for ties).
    from scipy.stats import rankdata  # local import keeps the dep optional

    spearman = float(np.corrcoef(rankdata(ek_v), rankdata(mk_v))[0, 1])

    sign_ek = np.sign(ek_v)
    sign_mk = np.sign(mk_v)
    agree = int(np.sum(sign_ek == sign_mk))
    direction_frac = agree / n_shared if n_shared else float("nan")

    summary = {
        "n_epykit_tested": n_ek,
        "n_methylkit_tested": n_mk,
        "n_shared": n_shared,
        "pearson_r_meth_diff": pearson,
        "spearman_rho_meth_diff": spearman,
        "direction_agree_count": agree,
        "direction_agree_frac": direction_frac,
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # 2-D histogram over the fractional effect-size plane for figure regen.
    edges = np.linspace(-1.0, 1.0, args.bins + 1)
    counts, _, _ = np.histogram2d(ek_v, mk_v, bins=[edges, edges])
    centers = 0.5 * (edges[:-1] + edges[1:])
    ek_c, mk_c = np.meshgrid(centers, centers, indexing="ij")
    hist = pl.DataFrame(
        {
            "epykit_meth_diff_bin": ek_c.ravel(),
            "methylkit_meth_diff_bin": mk_c.ravel(),
            "count": counts.ravel().astype(np.int64),
        }
    ).filter(pl.col("count") > 0)
    Path(args.out_hist).parent.mkdir(parents=True, exist_ok=True)
    hist.write_parquet(args.out_hist)
    print(f"wrote {args.out_hist} ({hist.height} occupied bins)")


if __name__ == "__main__":
    main()
