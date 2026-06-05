"""summarize_tabix_vs_epykit.py -- aggregate the tabix-vs-epykit per-run CSVs.

Reads the per-run CSV(s) produced by run_tabix_vs_epykit.py (simulator: many
seeds; gse263850: one whole-genome cell), reduces to a median per
(dataset, tool, cores), and derives the comparison ratios that answer the
question "does epykit beat methylKit's tabix backend?":
    - speedup_vs_epykit       : methylKit wall_s / epykit wall_s
    - rss_ratio_vs_epykit     : methylKit rss_peak / epykit rss_peak
    - uss_ratio_vs_epykit     : methylKit uss_peak / epykit uss_peak (fair for forks)
    - tabix_vs_ram_rss        : tabix rss_peak / ram rss_peak at the same cores
    - tabix_vs_ram_wall       : tabix wall_s  / ram wall_s  at the same cores

Writes <out-dir>/tabix_vs_epykit_summary.csv and tabix_vs_epykit_summary.md.

Usage:
    python summarize_tabix_vs_epykit.py \
        --per-run A.csv [B.csv ...] --out-dir <summaries dir>
"""

from __future__ import annotations

import argparse
import statistics as st
from pathlib import Path

import polars as pl

ORDER = {"epykit_lr": 0, "methylkit_ram": 1, "methylkit_tabix": 2}


def _median_rows(df: pl.DataFrame) -> pl.DataFrame:
    keep_ok = df.filter(pl.col("returncode") == 0)
    g = (
        keep_ok.group_by(["dataset", "tool", "backend", "cores"])
        .agg(
            n_runs=pl.len(),
            wall_s=pl.col("wall_s").median(),
            cpu_s=pl.col("cpu_s").median(),
            rss_peak_mb=pl.col("rss_peak_mb").median(),
            uss_peak_mb=pl.col("uss_peak_mb").median(),
            num_processes_peak=pl.col("num_processes_peak").max(),
            n_sites=pl.col("n_sites").median(),
        )
    )
    # any failed configs -> surface them with returncode
    failed = df.filter(pl.col("returncode") != 0).select(
        ["dataset", "tool", "backend", "cores", "returncode"]
    )
    return g, failed


def _derive(g: pl.DataFrame) -> pl.DataFrame:
    rows = g.to_dicts()
    # index epykit (single, cores=1) per dataset; and ram-by-cores per dataset
    ep = {}
    ram = {}
    for r in rows:
        if r["tool"] == "epykit_lr":
            ep[r["dataset"]] = r
        if r["tool"] == "methylkit_ram":
            ram[(r["dataset"], r["cores"])] = r
    out = []
    for r in rows:
        d = r["dataset"]
        e = ep.get(d)
        rr = {**r}
        if e and e["wall_s"]:
            rr["speedup_vs_epykit"] = round(r["wall_s"] / e["wall_s"], 1) if r["tool"] != "epykit_lr" else 1.0
        if e and e["rss_peak_mb"]:
            rr["rss_ratio_vs_epykit"] = round(r["rss_peak_mb"] / e["rss_peak_mb"], 2) if r["tool"] != "epykit_lr" else 1.0
        if e and e["uss_peak_mb"]:
            rr["uss_ratio_vs_epykit"] = round(r["uss_peak_mb"] / e["uss_peak_mb"], 2) if r["tool"] != "epykit_lr" else 1.0
        if r["tool"] == "methylkit_tabix":
            base = ram.get((d, r["cores"]))
            if base:
                if base["rss_peak_mb"]:
                    rr["tabix_vs_ram_rss"] = round(r["rss_peak_mb"] / base["rss_peak_mb"], 2)
                if base["uss_peak_mb"]:
                    rr["tabix_vs_ram_uss"] = round(r["uss_peak_mb"] / base["uss_peak_mb"], 2)
                if base["wall_s"]:
                    rr["tabix_vs_ram_wall"] = round(r["wall_s"] / base["wall_s"], 2)
        out.append(rr)
    out.sort(key=lambda r: (r["dataset"], r["cores"], ORDER.get(r["tool"], 9)))
    return pl.DataFrame(out)


def _md_table(df: pl.DataFrame) -> str:
    lines = []
    for dataset in df["dataset"].unique(maintain_order=True).to_list():
        sub = df.filter(pl.col("dataset") == dataset)
        lines.append(f"\n### {dataset}\n")
        lines.append("| tool | backend | cores | wall_s | RSS MB | USS MB | n_sites | "
                     "speedup vs epykit | RSS× vs epykit | USS× vs epykit | tabix/ram RSS | tabix/ram wall |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in sub.to_dicts():
            def g(k, suf=""):
                v = r.get(k)
                return f"{v}{suf}" if v is not None else "—"
            lines.append(
                f"| {r['tool']} | {r['backend']} | {r['cores']} | {g('wall_s')} | "
                f"{g('rss_peak_mb')} | {g('uss_peak_mb')} | {int(r['n_sites']) if r.get('n_sites') is not None else '—'} | "
                f"{g('speedup_vs_epykit','×')} | {g('rss_ratio_vs_epykit','×')} | {g('uss_ratio_vs_epykit','×')} | "
                f"{g('tabix_vs_ram_rss','×')} | {g('tabix_vs_ram_wall','×')} |"
            )
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-run", nargs="+", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    frames = [pl.read_csv(p) for p in args.per_run if p.exists()]
    if not frames:
        raise SystemExit("no per-run CSVs found")
    df = pl.concat(frames, how="vertical_relaxed")
    g, failed = _median_rows(df)
    summary = _derive(g)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "tabix_vs_epykit_summary.csv"
    summary.write_csv(csv_path)

    md = ["# methylKit tabix vs RAM vs epykit lr — summary",
          "",
          "Out-of-core head-to-head. tabix changes only storage/IO, not statistics "
          "(per-CpG p/q identical to RAM), so the axes are wall time and peak memory. "
          "USS (unique set size) is the fair memory metric under `mc.cores>1` because "
          "summed RSS double-counts copy-on-write shared pages across forked workers.",
          _md_table(summary)]
    if failed.height:
        md.append("\n### failed / OOM configs\n")
        for r in failed.to_dicts():
            md.append(f"- {r['dataset']} {r['tool']} cores={r['cores']}: returncode={r['returncode']}")
    md_path = args.out_dir / "tabix_vs_epykit_summary.md"
    md_path.write_text("\n".join(md) + "\n")

    print(f"wrote {csv_path}\nwrote {md_path}")
    print(summary.select([c for c in ("dataset","tool","cores","wall_s","rss_peak_mb",
          "uss_peak_mb","speedup_vs_epykit","uss_ratio_vs_epykit","tabix_vs_ram_uss")
          if c in summary.columns]))
    return 0


if __name__ == "__main__":
    sys.exit(main()) if (sys := __import__("sys")) else None
