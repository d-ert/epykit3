"""Profiled driver for run_dss_replication.R.

Launches Rscript with the paper-matched parameters, polls the R process
(and its children) every second via psutil, and writes:

    resources.csv      - per-1s time series: timestamp, elapsed_s, rss_mb,
                         uss_mb, cpu_percent, num_threads, num_children
    resources.json     - peak / mean / totals + system info
    parameters.json    - run config + DSS<->paper mapping
    summary.md         - headline numbers + paper-coord comparison
    README.md          - file index

Final outputs live in FINAL_REPORT/data/study3/dss/ alongside the
per-step R-side artifacts (dmr_dss.tsv/.csv, dmr_gene_links_100kb.csv,
step_timings.tsv, run_log.txt, dss_session_info.txt).
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import polars as pl
import psutil


# ---- Paths ------------------------------------------------------------------

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
RAW_DIR   = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

OUT_DIR     = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss"
R_SCRIPT    = REPO_ROOT / "FINAL_REPORT" / "scripts" / "run_dss_replication.R"
SAMPLESHEET = RAW_DIR / "samplesheet.csv"
REFGENE     = RAW_DIR / "refseq" / "refGene.txt.gz"
PAPER_DMR   = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Paper-matched parameters (passed to Rscript as CLI flags)
PARAMS = dict(
    coverage_min=5,
    p_threshold=1e-5,
    delta=0.0,
    minlen=50,
    minCG=3,
    dis_merge=100,
    pct_sig=0.5,
    gene_link_bp=100_000,
)

SAMPLE_INTERVAL_S = 1.0


# ---- Resource sampler ------------------------------------------------------

class ResourceSampler(threading.Thread):
    """Polls a psutil.Process tree for memory + CPU every interval seconds.

    `cpu_percent` is normalized to a single core (so a fully-utilized 8-core
    machine reports up to 800%). RSS is summed across the parent + all
    children; threads is summed too. USS (unique set size) is reported per
    process and summed where available (some processes deny access on
    Windows; falls back to NaN).
    """

    def __init__(self, root_proc: psutil.Process, interval: float = 1.0):
        super().__init__(daemon=True)
        self.root = root_proc
        self.interval = interval
        self.samples: list[dict] = []
        self.stop_flag = threading.Event()
        self.t0 = time.time()
        # Prime cpu_percent
        try:
            self.root.cpu_percent(interval=None)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

    def _walk(self):
        procs = [self.root]
        try:
            procs.extend(self.root.children(recursive=True))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        return procs

    def _sample(self) -> dict | None:
        procs = self._walk()
        if not procs:
            return None
        rss = 0
        uss = 0
        uss_avail = False
        threads = 0
        cpu = 0.0
        n_alive = 0
        for p in procs:
            try:
                with p.oneshot():
                    mi = p.memory_info()
                    rss += mi.rss
                    try:
                        u = p.memory_full_info().uss
                        uss += u
                        uss_avail = True
                    except (psutil.AccessDenied, AttributeError):
                        pass
                    threads += p.num_threads()
                    cpu += p.cpu_percent(interval=None)
                    n_alive += 1
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        if n_alive == 0:
            return None
        return dict(
            timestamp=datetime.now().isoformat(timespec="seconds"),
            elapsed_s=round(time.time() - self.t0, 2),
            rss_mb=round(rss / 1024**2, 2),
            uss_mb=round(uss / 1024**2, 2) if uss_avail else None,
            cpu_percent=round(cpu, 1),
            num_threads=threads,
            num_processes=n_alive,
        )

    def run(self) -> None:
        # First "warm-up" sample is discarded so cpu_percent is meaningful.
        time.sleep(self.interval)
        while not self.stop_flag.is_set():
            s = self._sample()
            if s is None:
                break
            self.samples.append(s)
            time.sleep(self.interval)


# ---- Main ------------------------------------------------------------------

def find_rscript() -> str:
    """Locate Rscript.exe; prefer the absolute path seen on this box."""
    candidates = [
        r"D:\Program Files\R\R-4.5.0\bin\Rscript.exe",
    ]
    for c in candidates:
        if Path(c).exists():
            return c
    # Fall back to PATH
    return "Rscript"


def run_r() -> tuple[int, list[dict], float, float]:
    rscript = find_rscript()
    cmd = [
        rscript,
        "--vanilla",
        str(R_SCRIPT),
        "--samplesheet",   str(SAMPLESHEET),
        "--out-dir",       str(OUT_DIR),
        "--refgene",       str(REFGENE),
        "--coverage-min",  str(PARAMS["coverage_min"]),
        "--p-threshold",   str(PARAMS["p_threshold"]),
        "--delta",         str(PARAMS["delta"]),
        "--minlen",        str(PARAMS["minlen"]),
        "--minCG",         str(PARAMS["minCG"]),
        "--dis-merge",     str(PARAMS["dis_merge"]),
        "--pct-sig",       str(PARAMS["pct_sig"]),
        "--gene-link-bp",  str(PARAMS["gene_link_bp"]),
    ]
    print("Launching:", " ".join(cmd))
    t0_wall = time.time()
    t0_proc = time.process_time()
    proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
    try:
        ps = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        proc.wait()
        return proc.returncode, [], time.time() - t0_wall, 0.0

    sampler = ResourceSampler(ps, SAMPLE_INTERVAL_S)
    sampler.start()
    rc = proc.wait()
    sampler.stop_flag.set()
    sampler.join(timeout=5)
    dt_wall = time.time() - t0_wall
    dt_cpu  = time.process_time() - t0_proc  # parent only; R CPU is in samples
    return rc, sampler.samples, dt_wall, dt_cpu


def summarize_samples(samples: list[dict], wall_seconds: float) -> dict:
    if not samples:
        return dict(samples_collected=0, wall_seconds=round(wall_seconds, 2))
    df = pd.DataFrame(samples)
    rss_peak = float(df["rss_mb"].max())
    rss_mean = float(df["rss_mb"].mean())
    uss_peak = (float(df["uss_mb"].max()) if df["uss_mb"].notna().any()
                else None)
    cpu_peak = float(df["cpu_percent"].max())
    cpu_mean = float(df["cpu_percent"].mean())
    # Effective core-seconds = mean cpu_percent / 100 * wall_seconds
    core_seconds = (df["cpu_percent"].mean() / 100.0) * wall_seconds
    # Average concurrency = core_seconds / wall_seconds
    avg_cores = core_seconds / wall_seconds if wall_seconds > 0 else 0.0
    threads_peak = int(df["num_threads"].max())
    return dict(
        samples_collected=int(len(df)),
        sample_interval_s=SAMPLE_INTERVAL_S,
        wall_seconds=round(wall_seconds, 2),
        rss_peak_mb=round(rss_peak, 1),
        rss_mean_mb=round(rss_mean, 1),
        uss_peak_mb=(round(uss_peak, 1) if uss_peak is not None else None),
        cpu_percent_peak=round(cpu_peak, 1),
        cpu_percent_mean=round(cpu_mean, 1),
        approx_core_seconds=round(core_seconds, 1),
        avg_concurrent_cores=round(avg_cores, 2),
        threads_peak=threads_peak,
    )


def write_resources_csv(samples: list[dict]) -> None:
    if not samples:
        return
    pd.DataFrame(samples).to_csv(OUT_DIR / "resources.csv", index=False)


def write_resources_json(summary: dict) -> None:
    sys_info = dict(
        platform=platform.platform(),
        python=platform.python_version(),
        psutil=psutil.__version__,
        cpu_count_logical=psutil.cpu_count(logical=True),
        cpu_count_physical=psutil.cpu_count(logical=False),
        cpu_freq_mhz=(round(psutil.cpu_freq().max)
                      if psutil.cpu_freq() is not None else None),
        total_ram_gb=round(psutil.virtual_memory().total / 1024**3, 1),
        hostname=platform.node(),
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    payload = dict(system=sys_info, resources=summary)
    (OUT_DIR / "resources.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


def write_parameters_json() -> None:
    payload = {
        "paper_reference": (
            "Farhangdoost et al., 2024. AKAP11 heterozygous LOF in human "
            "iPSC-derived cortical neurons (GSE263850), hg38."
        ),
        "dataset": "GSE263850 (n=6: 3 Het-AKAP11-KO + 3 WT)",
        "ingest": {
            "samplesheet": str(SAMPLESHEET),
            "assembly": "hg38",
            "coverage_min": PARAMS["coverage_min"],
            "unite_mode": "intersect (inner join on chrom+pos across samples)",
        },
        "dmltest": {
            "function": "DSS::DMLfit.multiFactor(smoothing = TRUE) -> "
                        "DMLtest.multiFactor(coef = <group level>)",
            "design": "~ group",
            "smoothing": True,
        },
        "callDMR": {
            "function": "DSS::callDMR",
            "p_threshold": PARAMS["p_threshold"],
            "delta":       PARAMS["delta"],
            "minlen":      PARAMS["minlen"],
            "minCG":       PARAMS["minCG"],
            "dis_merge":   PARAMS["dis_merge"],
            "pct_sig":     PARAMS["pct_sig"],
        },
        "annotation": {
            "tool": "HOMER-equivalent (UCSC refGene) reimplemented in R",
            "catalog": str(REFGENE),
            "paper_tool": "HOMER annotatePeaks.pl (hg38, RefSeq)",
        },
        "gene_linkage_100kb": {
            "rule": (
                "Every gene whose canonical TSS lies within 100 kb of the "
                "DMR midpoint. One row per (DMR, gene)."
            ),
            "max_bp": PARAMS["gene_link_bp"],
        },
        "outputs": [
            "dmr_dss_raw.tsv",
            "dmr_dss.csv", "dmr_dss.tsv",
            "dmr_gene_links_100kb.csv",
            "dmltest_per_cpg.tsv.gz",
            "step_timings.tsv",
            "resources.csv", "resources.json",
            "summary.md", "parameters.json",
            "run_log.txt", "dss_session_info.txt",
            "README.md",
        ],
    }
    (OUT_DIR / "parameters.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )


# ---- Paper comparison ------------------------------------------------------

def _interval_index(df: pd.DataFrame):
    by_ch: dict[str, list[tuple[int, int]]] = {}
    for _, r in df.iterrows():
        by_ch.setdefault(str(r["chrom"]), []).append(
            (int(r["start"]), int(r["end"]))
        )
    for ch in by_ch:
        by_ch[ch].sort()
    return by_ch


def _any_overlap(intervals, s, e) -> bool:
    for ps, pe in intervals:
        if pe < s:
            continue
        if ps > e:
            break
        return True
    return False


def compare_to_paper(dmr_csv: Path) -> dict:
    if not PAPER_DMR.exists():
        return {"paper_table_available": False}
    if not dmr_csv.exists():
        return {"paper_table_available": False, "note": "DMR csv missing"}
    paper = pd.read_excel(PAPER_DMR, sheet_name=0)
    paper_norm = paper[["chr", "start", "end"]].rename(columns={"chr": "chrom"}).copy()
    paper_norm["chrom"] = paper_norm["chrom"].astype(str)
    paper_norm["start"] = paper_norm["start"].astype(int)
    paper_norm["end"]   = paper_norm["end"].astype(int)
    ours = pd.read_csv(dmr_csv)[["chrom", "start", "end"]].copy()
    ours["chrom"] = ours["chrom"].astype(str)
    paper_by = _interval_index(paper_norm)
    ours_by  = _interval_index(ours)

    n_paper = len(paper_norm)
    n_ours  = len(ours)
    ours_overlap = sum(
        _any_overlap(paper_by.get(str(r["chrom"]), []),
                     int(r["start"]), int(r["end"]))
        for _, r in ours.iterrows()
    )
    paper_overlap = sum(
        _any_overlap(ours_by.get(str(r["chrom"]), []),
                     int(r["start"]), int(r["end"]))
        for _, r in paper_norm.iterrows()
    )
    return dict(
        paper_table_available=True,
        paper_n_dmr=n_paper,
        paper_n_hyper=int((paper["diff.meth_mean"] > 0).sum()),
        paper_n_hypo=int((paper["diff.meth_mean"] < 0).sum()),
        paper_median_length_bp=int((paper["end"] - paper["start"]).median()),
        our_n_dmr=n_ours,
        coord_overlap_count=int(ours_overlap),
        coord_recall_of_paper=round(paper_overlap / max(n_paper, 1), 4),
        coord_precision=round(ours_overlap / max(n_ours, 1), 4),
    )


def write_summary(res_summary: dict, paper_stats: dict) -> None:
    dmr_csv = OUT_DIR / "dmr_dss.csv"
    if not dmr_csv.exists():
        return
    dmr = pd.read_csv(dmr_csv)
    lengths = (dmr["end"] - dmr["start"])
    n_hyper = int((dmr["dmr_type"] == "hyper").sum())
    n_hypo  = int((dmr["dmr_type"] == "hypo").sum())
    links = pd.read_csv(OUT_DIR / "dmr_gene_links_100kb.csv")
    step = pd.read_csv(OUT_DIR / "step_timings.tsv", sep="\t")

    lines = []
    lines.append("# DSS replication of GSE263850 (paper-matched parameters)\n")
    lines.append("Per-CpG: `DSS::DMLfit.multiFactor(smoothing = TRUE)` -> "
                 "`DMLtest.multiFactor(coef = <group level>)`\n")
    lines.append("DMR calling: `DSS::callDMR(p.threshold = 1e-5, delta = 0, "
                 "minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5)`\n")
    lines.append("Annotation: HOMER-equivalent (UCSC refGene) — same catalog "
                 "HOMER ships, reimplemented in R for portability.\n")

    lines.append("\n## Pipeline\n")
    lines.append(f"- DMRs called: **{len(dmr):,}**")
    lines.append(f"- Hyper / hypo: **{n_hyper} / {n_hypo}** "
                 f"({100*n_hyper/max(n_hyper+n_hypo,1):.1f}% hyper)")
    lines.append(f"- Median length: **{int(lengths.median())} bp** "
                 f"(IQR {int(lengths.quantile(0.25))}-"
                 f"{int(lengths.quantile(0.75))} bp)")
    lines.append(f"- Range: {int(lengths.min())}-{int(lengths.max())} bp")
    lines.append("")

    if paper_stats.get("paper_table_available"):
        lines.append("\n## Coordinate-level comparison to paper (Supp Table 5)\n")
        lines.append(f"- Paper DMRs: **{paper_stats['paper_n_dmr']}** "
                     f"(hyper {paper_stats['paper_n_hyper']} / "
                     f"hypo {paper_stats['paper_n_hypo']}; "
                     f"median {paper_stats['paper_median_length_bp']} bp)")
        lines.append(f"- Our DMRs:   **{paper_stats['our_n_dmr']}**")
        lines.append(f"- Coord recall: **{paper_stats['coord_recall_of_paper']*100:.1f}%** "
                     f"({int(paper_stats['coord_recall_of_paper']*paper_stats['paper_n_dmr'])} "
                     f"/ {paper_stats['paper_n_dmr']})")
        lines.append(f"- Coord precision: **{paper_stats['coord_precision']*100:.1f}%** "
                     f"({paper_stats['coord_overlap_count']} / "
                     f"{paper_stats['our_n_dmr']})")
        lines.append("")

    n_links = len(links)
    n_genes = links["gene"].nunique() if n_links else 0
    n_dmrs_linked = links["dmr_index"].nunique() if n_links else 0
    lines.append("\n## 100 kb DMR-gene linkage\n")
    lines.append(f"- (DMR, gene) pairs: **{n_links:,}**")
    lines.append(f"- Unique genes within 100 kb: **{n_genes:,}**")
    lines.append(f"- DMRs with >= 1 linked gene: **{n_dmrs_linked:,} / {len(dmr):,}**")
    if n_links:
        per = links.groupby("dmr_index").size()
        lines.append(f"- Genes per DMR: median {int(per.median())}, "
                     f"mean {per.mean():.1f}, max {int(per.max())}")
    lines.append("")

    lines.append("\n## Resources\n")
    lines.append(f"- Wall time: **{res_summary.get('wall_seconds', 0):.1f} s**")
    lines.append(f"- Peak RSS (all R + child processes): "
                 f"**{res_summary.get('rss_peak_mb', 0):.0f} MB**")
    lines.append(f"- Mean RSS: {res_summary.get('rss_mean_mb', 0):.0f} MB")
    if res_summary.get("uss_peak_mb"):
        lines.append(f"- Peak USS (unique set): "
                     f"{res_summary.get('uss_peak_mb', 0):.0f} MB")
    lines.append(f"- Peak CPU%: {res_summary.get('cpu_percent_peak', 0):.1f}% "
                 f"(1 core = 100%; logical cores available: "
                 f"{psutil.cpu_count(logical=True)})")
    lines.append(f"- Mean CPU%: {res_summary.get('cpu_percent_mean', 0):.1f}%")
    lines.append(f"- Approx core-seconds consumed: "
                 f"**{res_summary.get('approx_core_seconds', 0):.0f}**  "
                 f"(avg {res_summary.get('avg_concurrent_cores', 0):.2f} "
                 f"cores throughout run)")
    lines.append(f"- Peak threads in R process tree: "
                 f"{res_summary.get('threads_peak', 0)}")
    lines.append(f"- Sample interval: {res_summary.get('sample_interval_s', 0)}s "
                 f"({res_summary.get('samples_collected', 0)} samples)")
    lines.append("")

    lines.append("\n## Per-step R-side timings (from step_timings.tsv)\n")
    lines.append("| step | wall s | total CPU s | CPU% of wall | R-mem peak MB |")
    lines.append("|---|---:|---:|---:|---:|")
    for _, r in step.iterrows():
        lines.append(
            f"| {r['step']} | {r['wall_seconds']:.2f} | "
            f"{r['total_cpu_seconds']:.2f} | "
            f"{r['cpu_pct_of_wall']:.1f}% | "
            f"{r['r_mem_peak_mb']:.0f} |"
        )
    lines.append("")
    (OUT_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_readme() -> None:
    readme = f"""# DSS replication of GSE263850 (paper-matched parameters)

A clean reproduction of the differential methylation analysis in
Farhangdoost et al. (2024) on GSE263850 using DSS, matching the paper's
exact `DMLfit.multiFactor(smoothing=TRUE)` + `callDMR` call with
parameters `p.threshold=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100,
pct.sig=0.5`. Annotation uses the same UCSC refGene catalog HOMER ships,
reimplemented in pure R for portability.

## Files

| File | Description |
|---|---|
| `dmr_dss.tsv` / `dmr_dss.csv` | Annotated DMR table (DSS callDMR output + HOMER refGene annotation). Same row count, columns: chrom, start, end, length, n_cpgs, meanMethy1, meanMethy2, diff.Methy, areaStat, dmr_type, feature_type, feature_gene, nearest_tss_gene, nearest_tss_distance. |
| `dmr_dss_raw.tsv` | Raw `callDMR` output before annotation (for cross-checking). |
| `dmr_gene_links_100kb.csv` | Long-form (DMR, gene) pairs where the gene's canonical TSS is within 100 kb of the DMR midpoint. Paper's exact gene-linkage rule. |
| `dmltest_per_cpg.tsv.gz` | DSS::DMLtest.multiFactor output for every CpG: chr, pos, stat, phi1, phi2, pvals, fdrs. |
| `step_timings.tsv` | Per-step R-side wall / CPU / R-mem-peak. |
| `resources.csv` | Per-1s OS-level samples: RSS, USS, CPU%, threads. |
| `resources.json` | Aggregated peak/mean/totals + host info. |
| `parameters.json` | Full DSS<->paper parameter mapping. |
| `summary.md` | Headline numbers + paper-coord comparison + resource summary. |
| `run_log.txt` | Full stdout/stderr of the R script. |
| `dss_session_info.txt` | `sessionInfo()` for full reproducibility. |

## Reproducing

```
py FINAL_REPORT\\scripts\\run_dss_replication.py
```

Single command. Wraps the R script with per-1s resource sampling.
Total wall time: depends on smoothing throughput; see `summary.md`.

## Three-way comparison

This folder is shaped identically to [../chain_merge/](../chain_merge/),
so downstream comparison scripts can iterate over both engines uniformly.
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")


# ---- Main ------------------------------------------------------------------

def main() -> None:
    print(f"=== DSS replication ({datetime.now().isoformat(timespec='seconds')}) ===")
    print(f"Output dir: {OUT_DIR}")

    rc, samples, wall_s, _ = run_r()
    print(f"\nRscript exited with code {rc}. wall={wall_s:.1f}s "
          f"samples={len(samples)}")

    write_resources_csv(samples)
    res_summary = summarize_samples(samples, wall_s)
    write_resources_json(res_summary)
    write_parameters_json()

    paper_stats = compare_to_paper(OUT_DIR / "dmr_dss.csv")
    write_summary(res_summary, paper_stats)
    write_readme()

    print("\n--- Resource summary ---")
    print(json.dumps(res_summary, indent=2))
    if paper_stats.get("paper_table_available"):
        print("\n--- Paper comparison ---")
        print(json.dumps(paper_stats, indent=2))

    sys.exit(rc)


if __name__ == "__main__":
    main()
