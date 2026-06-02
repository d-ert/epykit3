"""Resume DSS replication from cached dmltest_per_cpg.tsv.gz.

Wraps resume_dss_from_dmltest.R with reliable per-1s resource sampling
(RSS + cumulative CPU-time delta — fixes the cpu_percent zero-bug seen on
Windows when wrapping Rscript). After R exits, also:

  - merges resume_log.txt into run_log.txt
  - regenerates resources.csv / resources.json with corrected numbers
  - regenerates summary.md (paper-comparable headline) and README.md
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
import psutil


REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
RAW_DIR   = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

OUT_DIR     = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss"
R_SCRIPT    = REPO_ROOT / "FINAL_REPORT" / "scripts" / "resume_dss_from_dmltest.R"
SAMPLESHEET = RAW_DIR / "samplesheet.csv"
REFGENE     = RAW_DIR / "refseq" / "refGene.txt.gz"
PAPER_DMR   = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"

PARAMS = dict(
    coverage_min=5, p_threshold=1e-5, delta=0.0, minlen=50,
    minCG=3, dis_merge=100, pct_sig=0.5, gene_link_bp=100_000,
)

SAMPLE_INTERVAL_S = 1.0


class ResourceSampler(threading.Thread):
    """Polls a psutil.Process tree every interval seconds.

    Memory: rss + uss across parent + recursive children.
    CPU: derived from cumulative cpu_times() delta across consecutive
    samples. More reliable on Windows than cpu_percent(interval=None),
    which has been observed returning stale zeros when called from a
    background thread wrapping an Rscript subprocess.
    """

    def __init__(self, root_proc: psutil.Process, interval: float = 1.0):
        super().__init__(daemon=True)
        self.root = root_proc
        self.interval = interval
        self.samples: list[dict] = []
        self.stop_flag = threading.Event()
        self.t0 = time.time()
        self.prev_cpu_total: float | None = None
        self.prev_ts: float | None = None

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
        rss = uss = 0
        threads = 0
        cpu_total_now = 0.0
        n_alive = 0
        uss_avail = False
        for p in procs:
            try:
                with p.oneshot():
                    mi = p.memory_info()
                    rss += mi.rss
                    try:
                        u = p.memory_full_info().uss
                        uss += u; uss_avail = True
                    except (psutil.AccessDenied, AttributeError):
                        pass
                    threads += p.num_threads()
                    ct = p.cpu_times()
                    cpu_total_now += ct.user + ct.system
                    n_alive += 1
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        if n_alive == 0:
            return None
        now = time.time()
        if self.prev_cpu_total is None or self.prev_ts is None:
            cpu_pct = 0.0
        else:
            dt = now - self.prev_ts
            d_cpu = cpu_total_now - self.prev_cpu_total
            cpu_pct = 100.0 * d_cpu / max(dt, 1e-6)
        self.prev_cpu_total = cpu_total_now
        self.prev_ts = now
        return dict(
            timestamp=datetime.now().isoformat(timespec="seconds"),
            elapsed_s=round(now - self.t0, 2),
            rss_mb=round(rss / 1024**2, 2),
            uss_mb=round(uss / 1024**2, 2) if uss_avail else None,
            cpu_total_seconds=round(cpu_total_now, 2),
            cpu_percent=round(cpu_pct, 1),
            num_threads=threads,
            num_processes=n_alive,
        )

    def run(self) -> None:
        first = self._sample()
        if first is not None:
            self.samples.append(first)
        while not self.stop_flag.is_set():
            time.sleep(self.interval)
            s = self._sample()
            if s is None:
                break
            self.samples.append(s)


def run_r() -> tuple[int, list[dict], float]:
    cmd = [
        r"D:\Program Files\R\R-4.5.0\bin\Rscript.exe",
        "--vanilla", str(R_SCRIPT),
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
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
    try:
        ps = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        proc.wait()
        return proc.returncode, [], time.time() - t0
    sampler = ResourceSampler(ps, SAMPLE_INTERVAL_S)
    sampler.start()
    rc = proc.wait()
    sampler.stop_flag.set()
    sampler.join(timeout=5)
    return rc, sampler.samples, time.time() - t0


def summarize(samples: list[dict], wall_s: float, prior_samples: list[dict]) -> dict:
    all_samples = prior_samples + samples
    if not all_samples:
        return dict(samples_collected=0, wall_seconds=round(wall_s, 2))
    df = pd.DataFrame(all_samples)
    out = dict(
        samples_collected=int(len(df)),
        sample_interval_s=SAMPLE_INTERVAL_S,
        wall_seconds_resume_only=round(wall_s, 2),
        rss_peak_mb=round(float(df["rss_mb"].max()), 1),
        rss_mean_mb=round(float(df["rss_mb"].mean()), 1),
        uss_peak_mb=(round(float(df["uss_mb"].max()), 1)
                     if df["uss_mb"].notna().any() else None),
        cpu_percent_peak=round(float(df["cpu_percent"].max()), 1)
                         if "cpu_percent" in df.columns else None,
        cpu_percent_mean=round(float(df["cpu_percent"].mean()), 1)
                         if "cpu_percent" in df.columns else None,
        threads_peak=int(df["num_threads"].max()),
    )
    if "cpu_total_seconds" in df.columns and df["cpu_total_seconds"].notna().any():
        cpu_t = df["cpu_total_seconds"].dropna().astype(float)
        out["cpu_total_seconds_observed"] = round(float(cpu_t.max()), 1)
    return out


def write_resources(samples: list[dict], summary: dict) -> None:
    if samples:
        df = pd.DataFrame(samples)
        # Append-mode: keep prior resources.csv if present
        prior = OUT_DIR / "resources.csv"
        if prior.exists():
            # Add a "phase" column to distinguish
            df["phase"] = "resume"
            prior_df = pd.read_csv(prior)
            if "phase" not in prior_df.columns:
                prior_df["phase"] = "initial"
            df = pd.concat([prior_df, df], ignore_index=True, sort=False)
        df.to_csv(OUT_DIR / "resources.csv", index=False)

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
    (OUT_DIR / "resources.json").write_text(
        json.dumps(dict(system=sys_info, resources=summary), indent=2),
        encoding="utf-8",
    )


def merge_logs() -> None:
    """Concatenate the original run_log.txt + resume_log.txt into a single
    chronological log."""
    initial = OUT_DIR / "run_log.txt"
    resume  = OUT_DIR / "resume_log.txt"
    if not resume.exists():
        return
    parts = []
    if initial.exists():
        parts.append("# ---- Initial run (crashed at callDMR) ----\n")
        parts.append(initial.read_text(encoding="utf-8"))
    parts.append("\n\n# ---- Resume from cached DMLtest ----\n")
    parts.append(resume.read_text(encoding="utf-8"))
    initial.write_text("".join(parts), encoding="utf-8")


def _interval_index(df: pd.DataFrame):
    by_ch: dict[str, list] = {}
    for _, r in df.iterrows():
        by_ch.setdefault(str(r["chrom"]), []).append(
            (int(r["start"]), int(r["end"]))
        )
    for ch in by_ch:
        by_ch[ch].sort()
    return by_ch


def _any_overlap(ints, s, e) -> bool:
    for ps, pe in ints:
        if pe < s: continue
        if ps > e: break
        return True
    return False


def compare_to_paper(dmr_csv: Path) -> dict:
    if not PAPER_DMR.exists() or not dmr_csv.exists():
        return {"paper_table_available": False}
    paper = pd.read_excel(PAPER_DMR, sheet_name=0)
    paper_norm = paper[["chr", "start", "end"]].rename(
        columns={"chr": "chrom"}).copy()
    paper_norm["chrom"] = paper_norm["chrom"].astype(str)
    paper_norm["start"] = paper_norm["start"].astype(int)
    paper_norm["end"]   = paper_norm["end"].astype(int)
    ours = pd.read_csv(dmr_csv)
    ours["chrom"] = ours["chrom"].astype(str)
    ours = ours[["chrom", "start", "end"]].copy()
    pby = _interval_index(paper_norm)
    oby = _interval_index(ours)
    n_paper = len(paper_norm); n_ours = len(ours)
    o_hit = sum(_any_overlap(pby.get(str(r["chrom"]), []),
                             int(r["start"]), int(r["end"]))
                for _, r in ours.iterrows())
    p_hit = sum(_any_overlap(oby.get(str(r["chrom"]), []),
                             int(r["start"]), int(r["end"]))
                for _, r in paper_norm.iterrows())
    return dict(
        paper_table_available=True,
        paper_n_dmr=n_paper,
        paper_n_hyper=int((paper["diff.meth_mean"] > 0).sum()),
        paper_n_hypo=int((paper["diff.meth_mean"] < 0).sum()),
        paper_median_length_bp=int((paper["end"] - paper["start"]).median()),
        our_n_dmr=n_ours,
        coord_overlap_count=int(o_hit),
        coord_recall_of_paper=round(p_hit / max(n_paper, 1), 4),
        coord_precision=round(o_hit / max(n_ours, 1), 4),
    )


def write_summary(res: dict, paper_stats: dict) -> None:
    dmr_csv = OUT_DIR / "dmr_dss.csv"
    if not dmr_csv.exists():
        return
    dmr = pd.read_csv(dmr_csv)
    lengths = dmr["end"] - dmr["start"]
    n_hyper = int((dmr["dmr_type"] == "hyper").sum())
    n_hypo  = int((dmr["dmr_type"] == "hypo").sum())
    links = pd.read_csv(OUT_DIR / "dmr_gene_links_100kb.csv")

    L = []
    L.append("# DSS replication of GSE263850 (paper-matched parameters)\n")
    L.append("Per-CpG: `DSS::DMLfit.multiFactor(smoothing=TRUE)` -> "
             "`DMLtest.multiFactor(coef = groupWT)`\n")
    L.append("DMR calling: `DSS::callDMR(p.threshold=1e-5, delta=0, "
             "minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)`\n")
    L.append("Note: callDMR on multifactor output returns `chr, start, end, "
             "length, nCG, areaStat` — no `diff.Methy`. Per-DMR per-group "
             "mean methylation (and hence diff.Methy / dmr_type) is derived "
             "here directly from the 6 per-CpG BEDs.\n")
    L.append("\n## Pipeline\n")
    L.append(f"- DMRs called: **{len(dmr):,}**")
    L.append(f"- Hyper / hypo: **{n_hyper} / {n_hypo}** "
             f"({100*n_hyper/max(n_hyper+n_hypo,1):.1f}% hyper)")
    L.append(f"- Median length: **{int(lengths.median())} bp** "
             f"(IQR {int(lengths.quantile(0.25))}-{int(lengths.quantile(0.75))} bp)")
    L.append(f"- Range: {int(lengths.min())}-{int(lengths.max())} bp")
    L.append("")
    if paper_stats.get("paper_table_available"):
        L.append("\n## Coordinate-level comparison to paper (Supp Table 5)\n")
        L.append(f"- Paper DMRs: **{paper_stats['paper_n_dmr']}** "
                 f"(hyper {paper_stats['paper_n_hyper']} / "
                 f"hypo {paper_stats['paper_n_hypo']}; "
                 f"median {paper_stats['paper_median_length_bp']} bp)")
        L.append(f"- Our DMRs:   **{paper_stats['our_n_dmr']}**")
        L.append(f"- Coord recall: **{paper_stats['coord_recall_of_paper']*100:.1f}%** "
                 f"({int(paper_stats['coord_recall_of_paper']*paper_stats['paper_n_dmr'])} "
                 f"/ {paper_stats['paper_n_dmr']})")
        L.append(f"- Coord precision: **{paper_stats['coord_precision']*100:.1f}%** "
                 f"({paper_stats['coord_overlap_count']} / "
                 f"{paper_stats['our_n_dmr']})")
        L.append("")
    n_links = len(links); n_genes = links["gene"].nunique() if n_links else 0
    n_dmrs_linked = links["dmr_index"].nunique() if n_links else 0
    L.append("\n## 100 kb DMR-gene linkage\n")
    L.append(f"- (DMR, gene) pairs: **{n_links:,}**")
    L.append(f"- Unique genes within 100 kb: **{n_genes:,}**")
    L.append(f"- DMRs with >= 1 linked gene: **{n_dmrs_linked:,} / {len(dmr):,}**")
    if n_links:
        per = links.groupby("dmr_index").size()
        L.append(f"- Genes per DMR: median {int(per.median())}, "
                 f"mean {per.mean():.1f}, max {int(per.max())}")
    L.append("")
    L.append("\n## Resources\n")
    L.append(f"- Resume wall time: **{res.get('wall_seconds_resume_only',0):.1f} s**")
    L.append("- Initial run (DMLfit smoothing+test) wall: **~2,044 s** "
             "(see resume_log.txt header for details).")
    L.append(f"- Peak RSS across the whole run: "
             f"**{res.get('rss_peak_mb', 0):.0f} MB**")
    L.append(f"- Mean RSS: {res.get('rss_mean_mb', 0):.0f} MB")
    if res.get("uss_peak_mb"):
        L.append(f"- Peak USS: {res.get('uss_peak_mb', 0):.0f} MB")
    L.append(f"- Peak CPU% (1 core = 100%; logical cores: "
             f"{psutil.cpu_count(logical=True)}): "
             f"{res.get('cpu_percent_peak', 0):.1f}%")
    L.append(f"- Mean CPU%: {res.get('cpu_percent_mean', 0):.1f}%")
    if res.get("cpu_total_seconds_observed"):
        L.append(f"- Cumulative R-process CPU time observed: "
                 f"**{res.get('cpu_total_seconds_observed',0):.0f} s**")
    L.append(f"- Peak threads in R tree: {res.get('threads_peak', 0)}")
    L.append(f"- Total samples collected: {res.get('samples_collected', 0)}")
    L.append("")
    (OUT_DIR / "summary.md").write_text("\n".join(L), encoding="utf-8")


def write_readme() -> None:
    readme = """# DSS replication of GSE263850 (paper-matched parameters)

Clean reproduction of the differential methylation analysis in
Farhangdoost et al. (2024) on GSE263850 using DSS — matching the paper's
exact `DMLfit.multiFactor(smoothing=TRUE)` + `callDMR` call with
parameters `p.threshold=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100,
pct.sig=0.5`. Annotation uses the same UCSC refGene catalog HOMER ships,
reimplemented in R for portability.

Note: callDMR on multifactor DMLtest output returns
`[chr, start, end, length, nCG, areaStat]` — no `diff.Methy`. Per-DMR
per-group mean methylation (and hence `diff.Methy` / `dmr_type`) is
derived here by averaging per-sample (sumM / sumT) within each DMR and
then averaging within each group. This matches the paper's Table 5
`diff.meth_mean` column derivation.

## Files

| File | Description |
|---|---|
| `dmr_dss.tsv` / `dmr_dss.csv` | Annotated DMR table (chrom, start, end, length, n_cpgs, areaStat, dmr_id, meanMethy_treatment, meanMethy_control, diff.Methy, dmr_type, feature_type, feature_gene, nearest_tss_gene, nearest_tss_distance). |
| `dmr_dss_raw.tsv` | Raw `callDMR` output before any per-DMR mean-methylation derivation. |
| `dmr_gene_links_100kb.csv` | Long-form (DMR, gene) pairs where the gene's canonical TSS is within 100 kb of the DMR midpoint. |
| `dmltest_per_cpg.tsv.gz` | `DMLtest.multiFactor` output for every CpG (cached; used by the resume script). |
| `step_timings.tsv` / `step_timings_resume.tsv` | Per-step R-side wall / CPU / R-mem-peak. |
| `resources.csv` | Per-1s OS-level samples: RSS, USS, CPU%, threads. Two phases marked: `initial` (DMLfit+DMLtest) and `resume` (callDMR + downstream). |
| `resources.json` | Aggregated peak/mean/totals + host info. |
| `parameters.json` | DSS<->paper parameter mapping. |
| `summary.md` | Headline numbers + paper-coord comparison + resource summary. |
| `run_log.txt` | Concatenated log of both phases. |
| `dss_session_info.txt` | `sessionInfo()` for reproducibility. |
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    print(f"=== DSS resume ({datetime.now().isoformat(timespec='seconds')}) ===")

    # Load prior resource samples if any (initial DMLfit phase)
    prior_samples: list[dict] = []
    prior_csv = OUT_DIR / "resources.csv"
    if prior_csv.exists():
        prior_samples = pd.read_csv(prior_csv).to_dict(orient="records")

    rc, samples, wall_s = run_r()
    print(f"\nRscript exited with code {rc}. resume wall={wall_s:.1f}s "
          f"new samples={len(samples)}")

    summary = summarize(samples, wall_s, prior_samples)
    write_resources(samples, summary)
    merge_logs()

    paper_stats = compare_to_paper(OUT_DIR / "dmr_dss.csv")
    write_summary(summary, paper_stats)
    write_readme()

    print("\n--- Resources ---")
    print(json.dumps(summary, indent=2))
    if paper_stats.get("paper_table_available"):
        print("\n--- Paper comparison ---")
        print(json.dumps(paper_stats, indent=2))

    sys.exit(rc)


if __name__ == "__main__":
    main()
