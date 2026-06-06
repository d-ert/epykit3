# epykit Benchmark — Final Report

Consolidated outputs of three benchmark studies comparing **epykit** (a
Python-native pipeline for whole-genome bisulfite sequencing analysis) with
established R/CLI tools.

## What to read first

| If you have… | Read |
|---|---|
| ~5 minutes | [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md) (or `.docx`) |
| ~30 minutes | [paper/paper.md](paper/paper.md) (or `paper/paper.docx`) |
| Time to verify numbers | [report/REPORT.md](report/REPORT.md) + [report/methods_appendix.md](report/methods_appendix.md) |
| A specific figure | [figures/](figures/) (per-study + cross-benchmark summaries) |
| The underlying data | [data/](data/) + [data/README.md](data/README.md) |

## Folder layout

```
FINAL_REPORT/
├── README.md                       (this file)
├── EXECUTIVE_SUMMARY.md            (1-page TL;DR)
├── EXECUTIVE_SUMMARY.docx
├── paper/
│   ├── paper.md                    (unified scientific manuscript)
│   ├── paper.docx
│   └── refs.bib                    (Bibtex bibliography)
├── report/
│   ├── REPORT.md                   (full technical report, all tables)
│   ├── REPORT.docx
│   └── methods_appendix.md         (tool versions, parameters, schemas)
├── figures/
│   ├── study1_simulated_allPackages/  (Piao et al. 2021 panel comparison)
│   ├── study2_simulated_headToHead/   (epykit vs methylKit, simulated)
│   ├── study3_real_GSE263850/         (epykit vs methylKit, real WGBS)
│   └── summary/                       (3 new cross-benchmark figures)
├── data/
│   ├── study1/                     (eval_summary.parquet, ground truth, baselines)
│   ├── study2/                     (ground truth, methylKit per-scenario TSVs)
│   ├── study3/                     (benchmark CSVs, sig DMC/DMR tables, samplesheet)
│   └── README.md                   (data provenance and schemas)
└── scripts/
    └── make_summary_figures.py     (generates figures/summary/*.png)
```

## What the three studies are

| Study | What was compared | Data | Headline |
|---|---|---|---|
| **1. Simulated panel** | epykit vs 8 baselines (methylKit, DSS, RADMeth, BiSeq, methylSig, pooled Fisher, BSmooth, metilene) | Piao et al. 2021 simulated grid: 100K CpGs (DMC) + 3.97M CpGs (DMR), coverage 5×–25×, n = 2–10 | epykit `lr` competitive with methylKit/RADMeth/DSS across the grid; `lr+` reaches TPR ≥ 0.999 at every cell; FPR 100×–600× tighter at 5× |
| **2. Simulated head-to-head** | epykit vs methylKit, same machine, same harness | Same simulated grid | TPR/FPR/F1/AUROC identical to 3 dp at n ≥ 4; epykit 2× recall at n = 2; **43× faster** on full grid |
| **3. Real WGBS (3-way)** | epykit chain_merge vs methylKit-tile vs DSS-from-scratch, against the published call set in Farhangdoost et al. 2025 Supp Table 5 | GEO GSE263850, Het-AKAP11-KO vs SBP009 WT, hg38, 15.6 M CpGs after filter (22 M loaded) | Per-CpG Pearson r = 0.994 (methylKit), 94 % direction agreement. Paper-DMR coord recall: methylKit-tile 9 %, ek-chain_merge-100 53 %, ek-chain_merge-250 63 %, DSS-from-scratch 87.5 %. **100 % direction agreement on every matched DMR across every engine.** 12× faster than methylKit, 6× faster than DSS at matched-parameter chain_merge. |

## How to reproduce

Each of the three studies is independently reproducible from data in its
own source folder (sibling directories of `FINAL_REPORT/`). The
consolidated folder does not duplicate the heavy raw per-CpG outputs —
only the summary tables (`eval_summary.parquet`, ground truth, baseline
tables, methylKit per-scenario TSVs, Study 3 benchmark CSVs) are kept
under `data/`.

To reproduce a study end-to-end, see the recipes in
[report/methods_appendix.md](report/methods_appendix.md) §§ 7–9.

To regenerate the cross-benchmark summary figures (`figures/summary/`):

```bash
cd FINAL_REPORT/scripts
python make_summary_figures.py
```

This loads `data/study1/eval_summary.parquet`,
`data/study2/methylkit_results/timings.tsv`, and
`data/study3/benchmark/step_benchmarks.csv` to produce the three PNGs.

### Reproducing the Study 3 three-way analyses

The Study 3 chain_merge + DSS replication scripts and their outputs:

```bash
cd FINAL_REPORT

# 1. epykit chain_merge replication at paper-matched parameters
py scripts/run_chain_merge_replication.py
# -> data/study3/chain_merge/

# 2. dis.merge parameter sweep (re-uses chain_merge DMC cache)
py scripts/sweep_dis_merge.py
# -> data/study3/chain_merge_dis_merge_sweep/

# 3. DSS from-scratch (Rscript wrapped in psutil-sampled Python)
py scripts/run_dss_replication.py        # ~50 min wall (DMLfit dominates)
py scripts/resume_dss_from_dmltest.py    # 8 min finish after the column-bug fix
# -> data/study3/dss/

# 4. Pairwise + multi-way comparison scripts (~1 min each)
py scripts/compare_chain_merge_to_paper.py
py scripts/compare_epykit_to_dss.py
py scripts/per_dmr_stat_concordance.py
py scripts/annotation_distribution.py
py scripts/run_enrichment_three_way.py
# -> data/study3/comparisons/{chain_merge_vs_paper,epykit_vs_dss,...}/

# 5. Figure regeneration (F1-F10 + summary composite)
for f in FINAL_REPORT/scripts/figures_v2/f*.py; do py -X utf8 $f; done
# -> figures/study3_real_GSE263850/three_way/
```

Background on why Study 3 was restructured: see
[dmr_replication_investigation.md](dmr_replication_investigation.md)
for the narrative of how fixed-tile callers turned out to miss 90 %
of focused biological DMRs against the paper's DSS call set, and how
chain_merge + a DSS-from-scratch upper bound resolved the picture.

## Versions consolidated in this report

| Tool | Version |
|---|---|
| epykit | 1.0.0 (all studies; run at commit 60a71e0, engine tag v0.7.5-phase3-engines-frozen) |
| methylKit | 1.34.0 (Study 2); 1.36.0 (Study 3); 0.99.2 (Study 1 baseline) |
| R | 4.5.0 (Studies 2, 3) |
| Python | 3.12 |

### Reproducing the published numbers (pinned Python environment)

The exact resolved Python environment is committed as `uv.lock` at the repo
root. Reproduce the epykit numbers against that exact environment with:

```bash
uv sync --frozen --extra dev --extra all
```

All published epykit timings and per-CpG counts were produced **single-threaded**
so float-reduction order is deterministic across machines (parallel reductions
are not associative). Pin the thread pools before running:

```bash
export POLARS_MAX_THREADS=1
export OMP_NUM_THREADS=1
```

## Date

This consolidated report was assembled on **2026-05-22**. Study 1 was
completed on 2026-05-19; Study 2 on 2026-05-19; Study 3 on 2026-05-21.

## Limitations to be aware of

* **Simulator-real gap.** The Piao 2021 simulator is underdispersed
  (median φ ≈ 0.41 at 5×); real WGBS has heteroscedastic overdispersion
  (φ ≈ 1.5–5). Low-coverage TPR advantages observed in Studies 1 and 2
  may shrink or reverse on real data.
* **Single real dataset.** Study 3 is one tissue × one genome. Multi-dataset
  validation is future work.
* **Baseline software versions (Study 1).** Numbers come from 2021
  software releases via the published supplementary tables (audited
  cell-by-cell, 0 errors); relative ordering at low coverage / small n is
  robust across recent versions, but absolutes may have shifted.
* **Windows host (Study 2).** methylKit's `mc.cores` is a no-op on
  Windows; methylKit ran single-threaded by force. On Linux with parallel
  cores enabled, methylKit's DMR grid would drop from ~6 h to ~1–1.5 h
  but epykit would still be ~10× faster.

## Single-command Linux re-run (Track 1 Layer B)

The plan to re-run every benchmark from scratch on a Linux machine is
end-to-end driven by `scripts/regen_all.py --run-all`. The full sweep
covers:

* The Piao simulator at five `phi` values `{0.0, 0.01, 0.05, 0.1, 0.2}`
  spanning binomial through real-WGBS-like overdispersion.
* Twenty seeds (2026000..2026019) per phi cell.
* epykit's four DMC engines (lr, lr+, welch_t, fisher) on every cell.
* methylKit, DSS, dmrseq, BSmooth as locally-re-run R baselines.
* `k=1000`-shuffle null calibration on bare `lr` and `lr+` with Q-Q +
  KS plots.
* Both `truth_mode` labelings (intrinsic and threshold).
* DMR-caller parameter sensitivity sweep at ±50% of the default preset.

Setup once, on the Linux box:

```bash
git clone <repo> && cd epykit
docker build -t epykit-py -f Dockerfile.python .
docker build -t epykit-r  -f Dockerfile.r .

# One-shot R environment lockfile generation (writes benchmark/renv.lock).
docker run --rm -v "$PWD":/work -w /work epykit-r \
    Rscript benchmark/renv/install_packages.R

# Re-build the R image so it picks up the new lockfile.
docker build -t epykit-r -f Dockerfile.r .
```

Then run the full sweep:

```bash
# Sanity-check the dispatcher emits the plan you expect.
docker run --rm -v "$PWD":/work -w /work epykit-py \
    python benchmark/scripts/regen_all.py --run-all --dry-run --verbose

# Engine-regression slice (~30 s) -- catches any drift before the big run.
docker run --rm -v "$PWD":/work -w /work epykit-py \
    python benchmark/scripts/regen_small.py --update
git add benchmark/scripts/regen_small_hashes.json && git commit -m \
    "benchmark: snapshot Linux engine-regression hashes"

# The full pipeline.
docker run --rm -v "$PWD":/work -w /work epykit-py \
    python benchmark/scripts/regen_all.py --run-all --verbose
```

The dispatcher runs steps linearly; a failing step records the
failure and continues so a single broken comparator doesn't nuke
the whole run. `--only <step>` re-runs a single named step;
`--skip <step>` skips one or more steps. See `regen_all.py --help`
for the full step list.

After everything completes, commit the regenerated `benchmark/data/`
outputs:

```bash
git add benchmark/data && git commit -m "benchmark: full Linux re-run"
```

The TSV mirror under `benchmark/paper_data/` is re-derived from the
parquet sources via the converter (kept as a separate step so the
parquet sources are the single source of truth).

## Contact

Open an issue in the epykit repository for questions about specific
numbers or reproduction failures.
