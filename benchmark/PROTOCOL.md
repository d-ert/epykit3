# PROTOCOL — epykit Benchmark Re-do (Pre-registration)

**Status:** Locked at commit-of-record (to be tagged after first re-run begins).
**Companion:** [paper/paper.md](paper/paper.md), [report/REPORT.md](report/REPORT.md), [report/methods_appendix.md](report/methods_appendix.md).
**Purpose:** lock every methodological choice for the re-do that follows the
2026-05-22 integrity audit (see [is-this-a-good-cryptic-pike.md](../../is-this-a-good-cryptic-pike.md)).
Any deviation requires a written addendum at the bottom of this file — no
silent edits.

The single guiding principle: **whatever parameter exploration was done
on epykit's side must also be done on the baseline side, OR the headline
table must compare default-vs-default. Tuned-vs-default comparisons are
banned from the headline and may only appear as a clearly labelled
sensitivity panel.**

---

## 1. Tool versions (frozen)

| Tool | Version | Notes |
|---|---|---|
| Python | 3.12 | All studies |
| epykit | **1.0.0** | All studies, including the Study 3 re-run (executed at commit 60a71e0, engine tag v0.7.5-phase3-engines-frozen; version string 0.7.2 at run time, engine code unchanged through 1.0.0) |
| polars | 1.x (latest at time of run; record exact in `_seeds.json`) | I/O |
| NumPy, SciPy, statsmodels | latest at time of run | Recorded by `pip freeze` into `run_env.txt` per scenario |
| psutil | ≥ 5.9 | Per-scenario peak RSS, 50 ms sampling |
| R | **4.5.0** (2025-04-11, ucrt, x86_64-w64-mingw32) | Studies 2 and 3 |
| methylKit | **1.34.0** Study 2; **1.36.0** Study 3 | Two versions disclosed — Study 3 host's installed version differs from Study 2's |
| DSS | **2.56.0** | Study 3 (ceiling caller) |
| bsseq | **1.44.1** | Study 3 (DSS dependency) |
| BiocParallel | 1.42.1 | Study 3 |
| Study 1 panel baselines | methylKit 0.99.2, methylSig 0.4.4, DSS 2.12.0, RADMeth, BiSeq, Fisher (pooled), BSmooth, metilene | **Transcribed from Piao et al. 2021 supplementary tables — not re-run.** |

**Disclosure:** Study 1's eight-tool panel uses *published* TPR/FPR
values from Piao 2021. We do not claim to have re-run those tools at
their 2021 versions. Study 1 is therefore a panel cross-check against
the published literature, not a re-execution.

---

## 2. Hardware (frozen)

| Study | Host | CPU | Logical cores | OS | RAM |
|---|---|---|---|---|---|
| 1 & 2 | `DESKTOP-0GUMIA4` | AMD Ryzen 7 4800H | 8 physical / 16 logical | Windows 10 build 19045 | (record from `wmic computersystem get TotalPhysicalMemory`) |
| 3 (methylKit) | `pivoine` | (record `lscpu` snapshot in `data/study3/methylkit/host_info.txt`) | 24 logical | Linux (Debian-family) | (record) |
| 3 (epykit, DSS) | `DESKTOP-0GUMIA4` | same as Studies 1/2 | 16 logical | Windows 10 | same |

**Threading:** methylKit `mc.cores = 1` (no-op on Windows; explicit on
Linux for fair single-core comparison). epykit single-process for
timing; parallelism within epykit is internal and not toggled. DSS
`BPPARAM = SerialParam()`.

---

## 3. Datasets (frozen)

### 3.1 Studies 1 & 2 — Piao et al. 2021 simulated grid

- **Source:** IJERPH 18:7975 supplementary materials (or its companion repo).
- **Files:** `raw_sim_data/simulated_datasets/amp.coverage={5,10,15,20,25}.sampleN.txt` and `amp.replicate={2,4,6,8,10}.sampleN.txt`.
- **Status:** fixed published dataset — one instance, not a re-runnable simulator. *Variance over additional simulator draws is therefore not captured by this benchmark; see §6.*
- **Truth:** reconstructed from the 25× coverage samples via `scripts/_make_truth.py` at `TRUTH_THRESHOLD = 0.20`, `DMR_GAP_BP = 1000`, `DMR_MIN_CPGS = 5`. Recovers 19,999 DMCs and 35 reference DMRs (matches the paper's design).
- **Truth caveat (disclosed):** truth uses the same simulator output (high-coverage samples) that the tools read at lower coverage. This is mildly self-referential at low coverage; the alternative — using the simulator's internal `is_dmc` flag — is unavailable because the supplementary files don't carry it.

### 3.2 Study 3 — GSE263850 (Farhangdoost et al. 2025, *Genome Biol*)

- **Samples:** 6 (Clone16, Clone20, Clone21 vs SBP009 untreated 1/2/3), hg38.
- **Format:** 12-column strand-collapsed Bismark-style BEDs from the GEO submission. Cols 10–12 = M / T / pct after strand merging.
- **Both pipelines receive bit-identical per-CpG counts.** methylKit gets pre-converted 6-col `.cov.gz`; epykit reads cols 10–12 directly via `read_combined_strand_bed()`. Verified 0.02 % delta at edge cases only.
- **Reference call set for DMR reproduction:** the 813 DMRs from Farhangdoost et al. 2025 Supplementary Table 5 (DSS::callDMR output, `p.threshold=1e-5`, `minCG=3`, `minlen=50`, `dis.merge=100`).

---

## 4. Parameters (frozen)

All thresholds, recipes, and post-processing steps below are locked **before** the re-do runs. Any change later requires an addendum.

### 4.1 DMC pipelines — default recipes (headline row)

| Tool | Recipe |
|---|---|
| epykit | `ep.tl.dmc(test="lr", dispersion="site", fdr_method="fdr_bh", allow_n1=True)` |
| methylKit | `methRead(..., mincov=10) → normalizeCoverage(method="median") → unite(destrand=FALSE) → calculateDiffMeth(mc.cores=1)` |

**Cutoff (both):** `qvalue < 0.05` AND `|meth_diff| ≥ 0.25` (fractional; = 25 on methylKit's percent scale).

### 4.2 DMC pipelines — tuned recipes (sensitivity row, NEW in re-do)

| Tool | Recipe | Purpose |
|---|---|---|
| epykit | `lr+` = `lr` + neighbour_combine (Stouffer) + sep_fallback + `fdr_method="fdr_tsbh"` | already the existing opt-in recipe |
| methylKit | `calculateDiffMeth` + **post-hoc Stouffer combine across windows of 3 adjacent CpGs** (new R helper in `scripts/methylkit_stouffer_combine.R`) | analog of epykit's `lr+` neighbour_combine; same window of 3 |

**Cutoff (both):** unchanged (`q < 0.05`, `|Δ| ≥ 0.25`).

### 4.3 DMR pipelines

| Tool / engine | Recipe | Where it appears |
|---|---|---|
| epykit `dmr_chain_merge` | `alpha = 1e-5, delta = 0, minlen = 50, minCG = 3, pct.sig = 0.5, dis_merge_bp = 100` (paper-faithful match to DSS::callDMR) | **Headline** for Study 3 |
| epykit `dmr_chain_merge` `dis_merge_bp = 250` | same as above with merge gap relaxed | **Sensitivity** panel only — not headline |
| epykit `dmr_tile` | 1 kbp fixed tiles, ≥ 5 CpGs per tile | Studies 1, 2 (matches the Piao baseline framework) |
| methylKit `tileMethylCounts` | `win.size = 1000`, `step.size = 1000`, `cov.bases = 5` | Studies 1, 2 (default) |
| methylKit `tileMethylCounts` + Stouffer | + post-hoc Stouffer merge of adjacent significant tiles (gap ≤ 1000 bp) | Tuned-vs-tuned sensitivity panel |
| DSS `DMLfit.multiFactor` + `callDMR` | `p.threshold = 1e-5, minCG = 3, minlen = 50, dis.merge = 100` (paper) | Study 3 ceiling caller |

### 4.4 Ground-truth construction

Studies 1 & 2 only (Study 3 has no truth — see §5):

```python
# In data/study{1,2}/ground_truth/make_truth.py
TRUTH_THRESHOLD = 0.20         # effect-size threshold (fractional)
DMR_GAP_BP      = 1000         # max gap between same-direction true DMCs
DMR_MIN_CPGS    = 5            # min CpGs per reference DMR
```

Recovers 19,999 DMCs (20 %) and 35 reference DMRs. **Locked — do not retune.**

---

## 5. Decision rules (pre-committed)

These rules govern what gets reported in `EXECUTIVE_SUMMARY.md` and `paper/paper.md` after the re-do.

| Rule | What it means | Why |
|---|---|---|
| **R1. Default-vs-default in the headline.** | Headline tables compare epykit `lr` default to methylKit default (§4.1). Tuned-vs-tuned (§4.2) appears in a clearly labelled sensitivity panel below. | Tuned-vs-default is the single biggest defensibility risk in the current report. |
| **R2. CIs on every headline number.** | Wilson 95 % CI on TPR/FPR; bootstrap 95 % percentile CI on AUROC/F1 (B = 1000). No bare point estimates. | Single-run point estimates are the second biggest risk. |
| **R3. Same DMR semantics on both sides for each comparison.** | When comparing chain_merge to methylKit-tile, the report must explicitly note these are different region-definition semantics. When comparing chain_merge to DSS-callDMR, both share semantics — that's the appropriate apples-to-apples cell. | Comparing 102 chain_merge DMRs vs 37 1-kbp tiles by call count is methodologically meaningless. |
| **R4. Promote the discovery narrative.** | The "tile → chain_merge" pivot (8 % → 53 % paper-DMR recall) is described in the main text of paper.md as a methodological discovery, not buried in `benchmark/docs/historical/`. | A reviewer reading the GitHub history will find this anyway; better to own it. |
| **R5. Top-K agreement is part of the main report.** | Top-5 / 10 / 25 / 50 / 100 DMC and DMR agreement between epykit and methylKit goes into the body of paper.md, not as appendix-only. Top-5 = 0/5 must be stated honestly. | A reviewer who reads `top_k_report.md` and finds it isn't in the headline will read the paper's framing as evasive. |
| **R6. Version mismatch surfaced.** | The paper's Methods §1 explicitly names the epykit version per study; Study 3 was re-run, so all three studies are at 1.0.0 (engine-frozen at v0.7.5-phase3-engines-frozen). | Buried footnotes don't survive review. |
| **R7. Limitations section is non-trivial.** | `paper/paper.md` ends with a Limitations subsection that names: (a) one fixed simulator instance, (b) one tissue × one genome for real data, (c) Study 1 panel baselines transcribed not re-run, (d) Study 3's `mincov=10` deviates from the Piao recipe `mincov=0` for fairness, (e) bugs discovered during benchmarking. | Honest limitations sections protect against rejection; vague ones invite it. |
| **R8. Every claim is traceable to a script.** | Each numeric claim in EXECUTIVE_SUMMARY.md and paper.md gets a `<!-- source: scripts/X.py -->` comment immediately above. The script reads the same locked outputs and produces the same number. | Reproducibility is non-negotiable for peer review. |

---

## 6. Variance estimation

Because the Piao 2021 dataset is a fixed instance, we cannot estimate
simulator variance. We instead estimate the variance of the *evaluation*
on the data we have:

### 6.1 TPR / FPR — Wilson 95 % binomial CIs

TPR and FPR are proportions on a finite truth set:

- TPR = TP / (TP + FN), denominator = `len(truth_positive_set)` = 19 999.
- FPR = FP / (FP + TN), denominator = `len(truth_negative_set)` = 80 001.

For a binomial proportion `p̂ = k/n`, the Wilson CI is closed-form and is
the standard recommendation over the naive Wald CI for small or extreme
proportions. Implementation: `scipy.stats.binomtest(k, n).proportion_ci(method="wilson", confidence_level=0.95)`.

### 6.2 AUROC and F1 — bootstrap 95 % percentile CIs

These metrics depend on the joint ranking of all CpGs and don't have a
closed-form CI under the finite-truth-set design. We bootstrap CpGs with
replacement, B = 1000, and report the 2.5th and 97.5th percentiles of
each metric across bootstraps.

Implementation: `scripts/wilson_bootstrap_ci.py` (new, to be added in
Phase 2). Reads the existing single-run output tables; no methylKit or
epykit re-runs required for this variance story.

### 6.3 What is NOT captured

- Variance over independent simulator draws (we have one Piao instance).
- Variance over independent biological samples on real data (we have one
  GSE263850 cohort).
- Variance from tool-internal stochasticity (we set seeds for both tools
  where they accept them; methylKit and DSS are deterministic given input).

This is disclosed in §R7 of the paper's Limitations.

### 6.4 Future work — Piao simulator re-implementation (deferred)

Re-implementing Piao's binomial simulator (~200 lines of Python, ~1 day
of work) would let us draw N independent simulator instances and
estimate simulator variance properly. This is deferred to a robustness
appendix if reviewers demand it.

---

## 7. Metrics — exact formulas

For DMCs (per-CpG), with `truth` ∈ {DMC, non-DMC} from §4.4:

```
TP = | called ∩ truth_DMC |
FP = | called \ truth_DMC |
FN = | truth_DMC \ called |
TN = | (all CpGs \ called) \ truth_DMC |

TPR = TP / (TP + FN)
FPR = FP / (FP + TN)
Precision = TP / (TP + FP)  ; defined as 1.0 when TP + FP = 0
F1 = 2 · Precision · TPR / (Precision + TPR)
AUROC = sklearn.metrics.roc_auc_score(truth_labels, -log10(qvalue))
```

For DMRs, "match" between a called region and a reference region uses
**any-bp overlap by default** (lifted from the chain_merge replication
investigation). A second column reporting **strict 50%-reciprocal
overlap** appears in the same table as a stricter criterion. Both
columns are computed; neither is hidden.

For real-data DMR reproduction (Study 3), the truth is the 813-DMR
Farhangdoost et al. 2025 published set; match criterion is identical
(any-bp + strict 50 % reciprocal, both columns).

---

## 8. Output schema (locked)

Every scenario writes:

```
data/study<N>/<engine>/<scenario>/
  ├── dmc.parquet                # per-CpG: chrom, start, qvalue, meth_diff, called
  ├── dmr.parquet                # per-region: chrom, start, end, qvalue, dir, n_cpgs
  ├── metrics.json               # TP, FP, FN, TN, TPR, FPR, F1, AUROC, Wilson CIs, bootstrap CIs
  ├── timing.json                # wall_sec, cpu_sec, peak_mb (psutil-sampled)
  └── run_env.txt                # output of `pip freeze` or `sessionInfo()`
```

A driver script `scripts/regen_all.py` reads `PROTOCOL.md`'s parameter
section by reference and reproduces every metric in the headline tables
from these parquets. If `regen_all.py` produces a number that doesn't
match the paper, the paper is wrong, not the data.

### 8.1 Two parallel benchmark trees (data sources vs paper mirrors)

The benchmark directory has two structurally parallel data trees and
**both ship with the repo**. A reviewer that sees only one of them
will get a wrong picture of the reproducibility story.

| Tree | Format | What it is | Who reads it |
|---|---|---|---|
| `benchmark/data/` | Parquet (binary) | The frozen, byte-exact source-of-truth artefacts produced by `regen_all.py`. `.gitignore` whitelists each canonical artefact (e.g. `seeds.json`, every `MANIFEST.txt`, every `eval_*.parquet`); raw simulator outputs and intermediate run caches are ignored. | `regen_all.py` writes here; downstream scoring scripts read parquets directly. |
| `benchmark/paper_data/` | TSV + Markdown | Curated mirror of `benchmark/data/`, organised by paper section (`01_headline_piao` through `06_methodology`). Generated by a converter so reviewers can inspect every number in Excel / R / Python without a Parquet reader. | The paper text cites file paths in this tree. |

The migration to a TSV mirror happened at the 1.0 cut. A working-tree
view that shows `benchmark/data/*.parquet` as "deleted" is a
checkout staleness artefact, not a missing-artefact problem — the
parquets are committed in HEAD and `git restore -- benchmark/data`
brings them back. The Linux-side full regen (per the Track 1 plan)
rewrites the parquets and re-derives the TSVs in lockstep.

---

## 9. Phases (this protocol drives all phases)

This protocol is consumed by the four execution phases in
[../../is-this-a-good-cryptic-pike.md](../../is-this-a-good-cryptic-pike.md):

| Phase | Reads from this protocol | Writes |
|---|---|---|
| 1 | (this file is Phase 1) | this file |
| 2 — Study 2 re-runs + variance | §4.1, §4.2, §4.3 (tile + Stouffer tuned methylKit row), §6 | `scripts/wilson_bootstrap_ci.py`, `scripts/methylkit_stouffer_combine.R`, updated `_runs/`, updated `HEAD_TO_HEAD.md` |
| 3 — Study 3 real-data re-runs | §3.2, §4.3, §5 R3–R6 | re-run epykit (engine tag v0.7.5-phase3-engines-frozen); fix `convert_fwd.log`; verify DSS-from-scratch; top-K tables |
| 4 — Rewrite | §5 all rules; §6 variance; §7 metrics | new `EXECUTIVE_SUMMARY.md`, `paper/paper.md`, `report/REPORT.md` |
| 5 — Independent audit | (the audit reads this protocol to check the rewrites for compliance) | audit report |

---

## 10. Limitations explicitly disclosed (these go into paper.md verbatim or near-verbatim)

1. **Simulator instance.** Studies 1 and 2 use one fixed Piao 2021 simulated dataset. Simulator variance is not captured; Wilson + bootstrap CIs capture evaluation variance on the data we have.
2. **One cohort, one genome.** Study 3 is one tissue × one genome × six biological samples. We treat it as an existence proof of correctness on real biological data, not a generalisation.
3. **Study 1 panel baselines are transcribed.** The eight-tool panel uses Piao 2021's published TPR/FPR values, not local re-runs at current tool versions.
4. **`mincov` choice diverges from Piao.** Studies 2 and 3 use `mincov = 10` for both tools; the original Piao methylKit script used `mincov = 0`. This is a fairness-driven uplift; it favours stability over raw published-baseline-matching. Study 1 (transcribed) preserves the original `mincov = 0`.
5. **Bugs surfaced during benchmarking.** A `fisher` calibration bug and an F-distribution df-reference bug in pooled-dispersion modes were discovered while running this benchmark and fixed in epykit 0.7.2. The default `lr / site` mode was unaffected; Study 3's headline was re-verified to be unchanged after the fix.
6. **DMR-engine architectures differ.** Fixed-tile callers (methylKit-tile, epykit `dmr_tile`) and variable-width chain-merge callers (epykit `chain_merge`, DSS `callDMR`) solve subtly different problems. We disclose this in every DMR comparison cell.

---

## 11. Addenda

(Any deviation from this protocol after the locked commit is recorded below, with date, what changed, and why.)

— *Empty at lock time.*
