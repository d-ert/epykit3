# Handoff — epykit Bioinformatics / Genome Biology resubmission

**Branch:** `gb-resubmission-scaffolding` (not yet merged; stays open through resubmission).
**Reviewer report:** see `~/.claude/plans/here-is-a-complete-composed-mist.md` for the full plan with the 11 major (M1–M11) and 14 minor (m1–m14) issues. The reviewer recommended **Bioinformatics / Genome Biology (Software track)**, not Nature Methods. Strategy: hybrid — ship to Bioinformatics now, park Nature Methods novelty as a separate research thread.

## 1. Where we are right now (as of 2026-06-05)

Tree is clean. Branch carries 8 commits past A.1, all already on disk.

```
769c749 docs(paper): F.1-F.4 -- abstract + sections rewrite for Linux rerun
ceaae1b docs(paper): §2.5 -- promote DF_PHI_FLOOR=50 derivation to prose (M11a)
cd5da13 feat(benchmark): M3 dual-truth re-scoring of external baselines
6c9c359 feat(benchmark): integrate Linux rerun (M4/M5/M9/m13) + fix DMC resume-cache bug
6502cbb feat(benchmark): A.5 regen_all.py --run-all orchestrator + Linux runbook
12d5143 build(benchmark): A.4 Docker + renv containers + engine-regression CI lane
5a071b9 feat(benchmark): A.3 dmrseq + BSmooth baselines, phi-sweep simulator, truth_mode scoring
b62626a feat(benchmark): A.2 null calibration k=1000 + Q-Q + KS + EB prior validator
```

## 2. Reviewer scorecard (current state)

| ID | What | Status |
|---|---|---|
| **M1** | lr+ docs/impl mismatch | **Reviewer wrong** — all four knobs implemented; discoverability fix deferred (user chose to skip "B") |
| M2 | "~2×" claim ill-defined | ✅ Abstract now quotes absolute TPR 0.564 vs 0.302 |
| **M3** | Truth labelling asymmetry | ✅ Half-fix done — `M3_truth_mode_comparison.md`; epykit's threshold-truth row is a Linux-side TODO |
| **M4** | Linux 12–68× hides single-thread methylKit | ✅ Honest 33× on Linux + mc.cores=8 (methylKit 5.9× scaling) |
| **M5** | k=20 calibration insufficient | ✅ Exhaustive 10-partition calibration on GSE263850 (= complete null universe at n=6); mean p 0.506, frac<0.05 = 0.047, KS D = 0.051 |
| **M6** | 41 D-state benchmark artefacts | ✅ Closed in `6c9c359` |
| M7 | renv.lock / Dockerfile | ✅ Already done in `12d5143` |
| M8 | Tuned-vs-default asymmetry | ✅ lr+ empirically demoted by own rerun (14× FPR inflation for +7 pp TPR, lower F1+AUROC); framed as research knob in abstract + §4.2 |
| **M9** | DSS beats epykit by ~25 pp | ✅ Substantially defused via `dis_merge` sweep — now 10 pp (77.3 % vs 87.5 % any-bp), 100 % direction agreement on 713/713 overlapping DMRs |
| M10 | Simulator φ underdispersion | ✅ Already disclosed in §4 framing; own rerun empirically validates lr+ FPR inflation under realistic dispersion |
| M11(a) | DF_PHI_FLOOR=50 unjustified | ✅ §2.5 new paragraph (`ceaae1b`) |
| **M11(b)** | sep_threshold no ROC | ✅ Zero candidate sites on GSE263850 across thresholds 0.7–0.95; "rare-event safeguard" framing in §4.3. **Bonus: cache-key bug found and fixed.** |
| m13 | dis_merge sensitivity | ✅ §3.3.3 table |
| m10 | print() guard test | ✅ Pre-existing as `tests/test_no_print_outside_cli.py` |

## 3. What's still pending

### Linux-side (compute / methylstore-needed)

| # | Task | Effort | Priority |
|---|---|---|---|
| 1 | **TCGA tumour/normal WGBS cohort** (LUAD or COAD, n ≥ 6+6) | 12–24h end-to-end | high — single biggest remaining reviewer ask |
| 2 | **M3 full per-sample-noisy-count re-scoring** (if `amp.coverage*.txt` still on the Linux box; rerun-bundle README listed them as "Excluded") | 30 min | medium |
| 3 | **M3 epykit threshold-truth row** — re-run `_epykit_scoring` with `truth_mode="threshold"` on the simulator methylstore (Windows half-fix only covered methylkit/dss/dss_nosmooth) | 1h | medium |
| 4 | **M10 φ-sweep** at φ ∈ {0.5, 1.5, 3.0, 5.0} | ~8h | low (already empirically validated via lr+ demotion) |
| 5 | **External baselines** locally re-run (dmrseq, BSmooth, comb-p, RADMeth) on simulator + Study 3 — infra exists per commit `5a071b9` | days | defer to revision-2 |

### Windows-side / paper polish

Three deferred items, all flagged in `769c749` commit message:

1. **§3.3.2 Table 5a "headline coord-overlap"** — still has pre-rerun chain_merge n=702 / n=940 vs paper-813. Update to rerun n=852 / n=1139 vs DSS-922.
2. **§3.4 held-out simulator** — single-seed seed=2026000 numbers in the body text. The 20-seed median+IQR is in `eval_simulator_intrinsic_iqr.parquet`; switch the table to medians.
3. **`claim:` inline-comment validator pass** — run `regen_all.py` (or whatever the validator entrypoint is) to refresh the asserted values throughout the paper.

### Skipped per user instruction

- **B / A.2 power_stack discoverability fix** — user said "don't bother". The lr+ dispatch lives at `src/epykit/tl.py:494-532` while the engine code is in `dmc.py`; the reviewer reading dmc.py in isolation missed it. If a future reviewer hits the same confusion, consider relocating `_resolve_power_stack` into `dmc.py`.

## 4. Where things live

### Data
- **Linux rerun staging (gitignored, on disk):** `benchmark/rerun_outputs_2026-06-03/` — raw drop, ~3 GB with bulk artefacts
- **Committed headline tables:** `benchmark/data/` — the resubmission-bundle artefacts. Tracked subset only; `.gitignore` whitelists by name. Largest tracked file is `null_calibration/gse263850/lr_pvalues.parquet` at 7 MB.
- **Bulk excluded:** `dmltest_per_cpg.tsv.gz` (460 MB), `_ingest_store/` (2 GB), per-seed methylkit/dss TSVs, methylstore caches

### Key Linux rerun deliverables (now in `benchmark/data/`)
- `study1b_simulator/eval_seed_iqr.parquet` — 4-row median+IQR for 20 seeds
- `study1b_simulator/eval_simulator_intrinsic_truth_both_*` — M3 dual-truth result
- `study1b_simulator/M3_truth_mode_comparison.md` — paper-ready M3 framing
- `study3/chain_merge/dmr_chain_merge.parquet` — 852 DMRs at dis_merge=100
- `study3/dss/dmr_dss.csv` — 922 DSS-from-scratch DMRs
- `study3/comparisons/epykit_vs_dss/headline.json` — recall/precision/J≥0.5/direction
- `multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_vs_dss_sensitivity.csv` — m13 table
- `multi_thread_and_chain_sweep/methylkit_multicore/methylkit_multicore_timing.md` — M4 honest ratio
- `null_calibration/gse263850/lr_calibration_report.md` + `lr_qq.png` — M5 result
- `sep_threshold_roc/FINDINGS.md` + 2 CSVs — M11(b) result with cache-bug context

### Code
- **Manuscript:** `benchmark/paper/paper.md` (also `paper.docx` rendered)
- **Report:** `benchmark/report/REPORT.md` + `methods_appendix.md`
- **Scoring:** `benchmark/scripts/_epykit_scoring.py` (supports `truth_mode`); driver `eval_simulator_intrinsic.py` (now supports `--truth-mode {intrinsic,threshold,both}` + `--sim-root`)
- **Comparison:** `benchmark/scripts/compare_epykit_to_dss.py` (paths now relative to `epykit3` root via `Path(__file__).parents[1]`)
- **Null calibration:** `benchmark/scripts/run_null_calibration.py` (k=1000 + Q-Q + KS)
- **Engine source:** `src/epykit/dmc.py` (lr+ knobs at 947-983 / 775-820 / 2322-2500 / 2596-2640)
- **Cache fix:** `src/epykit/tl.py:551-575` (added 6 lr+ knobs to `resume_sig`)

### Reference docs
- **Full plan:** `~/.claude/plans/here-is-a-complete-composed-mist.md` (every reviewer item + execution status)
- **Project rules:** `CLAUDE.md` at repo root
- **Reviewer's report:** in the original conversation; not committed (the file `here-is-a-complete-composed-mist.md` summarises it with verdicts)

## 5. Important context (things easy to miss)

### The DMC resume-cache bug (now fixed)
Reviewing M11(b) surfaced that `tl.py:551-575`'s `resume_sig` did not include `sep_threshold`, `sep_fallback`, `neighbour_combine`, `neighbour_bp`, `fdr_method`, `power_stack`. Any user calling `ep.tl.dmc(..., resumable=True)` while varying these knobs silently hit a stale cache. Fixed in `6c9c359`. 22 existing resume/cache tests still pass. **No regression test was added per user instruction** — if you want one later, see plan §0.7 for the shape.

### lr+ has been honestly demoted
Your own 20-seed rerun shows lr+ trades 14× FPR for +7 pp TPR with *lower* F1 (0.746 vs 0.796) and *lower* AUROC (0.907 vs 0.928). The abstract, §4.2, and §4.3 now frame lr+ as a "research knob, not a recommended default" — bare lr is what the abstract numbers are reported against.

### DSS DMLfit.multiFactor has no multi-core option
Verified in DSS 2.58.0. So the 33× speedup vs DSS is not eroded by parallelizing DSS — single-thread by construction, not by choice. This is a clean defense, but **note the version**: PROTOCOL.md still pins 2.12.0; resubmission letter should reconcile.

### M5 exhaustive enumeration is *stronger* than the reviewer asked for
The reviewer asked for k ≥ 1000 random shuffles. For n=6 (3v3) the universe is C(6,3)/2 = 10 distinct partitions; we enumerated all 10. k=1000 random shuffles would draw the same 10 with replacement ~100× each, conveying no extra information. The abstract and §3.5 say this explicitly — pre-empt the "but the reviewer said 1000" comment.

### The simulator can't test sep_threshold
At cov ≥ 5 with 3v3, no quasi-separated CpGs exist. Same is true of real GSE263850 (0 candidate sites). sep_fallback fires only at pathologically low coverage that the standard coverage filter (≥10×/sample) removes upstream. Default 0.9 has zero effect on any reported number.

## 6. Open decisions waiting on you

1. **TCGA cohort go/no-go for the first submission.** Strongest single addition; 12–24h on Linux. Worth doing, but the manuscript is already credibly submittable without it.
2. **Polish the three §3.3.2 / §3.4 / claim-validator items before submitting**, or accept the current state? They're explicit in the latest commit message so a follow-up commit is cheap.
3. **DSS version reconciliation:** PROTOCOL.md 2.12.0 vs the multi-core check on 2.58.0. Either update PROTOCOL.md or add a one-line "verified no multi-core path in 2.12.0 → 2.58.0 release notes" check.
4. **Submission target — Bioinformatics vs Genome Biology Software track.** The plan was hybrid; pick one for the cover letter.

## 7. Fresh-session prompt template

If you're starting a new Claude session, paste this:

> I'm working on the epykit Bioinformatics resubmission on branch `gb-resubmission-scaffolding`. The full status is in `HANDOFF.md` at the repo root and the detailed plan at `~/.claude/plans/here-is-a-complete-composed-mist.md`. Read both before responding. The reviewer's 11 major issues are mostly closed (see scorecard in HANDOFF). Open items: [paste your specific ask].
