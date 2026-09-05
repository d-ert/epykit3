# epykit Pre-Submission Hardening — Design

**Date:** 2026-05-27
**Status:** Design, awaiting user review
**Companion files:** `D:\Coding\Projeler\methyl_lib\benchmarkin_merges\FINAL_REPORT\PROTOCOL.md`, `EXECUTIVE_SUMMARY.md`, `paper\paper.md`
**Trigger:** preparing epykit v0.7.2/0.7.3 + benchmark for journal submission. Three parallel audits (DMC engines, DMR/QC engines, benchmark methodology) identified findings that the current report and protocol do not fully address.

The user is in **aggressive cleanup** mode (drop/rename anything broken or misnamed) and wants benchmark hardening to come first, then bugs.

---

## 1. Scope and non-goals

**In scope**
- Amend `PROTOCOL.md` to close three gaps: method-shaped truth, no null calibration, headline numbers depend on buggy code paths.
- Catalogue every statistical/biological bug surfaced by the audits, with severity, file:line, and the exact patch shape.
- Decide which engines/knobs to remove or rename pre-submission, with deprecation pathway.
- Add the minimum set of new benchmark scripts needed to satisfy the amended protocol.
- Spell out the paper-level framing changes (abstract, intro, discussion) so the strength of the prose matches the strength of the evidence.

**Out of scope (this design)**
- Multi-dataset real-WGBS validation beyond GSE263850 (deferred per `EXECUTIVE_SUMMARY.md` "What's next").
- A second independent simulator beyond a Piao re-implementation (the re-implementation gives held-out + multi-seed, which is the marginal value reviewers will demand; an entirely different simulator like dmrseq is post-revision work).
- GPU backend (`_glm_gpu.py`) audit — gated behind extras, not benchmark-critical.
- DMR HMM rebrand as a fully fitted HMM (Baum-Welch) — too much code for the submission window; instead we rename it.

---

## 2. The three remaining gaps and how this design closes them

### 2.1 Gap A — Truth is method-shaped (the AUROC tautology)

**Current state.** `_make_truth.py` defines truth as `|mean_β_treat − mean_β_ctrl| ≥ 0.20` computed from the 25× simulator outputs. This is the same statistic that `lr`/`lr+` ranks by. AUROC of 0.9999 is the inevitable fixed-point of correlating something with itself plus a thin layer of noise. The protocol §3.1 discloses this but does not fix it.

**Fix.** Re-implement Piao et al.'s binomial DMC + DMR simulator in Python (≈ 200 lines, ~1 day per user confirmation). Add `benchmark/scripts/simulate_piao.py` that takes a seed and writes the same `amp.*.txt` schema. The simulator's intrinsic `is_dmc` flag (drawn from a `meth_diff ~ U(0.2, 1.0)` distribution per Piao §2) becomes the truth, replacing the threshold-reconstruction.

Once the simulator exists:
- Run N = 20 independent draws at the headline cell (cov=10×, n=3v3). Compute TPR/FPR/F1/AUROC across the 20 draws; report median + IQR. This gives **simulator variance**, which Wilson/bootstrap over CpGs does not.
- Lock the `lr+`, `chain_merge`, `power_stack="auto"` defaults using the Piao-as-distributed data ("development set"), then run the simulator-generated grid *once* with frozen defaults ("test set"). This kills the "you tuned on the test set" objection (audit findings F4/F5/F6).
- For Studies 1 and 2 headline tables: keep Piao-as-distributed as the comparator-rich panel (with the disclosed circularity), and add a parallel column with simulator-intrinsic truth on Piao-as-distributed reads, where available, to show the gap is small.

**Deliverable.** `benchmark/scripts/simulate_piao.py` + a new section §6.5 in PROTOCOL.md describing the simulator (seed handling, parameter mapping to Piao 2021, validation that the simulator reproduces Piao's marginal statistics).

### 2.2 Gap B — No null calibration run

**Current state.** Nothing in `FINAL_REPORT/` runs group-A-vs-group-A. The FPR claim ("100-600× lower than competitors") is not anchored to a calibration check.

**Fix.** Add `benchmark/scripts/run_null_calibration.py`. For each (engine, coverage) cell, randomly shuffle group labels (k=20 shuffles) and re-run the engine. Report:
- Empirical FDR at nominal `q < 0.05` per engine, with Wilson 95% CI.
- Inflation factor: empirical FDR / nominal threshold.
- QQ plot of observed p-values against Uniform(0,1) on the shuffled data.

Acceptance criteria for the paper: report observed FDR ± CI for every engine in Table S-Calib. If `lr` (default) lands at 0.05 ± CI on the simulator, the FPR-tightness story becomes "well-calibrated and conservative on real null cells"; if it lands at 10⁻⁵, the story becomes "the test is over-conservative on this simulator's underdispersed noise" (which is honest and weaker but defensible).

This same script can run on Study 3 (real GSE263850 data) by re-shuffling labels among the 6 samples (C(6,3)/2 = 10 unique splits). Calibration on real data is more persuasive than on simulated.

**Deliverable.** Script + Table S-Calib + paragraph in §3 of the paper.

### 2.3 Gap C — Headline numbers depend on buggy code paths

**Current state.** Several bugs surfaced in the audits affect numbers the headline reports. Re-running with the bugs unfixed reproduces the current numbers; fixing the bugs and re-running changes them. Either way, the paper has to ship the correct numbers and document the fix.

**Fix.** Section 3 of this design enumerates every bug with severity, file:line, and patch. Section 4 sequences bug fixes *before* the locked benchmark re-run. Section 5's new `regen_all.py` reads the patched code's outputs and produces every table from a single source of truth.

Per protocol R8 (every claim traceable to a script), each affected number gets a `<!-- source: scripts/X.py -->` comment in the paper. Per protocol §10.5, the bug-fix manifest (the table in §3 below) goes into the paper's Limitations section verbatim, owning the discovery rather than hiding it.

---

## 3. Bug-fix manifest (the §4.5 protocol addendum)

Severity scale: **P0** = changes headline numbers; **P1** = changes secondary/supplementary numbers or is publicly indefensible; **P2** = cosmetic, doc-only, or behind a flag.

### P0 — must land before the locked benchmark re-run

| # | Bug | File:line | Patch shape | Headline impact |
|---|---|---|---|---|
| P0-1 | `pvalue` silently overwritten by `pvalue_combined` when `lr+` engages, breaking downstream BH / empirical FDR / DMR engines | `tl.py:598-603`, `dmc.py:2685,2740,3042` | Keep `pvalue` as raw; expose `pvalue_combined` as a separate column throughout; update `apply_multiple_testing_correction` and `empirical_fdr_for_dmc` to take a `pvalue_col` arg with default `"pvalue"`; emit a warning if both columns are present and the caller did not specify | `lr+` q-values in Studies 1 & 2 will change. `empirical_fdr_for_dmc` results in DMR Study 1 will change. |
| P0-2 | DMR empirical-FDR denominator pools null DMRs across all permutations instead of using `n_perm` | `dmr.py:1418-1426` | Replace `(counts + 1) / (total_null + 1)` with rank-quantile comparison: per observed DMR with rank k in obs, compare to kth-best null DMR across perms; OR implement Storey-style π₀ · E[#null with p≤t] / #obs with p≤t. The simpler fix is the rank-quantile one. | All DMR empirical-FDR numbers in any cell where `empirical_fdr=True` was used. |
| P0-3 | `eb` is the actual default for `dispersion` in `tl.dmc`, but docstring (`tl.py:380-384`) and PROTOCOL.md §4.1 disagree (say `site`) | `tl.py:380-384`, `PROTOCOL.md` §4.1, `EXECUTIVE_SUMMARY.md`, `paper.md` | **User decision:** keep `eb` as default. Fix the disagreement by updating docstrings, PROTOCOL.md §4.1, EXECUTIVE_SUMMARY recommendation, and paper.md Methods to document `eb` as the default. Verify all headline runs were in fact `eb`; if any were `site` because the protocol said so, re-run with `eb`. | Documentation only if all runs used the code default; full headline regen if any cell used `site`. |
| P0-4 | Adaptive F(1, df_phi) at small `w_eb` is severely conservative | `dmc.py:967-974` | Floor `df_phi` at a minimum (e.g., 30) when `phi_eff > min_dispersion + ε`. Since `eb` is the default per P0-3, this fix is **mandatory** — we cannot defer by gating adaptive-F to non-`eb` modes. Validate the floor value with a calibration check (null simulation) before committing. | All headline DMC numbers in `eb` mode. This is the bug behind the unusually low FPR in `eb` mode; fixing it will likely raise `lr` FPR somewhat. Surface the per-cell delta in the bug-fix audit. |
| P0-5 | Tile-merge Stouffer combination has three math errors (one-sided isf on two-sided p; wrong denominator on chains > 2; no direction term) | `dmr.py:1053-1075`, `_merge_adjacent_tiles` | Maintain running `(sum_z, sum_w_sq, sign_majority)` across the chain; combine as `Σwᵢzᵢ / √Σwᵢ²` at the end; reject chains where direction is mixed | DMR tile results in Studies 1 & 2 when adjacent-tile merging is on. |
| P0-6 | Stouffer in `neighbour_combine` assumes independence of adjacent CpGs | `dmc.py:2451-2453, 2571` | Two options. (a) Document explicitly that `neighbour_combine` controls FDR only via the `min_sign_agreement` and `require_focal_signal` gates, *not* via the Stouffer null. (b) Implement Brown's correction using a lag-1 autocorrelation estimate of z-scores per chromosome. We pick (a) for the paper window (cheap, honest), and file (b) as a `v0.8` issue. | `lr+` FPR claims become more accurate; some headline FPR numbers may rise. |

### P1 — fix before publication; outside the locked re-run is acceptable if numbers don't change

| # | Bug | File:line | Patch shape |
|---|---|---|---|
| P1-1 | Fisher two-sided is doubled one-tail, not Fisher's exact two-sided | `dmc.py:253-256` | Route through `scipy.stats.fisher_exact(..., alternative="two-sided")` in vectorised form, OR document the doubled-mid-tail convention. Either is defensible; pick one and state it. |
| P1-2 | `bb_lr` uses `df_resid` instead of `df_phi` for F reference | `dmc.py:1656` | Use `df_phi` like the lr path. Re-run bb_lr cells if any are in the headline. |
| P1-3 | `meth_diff_ci_{lo,hi}` uses Welch-normal Wald CI; `newcombe_diff_ci` exists but is unwired | `dmc.py:1907-1912`, `_glm.py:923-963` | Wire `newcombe_diff_ci` into `_process_one_chromosome` for `lr`/`score` tests. Welch CI is fine for `welch_t`/`logit_t`. |
| P1-4 | GLM patsy contrast coding uses alphabetical reference silently; user has no way to set reference level explicitly | `_glm.py:42-197` | Log `design_info.column_names` and the chosen reference level; expose `reference_level` kwarg in `tl.dmc(formula=...)`. |
| P1-5 | IRLS non-convergence is silent (only separation is NaN'd) | `_glm.py:234, 297, 321-324` | Use the existing `converged` array to mask non-converged sites (NaN their Wald stats). Log fraction non-converged. |
| P1-6 | DMR empirical FDR ignores design (paired samples, covariates) and breaks at n_treat=1, n_ctrl=1 | `dmr.py:1370-1377` | Detect paired design (single covariate column with strict pairing) and shuffle within strata. Refuse to run with `n_treat=1, n_ctrl=1` and surface a `ValueError`. |
| P1-7 | DVC uses Bartlett on bounded U-shaped beta values | `dvc.py:58-86` | Implement two-pass Brown-Forsythe (one pass for medians, one for absolute deviations, F-test). The "Welford budget" excuse is moot — per-CpG sample counts are small. |
| P1-8 | HMM ships with hardcoded `state_means` and no fitting; rename to match | `dmr_hmm.py:33-41`, `_hmm.py` | Rename `dmr_hmm` → `dmr_segment` and the function `call_dmr_hmm` → `call_dmr_rule_segment`. Document explicitly that it is a rule-based segmenter with fixed priors, not a fitted HMM. Posterior-based region calls + per-segment Stouffer combine for region p-values (replacing `mean(qvalue) >= alpha`). |
| P1-9 | Sex check fails on single-sex cohorts (largest-gap 1D clustering on unimodal data) | `qc.py:474-495` | Compute Hartigan's dip test for bimodality first; if unimodal, fall back to the fixed 0.25 threshold and emit a warning. Document that Y-coverage is the preferred signal. |
| P1-10 | `_storey_pi0` hardcodes `lam=0.5` and can return 0 with all-significant input | `dmc.py:2419-2421` | Clamp `pi0 >= 1/n` (Storey's standard floor). Document that this is the plug-in estimator, not the smoother. |
| P1-11 | `log2_odds_ratio` column name is misleading for GLM/bb_lr backends (it's the logit coefficient, not log2 of OR) | `dmc.py:1673`, multiple sites | Rename per backend: `log2_odds_ratio_pooled` for count-pool tests; `coef_treatment_log2` for GLM tests. Document conventions. |

### P2 — doc / API hygiene

| # | Bug | File:line | Patch |
|---|---|---|---|
| P2-1 | `_BETA_EPSILON = 1e-6` is used for three different statistical purposes with no per-context justification | `dmc.py:208` | Split into `_LOGIT_EPS`, `_POOL_EPS`, `_OR_EPS` with comments naming each convention (Anscombe 3/8, Haldane 0.5, etc.). |
| P2-2 | `pct_sig` knob acknowledged dead at strict alpha | `dmr.py:717-721` | Drop from API (aggressive cleanup mode). Migration path: pass-through with deprecation warning for one release. |
| P2-3 | CGI shore/shelf widths hard-coded | `annotate.py:1327-1328` | Expose as kwargs (`shore_bp=2000`, `shelf_bp=4000`). Cite Irizarry 2009 / Bibikova 2011 in docstring. |
| P2-4 | `call_dmr_hmm` returns NaN p/q-values for every DMR; downstream consumers break | `dmr_hmm.py:135-136` | Compute per-segment Stouffer-combined p-values from constituent CpG p-values, apply BH per chromosome. |
| P2-5 | `n_perm=100` default for empirical FDR is too low | `dmr.py:1319` | Raise default to `1000`. |
| P2-6 | `_auroc` uses `1 - pvalue` as score; collides at numerical floor (`pvalue ≤ 1e-16`) | `benchmark/scripts/evaluate.py:line` | Use rank of `-log10(pvalue)` to break ties. AUROC reported to 4 decimal places needs this. |
| P2-7 | epykit version mismatch (paper.md says 0.7.1, REPORT.md says 0.7.2) | `paper.md:146`, `REPORT.md:5` | Tag a release `v0.7.3-paper`. Pin the paper to it. |
| P2-8 | Manuscript path inconsistency (`ground_truth/make_truth.py` vs actual `benchmark/scripts/_make_truth.py`) | `paper.md:141` | Reconcile. |

---

## 4. API cleanup (aggressive mode)

Per user direction. Each item below moves before submission, not after.

### Drop from public API

| Item | Reason | Migration |
|---|---|---|
| `tl.dmc(test="logit_t")` | Documented by code itself as miscalibrated near β=0/1; bound to be flagged | Map calls to `welch_t` with a `DeprecationWarning`. Internal helpers stay. |
| `tl.dmc(test="bb_lr")` | TPR < 8 % at n ≤ 4 in the benchmark; bug P1-2 affects what little signal it has | Keep code but move to `epykit.experimental.bb_lr`; emit `DeprecationWarning` on import. Remove from `tl.dmc` dispatcher. |
| `pct_sig` knob in `dmr_chain_merge` | Dead at strict alpha by your own admission (`dmr.py:717-721`) | Accept argument with `DeprecationWarning`, ignore. Remove in v0.8. |
| `tl.dmc(test="cmh")` if not in any headline table | Unusual stratification with unclear assumptions | If used: document properly. If not used: move to `epykit.experimental`. *Action: confirm headline usage before deciding.* |

### Rename

| Old | New | Reason |
|---|---|---|
| `dmr_hmm.call_dmr_hmm` | `dmr_segment.call_dmr_rule_segment` | It's a rule-based segmenter with fixed priors, not a fitted HMM. Calling it HMM invites the "did you Baum-Welch?" question you can't answer yes to. |
| `log2_odds_ratio` (in GLM/bb_lr output) | `coef_treatment_log2` | The number reported by GLM is the logit coefficient, not log₂(OR). Different scales. |
| `log2_odds_ratio` (in pooled-count output) | `log2_odds_ratio_pooled` | Distinguish from the GLM column. |

### Documentation overhaul

- Every test family (`lr`, `score`, `glm`, `welch_t`, `fisher`, `cmh`, `bb_lr`) gets a docstring section: **Assumptions** / **Boundary handling** / **Multiple-testing** / **When to use this engine instead of `lr`**.
- `_BETA_EPSILON` split into three named constants with comments naming each statistical convention.

---

## 5. New benchmark scripts (the minimum needed)

All paths relative to `D:\Coding\Projeler\methyl_lib\benchmarkin_merges\FINAL_REPORT\scripts\`.

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `simulate_piao.py` | Re-implementation of Piao et al.'s binomial simulator. Args: `--seed`, `--n_cpgs`, `--n_dmcs`, `--coverage`, `--replicates`, `--phi`. Output schema matches `amp.coverage=K.sampleN.txt`. | `_seeds.json` for reproducibility | `data/study1b_simulator/seed=K/amp.*.txt` and `is_dmc.parquet` (intrinsic truth) |
| `run_null_calibration.py` | For each engine in {lr, lr+, fisher, welch_t, chain_merge, tile, dmr_segment}, run k=20 label-shuffles per cell. Compute empirical FDR at nominal q<0.05 with Wilson CI. | Cached `_runs/` methylstores | `data/study1_null/calibration.parquet` and Table S-Calib |
| `wilson_bootstrap_ci.py` | Already promised by PROTOCOL §6 but listed as "to be added in Phase 2." Computes Wilson 95 % CI on TPR/FPR per cell; bootstrap percentile CI on AUROC/F1 (B=1000). | Existing single-run outputs | Adds CI columns to `eval_summary.parquet` |
| `regen_all.py` | Drives `simulate_piao.py` (N=20 seeds, headline cell only) + null calibration + the existing study runners. Verifies that every numeric claim in `paper.md` matches the parquets. | All locked outputs | Pass/fail per claim; non-zero exit if any mismatch |
| `methylkit_stouffer_combine.R` | Already promised by PROTOCOL §4.2 but not yet written. Adjacent-3-CpG Stouffer combine after `calculateDiffMeth`. | methylKit per-scenario TSVs | Tuned-methylKit per-scenario TSVs |
| `bug_fix_audit.py` | Diff every headline number between pre-fix and post-fix epykit runs. Reports per-cell delta with a sign so reviewers can see whether the fix helped or hurt epykit. | Pre-fix snapshot + post-fix run | `data/audit/bug_fix_deltas.parquet` and a paragraph in Limitations §10.5 |

---

## 6. Paper-level framing changes

### Abstract

Current abstract framing: "epykit is competitive with the strongest baselines and 10-1000× faster" (paraphrased from EXECUTIVE_SUMMARY). Three additions:

1. One sentence acknowledging the simulator-realism caveat **in the abstract**, not in §4 only: "Low-coverage TPR advantages observed on the underdispersed Piao 2021 simulator (φ ≈ 0.4) are not expected to transfer at the same magnitude to overdispersed real WGBS (φ ≈ 1.5–5); we report agreement with methylKit and DSS on real GSE263850 data as the principal real-data evidence." — already half-stated in `EXECUTIVE_SUMMARY.md` "What we do not claim", but the abstract is the only thing many reviewers read.
2. Replace "best-in-class" with "matches or exceeds the strongest baselines" (per audit finding F22).
3. Add the bug-fix manifest sentence ("Two calibration bugs were surfaced and fixed during this benchmark; numbers reported are post-fix and the per-cell deltas are documented in Supplementary Table S-Fix").

### Methods

- Adopt the PROTOCOL's §4 parameter freeze verbatim.
- Add §3.X for the Piao-simulator re-implementation (validation against Piao's marginals).
- Add §3.Y for the null calibration design and results.
- Move the "tile → chain_merge pivot" narrative from `idk_if_needed/` into Methods, per PROTOCOL R4. Frame it as a methodological observation: fixed-tile and chain-merge callers solve different problems, and on focused biological DMRs the chain-merge family dominates.

### Results

- Every TPR/FPR/F1/AUROC number gets a Wilson or bootstrap CI in the same cell. Bare point estimates are removed.
- Default-vs-default headline (PROTOCOL R1). Tuned-vs-tuned sensitivity panel below, clearly labelled.
- AUROC reported intra-epykit only (per audit F15). If we re-run methylKit and DSS to capture p-value vectors, then we can extend AUROC cross-tool — open question.
- Add Table S-Calib (null calibration), Table S-Sim (multi-seed simulator results), Table S-Fix (bug-fix per-cell deltas).

### Discussion

- The §4 underdispersion caveat moves up: it's the central methodological honesty of the paper and should appear by the second paragraph of the Discussion, not at the end.
- The bug-fix manifest is owned: "We found and fixed N bugs while running this benchmark. The Limitations section enumerates them; the comparison to the published Piao 2021 baselines is re-derived using the post-fix code." — this is more persuasive than burying it.

---

## 7. Execution sequence

Order matters because the locked benchmark re-run must use post-fix code.

1. **Apply all P0 fixes.** No headline re-run until P0 is clean. Commit each P0 fix with a CHANGELOG entry referencing this design.
2. **Write `simulate_piao.py`.** Validate marginals (Piao reproduces Piao within Monte Carlo noise; the implementation reproduces Piao within Monte Carlo noise).
3. **Write `run_null_calibration.py`.** Run on Piao-as-distributed first (no simulator yet needed for this part). Validate that `fisher`-with-the-old-bug gives observed FDR = 1.0 (sanity check).
4. **Run Phase 2 of the PROTOCOL** (Study 2 re-runs + variance) with post-fix code. Tag `v0.7.3-paper` after this.
5. **Run Phase 3 of the PROTOCOL** (Study 3 real-data re-runs at 0.7.3 to lift the 0.6.0 mismatch).
6. **Run the Piao-simulator multi-seed grid** (N=20 at headline cell + frozen-defaults sweep at all cells). Tag results as Study 1b.
7. **Apply P1 fixes** alongside the rewrite phase. Where P1 affects supplementary numbers, regenerate those tables.
8. **Apply P2 fixes** anytime.
9. **Rewrite the paper** per §6 of this design.
10. **Run `regen_all.py`** as the acceptance gate. Any mismatch between paper claims and parquet contents fails the gate.

Estimated wall time on a single laptop: P0 fixes ~2 days; simulator + null calibration ~2 days; locked re-runs ~1 day of compute; paper rewrite ~3 days. Roughly one calendar week to a submission-ready state, assuming no surprise findings during the bug fixes.

---

## 8. Open questions for the user

These are not blocking — I have a default for each — but they would tighten the design.

1. **Is the `epykit3/benchmark/` folder being kept in sync with `benchmarkin_merges/FINAL_REPORT/`, or is the latter authoritative and the former should be removed/archived?** The duplication is confusing for anyone navigating the repo. *Default: archive `epykit3/benchmark/` to `_legacy_benchmark/` and add a README redirect.*
2. **Is GSE263850 the only real dataset within budget for this submission?** A second real cohort would let us close the §6.1 "one tissue × one genome" caveat. *Default: yes, defer second cohort to revision response.*
3. **Are you willing to drop AUROC from cross-tool comparisons if we don't re-run baselines for p-value vectors?** It's currently used as evidence of epykit superiority despite baselines not having it. *Default: yes, drop from cross-tool, keep intra-epykit only.*
4. **What's the target journal?** Genome Biology, Bioinformatics, Briefings in Bioinformatics, and NAR have different expectations for reproducibility (Docker / Snakemake) and method novelty. *Default: design for Genome Biology / Bioinformatics tier (Docker recommended, code release required).*
5. **Should `cmh` and `score` engines be evaluated for keeping in the public API, or were they exploratory?** Per aggressive-cleanup mode. *Default: move to `epykit.experimental` if not used in the paper.*

---

## 9. Risks

- **Bug fixes may make the headline numbers worse.** Specifically, fixing P0-4 (eb conservative F, which now affects the *default* path per P0-3) will likely raise `lr` and `lr+` FPR — the eye-popping "100-600× tighter FPR" claim may shrink. Fixing P0-6 (Stouffer independence acknowledgement) doesn't change numbers but reframes the FPR-tightness story. The paper has to tell that story honestly. Mitigation: the bug-fix audit table (`bug_fix_audit.py`) makes the per-cell deltas first-class supplementary evidence; reviewers see the changes were intentional and quantified. The narrative pivot if needed: "after fixing the adaptive-F df bug, `lr` (eb) FPR is competitive with rather than tighter than baselines; TPR advantage at low coverage is unchanged."
- **Piao simulator re-implementation may diverge from Piao 2021's marginals.** If it does, we have to debug; the Piao-as-distributed grid is the fallback.
- **Null calibration may show `lr`'s FPR is mis-calibrated rather than just low.** This would invert the "FPR 100× tighter" story. Mitigation: report the calibration honestly; pivot the narrative to "lr is conservative on this simulator's noise regime; on real data, calibration agrees with methylKit to three decimal places."
- **API-removal of `logit_t` / `bb_lr` may break a user's existing notebook.** Mitigation: deprecation warning route; full removal scheduled for v0.8 (post-paper).

---

## 10. Acceptance criteria

The submission is ready when:

- [ ] All P0 fixes landed and verified by `tests/`.
- [ ] `simulate_piao.py` reproduces Piao 2021's marginal coverage and methylation-difference distributions within 2 % per coverage cell.
- [ ] `run_null_calibration.py` runs without crashing on Piao + simulator + Study 3 data.
- [ ] Every cell in every results table in `paper.md` carries a Wilson or bootstrap CI.
- [ ] `regen_all.py` exits 0 (all paper claims match parquets).
- [ ] The Limitations section enumerates: simulator instance, one cohort, transcribed panel, mincov choice, surfaced bugs, DMR engine architectural differences, and (new) any case where fixing a bug changed a headline number by > 1 pp.
- [ ] An external reader of the paper can `pip install epykit==0.7.3` and reproduce the headline cell from a single `docker run` command (or equivalent uv-based recipe).
