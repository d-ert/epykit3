# epykit Paper-Track Finish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the paper-side work for epykit 1.0 — commit the WIP simulator sweep, optionally run the 20-seed parallel column, generate the timings table, populate `claims.yaml`, rewrite `paper.md` to match spec §6 framing, gate via `regen_all.py --verify`, then merge to `main`. Runs on a `paper` branch off `p0-fixes` independently of the `1.0-prep` branch (which covers `src/epykit/**`).

**Architecture:** Single feature branch off `p0-fixes`. Compute-heavy steps (the simulator sweep) come first so writing can proceed while they run. Paper rewrites land as one commit per section so reviewers can read the diff section-by-section. Verification gate (`regen_all.py --verify`) runs after every claims-yaml change and after the final section commit. The companion `1.0-prep` branch is independent and merges separately.

**Tech Stack:** Python (polars, parquet) for data wrangling, R (methylKit, DSS) for the simulator sweep, markdown + Jinja-style claim markers for the paper, `regen_all.py` for the gate.

---

## Companion docs

- **Spec:** [`../specs/2026-06-01-library-1.0-and-paper-finish-design.md`](../specs/2026-06-01-library-1.0-and-paper-finish-design.md) — read §3 (paper track) and §1 (decisions) before starting.
- **Parent paper spec:** [`../specs/2026-05-27-paper-defendable-benchmark-design.md`](../specs/2026-05-27-paper-defendable-benchmark-design.md) — §6 has the per-section paper rewrite requirements (Abstract, Methods, Results, Discussion).
- **Phase 4 plan (predecessor):** [`2026-05-28-phase-4-locked-rerun.md`](2026-05-28-phase-4-locked-rerun.md) — Tasks 1-8 already shipped per commit log; Tasks 9-10 are what this plan picks up. Sub-checkboxes there are stale; Task 11 (the library 1.0-prep plan) is what marks them.
- **PROTOCOL:** [`/benchmark/PROTOCOL.md`](../../../benchmark/PROTOCOL.md) — tool versions, R rules, parameter freeze rules.

---

## Pre-flight verification

```powershell
# 1. Confirm branch state.
git status                                  # expect: 1 modified + 1 untracked, both paper-side
git branch --show-current                   # expect: p0-fixes

# 2. Confirm the WIP files match what we expect.
git diff --stat benchmark/scripts/eval_simulator_intrinsic.py
# expect: ~50 lines added (the --all-seeds mode + _score_one_seed helper)
ls benchmark/scripts/run_external_simulator_sweep.py
# expect: file exists, untracked

# 3. Confirm benchmark tests pass.
uv run pytest benchmark/scripts/tests/ -q

# 4. Confirm regen_all.py runs (even with empty claims.yaml -> exit 0).
uv run python benchmark/scripts/regen_all.py --verify
# expect: exit 0 (empty claims.yaml passes trivially)

# 5. Create the working branch.
git switch -c paper
```

---

## File map

| File | Role |
|---|---|
| `benchmark/scripts/run_external_simulator_sweep.py` | Currently untracked; gets committed in Task 1. Runs methylKit + DSS on every simulator seed. |
| `benchmark/scripts/eval_simulator_intrinsic.py` | `--all-seeds` mode change committed in Task 1. |
| `benchmark/idk_if_needed/*.md` | 8 historical exploration markdowns; renamed to `benchmark/docs/historical/` in Task 2. |
| `benchmark/data/study1b_simulator/seed=*/methylkit.tsv|dss.tsv|dss_nosmooth.tsv` | Outputs of the 20-seed sweep (Task 3). |
| `benchmark/data/study1b_simulator/eval_simulator_intrinsic_per_seed.parquet` | Output of `--all-seeds` scoring (Task 4). |
| `benchmark/data/study1/timings_post_phase3.parquet`, `benchmark/data/study1b_simulator/timings_simulator.parquet` | Existing input parquets for the timings table. |
| `benchmark/data/study1/timings_table.csv` | New summary CSV from Task 5. |
| `benchmark/scripts/claims.yaml` | Empty seed today; populated in Task 6. |
| `benchmark/paper/paper.md` | 789-line draft; rewritten section-by-section in Tasks 7-11. |
| `benchmark/data/audit/bug_fix_deltas.md` | Existing on disk; included verbatim into Limitations §10.5 in Task 11. |

---

## Task 1: Commit the WIP simulator sweep scripts

**Files:**
- Stage and commit: `benchmark/scripts/eval_simulator_intrinsic.py` (modified — `--all-seeds` mode)
- Stage and commit: `benchmark/scripts/run_external_simulator_sweep.py` (new — resumable methylKit+DSS sweep)

These are the Task-5 scope expansion: instead of methylKit + DSS on 1 seed, run all 20 simulator seeds. The pre-existing file (`run_external_simulator_sweep.py`) already has the resumable skip-if-exists logic and a header docstring documenting the ~3h walltime.

- [ ] **Step 1: Confirm both files exist and are in the working tree**

```powershell
git status --short benchmark/scripts/
# expect:
#  M benchmark/scripts/eval_simulator_intrinsic.py
# ?? benchmark/scripts/run_external_simulator_sweep.py
```

- [ ] **Step 2: Quickly review the diff to make sure nothing unintended is staged**

```powershell
git diff benchmark/scripts/eval_simulator_intrinsic.py
```

The diff should be only the additions: a new `_score_one_seed()` helper and a new `--all-seeds` branch in `main()`. No unrelated changes.

- [ ] **Step 3: Stage and commit**

```powershell
git add benchmark/scripts/eval_simulator_intrinsic.py benchmark/scripts/run_external_simulator_sweep.py
git commit -m "feat(benchmark): Task 5 extension -- 20-seed methylKit+DSS sweep + per-seed scoring"
```

- [ ] **Step 4: Confirm working tree is clean**

```powershell
git status
# expect: nothing to commit, working tree clean
```

---

## Task 2: Rename `benchmark/idk_if_needed/` → `benchmark/docs/historical/`

**Files:**
- Move: `benchmark/idk_if_needed/*.md` (8 files) → `benchmark/docs/historical/`

The directory is literally named "idk_if_needed" — leaving it that way at 1.0 is the visible-mess symptom the user flagged. The 8 markdown files inside are historical paper exploration notes worth keeping for provenance.

- [ ] **Step 1: Inspect the contents one more time to confirm nothing in there is in-flight work**

```powershell
ls benchmark/idk_if_needed/
# expect 8 .md files: dmr_replication_investigation, enrichment_vs_paper,
# panel_d_reproduction, paper_gene_check, paper_table_comparison,
# reactome_vs_paper, smoothed_dmr_vs_paper, top_k_report
```

- [ ] **Step 2: Verify the destination doesn't already exist**

```powershell
ls benchmark/docs/historical/ 2>$null
# expect: no such directory (it will be created by the move)
```

- [ ] **Step 3: Create destination and `git mv` the files**

```powershell
mkdir benchmark/docs/historical
git mv benchmark/idk_if_needed/dmr_replication_investigation.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/enrichment_vs_paper.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/panel_d_reproduction.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/paper_gene_check.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/paper_table_comparison.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/reactome_vs_paper.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/smoothed_dmr_vs_paper.md benchmark/docs/historical/
git mv benchmark/idk_if_needed/top_k_report.md benchmark/docs/historical/
```

- [ ] **Step 4: Remove the now-empty source directory**

```powershell
Remove-Item -Recurse benchmark/idk_if_needed
```

- [ ] **Step 5: Add a README.md to the new location explaining what these files are**

Create `benchmark/docs/historical/README.md`:

```markdown
# Historical paper exploration

These markdown files are paper-exploration scratch notes from the
pre-Phase-4 era (2026-05). They capture investigations that informed
the final paper (`benchmark/paper/paper.md`) but are not themselves
deliverables.

Files:
- `dmr_replication_investigation.md` — early DMR replication attempts vs. paper
- `enrichment_vs_paper.md` — enrichment-analysis comparison
- `panel_d_reproduction.md` — Figure 3D reproduction notes
- `paper_gene_check.md` — gene-list sanity checks
- `paper_table_comparison.md` — table-by-table comparison
- `reactome_vs_paper.md` — Reactome pathway comparison
- `smoothed_dmr_vs_paper.md` — smoothing-on-real-data investigation
- `top_k_report.md` — top-K gene-hit overlap analysis

Kept for provenance; not maintained.
```

- [ ] **Step 6: Commit**

```powershell
git add benchmark/docs/historical/
git commit -m "chore(benchmark): move idk_if_needed/ -> docs/historical/ with README"
```

- [ ] **Step 7: Verify no references in the codebase point at the old location**

```powershell
# Use Grep with pattern `idk_if_needed` across the whole repo
# (excluding .git and any committed CHANGELOG/plan files that reference it historically).
```

If any active script or paper text references `benchmark/idk_if_needed/`, update it to `benchmark/docs/historical/` and commit as a follow-up.

---

## Task 3: Run the 20-seed simulator sweep (~3h walltime, resumable)

**Files:**
- Reads: `benchmark/data/study1b_simulator/seed=*/` AMP files (existing)
- Writes: per-seed `methylkit.tsv`, `dss.tsv`, `dss_nosmooth.tsv` outputs

Resumable: re-running picks up where it left off. If a seed already has all three outputs, it's skipped.

- [ ] **Step 1: Pre-flight — confirm R is installed and `methylKit` + `DSS` are available**

```powershell
Rscript -e "library(methylKit); library(DSS); sessionInfo()"
```

Expected: R 4.5.0 (per `PROTOCOL.md` §1), methylKit 1.34.0+, DSS 2.56.0+ loaded without error.

If `Rscript` is not on PATH or the libraries are missing, this is the right time to stop and resolve. The sweep cannot run otherwise.

- [ ] **Step 2: Identify how many seeds need work vs. how many are already cached**

```powershell
$total = (Get-ChildItem benchmark/data/study1b_simulator -Directory -Filter "seed=*" | Measure-Object).Count
$done = (Get-ChildItem benchmark/data/study1b_simulator -Filter "methylkit.tsv" -Recurse | Measure-Object).Count
Write-Host "$done of $total seeds have methylkit.tsv"
```

This gives a rough estimate of remaining work. ~10 minutes per missing seed.

- [ ] **Step 3: Run the sweep (in the background; this is the long step)**

```powershell
# Foreground run is fine if you can spare ~3h of terminal. Otherwise
# start with run_in_background via your task runner.
uv run python benchmark/scripts/run_external_simulator_sweep.py
```

The script logs per-seed progress: which step ran, walltime, skip-if-exists hits. It writes `methylkit.tsv`, `dss.tsv`, `dss_nosmooth.tsv` per seed under `benchmark/data/study1b_simulator/seed=NNN/`.

- [ ] **Step 4: Verify all 20 seeds completed**

```powershell
Get-ChildItem benchmark/data/study1b_simulator -Directory -Filter "seed=*" | ForEach-Object {
    $seed = $_.Name
    $mk = Test-Path (Join-Path $_.FullName "methylkit.tsv")
    $dss = Test-Path (Join-Path $_.FullName "dss.tsv")
    $dssns = Test-Path (Join-Path $_.FullName "dss_nosmooth.tsv")
    Write-Host "$seed methylkit=$mk dss=$dss dss_nosmooth=$dssns"
}
```

Expected: all three columns `True` for all 20 seeds. If any are `False`, re-run Step 3 (resumable; only the missing seeds will execute).

- [ ] **Step 5: Skip the commit step**

These TSV outputs are data products, not source. They should be in `.gitignore` already (or covered by Task 12 of the 1.0-prep plan). Do NOT commit the raw TSVs.

---

## Task 4: Score all seeds intrinsically

**Files:**
- Writes: `benchmark/data/study1b_simulator/eval_simulator_intrinsic_per_seed.parquet`

- [ ] **Step 1: Run the new `--all-seeds` mode**

```powershell
uv run python benchmark/scripts/eval_simulator_intrinsic.py --all-seeds --coverage 10 -v
```

The script iterates every `seed=NNN/` with `methylkit.tsv + dss.tsv` present, scores each against `truth.parquet`, and writes:

- `eval_simulator_intrinsic_per_seed.parquet` (long-form, one row per (seed, tool, threshold, threshold_kind, meth_diff_bin))

It also prints an across-seed IQR summary at the headline cell (q < 0.05, all-bins) to the log.

- [ ] **Step 2: Sanity-check the output**

```powershell
uv run python -c "import polars as pl; df = pl.read_parquet('benchmark/data/study1b_simulator/eval_simulator_intrinsic_per_seed.parquet'); print(df.height, 'rows'); print(df.group_by('tool').len())"
```

Expected: `methylkit`, `dss`, and (if `dss_nosmooth.tsv` was present) `dss_nosmooth` rows; row count = N_seeds × N_thresholds × N_meth_diff_bins × N_tools.

- [ ] **Step 3: Skip the commit step**

Parquet output is a data product, not source.

---

## Task 5: Build the timings table

**Files:**
- Reads: `benchmark/data/study1/timings_post_phase3.parquet`, `benchmark/data/study1b_simulator/timings_simulator.parquet`
- Optionally reads: `benchmark/data/study1/timings.parquet` (pre-Phase-3) for delta
- Writes: `benchmark/data/study1/timings_table.csv` (paper-ready summary)

- [ ] **Step 1: Inspect what's in `timings_post_phase3.parquet`**

```powershell
uv run python -c "import polars as pl; df = pl.read_parquet('benchmark/data/study1/timings_post_phase3.parquet'); print(df.schema); print(df.head(20))"
```

Identify the columns. Typical: `tool`, `scenario`, `parameter_value` (or `coverage`), `wallclock_seconds`, `peak_rss_mb`, `machine_id` or similar.

- [ ] **Step 2: Decide what the paper's timings table will show**

Default recommendation (matches spec §3.3): one row per tool at the headline cell (coverage = 10, n = 5 vs 5, scenario = `dmc_coverage`), with columns:

| Tool | Wallclock (s) | Peak RSS (MB) | Notes |
|---|---|---|---|

If the pre/post-Phase-3 delta is interesting (Phase 3 changed engine math; some tools may have new runtime characteristics), add a second table comparing pre vs. post for epykit only. If pre/post deltas are noise (< 5%), drop the comparison and just report post-Phase-3 numbers in the main table.

- [ ] **Step 3: Generate the CSV**

Write a small one-off script — either inline at the powershell prompt or as `benchmark/scripts/build_timings_table.py` if you prefer it reusable. Inline version:

```powershell
uv run python -c @'
import polars as pl

post = pl.read_parquet("benchmark/data/study1/timings_post_phase3.parquet")
# Filter to the headline cell. Adjust column names to match the actual schema.
headline = post.filter(
    (pl.col("scenario") == "dmc_coverage") & (pl.col("parameter_value") == 10)
)
summary = headline.group_by("tool").agg(
    pl.col("wallclock_seconds").median().alias("median_wallclock_s"),
    pl.col("peak_rss_mb").median().alias("median_peak_rss_mb"),
).sort("median_wallclock_s")
summary.write_csv("benchmark/data/study1/timings_table.csv")
print(summary)
'@
```

Adjust the column names if the schema differs.

- [ ] **Step 4: Verify the CSV**

```powershell
cat benchmark/data/study1/timings_table.csv
```

Expected: one row per tool at the headline cell with wallclock + RSS columns. The slowest-tool row is what reviewers care about; make sure it's present.

- [ ] **Step 5: Decide whether to compute the pre/post delta**

```powershell
uv run python -c @'
import polars as pl
pre = pl.read_parquet("benchmark/data/study1/timings.parquet")
post = pl.read_parquet("benchmark/data/study1/timings_post_phase3.parquet")
# Filter both to epykit at the headline cell and compare.
def hl(df):
    return df.filter(
        (pl.col("tool").str.starts_with("epykit"))
        & (pl.col("scenario") == "dmc_coverage")
        & (pl.col("parameter_value") == 10)
    ).group_by("tool").agg(pl.col("wallclock_seconds").median())
print("Pre-Phase-3:"); print(hl(pre))
print("Post-Phase-3:"); print(hl(post))
'@
```

If the median delta is > 5% for any epykit variant, add a one-paragraph note in the paper Methods (Task 8) about what changed. If < 5%, skip the pre/post comparison and use only `timings_post_phase3.parquet`.

- [ ] **Step 6: Commit (csv + any new script)**

```powershell
git add benchmark/data/study1/timings_table.csv benchmark/scripts/build_timings_table.py 2>$null
git commit -m "feat(benchmark): timings table from timings_post_phase3.parquet for paper"
```

If you didn't create a standalone script (used the inline one-liner), only commit the CSV.

---

## Task 6: Populate `claims.yaml` with headline cells

**Files:**
- Modify: `benchmark/scripts/claims.yaml`

Per Phase 4 plan §9 and spec §3.1: pick the ~6-10 headline cells the paper cites, add entries so `regen_all.py --verify` can gate the paper against parquet drift.

- [ ] **Step 1: List the cells the paper will cite as headline numbers**

Look at `benchmark/paper/paper.md` and `benchmark/report/REPORT.md` for the numbers that currently appear in the abstract, headline result tables, and figure captions. Typical headline cells:

| Cell ID | Source parquet | Filter | Expected value |
|---|---|---|---|
| `study1_lr_auroc_cov10` | `eval_summary_post_phase3.parquet` | tool=`epykit_lr`, scenario=`dmc_coverage`, parameter_value=10 | AUROC |
| `study1_lrplus_auroc_cov10` | `eval_summary_post_phase3.parquet` | tool=`epykit_lrplus`, scenario=`dmc_coverage`, parameter_value=10 | AUROC |
| `study1_methylkit_auroc_cov10` | `eval_summary_post_phase3.parquet` | tool=`methylkit`, scenario=`dmc_coverage`, parameter_value=10 | AUROC |
| `study1_methylkit_tuned_auroc_cov10` | `eval_summary_post_phase3.parquet` | tool=`methylkit_tuned`, scenario=`dmc_coverage`, parameter_value=10 | AUROC |
| `study2_lr_f1_n2` | study2 eval parquet | tool=`epykit_lr`, n=2 cell | F1 |
| `simulator_lrplus_f1_cov10` | `eval_simulator_intrinsic_per_seed.parquet` | tool=`epykit_lrplus`, median across 20 seeds | F1 |
| `null_calib_lr_piao` | `benchmark/data/null_calibration/summary.parquet` | engine=`lr`, dataset=`piao` | observed FDR @ q=0.05 |
| `study3_dmr_iou_epykit_vs_dss` | `study3/comparisons/dmr_iou.parquet` | tools=`epykit`,`dss` | Intersection-over-Union |
| `headline_wallclock_epykit_lr` | `study1/timings_post_phase3.parquet` | tool=`epykit_lr`, cov=10 | median wallclock seconds |

Adjust this list to match the actual numbers the paper cites. The principle is: every distinct float that appears in the paper Abstract, Results headline tables, or figure captions needs a claim entry.

- [ ] **Step 2: Look up the actual value for each cell**

For each cell, query the parquet to get the real number. Example for `study1_lr_auroc_cov10`:

```powershell
uv run python -c @'
import polars as pl
df = pl.read_parquet("benchmark/data/study1/eval_summary_post_phase3.parquet")
hit = df.filter(
    (pl.col("tool") == "epykit_lr")
    & (pl.col("scenario") == "dmc_coverage")
    & (pl.col("parameter_value") == 10)
).select("auroc")
print(hit)
'@
```

Write each lookup down (in a temporary scratch file is fine).

- [ ] **Step 3: Populate `claims.yaml`**

Replace the empty `[]` in `benchmark/scripts/claims.yaml` with entries. Use the existing schema (shown in the file's header comment):

```yaml
# Phase 3 seed: empty. Phase 4 populates during the locked re-run.
# Schema per entry:
#   - claim_id: <id>
#     parquet:  <path>
#     column:   <col>
#     filter:   {col: val, ...}
#     expected: <float>
#     precision: <float>

- claim_id: study1_lr_auroc_cov10
  parquet: benchmark/data/study1/eval_summary_post_phase3.parquet
  column: auroc
  filter:
    tool: epykit_lr
    scenario: dmc_coverage
    parameter_value: 10
  expected: 0.987       # value from Step 2
  precision: 0.005

- claim_id: study1_lrplus_auroc_cov10
  parquet: benchmark/data/study1/eval_summary_post_phase3.parquet
  column: auroc
  filter:
    tool: epykit_lrplus
    scenario: dmc_coverage
    parameter_value: 10
  expected: 0.992
  precision: 0.005

# ... continue for all headline cells.
```

`precision` is the tolerance used by `regen_all.py --verify`. Set it to roughly the bootstrap CI half-width or 0.005, whichever is larger.

- [ ] **Step 4: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If any expected value mismatches the parquet by more than `precision`, the gate fails with a clear "claim X expected Y, got Z" message. Fix the `expected:` value and re-run.

- [ ] **Step 5: Commit**

```powershell
git add benchmark/scripts/claims.yaml
git commit -m "feat(benchmark): populate claims.yaml with paper headline cells"
```

---

## Task 7: Paper rewrite — Abstract

**Files:**
- Modify: `benchmark/paper/paper.md` (the abstract section at the top, before `# 1. Introduction` at line 42)

Per parent paper spec §6 Abstract requirements and Phase 4 plan §10 Step 1:

- Add the simulator-realism caveat: "Low-coverage TPR advantages observed on the underdispersed Piao 2021 simulator (φ ≈ 0.4) are not expected to transfer at the same magnitude to overdispersed real WGBS (φ ≈ 1.5–5)..."
- Replace "best-in-class" (if present) with "matches or exceeds the strongest baselines."
- Add the bug-fix manifest sentence: "We discovered and fixed N statistical bugs while running this benchmark; the post-fix numbers reported here..."

Plus per spec §1 decision (lr+ defaults stay OFF at 1.0): when the abstract cites headline numbers, report both bare `lr` and `lr+` transparently.

- [ ] **Step 1: Read the current abstract**

```powershell
sed -n '1,40p' benchmark/paper/paper.md
```

- [ ] **Step 2: Identify what to change**

Look for:
- Any phrasing claiming "best-in-class" / "outperforms all" / "state-of-the-art" — soften to "matches or exceeds the strongest baselines."
- Whether a single headline number is reported (e.g., "AUROC = 0.99 at coverage 10"). If so, replace with both `lr` and `lr+` numbers (e.g., "AUROC = 0.987 (bare `lr`) / 0.992 (`lr+` power stack) at coverage 10").
- Whether the bug-fix manifest is mentioned. If absent, add a sentence.
- Whether the simulator-realism caveat is present. If absent, add it.

- [ ] **Step 3: Apply the edits in `paper.md`**

Make focused replacements. Each edit should be a discrete diff hunk.

For each headline number cited, embed a claim marker right next to it in the markdown:

```markdown
The `lr+` power stack achieves AUROC = 0.992 <!-- claim: study1_lrplus_auroc_cov10 --> at coverage 10 vs. AUROC = 0.987 <!-- claim: study1_lr_auroc_cov10 --> for the bare `lr` engine.
```

`regen_all.py --verify` will find these markers and gate the cited numbers against `claims.yaml`.

- [ ] **Step 4: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If the gate fails, the cited number disagrees with `claims.yaml` (or the marker references a claim_id that doesn't exist in the yaml).

- [ ] **Step 5: Commit**

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): rewrite abstract per spec section 6"
```

---

## Task 8: Paper rewrite — Methods

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 2. Materials and Methods`, lines 94-228)

Per parent paper spec §6 Methods and Phase 4 plan §10 Step 2:

- Adopt PROTOCOL §4 parameter freeze verbatim
- Add a new subsection documenting the `simulate_piao.py` re-implementation (validation against Piao 2021 marginals — tested via `test_simulate_piao.py`)
- Add a new subsection on null calibration design
- Move the tile→chain_merge pivot narrative from elsewhere into Methods (per PROTOCOL R4)
- Document the 4 surviving engines (`lr`, `welch_t`, `fisher`, `glm`) + the 4 removed engines (`logit_t`, `bb_lr`, `score`, `cmh`) with migration hints
- Update §2.5 ("The `lr` and `lr+` engines") to match the 1.0 framing: bare `lr` is the default; `lr+` is the opt-in power stack engaged via `power_stack="lr+"`

Plus: cite the timings table generated in Task 5.

- [ ] **Step 1: Read the existing Methods section**

```powershell
sed -n '94,228p' benchmark/paper/paper.md
```

- [ ] **Step 2: Identify sections to add or modify**

Existing subsections (line numbers):
- 2.1 Simulated data (96)
- 2.2 Real data (115)
- 2.3 Ground truth (132)
- 2.4 Tools, versions, and parameters (150)
- 2.5 The `lr` and `lr+` engines (178)
- 2.6 Evaluation metrics (208)

Per spec, add:
- 2.4.X: Parameter freeze (verbatim from PROTOCOL §4 — copy the table into the paper)
- 2.5 update: lr is the default in 1.0; lr+ is opt-in via `power_stack`. Describe both.
- New 2.5.1 or sidebar: removed engines (`logit_t`, `bb_lr`, `score`, `cmh`) and migration hints (see CLAUDE.md or the 0.7.5 commit messages)
- New 2.7: Simulator re-implementation (`simulate_piao.py`) — describe validation against original marginals
- New 2.8: Null calibration design — describe `run_phase4_null_calibration.py` and the 12-cell sweep

- [ ] **Step 3: Apply the edits**

Add the new subsections. Update existing §2.5 to clarify lr (default) vs. lr+ (opt-in). Pull the PROTOCOL §4 table verbatim. Cite the timings table in §2.4 (e.g., "Wallclock and peak RSS per tool are reported in `benchmark/data/study1/timings_table.csv`; the headline values appear in §3.X").

Embed `<!-- claim: -->` markers next to any specific numbers cited (e.g., simulator dispersion φ ≈ 0.4 if that's a verifiable claim).

- [ ] **Step 4: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0.

- [ ] **Step 5: Commit**

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): rewrite methods per spec section 6 (parameter freeze, simulator, null calibration)"
```

---

## Task 9: Paper rewrite — Results

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 3. Results`, lines 230-641)

Per parent paper spec §6 Results and Phase 4 plan §10 Step 3:

- Every cell gets a Wilson or bootstrap CI in the same row
- Default-vs-default headline; tuned-vs-tuned panel clearly labelled
- AUROC reported intra-epykit only (per audit F15) — or extend cross-tool if you have p-value vectors from methylKit/DSS
- Add Table S-Calib (null calibration), Table S-Sim (multi-seed simulator from Task 4), Table S-Fix (bug-fix audit from `bug_fix_deltas.md`)
- Per spec §1 decision: lr and lr+ reported side-by-side

- [ ] **Step 1: Read the existing Results section**

```powershell
sed -n '230,420p' benchmark/paper/paper.md   # studies 1 and 2
sed -n '420,640p' benchmark/paper/paper.md   # study 3
```

- [ ] **Step 2: Add CI columns to result tables**

For each headline result table in §3.1, §3.2, §3.3: every metric cell needs a CI. Wilson CIs are appropriate for binary metrics (TPR, FPR); bootstrap CIs for AUROC and F1. Run `wilson_bootstrap_ci.py` (already on disk) if any CI data is missing.

If CI data already exists in `eval_summary_post_phase3.parquet` (per earlier work, columns like `auroc_ci_lo`, `auroc_ci_hi` should be present), reference them directly in the paper table.

- [ ] **Step 3: Add the lr / lr+ split to the headline tables**

For Study 1 headline (§3.1) and Study 2 headline (§3.2): show both `epykit_lr` (bare) and `epykit_lrplus` (power stack) as separate rows. Per spec §1 decision, this is transparent — the paper doesn't pick one as "the" epykit; it shows both and lets the reader decide.

- [ ] **Step 4: Add Table S-Sim (multi-seed simulator)**

From `eval_simulator_intrinsic_per_seed.parquet` (Task 4 output), build a table showing TPR/FPR/F1 median + IQR across 20 simulator seeds for `epykit_lrplus`, `methylkit_tuned`, and `dss` at the headline cell (cov=10, q<0.05, all-bins). Embed in a Results subsection or in the supplement (the spec's §6 leaves placement open).

- [ ] **Step 5: Add Table S-Calib (null calibration)**

From `benchmark/data/null_calibration/summary.parquet`, build a table with columns `(engine, dataset, scenario, observed_fdr_median, observed_fdr_q1, observed_fdr_q3)` — the 12 cells that ran. Add a note explaining the deferred `fisher@gse263850` cell (Task 11 covers documenting this in Limitations).

- [ ] **Step 6: Add Table S-Fix (bug-fix audit)**

Insert the existing `benchmark/data/audit/bug_fix_deltas.md` table verbatim into the supplement, with a one-paragraph framing: "We discovered and fixed N statistical bugs over the course of this benchmark; each row shows the per-cell delta between pre-fix and post-fix numbers, attributed to a specific commit."

- [ ] **Step 7: Embed claim markers next to every cited number**

Every floating-point number in the Results that appears as a paper-cited claim gets a `<!-- claim: cell_id -->` marker matching `claims.yaml`. If the marker references a claim_id not in the yaml, add the entry to `claims.yaml` first (Task 6 already covered the headline cells; new cells discovered while writing get added here).

- [ ] **Step 8: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If new claim markers were added without corresponding yaml entries, the gate fails — add them and re-verify.

- [ ] **Step 9: Commit**

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): rewrite results per spec section 6 (CIs, lr/lr+ split, supplementary tables)"
```

---

## Task 10: Paper rewrite — Discussion

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 4. Discussion`, lines 643-751)

Per parent paper spec §6 Discussion and Phase 4 plan §10 Step 4:

- Move the §4 underdispersion caveat to paragraph 2 of Discussion (currently at §4.2 or similar — promote it)
- Own the bug-fix manifest: "We found and fixed N bugs while running this benchmark..."

- [ ] **Step 1: Read the existing Discussion section**

```powershell
sed -n '643,751p' benchmark/paper/paper.md
```

Current structure:
- 4.1 What the three studies establish (645)
- 4.2 The calibration–sensitivity trade-off (667)
- 4.3 Two bugs we found and fixed (684) — may need to be renumbered/expanded
- 4.4 Limitations (705)

- [ ] **Step 2: Promote the underdispersion caveat**

Move the simulator-realism caveat ("the Piao 2021 simulator has φ ≈ 0.4, real WGBS has φ ≈ 1.5–5...") to paragraph 2 of §4.1 or as its own §4.1.1, so it's read before the calibration–sensitivity trade-off.

- [ ] **Step 3: Update the bug-fix manifest**

§4.3 says "Two bugs we found and fixed" — but per `bug_fix_deltas.md` and the audit, multiple bugs (P0-1 through P1-11) shipped. Update the section title to "The bug-fix manifest" or similar, and rewrite to own the full scope: "Over the course of the Phase 1-3 work, we discovered and fixed N statistical bugs in epykit. The post-fix numbers reported in §3 reflect..." Add a short paragraph for each significant fix or reference Table S-Fix.

- [ ] **Step 4: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0.

- [ ] **Step 5: Commit**

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): rewrite discussion per spec section 6 (underdispersion caveat, bug-fix manifest)"
```

---

## Task 11: Paper rewrite — Limitations §10.5 (fisher cell + bug-fix table)

**Files:**
- Modify: `benchmark/paper/paper.md` (the §4.4 Limitations section at line 705, specifically §10.5)

Per spec §3.4: the `fisher@gse263850` null-calibration cell did not run (deferred pending parallel backend). Document it as a limitation rather than re-run.

Also per spec §3.2 / Phase 4 plan §10 Step 5: insert the `bug_fix_deltas.md` table verbatim under Limitations §10.5.

- [ ] **Step 1: Read the existing Limitations section**

```powershell
sed -n '705,755p' benchmark/paper/paper.md
```

- [ ] **Step 2: Add the fisher cell limitation**

Add a new bullet (or subsection) in Limitations §10.5:

```markdown
**Null calibration coverage.** The null-calibration sweep (§3.X, Table S-Calib)
covers 12 of 13 engine × dataset cells. The `fisher@gse263850` cell — the
small-n Fisher engine on the 12-sample real GSE263850 dataset — was deferred
because the closure requires a parallel backend to be tractable. The
fisher engine is documented as a small-n fallback (`n < 2` per group);
its null behavior on a 12-sample real cohort is the least informative
cell in the sweep. The remaining 12 cells confirm engine calibration at
nominal q ∈ {0.01, 0.05, 0.10}.
```

- [ ] **Step 3: Inline the bug-fix delta table**

The existing markdown at `benchmark/data/audit/bug_fix_deltas.md` contains a per-cell delta table. Copy it verbatim into the paper Limitations §10.5 (or reference it as Supplementary Table S-Fix if it's already in the supplement per Task 9).

```powershell
cat benchmark/data/audit/bug_fix_deltas.md
```

- [ ] **Step 4: Run the verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0.

- [ ] **Step 5: Commit**

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): limitations 10.5 -- fisher@gse263850 deferred cell + bug-fix table"
```

---

## Task 12: Full claims gate + paper proofread pass

**Files:**
- Read: `benchmark/paper/paper.md`, `benchmark/scripts/claims.yaml`
- Modify (as needed): `benchmark/paper/paper.md`, `benchmark/scripts/claims.yaml`

- [ ] **Step 1: Re-run the full verification gate**

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If non-zero, fix the cited number / claim entry / precision and re-run.

- [ ] **Step 2: Verify every claim marker references a yaml entry**

```powershell
# Extract claim IDs from paper.md and from claims.yaml; diff.
$paper_ids = (Select-String -Path benchmark/paper/paper.md -Pattern "<!-- claim: ([\w_]+) -->" -AllMatches | ForEach-Object { $_.Matches } | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
$yaml_ids = (Select-String -Path benchmark/scripts/claims.yaml -Pattern "^- claim_id: (\w+)" | ForEach-Object { $_.Matches.Groups[1].Value } | Sort-Object -Unique)
$missing_in_yaml = Compare-Object $paper_ids $yaml_ids -PassThru | Where-Object { $_ -in $paper_ids }
$missing_in_paper = Compare-Object $paper_ids $yaml_ids -PassThru | Where-Object { $_ -in $yaml_ids }
Write-Host "Claim IDs in paper but not yaml: $missing_in_yaml"
Write-Host "Claim IDs in yaml but not paper: $missing_in_paper"
```

Expected: both lists empty. Add missing yaml entries OR remove orphan yaml entries until clean.

- [ ] **Step 3: Proofread for the spec §1 decision**

Search the paper for any remaining "recommended" / "default" language that describes `lr+` as the default or recommended setup. Per spec §1, `lr+` is opt-in; bare `lr` is the default. Adjust any leftover phrasing.

```powershell
Select-String -Path benchmark/paper/paper.md -Pattern "lr\+.*default|default.*lr\+|recommended.*lr\+|lr\+.*recommended" -Context 1
```

If any matches return, edit them to match the new framing (`lr+` is an opt-in power stack).

- [ ] **Step 4: Verify all section anchors resolve**

If the paper uses cross-references like `(§3.X)`, `(Table S-Calib)`, `(Figure 2)` — confirm the referenced section/table/figure actually exists. The paper has been edited section-by-section; a stale reference from an earlier draft is a likely artifact.

- [ ] **Step 5: If any edits were needed, commit**

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): final proofread pass -- claim coverage, lr+ framing, cross-references"
```

If no edits were needed, skip the commit.

---

## Task 13: Open PR and merge to main

**Files:** none — git/GH workflow.

- [ ] **Step 1: Push the branch**

```powershell
git push -u origin paper
```

- [ ] **Step 2: Open the PR**

```powershell
gh pr create --title "epykit 1.0 -- paper finish (claims.yaml, paper.md rewrite, 20-seed sweep)" --body "$(cat <<'EOF'
## Summary
- Implements §3 (paper track) + paper-side §5 hygiene of the [1.0 + paper-finish design spec](docs/superpowers/specs/2026-06-01-library-1.0-and-paper-finish-design.md)
- 20-seed simulator parallel column (methylKit + DSS); scored intrinsically against truth.parquet
- Timings table built from timings_post_phase3.parquet
- claims.yaml populated with paper headline cells; regen_all.py --verify gating in place
- paper.md rewritten per parent spec section 6: abstract, methods, results, discussion, limitations
- fisher@gse263850 null-calibration cell documented as limitation 10.5
- bug_fix_deltas table inlined into paper limitations
- benchmark/idk_if_needed/ renamed to benchmark/docs/historical/ with README

## Test plan
- [ ] uv run pytest benchmark/scripts/tests/ -q passes
- [ ] uv run python benchmark/scripts/regen_all.py --verify exits 0
- [ ] All 20 simulator seeds have methylkit.tsv + dss.tsv + dss_nosmooth.tsv
- [ ] eval_simulator_intrinsic_per_seed.parquet exists with expected schema
- [ ] All claim markers in paper.md resolve to claims.yaml entries
- [ ] Paper renders cleanly (visual review of rendered markdown)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 3: Wait for CI**

```powershell
gh pr checks --watch
```

If any check fails, investigate and fix on this branch.

- [ ] **Step 4: Merge to main**

```powershell
gh pr merge --merge   # preserve per-task commits for bisectability
```

- [ ] **Step 5: Switch back to main and pull**

```powershell
git switch main
git pull --ff-only origin main
```

---

## Definition of Done

The paper-track branch is complete when ALL of the following are true:

- [ ] `git log --oneline main..origin/main` shows the paper-track merge commit on main
- [ ] `benchmark/scripts/claims.yaml` has ≥ 6 entries covering every floating-point number cited in the paper Abstract and Results headline tables
- [ ] `uv run python benchmark/scripts/regen_all.py --verify` exits 0
- [ ] `benchmark/paper/paper.md` has been rewritten per parent spec §6 — abstract, methods, results, discussion, limitations all match the spec's framing
- [ ] Paper Methods section reports bare `lr` and `lr+` numbers side-by-side per spec §1 decision
- [ ] `fisher@gse263850` null-calibration cell is documented as a limitation in §10.5
- [ ] `bug_fix_deltas.md` table is inlined into Limitations (or referenced as a supplementary table)
- [ ] `benchmark/data/study1b_simulator/eval_simulator_intrinsic_per_seed.parquet` exists (20 seeds × 3 tools × thresholds)
- [ ] `benchmark/data/study1/timings_table.csv` exists with one row per tool at the headline cell
- [ ] `benchmark/idk_if_needed/` no longer exists; contents are at `benchmark/docs/historical/` with a README
- [ ] Working tree is clean on main after merge

---

## Notes for the implementer

1. **Task 3 is the long pole** (~3 hours of methylKit + DSS sweep). Start it in the background; Tasks 1, 2, 5, 6 can run in parallel while it executes.

2. **`regen_all.py --verify` is the safety net.** Run it after every paper edit. If it fails, the cited number drifted from the parquet — don't push past it.

3. **Resist editing the paper out of order.** Each section commit (Tasks 7-11) lands one section. Reviewing the diff section-by-section is what makes a section-by-section commit worthwhile.

4. **The Phase 4 plan has additional P2 hygiene items** ([Task 11 in `2026-05-28-phase-4-locked-rerun.md`](2026-05-28-phase-4-locked-rerun.md)) — `_BETA_EPSILON` split, drop `pct_sig`, CGI shore/shelf width kwargs, raise empirical FDR `n_perm` default, rank `-log10(pvalue)` in `_auroc`, manuscript path inconsistency. They are NOT in scope for this plan. If the user wants them, they can be picked up as a separate sweep after the paper merges.

5. **The companion `1.0-prep` branch is independent.** It can merge before or after this branch. The CHANGELOG entry for 1.0 lives on `1.0-prep`; this branch does not need to add to it.

6. **If Task 3 (the 20-seed sweep) cannot run** (Rscript not installed, methylKit unavailable, walltime exhausted) — fall back to the 1-seed parallel column already on disk and document the multi-seed plan as future work in Limitations. The paper is still defensible with 1 seed; 20 seeds is a strengthening, not a requirement.

7. **Paper word count and figure count are not bounded by this plan.** The spec leaves those open. Use judgment — match the length and figure count of comparable methods-paper preprints in the field.
