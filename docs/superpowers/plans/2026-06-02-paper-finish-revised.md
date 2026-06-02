# epykit Paper Finish — Revised Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Status:** This plan supersedes [`2026-06-01-paper-finish.md`](2026-06-01-paper-finish.md). The original assumed the 20-seed simulator sweep and per-seed scoring still needed to run; they already shipped. This revision drops the done items and reframes the remaining work as a paper-rewrite-focused sprint.

**Goal:** Take the existing benchmark deliverables (eval parquets, null calibration summary, timings, parallel-column prose) and produce a paper-ready `paper.md` + populated `claims.yaml` gated by `regen_all.py --verify`. Then merge to `main`.

**Architecture:** Single feature branch `paper` off the new `main` (post-1.0 merge). Paper sections land as one commit per section. Verification gate (`regen_all.py --verify`) runs after every `claims.yaml` change and after the final section commit.

**Tech Stack:** Polars for parquet inspection, R subprocesses already ran (nothing new needs to compute), markdown + Jinja-style claim markers for the paper, `regen_all.py` for the gate.

---

## What's already done (confirmed on main as of 2026-06-02)

- **20-seed simulator sweep** — all `benchmark/data/study1b_simulator/seed=2026000` through `seed=2026019` have `methylkit.tsv` / `dss.tsv` / `dss_nosmooth.tsv` outputs (Task 3 of the old plan)
- **Per-seed intrinsic scoring** — `eval_simulator_intrinsic_per_seed.parquet` (540 rows: 20 seeds × tools × thresholds) and `_iqr.parquet` summary exist (Task 4 of the old plan)
- **WIP scripts committed** — `run_external_simulator_sweep.py` is on main; `eval_simulator_intrinsic.py --all-seeds` mode is on main (Task 1 of the old plan)
- **Parallel-column prose** — `benchmark/data/study1b_simulator/parallel_column_summary.md` has well-drafted narrative for the 7-tool comparison; can be pulled into the paper Results section
- **Phase 4 deliverables** — `eval_summary_post_phase3.parquet` (15 tools, 3 scenarios, with CI columns), `null_calibration/summary.parquet` (12 cells), `data/audit/bug_fix_deltas.md/.parquet`, `study3/comparisons_post_phase3/` all on main
- **Timings parquets** — `timings_post_phase3.parquet`, `timings.parquet`, `timings_simulator.parquet`, `eval_external_timings_per_seed.parquet` all on main
- **Stashed draft** — `stash@{0}` carries a 180-line `benchmark/docs/timing-comparison.md` rewrite from an earlier session; useful input for Task 4 (paper Methods)

## What's still pending

- **`benchmark/idk_if_needed/`** — still exists with 8 historical exploration markdowns; needs rename to `benchmark/docs/historical/`
- **`benchmark/scripts/claims.yaml`** — still 9 lines (empty seed); needs population
- **`benchmark/paper/paper.md`** — 789 lines, last touched in commit `1e94b9d` (pre-Phase-4 consolidation). Needs the full §6 rewrite per parent paper spec
- **Stash decision** — pop `stash@{0}` and integrate, or discard if obsolete
- **Optional: `benchmark/data/study1/timings_table.csv`** — paper-ready summary CSV from the three timings parquets

---

## Companion docs

- **Spec:** [`../specs/2026-06-01-library-1.0-and-paper-finish-design.md`](../specs/2026-06-01-library-1.0-and-paper-finish-design.md) — read §3 (paper track) and §1 (decisions) before starting
- **Parent paper spec:** [`../specs/2026-05-27-paper-defendable-benchmark-design.md`](../specs/2026-05-27-paper-defendable-benchmark-design.md) — §6 has per-section paper rewrite requirements
- **Original (superseded) paper plan:** [`2026-06-01-paper-finish.md`](2026-06-01-paper-finish.md) — kept for historical reference

---

## Pre-flight verification

```powershell
# 1. Confirm on main with clean tree.
git status                                  # expect: clean (stash separately)
git branch --show-current                   # expect: main
git log --oneline -1                        # expect: tip is "Merge pull request #1" (d7c588f or later)

# 2. Confirm the existing benchmark tests pass.
uv run pytest benchmark/scripts/tests/ -q   # expect: green

# 3. Confirm regen_all.py runs (currently exits 0 because claims.yaml is empty).
uv run python benchmark/scripts/regen_all.py --verify

# 4. Inspect the stash; decide whether to pop, drop, or hold.
git stash show -p stash@{0}                 # read the 180-line diff
# Decision lives in Task 1 of this plan.

# 5. Create the working branch.
git switch -c paper
```

---

## File map

| File | Role |
|---|---|
| `benchmark/idk_if_needed/*.md` | Renamed to `benchmark/docs/historical/` in Task 1 |
| `benchmark/docs/timing-comparison.md` | Stashed rewrite — integrate or discard in Task 1 |
| `benchmark/data/study1/timings_table.csv` | Optional: generated in Task 2 from existing parquets |
| `benchmark/scripts/claims.yaml` | Populated in Task 3 with headline cells |
| `benchmark/paper/paper.md` | Section-by-section rewrite in Tasks 4-8 |
| `benchmark/data/audit/bug_fix_deltas.md` | Inlined verbatim into Limitations §10.5 in Task 8 |

---

## Task 1: Hygiene — `idk_if_needed/` rename + stash decision

**Files:**
- Move: `benchmark/idk_if_needed/*.md` (8 files) → `benchmark/docs/historical/`
- Create: `benchmark/docs/historical/README.md`
- Update: `.gitignore` (line 19 references `benchmark/idk_if_needed/`)
- Optionally apply or discard: `stash@{0}`

### Step 1: Move the 8 markdown files

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
Remove-Item -Recurse benchmark/idk_if_needed
```

### Step 2: Add a README at the new location

Create `benchmark/docs/historical/README.md`:

```markdown
# Historical paper exploration

These markdown files are paper-exploration scratch notes from the
pre-Phase-4 era (2026-05). They captured investigations that informed
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

### Step 3: Update `.gitignore`

The current `.gitignore` line 19 reads `benchmark/idk_if_needed/`. The directory no longer exists, but the line should be removed (or commented out) for cleanliness. Edit `.gitignore`: delete line 19.

### Step 4: Verify no active references

```powershell
# Use Grep with pattern `idk_if_needed` across the whole repo
# (excluding committed CHANGELOG/plan files that may reference it historically).
```

If any active script references `benchmark/idk_if_needed/`, update it to `benchmark/docs/historical/` and stage.

### Step 5: Decide on the stash

```powershell
git stash show -p stash@{0}
```

The stash contains a 180-line rewrite of `benchmark/docs/timing-comparison.md`. Inspect it:
- If it's a useful narrative (likely — it's a comprehensive 7-tool timing + accuracy comparison), **apply and commit** as part of this task: `git stash pop` then `git add benchmark/docs/timing-comparison.md`.
- If it's stale or duplicates `parallel_column_summary.md`, **drop**: `git stash drop stash@{0}`.

Recommended: pop and commit. The timing-comparison.md will likely be cited in the paper Methods §2.4 (performance reporting).

### Step 6: Commit

```powershell
git add benchmark/docs/historical/ benchmark/idk_if_needed/ .gitignore benchmark/docs/timing-comparison.md
git commit -m "chore(benchmark): rename idk_if_needed/ -> docs/historical/ + integrate timing-comparison rewrite"
```

---

## Task 2: Generate `timings_table.csv` (optional but recommended)

**Files:**
- Read: `benchmark/data/study1/timings_post_phase3.parquet`, `benchmark/data/study1b_simulator/timings_simulator.parquet`
- Write: `benchmark/data/study1/timings_table.csv`

A paper-ready summary CSV with one row per tool at the headline cell (coverage=10, 3v3 design). The paper Methods §2.4 / Results performance table will cite numbers from this CSV.

### Step 1: Inspect schemas

```powershell
uv run python -c "import polars as pl; df = pl.read_parquet('benchmark/data/study1/timings_post_phase3.parquet'); print('schema:', df.schema); print(df.head(20))"
```

Note the actual column names (tool, wallclock_seconds, peak_rss_mb, parameter_value, scenario, etc.).

### Step 2: Build the headline summary

Write an inline script (no need for a standalone .py file):

```powershell
uv run python -c @'
import polars as pl

post = pl.read_parquet("benchmark/data/study1/timings_post_phase3.parquet")
# Adjust filter columns to match actual schema from Step 1.
headline = post.filter(
    (pl.col("scenario") == "dmc_coverage") & (pl.col("parameter_value") == 10)
)
# Pick the right column names from Step 1 - "wallclock_seconds" is the typical convention.
summary = headline.group_by("tool").agg(
    pl.col("wallclock_seconds").median().alias("median_wallclock_s"),
    pl.col("peak_rss_mb").median().alias("median_peak_rss_mb"),
).sort("median_wallclock_s")
summary.write_csv("benchmark/data/study1/timings_table.csv")
print(summary)
'@
```

If schema columns differ, adapt the filter and aggregation accordingly.

### Step 3: Optional pre/post delta

If `timings.parquet` (pre-Phase-3) shows materially different numbers from `timings_post_phase3.parquet` at the headline cell, compute the delta and decide whether to mention in the paper Methods. If deltas are < 5%, drop the comparison and just report post-Phase-3 numbers.

### Step 4: Commit

```powershell
git add benchmark/data/study1/timings_table.csv
git commit -m "feat(benchmark): timings table CSV for paper headline performance row"
```

---

## Task 3: Populate `claims.yaml` with headline cells

**Files:**
- Modify: `benchmark/scripts/claims.yaml`

Per Phase 4 plan §9 and parent spec §6: pick the ~8-12 headline cells the paper will cite, look up actual values from parquets, populate the yaml. `regen_all.py --verify` then gates against drift.

### Step 1: List the headline cells

Cells the paper will likely cite (adjust based on what the paper rewrite ends up emphasizing):

| Cell ID | Source parquet | Filter | Column |
|---|---|---|---|
| `study1_lr_auroc_cov10` | `study1/eval_summary_post_phase3.parquet` | `tool=epykit_lr, scenario=dmc_coverage, parameter_value=10, meth_diff_bin=all` | `auroc` |
| `study1_lrplus_auroc_cov10` | same | `tool=epykit_lrplus, ...` | `auroc` |
| `study1_methylkit_auroc_cov10` | same | `tool=methylkit, ...` | `auroc` |
| `study1_methylkit_tuned_auroc_cov10` | same | `tool=methylkit_tuned, ...` | `auroc` |
| `study1_lr_f1_cov10` | same | `tool=epykit_lr, ...` | `f1` |
| `simulator_lrplus_auroc_seed0` | `study1b_simulator/eval_simulator_intrinsic_per_seed.parquet` | `tool=epykit_lrplus, seed=2026000` | `auroc` |
| `simulator_methylkit_auroc_seed0` | same | `tool=methylkit, seed=2026000` | `auroc` |
| `null_calib_lr_piao_q05` | `null_calibration/summary.parquet` | `engine=lr, dataset=piao, scenario=dmc_coverage` | `observed_fdr_median` |
| `null_calib_lrplus_piao_q05` | same | `engine=lrplus, dataset=piao, ...` | `observed_fdr_median` |
| `headline_wallclock_epykit_lr` | `study1/timings_table.csv` | `tool=epykit_lr` | `median_wallclock_s` |
| `headline_wallclock_methylkit` | same | `tool=methylkit` | `median_wallclock_s` |

Adjust if the paper cites different numbers. Add cells as needed during the paper rewrite (Tasks 4-8) — `claims.yaml` is a living file across this branch.

### Step 2: Look up actual values

For each cell, query the parquet (or CSV for timings) and record the value:

```powershell
uv run python -c @'
import polars as pl
df = pl.read_parquet("benchmark/data/study1/eval_summary_post_phase3.parquet")
hit = df.filter(
    (pl.col("tool") == "epykit_lr")
    & (pl.col("scenario") == "dmc_coverage")
    & (pl.col("parameter_value") == 10)
    & (pl.col("meth_diff_bin") == "all")
).select("auroc", "tpr", "fpr", "f1")
print(hit)
'@
```

Repeat for every cell. The parallel-column summary already has many of the seed=0 simulator numbers in prose form (see `study1b_simulator/parallel_column_summary.md`) — cross-check.

### Step 3: Write entries to `benchmark/scripts/claims.yaml`

Replace the empty `[]` with entries:

```yaml
- claim_id: study1_lr_auroc_cov10
  parquet: benchmark/data/study1/eval_summary_post_phase3.parquet
  column: auroc
  filter:
    tool: epykit_lr
    scenario: dmc_coverage
    parameter_value: 10
    meth_diff_bin: all
  expected: 0.9999       # value from Step 2 (per parallel_column_summary.md the Piao-as-distributed AUROC is ~0.9999)
  precision: 0.0005

- claim_id: simulator_methylkit_auroc_seed0
  parquet: benchmark/data/study1b_simulator/eval_simulator_intrinsic_per_seed.parquet
  column: auroc
  filter:
    tool: methylkit
    seed: 2026000
    threshold_kind: qvalue
    threshold: 0.05
    meth_diff_bin: all
  expected: 0.9246
  precision: 0.001
```

`precision` is the tolerance used by `regen_all.py --verify`. Default to 0.001 for AUROC, 0.005 for TPR/F1, 0.5 for wallclock seconds.

### Step 4: Run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If any expected value mismatches the parquet by more than `precision`, the gate fails with a "claim X expected Y, got Z" message. Fix the `expected:` value and re-run.

### Step 5: Commit

```powershell
git add benchmark/scripts/claims.yaml
git commit -m "feat(benchmark): populate claims.yaml with paper headline cells"
```

---

## Task 4: Paper rewrite — Abstract

**Files:**
- Modify: `benchmark/paper/paper.md` (lines 1-41, before `# 1. Introduction`)

Per parent paper spec §6 Abstract requirements:

- Add the simulator-realism caveat: "Low-coverage TPR advantages observed on the underdispersed Piao 2021 simulator (φ ≈ 0.4) are not expected to transfer at the same magnitude to overdispersed real WGBS (φ ≈ 1.5–5)..."
- Replace "best-in-class" with "matches or exceeds the strongest baselines"
- Add the bug-fix manifest sentence: "We discovered and fixed N statistical bugs while running this benchmark; the post-fix numbers reported here..." (N comes from `data/audit/bug_fix_deltas.md` — count the rows)
- Per spec §1 decision (lr+ defaults stay OFF): when the abstract cites headline numbers, report both bare `lr` and `lr+` transparently

### Step 1: Read the current abstract

```powershell
Get-Content benchmark/paper/paper.md -TotalCount 41
```

### Step 2: Apply targeted edits

For each cited number, embed a claim marker:

```markdown
The `lr+` power stack achieves AUROC = 0.992 <!-- claim: study1_lrplus_auroc_cov10 --> at coverage 10 vs. AUROC = 0.987 <!-- claim: study1_lr_auroc_cov10 --> for the bare `lr` engine on Piao-as-distributed.
```

### Step 3: Run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0. If a claim marker references an id not in `claims.yaml`, the gate fails — add the entry to `claims.yaml` first.

### Step 4: Commit

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): rewrite abstract per spec section 6 (simulator caveat + bug-fix manifest)"
```

---

## Task 5: Paper rewrite — Methods

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 2. Materials and Methods`, lines 94-228)

Per parent paper spec §6 Methods:

- Adopt PROTOCOL §4 parameter freeze verbatim
- Add new subsection §2.X documenting the `simulate_piao.py` re-implementation (validate against Piao 2021 marginals)
- Add new §2.Y on null calibration design (the 12-cell `run_phase4_null_calibration.py` sweep)
- Move the tile→chain_merge pivot narrative into Methods (per PROTOCOL R4)
- Update §2.5 ("The `lr` and `lr+` engines"): bare `lr` is the default; `lr+` opt-in via `power_stack="lr+"`. Describe both.
- Document the 4 surviving engines + 4 removed engines with migration hints
- Cite the timings table from Task 2 in §2.4

The existing §2 subsection structure (lines 96-228 per the original paper outline):
- 2.1 Simulated data (96)
- 2.2 Real data (115)
- 2.3 Ground truth (132)
- 2.4 Tools, versions, and parameters (150)
- 2.5 The `lr` and `lr+` engines (178)
- 2.6 Evaluation metrics (208)

### Step 1: Read the current Methods section

```powershell
Get-Content benchmark/paper/paper.md | Select-Object -Skip 93 -First 135
```

### Step 2: Apply edits

For each cited number, embed a claim marker. Pull data citations from `PROTOCOL.md`, `data/audit/commits.json`, and the existing parquets.

If `benchmark/docs/timing-comparison.md` (from the stash) has good prose for §2.4 performance reporting, paraphrase / link into Methods rather than duplicating.

### Step 3: Run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

### Step 4: Commit

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): rewrite methods (parameter freeze, simulator, null calibration)"
```

---

## Task 6: Paper rewrite — Results

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 3. Results`, lines 230-641)

Per parent paper spec §6 Results:

- Every cited cell gets a Wilson or bootstrap CI in the same row (data already in `eval_summary_post_phase3.parquet` via `*_ci_lo` / `*_ci_hi` columns from Phase 4 work)
- Default-vs-default headline; tuned-vs-tuned panel clearly labelled
- AUROC reported intra-epykit only OR extended cross-tool (use the parallel-column data)
- Add Table S-Calib (null calibration from `null_calibration/summary.parquet`)
- Add Table S-Sim (multi-seed simulator IQR from `eval_simulator_intrinsic_iqr.parquet`)
- Add Table S-Fix (bug-fix audit from `data/audit/bug_fix_deltas.md`)
- Per spec §1 decision: lr and lr+ reported side-by-side

The 7-tool parallel-column comparison in `study1b_simulator/parallel_column_summary.md` is well-drafted prose — pull it into §3 as a Results subsection.

Existing §3 structure (line 230+):
- 3.1 Study 1 — Panel comparison on simulated data
- 3.2 Study 2 — Head-to-head with methylKit
- 3.3 Study 3 — Real WGBS data (GSE263850), 3-way DMR-caller comparison

### Step 1: Read the current Results section

```powershell
Get-Content benchmark/paper/paper.md | Select-Object -Skip 229 -First 410
```

### Step 2: Apply edits

For each table:
- Add CI columns reading from `eval_summary_post_phase3.parquet`'s `*_ci_lo`/`*_ci_hi` columns
- Embed `<!-- claim: -->` markers next to each cited number

For the parallel-column subsection, pull prose verbatim from `parallel_column_summary.md` and adapt to paper voice.

For supplementary Tables S-Calib, S-Sim, S-Fix: either inline at the end of §3 or in an appendix — pick consistent placement.

### Step 3: Run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

### Step 4: Commit

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): rewrite results (CIs, lr/lr+ split, supplementary tables S-Calib/S-Sim/S-Fix)"
```

---

## Task 7: Paper rewrite — Discussion

**Files:**
- Modify: `benchmark/paper/paper.md` (section `# 4. Discussion`, lines 643-751)

Per parent paper spec §6 Discussion:

- Move the underdispersion caveat to paragraph 2 of Discussion (currently in §4.2; promote it)
- Own the bug-fix manifest: "We found and fixed N bugs while running this benchmark..." with N from `bug_fix_deltas.md`
- Update §4.3 from "Two bugs we found and fixed" to reflect actual count

Existing §4 structure:
- 4.1 What the three studies establish (645)
- 4.2 Calibration–sensitivity trade-off (667)
- 4.3 Two bugs we found and fixed (684) — needs rename + expansion
- 4.4 Limitations (705)

### Step 1-3: Read, edit, gate (same pattern as Tasks 4-6)

### Step 4: Commit

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): rewrite discussion (underdispersion caveat promotion + bug-fix manifest)"
```

---

## Task 8: Paper rewrite — Limitations §10.5 + bug-fix table

**Files:**
- Modify: `benchmark/paper/paper.md` (section §4.4 Limitations, line 705+)

Per spec §3.4 and §3.2:

- Add the fisher@gse263850 deferred cell as a limitation bullet
- Inline `benchmark/data/audit/bug_fix_deltas.md` table verbatim under Limitations §10.5

### Step 1: Add the fisher cell limitation

Append a new bullet/subsection in §4.4:

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

### Step 2: Inline the bug-fix delta table

Read `benchmark/data/audit/bug_fix_deltas.md` and append it verbatim (or as Supplementary Table S-Fix) under Limitations §10.5.

### Step 3: Run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

### Step 4: Commit

```powershell
git add benchmark/paper/paper.md
git commit -m "docs(paper): limitations 10.5 -- fisher@gse263850 deferred + bug-fix audit table"
```

---

## Task 9: Final claims gate + proofread

**Files:**
- Read: `benchmark/paper/paper.md`, `benchmark/scripts/claims.yaml`
- Modify (as needed): both

### Step 1: Re-run the gate

```powershell
uv run python benchmark/scripts/regen_all.py --verify
```

Expected: exit 0.

### Step 2: Verify every claim marker resolves

```powershell
$paper_ids = (Select-String -Path benchmark/paper/paper.md -Pattern "<!-- claim: ([\w_]+) -->" -AllMatches | ForEach-Object { $_.Matches } | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
$yaml_ids = (Select-String -Path benchmark/scripts/claims.yaml -Pattern "^- claim_id: (\w+)" | ForEach-Object { $_.Matches.Groups[1].Value } | Sort-Object -Unique)
Compare-Object $paper_ids $yaml_ids
```

Expected: empty Compare-Object output (both lists identical).

### Step 3: Spec §1 framing sweep

```powershell
Select-String -Path benchmark/paper/paper.md -Pattern "lr\+.*default|default.*lr\+|recommended.*lr\+|lr\+.*recommended" -Context 1
```

If any matches return, edit to the new framing (`lr+` is an opt-in power stack, not a default).

### Step 4: Cross-reference sweep

If the paper uses `§3.X`, `Table S-Calib`, `Figure 2`, confirm each reference resolves. Stale references from earlier drafts get cleaned up.

### Step 5: Commit if any edits

```powershell
git add benchmark/paper/paper.md benchmark/scripts/claims.yaml
git commit -m "docs(paper): final proofread (claim coverage, lr+ framing, cross-references)"
```

If no edits were needed, skip.

---

## Task 10: PR + merge to main

**Files:** none — git/GH workflow.

### Step 1: Push

```powershell
git push -u origin paper
```

### Step 2: Open the PR

```powershell
gh pr create --title "epykit benchmark paper finish (claims.yaml + paper.md rewrite)" --body "$(cat <<'EOF'
## Summary
- Implements paper-track work from [1.0 + paper-finish design spec](docs/superpowers/specs/2026-06-01-library-1.0-and-paper-finish-design.md) section 3
- paper.md rewritten per parent spec section 6: abstract, methods, results, discussion, limitations
- claims.yaml populated with paper headline cells; regen_all.py --verify gating in place
- benchmark/idk_if_needed/ renamed to benchmark/docs/historical/ with README
- timings_table.csv generated for paper Methods 2.4
- fisher@gse263850 null-calibration cell documented as limitation
- bug_fix_deltas table inlined into limitations

## Test plan
- [ ] uv run pytest benchmark/scripts/tests/ -q passes
- [ ] uv run python benchmark/scripts/regen_all.py --verify exits 0
- [ ] All claim markers in paper.md resolve to claims.yaml entries
- [ ] Paper renders cleanly (visual review)

EOF
)"
```

### Step 3: Wait for CI

```powershell
gh pr checks --watch
```

### Step 4: Merge

```powershell
gh pr merge --merge   # preserve per-task commits
```

### Step 5: Cleanup

```powershell
git switch main
git pull --ff-only origin main
```

---

## Definition of Done

- `benchmark/scripts/claims.yaml` has ≥ 8 entries, all gating real paper citations
- `uv run python benchmark/scripts/regen_all.py --verify` exits 0
- `benchmark/paper/paper.md` is rewritten per parent spec §6
- Paper Methods reports both `lr` and `lr+` per spec §1
- `fisher@gse263850` documented as limitation
- `bug_fix_deltas.md` inlined into Limitations §10.5
- `benchmark/idk_if_needed/` no longer exists; contents at `benchmark/docs/historical/` with README
- `benchmark/data/study1/timings_table.csv` exists with headline cell rows
- Working tree clean on main after merge

---

## Notes for the implementer

1. **`regen_all.py --verify` is the safety net.** Run it after every paper edit. If it fails, the cited number drifted from the parquet — don't push past it.

2. **The parallel-column prose at `benchmark/data/study1b_simulator/parallel_column_summary.md` is gold.** It's a well-drafted 7-tool comparison narrative — pull into the paper Results.

3. **The stash decision in Task 1 affects later tasks.** If you pop the stash, `benchmark/docs/timing-comparison.md` becomes a source for Methods §2.4. If you drop it, Methods §2.4 cites `timings_table.csv` only.

4. **The original 13-task plan at `2026-06-01-paper-finish.md` is superseded by this revision** — kept for historical reference. Do not execute it.

5. **20-seed sweep ran with `power_stack="off"` and `power_stack="lr+"` cases?** Worth verifying. If the eval parquet only has bare-lr results, an lr+ pass on the simulator data may be a stretch goal. Check `eval_simulator_intrinsic_per_seed.parquet`'s `test` column for what was actually scored. If lr+ is missing and the paper cites it, regenerate before Task 6.

6. **The 1.1 P2 hygiene items** (from the Phase 4 plan §11) are NOT in scope. They include: `_BETA_EPSILON` constants split, `pct_sig` knob removal, CGI shore/shelf width kwargs, etc. Pick up post-paper if desired.
