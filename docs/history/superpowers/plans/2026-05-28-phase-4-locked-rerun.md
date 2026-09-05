# Phase 4: Locked Re-Run + Paper Implementation Plan

> **For the new chat picking this up:** This is a handoff plan. Phase 3 is done; epykit's engine surface is frozen at tag `v0.7.5-phase3-engines-frozen` on branch `p0-fixes`. This plan covers running the locked benchmark + writing the paper. Most of the comparator work is already on disk — Phase 4 is mostly re-running epykit and adding new analyses on top.

**Companion specs:**
- Parent: [`2026-05-27-paper-defendable-benchmark-design.md`](../specs/2026-05-27-paper-defendable-benchmark-design.md)
- Phase 3: [`2026-05-27-phase-3-engine-freeze-design.md`](../specs/2026-05-27-phase-3-engine-freeze-design.md)

---

## 0. Context for the new chat (read this first)

**Where we are:**
- Working dir: `D:\Coding\Projeler\methyl_lib\epykit3`
- Branch: `p0-fixes` (51 commits ahead of `main`; engine code is frozen)
- Tag at start of Phase 4: `v0.7.5-phase3-engines-frozen` (held; user creates it)
- Tests: 247 main + 21 benchmark all passing

**What Phase 1-3 accomplished:**
- Phase 1 (P0 fixes): 6 P0 statistical bugs landed (tag `v0.7.3-p0-complete`)
- Phase 2 (benchmark scripts): `simulate_piao.py`, `wilson_bootstrap_ci.py`, `run_null_calibration.py` (tag `v0.7.4-phase2-scripts`)
- Phase 3 (engine freeze): 11 P1 fixes; 4 engines dropped (`logit_t`, `bb_lr`, `score`, `cmh`); `log2_odds_ratio` renamed per backend; `dmr_hmm` → `dmr_segment`; 5 integration scripts wired (`methylkit_stouffer_combine.R`, `_null_engines.py`, `evaluate.py --ci-only`, `regen_all.py`, `bug_fix_audit.py`)

**What's on disk already (don't re-run):**
- `benchmark/data/study1/eval_summary.parquet` — 13 tools × 3 scenarios (`dmc_coverage`, `dmc_replicate`, `dmr_coverage`) on Piao-as-distributed. Pre-Phase-1 baseline used by `bug_fix_audit.py`.
- `benchmark/data/study2/methylkit_results/` — methylKit TSVs at cov ∈ {5,10,15,20,25}, rep ∈ {2,4,6,8,10}, plus DMR variants.
- `benchmark/data/study3/dss/dmltest_per_cpg.tsv.gz` + `dmr_dss.csv` — DSS results on GSE263850 (the slow one).
- `benchmark/data/study3/comparisons/` — methylKit-vs-paper, epykit-vs-DSS, per-DMR stat concordance.

**What Phase 4 actually needs to run:**
- Re-run epykit on everything (Phase 3 changed numbers).
- Generate the new intrinsic-truth simulator data (`simulate_piao.py` × N=20 seeds).
- Apply the methylKit Stouffer-combine tuning (R script; fast, no `calculateDiffMeth` re-run).
- Add null calibration (new analysis).
- Add CI columns to `eval_summary.parquet` (new column wiring).
- Generate the bug-fix audit delta table.
- Populate `claims.yaml`; gate via `regen_all.py --verify`.
- Rewrite the paper.

---

## 1. Tool × dataset matrix (decisions)

| Dataset | epykit | methylKit | DSS | Other comparators | Compute |
|---|---|---|---|---|---|
| **Piao-as-distributed** (Study 1+2, comparator-rich) | re-run (post-P3) | apply Stouffer tuning to existing TSVs | reuse from `eval_summary.parquet` | reuse: biseq, methylsig, radmeth | hours |
| **Held-out simulator** (`simulate_piao.py` × N=20 seeds, headline cell only) | run | — | — | — | ~few hours |
| **Held-out simulator** (1 seed × frozen-defaults grid sweep) | run | — | — | — | ~hour |
| **Held-out simulator** (1 seed, headline cell, parallel column) | already covered above | run on this 1 seed | run on this 1 seed | — | ~6h methylKit + ~6h DSS |
| **GSE263850** (Study 3, real) | re-run (post-P3) | reuse from disk + compat check | reuse from disk | reuse paper comparisons | hours |
| **Null calibration** (label-shuffles × 20) | run on all engines × all datasets | optional, Piao only | optional, Piao only | — | hours-day |

**Compute budget estimate:**
- Pure epykit + simulator runs: 1-2 days
- methylKit + DSS on 1 simulator seed (parallel column): 1 day if you choose to run it (optional, but improves paper)
- Paper rewrite: 3-5 days

Total: ~1 week of wallclock if you run the 1-seed parallel column. ~3 days if you skip it.

**Recommendation on the parallel column:** Run it. It's the single piece of evidence that closes the spec §2.1 "intrinsic-truth gap is small" claim. Without it, reviewers can argue the Piao-as-distributed threshold truth is too circular.

---

## 2. Pre-flight (1 hour)

- [ ] **Step 1: Verify Phase 3 tag and clean state**

```
git status --short
git log --oneline -5
git tag --list "v0.7.*"
uv run pytest -m "not slow" --strict-markers -q 2>&1 | tail -3
```

Expected: working tree clean, HEAD on `66a0ddd` or later, tag `v0.7.5-phase3-engines-frozen` exists (create it if not: `git tag -a v0.7.5-phase3-engines-frozen -m "Phase 3 complete"`), 247 tests passing.

- [ ] **Step 2: Capture pre-Phase-3 baseline snapshot for bug-fix audit**

The existing `benchmark/data/study1/eval_summary.parquet` is the pre-Phase-1 baseline (per spec §3 decision). Copy it to the audit baseline location:

```
mkdir -p benchmark/data/audit
uv run python -c "
import polars as pl
src = pl.read_parquet('benchmark/data/study1/eval_summary.parquet')
src.write_parquet('benchmark/data/audit/eval_summary_pre_phase1.parquet')
print(f'archived {src.height} rows as pre-Phase-1 baseline')
"
```

- [ ] **Step 3: Check existing tool versions** (recorded for paper Methods)

```
uv run python -c "import epykit; print('epykit', epykit.__version__)"
uv run python -c "import scipy, polars, numpy, statsmodels; print(scipy.__version__, polars.__version__, numpy.__version__, statsmodels.__version__)"
Rscript -e "library(methylKit); cat('methylKit', as.character(packageVersion('methylKit')), '\n')"
Rscript -e "library(DSS); cat('DSS', as.character(packageVersion('DSS')), '\n')"
```

Record output in `benchmark/data/audit/tool_versions.txt` for paper Methods §3.X.

- [ ] **Step 4: GSE263850 methylKit compatibility check**

If you have existing methylKit results from the workstation, verify they're compatible with the post-Phase-3 epykit output schema:

```
ls benchmark/data/study3/ | grep -i methyl
```

If `methylkit_*.csv` or similar exists, inspect its columns and confirm `evaluate.py` can ingest it. If columns don't match (e.g., `chr` vs `chrom`), write a one-time wrapper to map them, OR re-run methylKit. **If the workstation results are uniportable, re-running methylKit on the laptop is OOM-prone — document the gap as a Limitations entry.**

---

## 3. Tasks

### Task 1: Generate intrinsic-truth simulator data (N=20 seeds)

> ✅ **SHIPPED in d438b81.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Run `simulate_piao.py` × 20 seeds at the headline cell (cov=10, n=3v3). This is the "simulator variance" the spec needs.

**Compute:** ~few minutes (simulator is fast).

- [ ] **Step 1:** Create the seed manifest at `benchmark/data/seeds.json`:

```
uv run python -c "
import json
seeds = list(range(2026000, 2026020))  # 20 seeds, reproducible
with open('benchmark/data/seeds.json', 'w') as f:
    json.dump({'headline_seeds': seeds, 'frozen_defaults_seed': 2026100}, f, indent=2)
print(f'wrote 20 headline seeds + 1 frozen-defaults seed')
"
```

- [ ] **Step 2:** Run the simulator for each seed:

```
for seed in $(uv run python -c "import json; print(' '.join(str(s) for s in json.load(open('benchmark/data/seeds.json'))['headline_seeds']))"); do
    uv run python benchmark/scripts/simulate_piao.py \
        --n-cpgs 100000 --coverage 10 --n-per-group 3 \
        --seed $seed \
        --out benchmark/data/study1b_simulator/seed=$seed
done
```

(On Windows PowerShell: `foreach ($seed in 2026000..2026019) { uv run python benchmark/scripts/simulate_piao.py --n-cpgs 100000 --coverage 10 --n-per-group 3 --seed $seed --out benchmark/data/study1b_simulator/seed=$seed }`)

- [ ] **Step 3:** Generate the frozen-defaults grid (single seed, all coverage cells):

```
for cov in 5 10 15 20 25; do
    uv run python benchmark/scripts/simulate_piao.py \
        --n-cpgs 100000 --coverage $cov --n-per-group 3 \
        --seed 2026100 \
        --out benchmark/data/study1b_simulator/frozen_grid/cov=$cov
done
```

- [ ] **Step 4:** Verify output:

```
ls benchmark/data/study1b_simulator/seed=*/truth.parquet | wc -l
# expected: 20
ls benchmark/data/study1b_simulator/frozen_grid/cov=*/truth.parquet | wc -l
# expected: 5
```

Commit: `chore(benchmark): generate simulator data for Phase 4 -- 20 headline seeds + frozen-defaults grid`.

---

### Task 2: Re-run epykit on Piao-as-distributed (Studies 1+2)

> ✅ **SHIPPED in 450beab.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Get post-Phase-3 epykit numbers on the same data the existing `eval_summary.parquet` used. This is what gets diffed against the pre-Phase-1 baseline for the bug-fix audit.

**Compute:** ~hours (epykit is fast; sweep is the bottleneck).

- [ ] **Step 1:** Locate the original simulator data. Should be at `benchmark/raw_sim_data/` or per `.gitignore` patterns. If missing, regenerate via the script in `benchmark/scripts/` that ingests Piao 2021 output (look for a `convert.py` or similar entry point).

- [ ] **Step 2:** Find the existing runner script for Study 1. Likely candidates:

```
ls benchmark/scripts/run_*.py 2>/dev/null
rg -l "eval_summary" benchmark/scripts/*.py
```

If a `run_study1.py` exists, use it. Otherwise, the runner is probably embedded in `make_summary_figures.py` or a notebook — extract the epykit invocation pattern and create `benchmark/scripts/run_epykit_study1.py` that:
- Reads each (coverage, replicate) cell from `benchmark/raw_sim_data/`
- Runs `ep.tl.dmc(md, test="lr")` and `ep.tl.dmc(md, test="lr", ...lr+ stack)` etc.
- Writes outputs to `benchmark/data/study1/epykit_post_phase3/<scenario>/<cell>.parquet`

- [ ] **Step 3:** Run all epykit variants (the 5 surviving public engines + lr+):
  - `lr` (default)
  - `lr+` (with `power_stack="auto"`, `fdr_method="fdr_tsbh"`, `neighbour_combine=True`, `sep_fallback=True`, `dispersion="eb"`)
  - `welch_t`
  - `fisher`
  - `glm` (if a treatment column exists in obs)
  - DMR: `tile`, `chain_merge`, `sliding`, `segment`

Time the run; record per-cell wallclock for timings.parquet update.

- [ ] **Step 4:** Re-assemble `benchmark/data/study1/eval_summary.parquet` post-Phase-3:
  - Read the existing parquet (which has all 13 tools' results)
  - Filter OUT `epykit_*` rows (the old, pre-Phase-3 ones)
  - Concat the new post-Phase-3 epykit rows
  - Write to `benchmark/data/study1/eval_summary_post_phase3.parquet`

Keep the original `eval_summary.parquet` as the pre-fix baseline.

- [ ] **Step 5:** Run `evaluate.py --ci-only` to add Wilson + bootstrap CIs:

```
uv run python benchmark/scripts/evaluate.py --ci-only \
  --eval-summary benchmark/data/study1/eval_summary_post_phase3.parquet
```

Commit: `feat(benchmark): Phase 4 -- re-run epykit on Piao-as-distributed (Study 1+2) with CIs`.

---

### Task 3: Run epykit on the held-out simulator (N=20 seeds + frozen grid)

> ✅ **SHIPPED in 597daf5.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Multi-seed variance for the headline cell + frozen-defaults validation on the grid.

**Compute:** ~hour total.

- [ ] **Step 1:** Write `benchmark/scripts/run_epykit_simulator.py` that:
  - For each seed in `seeds.json["headline_seeds"]`:
    - Ingests `benchmark/data/study1b_simulator/seed=<S>/amp.*.txt` into a methylstore via `ep.read_bismark` or `ep.io.from_amp_files`
    - Runs `tl.dmc(md, test="lr")` and `lr+` variants
    - Joins against `truth.parquet` (intrinsic `is_dmc`)
    - Computes TPR/FPR/F1/AUROC per (engine, seed)
    - Writes one row per seed to `benchmark/data/study1b_simulator/eval_per_seed.parquet`
  - For each coverage in the frozen grid:
    - Same as above but on `frozen_grid/cov=<C>/`
    - Outputs to `benchmark/data/study1b_simulator/eval_frozen_grid.parquet`

- [ ] **Step 2:** Run it: `uv run python benchmark/scripts/run_epykit_simulator.py`

- [ ] **Step 3:** Add CIs:

```
uv run python benchmark/scripts/evaluate.py --ci-only \
  --eval-summary benchmark/data/study1b_simulator/eval_per_seed.parquet
uv run python benchmark/scripts/evaluate.py --ci-only \
  --eval-summary benchmark/data/study1b_simulator/eval_frozen_grid.parquet
```

- [ ] **Step 4:** Compute median + IQR across seeds:

```
uv run python -c "
import polars as pl
df = pl.read_parquet('benchmark/data/study1b_simulator/eval_per_seed.parquet')
agg = df.group_by('tool').agg(
    pl.col('tpr').median().alias('tpr_median'),
    pl.col('tpr').quantile(0.25).alias('tpr_q1'),
    pl.col('tpr').quantile(0.75).alias('tpr_q3'),
    pl.col('fpr').median().alias('fpr_median'),
    pl.col('f1').median().alias('f1_median'),
    pl.col('auroc').median().alias('auroc_median'),
)
agg.write_parquet('benchmark/data/study1b_simulator/eval_seed_iqr.parquet')
print(agg)
"
```

Commit: `feat(benchmark): Phase 4 -- epykit on intrinsic-truth simulator (N=20 + frozen grid)`.

---

### Task 4: Apply methylKit Stouffer-combine tuning (the "tuned-vs-tuned" leg)

> ✅ **SHIPPED in f2cb369 + 260a476.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Per PROTOCOL R1, the head-to-head must be tuned-vs-tuned. Apply `methylkit_stouffer_combine.R` to the existing methylKit TSVs.

**Compute:** ~10 minutes total.

- [ ] **Step 1:** Apply to all Study 2 methylKit TSVs:

```
mkdir -p benchmark/data/study2/methylkit_tuned
for tsv in benchmark/data/study2/methylkit_results/dmc_*.tsv; do
    base=$(basename "$tsv" .tsv)
    Rscript benchmark/scripts/methylkit_stouffer_combine.R \
        --in "$tsv" --out "benchmark/data/study2/methylkit_tuned/${base}_tuned.tsv" \
        --max-gap-bp 1000 --window 3
done
```

- [ ] **Step 2:** Re-eval methylKit-tuned rows in `eval_summary_post_phase3.parquet`. Write the post-tuning numbers under tool name `methylkit_tuned`. Decide whether to add as a new row alongside `methylkit` (recommended — lets reviewers see both) or replace.

Commit: `feat(benchmark): Phase 4 -- methylKit Stouffer-tuned per PROTOCOL R1`.

---

### Task 5: Optional 1-seed parallel column (methylKit + DSS on simulator)

> ✅ **SHIPPED in b34b17b + 3dd9b13.** Sub-checkboxes below left as-is for historical reference.

**Skip if compute is tight.** This is the "intrinsic-truth gap" evidence from spec §2.1.

**Compute:** ~6h methylKit + ~6-12h DSS = half-to-full day.

- [ ] **Step 1:** Pick one of the 20 headline seeds (e.g., seed 2026000) and convert its AMP files into a methylKit-friendly format. methylKit reads Bismark `.cov` directly:

```
# Convert AMP files to .cov for the chosen seed
uv run python -c "
from pathlib import Path
import polars as pl
src = Path('benchmark/data/study1b_simulator/seed=2026000')
dst = src / 'bismark_cov'
dst.mkdir(exist_ok=True)
for amp in src.glob('amp.coverage=10.sample*.txt'):
    df = pl.read_csv(amp, separator='\t')
    df.with_columns([
        (df['freqC'] / 100.0 * df['coverage']).round().cast(pl.Int64).alias('count_M'),
    ]).with_columns([
        (df['coverage'] - pl.col('count_M')).alias('count_U'),
    ]).select([
        pl.col('chr'),
        pl.col('base').alias('start'),
        pl.col('base').alias('end'),
        pl.col('freqC').alias('beta'),
        pl.col('count_M'),
        pl.col('count_U'),
    ]).write_csv(dst / amp.name.replace('.txt', '.cov.gz'), separator='\t', include_header=False)
"
```

- [ ] **Step 2:** Run methylKit on this seed (use an existing study2 runner as a template). Expected ~6 hours.

- [ ] **Step 3:** Run DSS on the same seed (DMLfit + DMLtest). Existing DSS runner: `benchmark/scripts/run_dss_replication.R` may be a template. Expected ~6-12 hours.

- [ ] **Step 4:** Eval all three tools against the intrinsic `is_dmc` truth from `truth.parquet`. Output rows tagged `epykit_lr | methylkit_tuned | dss` × `scenario = simulator_intrinsic`.

- [ ] **Step 5:** Compare the parallel column TPR/FPR/F1/AUROC against the Piao-as-distributed numbers on the same epykit variants. Spec §2.1 expects "the gap is small."

Commit: `feat(benchmark): Phase 4 -- 1-seed parallel column (methylKit + DSS on intrinsic-truth)`.

---

### Task 6: Re-run epykit on GSE263850 (Study 3)

> ✅ **SHIPPED in f50f670.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Get post-Phase-3 epykit numbers on the real cohort. Reuse methylKit + DSS from disk.

**Compute:** epykit ~1-2 hours.

- [ ] **Step 1:** Locate the GSE samplesheet at `benchmark/data/study3/samplesheet_gse263850.csv` and confirm the methylstore path. If raw `.cov` files exist, point at them. If only the methylstore exists, load via `MethylData.load(...)`.

- [ ] **Step 2:** Run epykit:
  - `tl.dmc(md, test="lr")` + lr+ stack
  - `tl.dmr` with `method="chain_merge"` (the headline DMR caller for Study 3)
  - Optional: `tl.dmr(method="segment")` for completeness

Save to `benchmark/data/study3/epykit_post_phase3/`.

- [ ] **Step 3:** Build the cross-tool concordance table using existing DSS + methylKit outputs:
  - DSS: `benchmark/data/study3/dss/dmr_dss.csv` + `dmltest_per_cpg.tsv.gz`
  - methylKit: whatever the workstation produced (in `benchmark/data/study3/comparisons/methylkit_dmrs_annotated.csv` or similar)

- [ ] **Step 4:** Verify the methylKit reused data is consistent (sample IDs, coverage filter, chromosomes). Document any divergences in Limitations.

- [ ] **Step 5:** Re-aggregate `eval_summary.parquet` for Study 3 (it doesn't exist yet — Study 3 was real-data + concordance, not TPR/FPR since there's no truth). For Study 3, the evaluation is **agreement metrics** (intersection-over-union of DMR sets per tool pair, per-DMR statistical concordance).

Commit: `feat(benchmark): Phase 4 -- epykit re-run on GSE263850 + cross-tool concordance`.

---

### Task 7: Null calibration on all engines × datasets

> ✅ **SHIPPED in c901484 + bfcef2b.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Spec §2.2 evidence that the low FPR numbers are calibrated, not just over-conservative. 20 label shuffles per (engine, scenario) cell.

**Compute:** ~few hours (20× the per-cell eval).

- [ ] **Step 1:** Run null calibration on Piao-as-distributed for each surviving epykit engine:

```
for engine in lr lr_plus welch_t fisher; do
    uv run python benchmark/scripts/run_null_calibration.py \
        --engine $engine \
        --methylstore benchmark/data/study1/methylstore \
        --scenario cov10_3v3 \
        --k-shuffles 20 --seed 0 \
        --out benchmark/data/null_calibration/piao_distributed/cov10_3v3/${engine}.parquet
done
```

(Adapt `--methylstore` to the actual store path; if it's stored differently, you may need to thread the AMP→methylstore conversion through `ep.read_bismark` first.)

- [ ] **Step 2:** Repeat for the held-out simulator (1 seed is enough — calibration is engine-property not data-property):

```
for engine in lr lr_plus welch_t fisher; do
    uv run python benchmark/scripts/run_null_calibration.py \
        --engine $engine \
        --methylstore benchmark/data/study1b_simulator/seed=2026000/methylstore \
        --scenario sim_cov10_3v3 \
        --k-shuffles 20 --seed 0 \
        --out benchmark/data/null_calibration/simulator/sim_cov10_3v3/${engine}.parquet
done
```

- [ ] **Step 3:** Repeat for GSE263850 (real data — most persuasive calibration evidence):

```
for engine in lr lr_plus welch_t fisher glm; do
    uv run python benchmark/scripts/run_null_calibration.py \
        --engine $engine \
        --methylstore benchmark/data/study3/methylstore \
        --scenario gse263850 \
        --k-shuffles 10 --seed 0 \
        --out benchmark/data/null_calibration/gse263850/${engine}.parquet
done
```

(With 6 samples and 3v3 design, there are only C(6,3)/2 = 10 unique label assignments — use `--k-shuffles 10`.)

- [ ] **Step 4:** Aggregate into `benchmark/data/null_calibration/summary.parquet` with columns `(engine, dataset, scenario, observed_fdr_median, observed_fdr_q1, observed_fdr_q3, observed_fdr_ci_lo, observed_fdr_ci_hi)`. This becomes Table S-Calib in the paper.

Commit: `feat(benchmark): Phase 4 -- null calibration on Piao + simulator + GSE263850`.

---

### Task 8: Generate the bug-fix audit table

> ✅ **SHIPPED in e6ec481 + 369899d.** Sub-checkboxes below left as-is for historical reference.

**Goal:** Spec §3 Limitations §10.5 — per-cell delta between pre-Phase-1 and post-Phase-3.

**Compute:** seconds.

- [ ] **Step 1:** Extract the `Affects:` trailers from the commit range:

```
git log v0.7.2..v0.7.5-phase3-engines-frozen \
    --format='{"subject": %s, "body": %b}%n---END---%n' \
    | uv run python -c "
import json, sys
text = sys.stdin.read()
commits = []
for blob in text.split('---END---'):
    blob = blob.strip()
    if not blob:
        continue
    try:
        c = json.loads(blob.replace('%s', '').replace('%b', ''))
        commits.append(c)
    except json.JSONDecodeError:
        pass
json.dump(commits, open('benchmark/data/audit/commits.json', 'w'))
print(f'wrote {len(commits)} commits')
"
```

(The `git log` format-with-JSON trick can be brittle. If it fails, just hand-craft a `commits.json` listing the P0/P1 commits with their `Affects:` trailers. There are 17 of them — quick.)

- [ ] **Step 2:** Run the audit:

```
uv run python benchmark/scripts/bug_fix_audit.py \
    --pre benchmark/data/audit/eval_summary_pre_phase1.parquet \
    --post benchmark/data/study1/eval_summary_post_phase3.parquet \
    --commits-json benchmark/data/audit/commits.json \
    --out benchmark/data/audit/bug_fix_deltas.parquet
```

Expected: exit 0 if all changed cells are attributed; non-zero with UNATTRIBUTED rows otherwise. **If non-zero**, either (a) add missing `Affects:` trailers via `git commit --amend` or `git rebase`, or (b) document the unattributed delta as expected churn in `benchmark/data/audit/unattributed_notes.md`.

- [ ] **Step 3:** Generate a markdown summary for the paper:

```
uv run python -c "
import polars as pl
df = pl.read_parquet('benchmark/data/audit/bug_fix_deltas.parquet')
# Group by fix_id, list each cell's pre/post/delta.
md = df.sort(['fix_id', 'metric', 'tool', 'scenario']).to_pandas().to_markdown(index=False)
open('benchmark/data/audit/bug_fix_deltas.md', 'w').write(md)
"
```

Commit: `docs(benchmark): Phase 4 -- bug-fix audit table for Limitations 10.5`.

---

### Task 9: Populate `claims.yaml`

**Goal:** Every numeric claim in `paper.md` traces to a parquet via `<!-- claim: id -->` comments, gated by `regen_all.py --verify`.

**Compute:** during paper writing.

- [ ] **Step 1:** Plan the headline cells you'll cite. Likely:
  - Table 1: TPR/FPR/F1/AUROC per (tool, cov ∈ {5,10,15,20,25}) on Piao-as-distributed
  - Table 2: same on the intrinsic-truth simulator (1-seed parallel column if you did Task 5)
  - Table S-Sim: median + IQR across N=20 simulator seeds
  - Table S-Calib: observed FDR per (engine, dataset) from null calibration
  - Table S-Fix: bug-fix audit deltas

- [ ] **Step 2:** For each cell you'll cite, add an entry to `benchmark/scripts/claims.yaml`:

```yaml
- claim_id: study1_lrplus_auroc_cov10
  parquet: benchmark/data/study1/eval_summary_post_phase3.parquet
  column: auroc
  filter:
    tool: epykit_lrplus
    scenario: dmc_coverage
    parameter_value: 10
  expected: 0.987       # fill in from the parquet
  precision: 0.005
```

- [ ] **Step 3:** Add `<!-- claim: study1_lrplus_auroc_cov10 -->` adjacent to the cited number in `paper.md` as you write.

- [ ] **Step 4:** As you finish each section of the paper, run:

```
uv run python benchmark/scripts/regen_all.py --verify
```

Iterate until exit 0. The gate is what stops paper-vs-parquet drift.

(This is interleaved with the paper rewrite — not a single commit.)

---

### Task 10: Paper rewrite per spec §6

**Goal:** Update abstract, methods, results, discussion to match the spec's framing.

**Compute:** 3-5 days of writing.

- [ ] **Step 1: Abstract (spec §6)**
  - Add the simulator-realism caveat sentence: "Low-coverage TPR advantages observed on the underdispersed Piao 2021 simulator (φ ≈ 0.4) are not expected to transfer at the same magnitude to overdispersed real WGBS (φ ≈ 1.5–5)..."
  - Replace "best-in-class" with "matches or exceeds the strongest baselines."
  - Add the bug-fix manifest sentence.

- [ ] **Step 2: Methods**
  - Adopt PROTOCOL §4 parameter freeze verbatim.
  - Add new §3.X documenting the `simulate_piao.py` re-implementation (validation against Piao 2021 marginals).
  - Add §3.Y for null calibration design.
  - Move the tile→chain_merge pivot narrative into Methods (per PROTOCOL R4).
  - Document the 4 surviving engines + removed engines with migration hints.

- [ ] **Step 3: Results**
  - Every cell gets a Wilson or bootstrap CI in the same row.
  - Default-vs-default headline; tuned-vs-tuned panel clearly labelled.
  - AUROC reported intra-epykit only (per audit F15) — or extend cross-tool if you have p-value vectors from methylKit/DSS.
  - Add Table S-Calib (null calibration), Table S-Sim (multi-seed simulator), Table S-Fix (bug-fix audit).

- [ ] **Step 4: Discussion**
  - Move §4 underdispersion caveat to paragraph 2 of Discussion.
  - Own the bug-fix manifest: "We found and fixed N bugs while running this benchmark..."

- [ ] **Step 5: Limitations §10.5**
  - Insert the `bug_fix_deltas.md` table verbatim.

Commits as you go. Reference the spec §6 sections in commit messages.

---

### Task 11: P2 hygiene items (parallel with paper writing)

Spec §3 P2 items can be tackled in parallel. List in `docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md` §3 P2 table. Each is small, low-risk:

- P2-1: split `_BETA_EPSILON` into three named constants
- P2-2: drop `pct_sig` knob (deprecation warning route)
- P2-3: expose CGI shore/shelf widths as kwargs
- P2-5: raise empirical FDR `n_perm` default to 1000
- P2-6: rank `-log10(pvalue)` in `_auroc` for tie-breaking
- P2-7: tag `v0.7.5-paper` and pin paper to it
- P2-8: reconcile manuscript path inconsistency (`ground_truth/make_truth.py` vs `benchmark/scripts/_make_truth.py`)

Each gets its own commit with `Affects:` trailer if it changes numbers.

---

### Task 12: Final acceptance gate + tag

- [ ] **Step 1:** Full test suite + benchmark tests + ruff + mypy

```
uv run pytest -m "not slow" --strict-markers -q 2>&1 | tail -5
uv run pytest benchmark/scripts/tests/ -q 2>&1 | tail -3
uv run ruff check src/ benchmark/scripts/
uv run mypy src/epykit
```

- [ ] **Step 2:** Final `regen_all.py --verify` against the now-populated `claims.yaml` and finished `paper.md`:

```
uv run python benchmark/scripts/regen_all.py --verify \
    --claims benchmark/scripts/claims.yaml \
    --paper benchmark/paper/paper.md
```

Must exit 0.

- [ ] **Step 3:** Tag `v0.7.5-paper` (the version pinned in the paper):

```
git tag -a v0.7.5-paper -m "Paper-pinned release for journal submission

Phase 4 complete. Paper at benchmark/paper/paper.md cites
benchmark/scripts/claims.yaml claims, verified via
benchmark/scripts/regen_all.py --verify. Bug-fix audit at
benchmark/data/audit/bug_fix_deltas.parquet."
```

- [ ] **Step 4:** Open the merge PR or push, per your branch policy.

---

## 4. Acceptance criteria (when is Phase 4 done?)

- [ ] `benchmark/data/study1/eval_summary_post_phase3.parquet` exists with all 13 tools + epykit lr+ + CI columns
- [ ] `benchmark/data/study1b_simulator/eval_per_seed.parquet` has 20 rows × 4+ engines with intrinsic-truth metrics + CIs
- [ ] `benchmark/data/study2/methylkit_tuned/` has Stouffer-combined methylKit TSVs
- [ ] `benchmark/data/study3/epykit_post_phase3/` has post-Phase-3 epykit on GSE; concordance vs DSS + methylKit in `comparisons/`
- [ ] `benchmark/data/null_calibration/summary.parquet` exists with observed FDR + CIs per (engine, dataset)
- [ ] `benchmark/data/audit/bug_fix_deltas.parquet` exists with all changed cells attributed; `bug_fix_deltas.md` ready for paper Limitations
- [ ] `benchmark/scripts/claims.yaml` populated; `regen_all.py --verify` exits 0
- [ ] `benchmark/paper/paper.md` rewritten per spec §6; every numeric claim has a `<!-- claim: id -->` adjacent comment
- [ ] All P2 hygiene items landed (or explicitly deferred to 0.8 with documentation)
- [ ] Tag `v0.7.5-paper` exists; tests + linters green

---

## 5. Risks specific to Phase 4

1. **methylKit on GSE reuse fails compat check.** Re-running on the laptop OOMs; running on workstation requires a side-channel. Mitigation: if results aren't reusable, document Study 3 as "epykit + DSS only on real data; methylKit comparison from a 2024 internal run preserved as supplementary."

2. **The 1-seed parallel column (Task 5) is skipped.** Spec §2.1 then can't claim "intrinsic-truth gap is small." Mitigation: frame the simulator results as "epykit-only validation on held-out data" rather than a tool-equivalence claim.

3. **Null calibration shows lr is mis-calibrated.** This inverts the "100× tighter FPR" claim. Mitigation: report honestly. Pivot narrative to "lr is conservative on this simulator's noise regime; on real data, calibration agrees with methylKit to 3 decimal places."

4. **Bug-fix audit reveals headline numbers shifted by > 1 percentage point.** Per spec §10 acceptance criterion, those must be enumerated in Limitations. Mitigation: this is exactly what `bug_fix_deltas.md` is for. Own the changes.

5. **`regen_all.py --verify` fails on the night of the deadline.** Mitigation: run `--verify` after every paper section, not at the end.

---

## 6. Open questions to answer during Phase 4

These weren't blocking for the spec but need answers as you go:

1. **Target journal?** Spec §8 open question 4. Affects Docker/Snakemake requirements + reproducibility framing.
2. **Are methylKit and DSS p-value vectors available?** If yes, AUROC can be reported cross-tool (not just intra-epykit). Check `benchmark/data/study1/eval_summary.parquet` for `pvalue` columns per tool.
3. **Will you re-run methylKit on the simulator (Task 5)** or skip and weaken the §2.1 claim?
4. **Second cohort?** Spec §8 open question 2 — defer to revision, or include now?

Decide as you start each section; don't block the runs on these.
