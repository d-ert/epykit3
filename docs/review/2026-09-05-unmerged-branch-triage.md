# Triage of the four unmerged branches from June and August

Date: 2026-09-05. Comparison base: `origin/main` at `6cfd610` (merged 2026-09-05).
This is a read-only review. No branch was checked out, merged or rebased; every
claim below comes from `git log`, `git diff`, `git show` and `git merge-tree`
against the remote refs, plus reading the current tree on `main`.

## Summary

| Branch | Tip | Size vs main | Results change | Rebase cost | Recommendation |
|---|---|---|---|---|---|
| `feat/canonical-chrom-filter` | `cd9f89b`, 2026-06-08, 4 commits | 27 files, +1870 / -93 | Yes. q-values and call sets move under the new `canonical_only=True` default. `empirical_qvalue` becomes a different statistic under the new `fdr_method="region"` default. | Medium. The four Python conflicts are formatting-only. Real work: one deleted test file, two doc conflicts, one moved directory, the lint ratchet. | Keep, in two parts. Re-apply the count-ratio FDR and chain_merge permutation commits on main. Land the canonical filter opt-in first and flip the default in its own PR after a real-data re-run. |
| `fine-tune` | `47d7878`, 2026-06-18, 1 commit | 6 files, +900 / -579 (uv.lock is -568) | Yes. Every `lr` p-value (`smoothing=True`), every q-value (`fdr_tsbh`), chain_merge DMR sets (`min_cpgs` 3), the permissive preset. | Low. Three formatting-only Python conflicts. uv.lock conflict resolves by taking main's. | Drop. `bench-tune` contains it and reverts its DMC part. |
| `bench-tune` | `66ee2ef`, 2026-06-19, 2 commits (includes `fine-tune`) | 6 files, +895 / -574 | Yes, DMR only. `tl.dmr` and CLI chain_merge `min_cpgs` 5 to 3, permissive `pct_sig` 0.5 to 0.4. | Same as `fine-tune`. | Re-implement the additive parts by hand (CLI smoothing flags, recorded smoothing keys, docstring fix). Treat `min_cpgs=3` as its own decision. Drop AGENTS.md and the lock. |
| `optimize/autonomous-v1` | `99b2b4f`, 2026-08-11, 1 commit | 2 files, +67 / -11 | Yes. `pvalue_combined` and `qvalue_combined` for every `neighbour_combine` run, and chain_merge silently starts consuming them. | Low. One formatting-only conflict. | Drop. If correlation-aware combining is wanted, re-implement it as an opt-in knob with calibration evidence. |

## What moved on main since June

Three merged PRs separate every branch from `main`. PR #3 (`b6757c2`) ran a
one-time `ruff format` over the tree and widened the lint selection to
F, I, W, UP, E and B. It rewrote `cli.py`, `convert.py`, `dmc.py`, `dmr.py`,
`pl/_compute.py` and `tl.py`, which is why every branch conflicts. PR #2
(`8d07a40`) moved `docs/superpowers/` to `docs/history/superpowers/` and folded
the ticket-named tests. `tests/test_p0_dmr_empirical_fdr_denominator.py` no
longer exists; its two tests now live in `tests/test_dmr_empirical_fdr.py`.
PR #4 (`6cfd610`) upgraded dependencies, set the Python floor to 3.10 and made CI
check `uv.lock`.

Two fixes relevant to the deferred list landed on 2026-06-07, before any of the
four branches forked. `298f258` made stratified permutation shuffle labels
within each stratum (M-DMR3). `7484fdf` switched `dmr_segment` to the signed
Stouffer combine (M-DMR6). None of the four branches touches either item.

## feat/canonical-chrom-filter

Four commits by Deniz, all dated 2026-06-08, forked from `59f7f58`:
`b36e7ee` (count-ratio FDR as the `empirical_fdr` default), `91ccf16` (design
note), `d5cb254` (permutation FDR for chain_merge) and `cd9f89b`
(`canonical_only`). The commit messages match the diff.

### What it changes

The branch adds `src/epykit/_chroms.py` with `is_canonical_chrom`,
`filter_canonical`, `filter_canonical_logged`, `CANONICAL_CHROM_CORES` and
`CANONICAL_CHROMS_UCSC`. The predicate accepts chr1 to chr22, X, Y, M and MT
under UCSC or Ensembl naming and rejects everything else.

In `dmc.py`, `process_chromosomes_dmc` gains `canonical_only: bool = True`.
When `chromosomes` is auto-detected the list is filtered and one INFO line names
the dropped contigs. An explicit `chromosomes=` list is honoured verbatim. The
`empirical_fdr_for_dmc` docstring is corrected to describe the new DMR default.

In `dmr.py`, `call_dmr_tile_based` gains the same `canonical_only=True`
parameter. `empirical_fdr_for_dmr` gains `fdr_method` with default `"region"`
and legacy `"max_t"`. New helpers: `_is_self_or_mirror_perm`,
`_region_count_ratio_fdr`, `_aggregate_region_perm_results`,
`_chain_merge_perm_survivors` and the public `empirical_fdr_for_chain_merge`.
The output gains a constant `empirical_fdr_set` column. A `UserWarning` fires
when either group has fewer than four samples, and another fires on every
chain_merge permutation run because each shuffle recomputes the genome-wide DMC.

In `tl.py`, `dmc`, `_run_dmc_contrast` and `dmr` gain `canonical_only=True`.
`dmr` also gains `fdr_method="region"` and its `NotImplementedError` gate now
admits `chain_merge`. `md.uns["dmr_params"]` gains `fdr_method` and
`empirical_fdr_set` on the tile path and `empirical_fdr`, `n_perm`, `perm_seed`,
`fdr_method` and `empirical_fdr_set` on the chain_merge path. The chain_merge
permutation path reads `test_used`, `unite`, `min_samples_*`, `dispersion` and
`reference` back from `md.uns["dmc"]` and raises `NotImplementedError` unless
the DMC test was `lr`, `welch_t` or `fisher`.

In `cli.py`, `convert` gains `--canonical-only` (default off), and `dmc` and
`dmr` gain `--all-contigs` (default off, so the filter is on). The CLI does not
expose chain_merge empirical FDR. In `convert.py`, `convert_sample`,
`ensure_converted_sample` and `_can_reuse_sample` gain `canonical_only=False`
and the per-sample manifest records it, so a changed value invalidates the
conversion cache. In `io.py`, `read_bismark`, `read_methyldackel` and
`read_combined_strand_bed` forward `canonical_only=False`. In
`pl/_compute.py`, `compute_manhattan_data` imports the shared list instead of
its own hardcoded copy; behaviour is unchanged.

Docs: CHANGELOG, README, `docs/advanced/architecture.md`, the dmc, dmr,
annotate, cli and read-bismark pages, and two design notes,
`docs/review/2026-06-08-region-empirical-fdr-design.md` and
`docs/superpowers/specs/2026-06-08-canonical-chrom-filter-design.md`.

### Whether it changes results

Yes, in two independent ways.

The `canonical_only=True` default in `process_chromosomes_dmc`, `tl.dmc`,
`_run_dmc_contrast`, `call_dmr_tile_based` and `tl.dmr` removes non-canonical
contigs from every default run on real data. Per-CpG `pvalue` and `meth_diff`
on the retained chromosomes do not move, because the engines work one
chromosome at a time and the default `eb` dispersion prior is fitted from the
chromosome-wide distribution of per-site phi inside `_score_finalize` (the
`process_chromosomes_dmc` docstring says "across chromosomes", but the code
pools per chromosome). Every
`qvalue` can move, because `apply_multiple_testing_correction` runs over the
whole store and the BH family shrinks. The DMR callers that read the DMC table
(`call_dmr_chain_merge`, `call_dmr_sliding_window`, the segment caller) inherit
the smaller set. On the branch's own GSE263850 run, 156 of 181 contigs were
dropped and 56 significant scaffold rows vanished; I did not reproduce this.

The `fdr_method="region"` default in `empirical_fdr_for_dmr` and `tl.dmr`
replaces the Westfall-Young min-P statistic with a count-ratio target-decoy
estimate. `empirical_pvalue` becomes the pooled-null tail fraction and
`empirical_qvalue` becomes the monotone suffix-min of decoy count over observed
count. These are different quantities from the old columns, not a recalibration
of them. The `"max_t"` path keeps the old numbers; I checked that its
aggregation in `_aggregate_region_perm_results` still counts only permutations
with at least one null region and does not apply the new self-or-mirror
exclusion.

The ingestion filter defaults to off and changes nothing unless requested.

### Test coverage

Five new test files with 30 test functions, one new slow test in an edited
file, and one edited legacy file. `tests/test_chroms.py` covers the
predicate and the order-preserving filter. `tests/test_canonical_chrom_filter.py`
covers the DMC default, the `canonical_only=False` override, the explicit
`chromosomes=` override, the audit log line, the tile DMR default and override,
and the ingestion opt-in. `tests/test_region_count_ratio_fdr.py` has nine unit
tests on `_is_self_or_mirror_perm` and `_region_count_ratio_fdr`, including a
hand-computed case, order preservation, capping at one and monotonicity.
`tests/test_dmr_region_fdr_mode.py` monkeypatches the tile caller and checks
that region mode finds signal where max-T saturates, that region is the default,
that max-T still saturates, the small-n warning and the unknown-method error.
`tests/test_chain_merge_empirical_fdr.py` monkeypatches
`_chain_merge_perm_survivors` for four integration checks.
`tests/test_empirical_fdr_method_coverage.py` drops `chain_merge` from the
raise list and adds a slow end-to-end test on the synth bundle. The branch
commit says the full not-slow suite passed with 497 tests on the June tree. I
did not run it.

Against current main, three things need attention. First,
`tests/test_p0_dmr_empirical_fdr_denominator.py` was edited on the branch to
pass `fdr_method="max_t"` but was deleted on main; the same two tests,
`test_emp_p_uses_n_perm_denominator_not_pooled` and
`test_emp_p_correct_on_pure_null`, now sit in `tests/test_dmr_empirical_fdr.py`
and call `empirical_fdr_for_dmr` without `fdr_method`. Under the branch's new
default they would run in region mode and fail, so the re-pointing has to be
redone there. Second, the regex in the coverage test,
`match="empirical_fdr.*tile"`, still matches the branch's new error text, so the
`sliding_window` and `segment` cases keep passing. The CLI gate test keeps
`chain_merge` in its raise list, which is correct because the CLI is unchanged.
Third, other main tests that exercise `empirical_fdr_for_dmr`
(`tests/test_empirical_fdr_stratified.py`, the rest of
`tests/test_dmr_empirical_fdr.py`) need a read to see which of them assert
min-P semantics. The synth fixtures use chr1 to chr5, so `canonical_only` does
not remove anything in the existing suite.

### Rebase cost

`git merge-tree --write-tree origin/main origin/feat/canonical-chrom-filter`
reports eight conflicting paths. `cli.py`, `io.py`, `tl.py`, CHANGELOG.md and
CLAUDE.md auto-merge.

- `src/epykit/convert.py`, `src/epykit/dmc.py`, `src/epykit/dmr.py`,
  `src/epykit/pl/_compute.py`: content conflicts. I re-ran main's
  `ruff check --fix` and `ruff format` on the merge-base and branch versions
  of each file and repeated the three-way merge; all four merge cleanly and the
  results compile. These are formatting-only conflicts.
- `tests/test_p0_dmr_empirical_fdr_denominator.py`: modify/delete. Real work,
  small (see above).
- `README.md`, `docs/advanced/architecture.md`: content conflicts in prose that
  PR #2 also rewrote. Real work, small.
- `docs/superpowers/specs/2026-06-08-canonical-chrom-filter-design.md`: added
  inside a directory main renamed. Move it to
  `docs/history/superpowers/specs/`.

Beyond conflicts, the new code has to pass the widened lint set. I saw at least
`from typing import Iterable` in `_chroms.py`, which the UP rules will want from
`collections.abc`, and quoted forward-reference annotations that the tree no
longer uses.

### Relation to the documented plan

The branch closes the "permutation FDR for non-tile callers" follow-up named in
`tl.dmr`'s `NotImplementedError` for `chain_merge` only; `sliding_window` and
`segment` still raise. It addresses part of M-DMR5 from the remediation
summary: the count-ratio construction applies the same post-filter to the
decoys, so the selective-inference objection no longer applies to
`empirical_qvalue`. The asymptotic `combined_qvalue` path that M-DMR5 also names
is untouched, and the fix route differs from the one the summary proposed
(correct over all candidates, drop the second BH). This is my reading of the
design note; the branch does not cite M-DMR5. M-DMR3 and M-DMR6 were already
fixed on main before the fork. The canonical filter is not in CONCERNS.md or the
remediation summary; it is new scope motivated by the GSE263850 annotation
warning. The CONCERNS.md item "Permutation FDR tested sparsely" asks for a
real-data slow-tier test; the branch adds a slow test on the synth bundle, not
real data.

### Recommendation

Keep it, as two separate PRs re-applied on current main rather than a rebase of
the whole branch.

The first PR takes `b36e7ee` and `d5cb254`: the count-ratio FDR core, the
`fdr_method` switch, the chain_merge harness and their five test files. The
code is self-contained, the helpers have hand-computed unit tests, and the
conflicts are formatting-only. Making `"region"` the default is a results
change on the empirical path only. `benchmark/scripts/regen_small.py` hashes
DMC output from `tl.dmc` and does not touch `tl.dmr`, so the hash slice should
not move; run it to confirm. The paper's DMR tables in
`benchmark/paper/report/REPORT.md` do not use `empirical_fdr`, so they are
unaffected. Re-point the two folded tests in `tests/test_dmr_empirical_fdr.py`
to `fdr_method="max_t"`.

The second PR takes `cd9f89b` with the default set to `canonical_only=False` on
every function, so it adds the helper, the flags and the audit log without
moving results. Flip the default in a third PR that also re-runs the real-data
cohort. That flip changes the default output of `tl.dmc`, `tl.dmr` and the CLI
on any assembly with scaffolds. The Piao simulator writes a single `chr1`
contig, so `regen_small.py` should not change; re-run it and commit the hash
file if it does. The GSE263850 numbers in REPORT.md section 3 (for example the
30,957 default `lr` calls in the dispersion-sensitivity table) were computed
over all 181 contigs and will change.

## fine-tune

One commit by Deniz on 2026-06-18, forked from `1f501e0`. The commit subject
ends with "Changed:" and the body is empty, so the diff is the only record.

### What it changes

In `tl.py`, `dmc` changes two defaults: `smoothing` from `False` to `True` and
`fdr_method` from `"fdr_bh"` to `"fdr_tsbh"`. `_run_dmc_contrast` changes its
`fdr_method` default the same way. Two further `md.uns["dmc"]` records gain
`smoothing` and `smoothing_span_bp` keys; main already records them in the
primary record. The `dmr` docstring corrects the permissive preset text
(`dis_merge_bp` 1000, not 200) and describes the new `min_cpgs` default.

In `cli.py`, `dmc --fdr-method` defaults to `fdr_tsbh`, new `--smoothing`
(default on) and `--no-smoothing` flags are added, and `smoothing` is forwarded
to `tl.dmc` and `process_chromosomes_dmc`.

In `dmr.py`, `DMR_PRESETS["permissive"]["pct_sig"]` drops from 0.5 to 0.4. A
new `_DMR_DEFAULT_CHAIN_MERGE_MIN_CPGS = 3` replaces `_DMR_DEFAULT_MIN_CPGS`
inside `resolve_layer_min_cpgs`, so `tl.dmr(method="chain_merge")` and
`epykit dmr` without a preset now use 3 CpGs instead of 5. The docstring that
warned against exactly this change is deleted.

`AGENTS.md` is a copy of CLAUDE.md addressed to Codex; it already lags main
(it names a py3.9 matrix and `--extra dev`).
`docs/superpowers/plans/2026-06-18-dmc-dmr-power-calibration.md` is a 740-line,
eight-task plan whose first task is to add calibration tests before any default
changes. `uv.lock` is regenerated against the June pyproject; it removes the
`gpu-jax` extra and adds `parallel`, both of which main's lock from PR #4
already reflects.

### Whether it changes results

Yes, on every default DMC run. `smoothing=True` changes every `lr` p-value; the
REPORT.md dispersion-sensitivity table records `lr / site + smoothing` at
10,691 calls and 5.0 percent recall of methylKit against 30,957 and 31.7
percent for the default. `fdr_tsbh` changes every q-value and is one of the four
`lr+` components that CLAUDE.md says must not be promoted to a default without
re-running the ablations. The `min_cpgs` and `pct_sig` changes alter the
chain_merge DMR set from `tl.dmr` and the CLI, and the permissive preset from
`call_dmr_chain_merge`.

### Test coverage

None added. On main, `test_contrast_path_writes_canonical_record` in
`tests/test_dmc_metadata.py` calls `tl.dmc` without `fdr_method` and asserts
`"fdr_bh"`. `tests/test_dmr_min_cpgs_parity.py` asserts
`resolve_layer_min_cpgs(None, None) == 5`, and `test_tl_dmr_chain_merge_bare_default_preserved`
and `test_cli_chain_merge_bare_default_is_five` pin the 5 at the API and CLI
layers. All of these fail on this branch. The plan's own first task,
`tests/test_power_calibration_defaults.py`, was never written.

### Rebase cost

`git merge-tree` reports `cli.py`, `dmr.py`, `tl.py` and `uv.lock` as content
conflicts and the plan document as a moved-directory conflict. The three Python
conflicts vanish after normalising both sides with main's ruff, so they are
formatting-only. `uv.lock` is a real conflict but the resolution is to discard
the branch version, because CI now checks the lock against `pyproject.toml`.

### Relation to the documented plan

None of M-DMR3, M-DMR5 or M-DMR6. It is a defaults calibration effort that
skips its own plan's test-first step and contradicts two documented positions:
the CLAUDE.md note that `lr+` components are research knobs, and the REPORT.md
ablation showing smoothing dilutes recall on real data.

### Recommendation

Drop the branch. `bench-tune` contains this commit and reverts its DMC part, so
anything worth keeping is assessed there. If the branch's defaults were adopted,
`regen_small.py` would change on both the `lr` and `lr_plus` hashes, and every
epykit row in REPORT.md and `benchmark/paper/paper.md` is claimed around
`smoothing=False` and `fdr_bh`.

## bench-tune

Two commits by Deniz: `fine-tune`'s `47d7878` and `66ee2ef` on 2026-06-19. The
second commit's subject says it "changed and then reverted some dmc setting",
and the diff confirms it: `fdr_method` returns to `fdr_bh` in `tl.dmc`,
`_run_dmc_contrast` and the CLI, and `smoothing` returns to `False` in `tl.dmc`
and the CLI.

### What it changes

The net diff against main is `fine-tune` minus the DMC default flips. What
remains: the `--smoothing` (now default off) and `--no-smoothing` CLI flags
with `smoothing` forwarded to `tl.dmc` and `process_chromosomes_dmc`; the
`smoothing` and `smoothing_span_bp` keys added to two `md.uns["dmc"]` records;
`_DMR_DEFAULT_CHAIN_MERGE_MIN_CPGS = 3` in `resolve_layer_min_cpgs`; the
permissive preset `pct_sig` 0.4; the `tl.dmr` docstring fixes; `AGENTS.md`; the
calibration plan; the stale `uv.lock`.

### Whether it changes results

Yes, on DMR only. `tl.dmr(method="chain_merge")` and `epykit dmr` with no
preset and no `--min-cpgs` accept 3-CpG regions that main rejects, through
`resolve_layer_min_cpgs`. `preset="permissive"` in `call_dmr_chain_merge` and
the high-level layers accepts chains where 40 percent of CpGs are significant
instead of 50. DMC output is identical to main.

Main's docstring on `_DMR_DEFAULT_MIN_CPGS` says the paper's chain_merge numbers
depend on the 5. I could not confirm that from the benchmark scripts on main:
`benchmark/scripts/run_chain_merge_replication.py` passes `min_cpgs=3`,
`dis_merge_bp=100` and `pct_sig=0.5` explicitly, and
`benchmark/scripts/sensitivity_sweep.py` calls `call_dmr_chain_merge` directly,
whose engine default is already 3. No benchmark script calls `tl.dmr` with the
bare default. So the 3 versus 5 question is a product-default decision, not a
paper-reproducibility one, unless a run cache not in git depends on it.

### Test coverage

None added. The three `min_cpgs` parity tests named under `fine-tune` fail. The
`fdr_bh` metadata test passes again. The preset-ordering test in
`tests/test_dmr_presets_and_diagnose.py` does not compare `pct_sig`, so the 0.4
passes it.

### Rebase cost

Identical to `fine-tune`: three formatting-only Python conflicts, a discardable
`uv.lock`, and the plan document to move under `docs/history/superpowers/plans/`.

### Relation to the documented plan

Same as `fine-tune`: none of the deferred M-DMR items. The `min_cpgs=3` change
is Task 2 of the branch's own calibration plan, applied without Task 1's tests.

### Recommendation

Do not rebase the branch; re-implement the additive parts by hand on main in
one small PR: the `--smoothing` and `--no-smoothing` CLI flags with the
forwarding, the two extra `uns["dmc"]` records, and the `tl.dmr` docstring
correction for the permissive preset. These change no output.

Open `min_cpgs=3` and `pct_sig=0.4` as a separate decision. If adopted, update
the three parity tests in `tests/test_dmr_min_cpgs_parity.py`, rewrite the
`_DMR_DEFAULT_MIN_CPGS` docstring, and re-run the chain_merge DMR benchmark.
`regen_small.py` hashes DMC only, so it will not move, but the DMR rows in
REPORT.md are claimed around the current presets. Drop `AGENTS.md`; CLAUDE.md
is the maintained copy and the branch's version is already stale. Drop the lock.
The calibration plan can be archived under `docs/history/superpowers/plans/`
if someone wants the eight-task outline; it is not a record of work done.

## optimize/autonomous-v1

One commit on 2026-08-11 by "EC2 Default User" from an EC2 hostname, forked
from `1f501e0`. No author identity, no PR, no linked plan. The message claims
it resolves severe FPR inflation and "fixes regional sensitivity metrics when
tuned with Optuna thresholds". There is no Optuna reference anywhere in the
repository.

### What it changes

In `dmc.py`, inside the per-chromosome loop of `combine_neighbour_pvalues`,
the branch inserts a correlation estimate before the sliding window. It marks
sites with finite `z` and `|z| < 2` as null, then for lags 1 to
`min(100, n - 1)` in index space accumulates `z[i] * z[i + lag]` into 50-bp
distance bins up to `neighbour_bp`. Bins with more than 50 pairs get the mean
product; other bins get 0. It then clips negatives to 0, forces bin 0 to 1.0
and enforces a non-increasing profile. In the window loop it builds the full
pairwise distance matrix of the window's valid sites, looks up the binned
profile for every pair, sums it into `z_var[i]`, and divides `z_sum` by
`sqrt(z_var)` instead of `sqrt(n_in)`. The sign-agreement gate
(`min_sign_agreement`, `n_agree`), `require_focal_signal`, `focal_p_thresh` and
the NaN handling are otherwise unchanged. The comment block explaining the D9
NaN fix is deleted but the behaviour it described is kept. No new parameter is
added, so there is no way to turn the change off, and the
`pvalue_combined_n_neighbours` audit column still reports `n_in` rather than
the effective count that now sits in the denominator.

In `dmr.py`, all three input branches of `call_dmr_chain_merge` change the
`sig_col` choice. With `use_q_for_sig=True` the column is `qvalue_combined` if
present, else `qvalue`, else `pvalue`. With `use_q_for_sig=False` it is
`pvalue_combined` if present, else `pvalue`. The cache key already includes the
column name, so a stale cache is not a risk, but the selection is silent.

### Whether it changes results

Yes. `pvalue_combined`, `qvalue_combined` and `qvalue_combined_reject` change
for every run with `neighbour_combine=True`, which includes
`power_stack="lr+"`, `"auto"` and `"conservative"` at small n. Raw `pvalue` and
`qvalue` do not move. `regen_small.py` hashes the `lr_plus` output, so its
`lr_plus` hash changes and CI fails until the hash file is regenerated. The
REPORT.md `lr+` row (406,515 calls, 92.9 percent recall) is claimed around the
Stouffer combine. The `dmr.py` change additionally makes `call_dmr_chain_merge`
consume the combined columns whenever they exist, which reverses the contract
in CLAUDE.md that downstream code must opt in to `_combined` columns
explicitly; every `lr+` run followed by chain_merge would produce a different
DMR set.

### The statistics in plain terms

Brown's method, as published, is for combining dependent p-values through
Fisher's statistic, with the chi-square scale and degrees of freedom adjusted
by the sum of pairwise covariances. What this commit implements is the Stouffer
analogue: keep the sum of signed z-scores but set its variance to the sum of
all pairwise correlations in the window instead of the count. Under
independence the two coincide, so the change only matters when the estimated
correlations are non-zero.

The current `combine_neighbour_pvalues` assumes independence inside the window
and compensates with two gates: the focal site must show its own signal and a
majority of neighbours must agree in sign. The docstring on main states the
independence violation openly, and `tests/test_neighbour_combine.py` pins that
statement. The commit keeps both gates and replaces only the denominator.

What the estimate assumes: that correlation depends on genomic distance only,
is the same everywhere on the chromosome, is non-negative, and never increases
with distance. Those are plausible for background CpG dependence. Where the
implementation departs from them: any two sites closer than 50 bp are treated
as perfectly correlated, because bin 0 is overwritten with 1.0 after the
estimate; that is the diagonal's value, not a between-site correlation, and it
removes power inside dense CpG islands. The null set is chosen by thresholding
the very statistic being combined, so true signal with small `z` leaks into the
estimate and sites near real DMCs are counted as null. The product moment is
computed over truncated `z` without dividing by the truncated variance, so it is
not a correlation coefficient on the usual scale. Chromosomes with fewer than
51 pairs in a bin silently fall back to independence for that bin. The
per-site pairwise matrix makes the loop quadratic in the window size, which
the original two-pointer design avoided.

The commit carries no test, no simulation, no FPR or TPR numbers, no ablation
against the existing sign-agreement gate, and no calibration on either the Piao
simulator or GSE263850. The claim that it resolves FPR inflation is therefore
unsupported. None of the existing tests on main would fail, because
single-site and NaN-neighbour cases reduce to a variance of 1.0 and no test
pins the multi-site Stouffer value; only the `lr_plus` hash in CI would notice.

### Test coverage

None.

### Rebase cost

`src/epykit/dmc.py` is a content conflict that disappears after normalising
both sides with main's ruff; `src/epykit/dmr.py` auto-merges. Low.

### Relation to the documented plan

It targets the independence violation acknowledged in the
`combine_neighbour_pvalues` docstring and the CONCERNS.md item on `lr+` FPR
drift. It is not one of M-DMR3, M-DMR5 or M-DMR6.

### Recommendation

Drop the branch. Do not cherry-pick either hunk. If a correlation-aware
denominator is wanted, re-implement it on main as an opt-in argument to
`combine_neighbour_pvalues`, keep bin 0 as an estimate rather than a constant,
estimate the profile from a held-out null rather than from `|z| < 2`, report the
effective count in the audit column, and ship it with a dispersion sweep on the
Piao simulator and a unit test with a hand-computed variance. Never let
`call_dmr_chain_merge` pick `_combined` columns implicitly; if that behaviour is
wanted it belongs behind an explicit `pvalue_col` argument.

## Suggested order

Three pieces are worth keeping, in this order. First, the count-ratio FDR and
the chain_merge permutation harness from `feat/canonical-chrom-filter`
(`b36e7ee`, `d5cb254`), because they close a documented follow-up, carry their
own tests and do not move any benchmarked number. Second, the canonical
chromosome helper and flags from `cd9f89b` with the default left off, followed
by a separate default-flip PR paired with a real-data re-run. Third, the
additive CLI and docstring pieces of `bench-tune`, re-typed on main, with the
`min_cpgs=3` question raised as its own decision. `fine-tune` and
`optimize/autonomous-v1` are dropped.

## How this was checked

All commands ran from the review worktree on 2026-09-05 with `/usr/bin/git`.
`git fetch origin` first. For each branch:
`git log --format='%H %an %ad %s%n%b' origin/main..origin/<branch>`,
`git merge-base origin/main origin/<branch>`,
`git rev-list --count origin/main..origin/<branch>`,
`git diff --stat origin/main...origin/<branch>`,
`git diff origin/main...origin/<branch> -- src/ tests/ CHANGELOG.md uv.lock`,
`git show 66ee2ef -- src/`, and
`git merge-tree --write-tree --name-only origin/main origin/<branch>`.
To separate formatting-only conflicts from real ones I ran, for every
conflicting Python file, `git show <merge-base>:<file>` and
`git show origin/<branch>:<file>` through main's
`uv run ruff check --fix --exit-zero --stdin-filename <file> -` and
`uv run ruff format --stdin-filename <file> -` (ruff 0.16.6), then
`git merge-file -p <main> <base-formatted> <branch-formatted>` and
`python -m py_compile` on the result; zero conflict markers and a clean compile
were the criterion. I read `git show origin/<branch>:<path>` for the new test
files, `AGENTS.md`, the calibration plan and the two design notes. On main I
read `docs/review/2026-06-06-remediation-summary.md`,
`docs/review/2026-06-06-epykit-peer-review.md`,
`docs/history/planning/codebase/CONCERNS.md`, `benchmark/scripts/regen_small.py`,
`benchmark/scripts/simulate_piao.py`,
`benchmark/scripts/run_chain_merge_replication.py`,
`benchmark/scripts/sensitivity_sweep.py`, `benchmark/paper/report/REPORT.md`,
and the relevant functions in `src/epykit/dmc.py`, `dmr.py`, `tl.py` and the
tests named above, using `grep -n` and `sed -n`. `git show --stat b6757c2`,
`8d07a40` and `6cfd610` established which main commits touched the conflicting
files, and `git log -S` located the M-DMR3 fix. No test suite was run and no
branch was checked out.
