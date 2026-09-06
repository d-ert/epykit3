# S1: count-ratio empirical FDR and the chain_merge permutation harness

Branch `salvage-1-region-fdr` from `main` after L merged. Bottom of the salvage stack; PR against `main`. Decision: map ticket 17. Read `README.md` here, `CLAUDE.md`, and `docs/review/2026-09-05-unmerged-branch-triage.md` (the `feat/canonical-chrom-filter` section) first.

You own `src/epykit/dmr.py`, the `tl.dmr` function in `src/epykit/tl.py`, their tests, `docs/` pages for dmr, and `CHANGELOG.md`. You do not touch `tl.dmc`, `_run_dmc_contrast`, `dmc.py`, `cli.py`.

## Commits, in this order

1. `feat(dmr): count-ratio empirical FDR as an opt-in fdr_method`
   - `git cherry-pick -n b36e7ee`, then `ruff format` and `ruff check --fix` on the touched files, then resolve by hand. The triage report says the `dmr.py` conflict is formatting-only.
   - `empirical_fdr_for_dmr(..., fdr_method: Literal["max_t", "region"] = "max_t")`. The default is `max_t`, not `region` as on the branch: today's `empirical_pvalue` and `empirical_qvalue` must be byte-identical. `region` is the count-ratio target-decoy estimate. Keep the constant `empirical_fdr_set` column on both methods. Keep the small-n `UserWarning` (fewer than four samples per group).
   - `tl.dmr` gains `fdr_method` with the same default and records it in `md.uns["dmr_params"]`.
   - Drop the branch's `dmc.py` docstring edit; it is handed to the refactor stack (name it in `worker_done`).
2. `feat(dmr): permutation empirical FDR for chain_merge`
   - `git cherry-pick -n d5cb254`, same treatment. `empirical_fdr_for_chain_merge` becomes public; helpers `_is_self_or_mirror_perm`, `_region_count_ratio_fdr`, `_aggregate_region_perm_results`, `_chain_merge_perm_survivors` come with it.
   - `tl.dmr`: the `NotImplementedError` gate admits `chain_merge`; `sliding_window` and `segment` still raise, with the message updated. The chain_merge permutation reads `test_used`, `unite`, `min_samples_treatment`, `min_samples_control`, `dispersion` and `reference` from `md.uns["dmc"]` (the canonical record) and raises `NotImplementedError` for `glm` and `glm_contrast` as on the branch. `md.uns["dmr_params"]` gains `empirical_fdr`, `n_perm`, `perm_seed`, `fdr_method`, `empirical_fdr_set` on this path.
   - The per-run "each shuffle recomputes the genome-wide DMC" `UserWarning` becomes one `logger.info` line stating `n_perm` and that each permutation reruns the DMC engine.
3. `test(dmr): carry the region-FDR and chain_merge permutation tests`
   - From the branch: `tests/test_region_count_ratio_fdr.py` (unit tests on the two helpers, unchanged), `tests/test_dmr_region_fdr_mode.py` (invert the "region is the default" assertion; keep the rest), `tests/test_chain_merge_empirical_fdr.py`, and the edits to `tests/test_empirical_fdr_method_coverage.py` (chain_merge leaves the raise list; the slow end-to-end test on the synth bundle stays slow-marked). The two tests in `tests/test_dmr_empirical_fdr.py` need no change. Do not carry `test_chroms.py` or `test_canonical_chrom_filter.py`; they belong to S2 and L.
4. `docs(dmr): region FDR design note and changelog`
   - Re-add `docs/review/2026-06-08-region-empirical-fdr-design.md` from the branch, with a one-line header note that the default stayed `max_t`. Update the dmr docs page for `fdr_method` and chain_merge empirical FDR. CHANGELOG under `[Unreleased]` / `### Added`.

## Contract

Every existing number is unchanged: `max_t` stays default and the benchmark slice does not run empirical FDR. Only opt-in calls produce new columns or values. Regen hashes unchanged.

## Deliver

PR title: `Empirical FDR for DMRs: count-ratio estimator as an opt-in method, permutation harness for chain_merge`. Body per the run README, with the exact semantics of `empirical_pvalue` and `empirical_qvalue` under each method. Then `worker_done`, naming the `dmc.py` docstring hand-off.
