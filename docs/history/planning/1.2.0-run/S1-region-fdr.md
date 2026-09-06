# S1: add opt-in region FDR and chain_merge permutations

Start from main after L merges. Use `salvage-1-region-fdr` and target main.
Read the [run rules](README.md), issue [17](https://github.com/d-ert/epykit3/issues/17), and `docs/review/2026-09-05-unmerged-branch-triage.md`.

Own `src/epykit/dmr.py`, `tl.dmr` and its required imports in `src/epykit/tl.py`, the DMR tests, DMR and architecture docs, and this layer's changelog entry.
The refactor stack owns `dmc.py` and DMC orchestration. Keep the CLI tile-only empirical FDR behavior.

## Port the behavior, not the whole commits

Read commits `b36e7eea9ac5aec099fe2cac5fb257f19e0e574d` and `d5cb25448e5067ac3256abb9ce1c542a20ee68ba`.
Port the owned changes into current main. Do not cherry-pick a commit that also edits another stack's files, restores folded tests, or copies obsolete docs.

1. Add `fdr_method: Literal["max_t", "region"] = "max_t"` to the tile empirical FDR function and as a keyword-only option to `tl.dmr`. Preserve existing positional argument meaning.
2. Port `_region_count_ratio_fdr`, `_is_self_or_mirror_perm`, and the shared aggregation helper. Validate the method and positive permutation count before expensive work.
3. Add `empirical_fdr_for_chain_merge` with the same default. Admit `tl.dmr(method="chain_merge", empirical_fdr=True)`. Sliding-window and segment methods still raise.
4. Record the chosen method, requested permutation count, seed, and set estimate on the empirical paths in `md.uns["dmr_params"]`.
5. Keep the region-mode small-sample warning. Log the genome-wide recomputation cost once per run at INFO.

## Keep the statistics explicit

For `max_t`, preserve the current random assignments, denominator, NaN handling, empirical p-values, and their BH transform.
Do not exclude self or mirror assignments from this legacy method.
A new constant `empirical_fdr_set` column is NaN for `max_t`. This is an additive schema change on an existing opt-in path.

For `region`, use mean null-survivor counts divided by observed-survivor counts at each threshold, followed by the donor's monotone adjustment.
Document that `empirical_pvalue` is the pooled-null tail fraction in this mode, a diagnostic rather than a calibrated individual-region p-value.
Exclude self or mirror assignments and failed runs. Count a successful run with zero survivors as a zero contribution.
If no usable assignments remain, return NaN estimates with a warning. Preserve nonfinite observed values as NaN.
Describe these limits in the docs. The LR hash gate alone does not establish FDR calibration.

## Replay the observed chain_merge analysis

The donor helper calls the engine eagerly and omits the observed DMC correction step. Correct both during the port.

- Obtain the observed chromosome set from `DMCStore.chroms()` or the full materialized DMC table. Do not derive it from surviving DMRs.
- Replay the canonical DMC metadata: `test_used`, `unite`, both minimum-sample counts, `dispersion`, `reference`, `smoothing`, `smoothing_span_bp`, `sep_fallback`, `sep_threshold`, and the DMC multiple-testing method.
- Stream each permutation through `process_chromosomes_dmc(return_store=True)`, apply the observed DMC multiple-testing correction, then chain-merge and apply the same region filters. Do not pass DMR `fdr_method="region"` as the DMC correction method.
- Keep the observed store intact. Per-permutation caches must not overwrite its manifest or parquet files, including with `perm_n_jobs>1`.
- Reject GLM, contrast, unavailable metadata, and `use_smoothed=True` before work. Reconstructing a missing temporary pseudo-count store is outside this layer.
- Raw p-value and q-value columns remain authoritative when neighbour-combined columns also exist. Do not substitute combined columns during a permutation.
- Preserve within-stratum assignments. If a requested strata column is missing or does not cover the samples, raise instead of silently using unrestricted shuffles.
- If an explicit chromosome restriction differs from the observed DMC universe, reject it and tell the caller to rerun DMC with that restriction.

## Accept when

Port the donor's focused tests into `test_region_count_ratio_fdr.py`, `test_dmr_region_fdr_mode.py`, and `test_chain_merge_empirical_fdr.py`.
Update the API cases in `tests/test_empirical_fdr_method_coverage.py`; keep its CLI refusal cases.
Preserve the existing tests in `tests/test_dmr_empirical_fdr.py`. The donor's old P0 test file was folded into that file.

Prove these behaviors with small deterministic cases:

- Explicit `max_t` and omitted method match main, including failed and empty permutations.
- Region mode distinguishes failed, empty, and excluded assignments and handles zero usable assignments.
- Chain-merge permutations receive the observed knobs, chromosome universe, DMC correction, and region filters.
- An observed store remains readable and unchanged after serial and parallel permutations.
- A slow-marked synthetic end-to-end case exercises chain_merge FDR. The seed gives reproducible output.
- Unsupported modes and missing strata fail before a permutation starts.

Update the DMR and architecture docs. Port the donor design note with its defaults and paths corrected throughout.
R2 owns the cross-reference in `empirical_fdr_for_dmc`; specify that only the DMR API has region mode.
Run all code-layer gates.

PR title: `Add opt-in count-ratio FDR and chain_merge permutations`.
