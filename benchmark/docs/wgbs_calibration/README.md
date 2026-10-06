# First WGBS calibration development results

Source: v4 frozen whole-genome benchmark; development branch `codex/wgbs-calibration-20261006`. Detailed rationale and reproduction: `docs/wgbs_calibration_development.md`.

## Chr22 development pilot

Five samples per group; 562,970 calibrated CpGs; 63 true signal regions; 99 complete label-shuffled scans, seed 20261006. One self/mirror assignment occurred and was retained in both modes. All sample counts were copied into private pilot stores. Truth was used only for scoring.

| Inference | Signal calls | Matched true regions | Region precision | Region recall | Null calls |
|---|---:|---:|---:|---:|---:|
| Existing asymptotic ranking / BH | 98 | 28 | 28.6% | 44.4% | 47 |
| Complete-pipeline permutation, count ratio | 33 | 25 | 75.8% | 39.7% | 0 |
| Complete-pipeline permutation, max-T | 8 | 8 | 100% | 12.7% | 0 |

Matches use the v4 scorer's 50% reciprocal overlap in eligible CpGs and maximum-cardinality one-to-one matching. The count-ratio path retains 25/28 of the original matched regions. Max-T is much more conservative on this input. The former remains an approximate FDR estimator; these observed precision values include localization mismatches and do not themselves equal a calibrated biological region FDR.

These are method comparisons on one already-inspected chromosome/seed. Permutation support existed before the new patch: the gain from activating that path is not attributable solely to the new handling of identity/mirror assignments and failed scans. No defaults were promoted based on this pilot. New-seed validation remains necessary before a performance or calibration claim.

## Whole-genome diagnosis

- Epykit: 317/323 false-positive CpGs are >1 kb from a true region.
- DSS: 103,628/112,195 false-positive CpGs are >1 kb from a true region.
- Null CpG false positives increase with observed mean coverage. The null raw-p rejection rate below .05 is 4.25% at 10–20×, 6.96% at 20–40×, 9.57% at 40–80×, and 16.78% above 80×.
- The previously available `bb_eb` option substantially suppresses the chr22 null tail; it also becomes strongly conservative. Its signal power must be investigated before any default change.

## Code and validation

The package patch includes all valid label assignments in count-ratio inference, refuses any failed/incomplete scan, validates scores in both modes, and corrects max-T documentation. The relevant tests, including engine replay and stratified assignment checks, passed: **84 passed**. `tests.log` contains the executed result. All CSVs and the pilot manifest are produced from the retained script outputs; bulk counts/caches remain under ignored `benchmark/results/wgbs_calibration/`.

The old whole-genome all-tool run remains incomplete on real methylKit. It is retained as the simulation reference; the branch has not declared it a completed six-tool real-data comparison.
