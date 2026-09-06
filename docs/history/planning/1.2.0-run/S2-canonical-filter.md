# S2: canonical chromosome filter, opt-in

Branch `salvage-2-canonical-filter` from `salvage-1-region-fdr`. PR against `salvage-1-region-fdr`. Decision: map ticket 27. The helper `src/epykit/_chroms.py` is already on `main` from L; use it, do not re-add it.

You own `dmr.py`, `tl.dmr`, `convert.py`, `io.py`, `pl/_compute.py`, their tests and docs. The `tl.dmc`, `_run_dmc_contrast`, `process_chromosomes_dmc` and CLI parts of this feature are R6 in the refactor stack; do not touch them.

## Commits, in this order

1. `feat(dmr): canonical_only on the tile caller and tl.dmr, default off`
   - `call_dmr_tile_based(canonical_only: bool = False)` and `tl.dmr(canonical_only: bool = False)`. When on and `chromosomes` is auto-detected, filter with `filter_canonical_logged` and log one INFO line naming the dropped contigs; an explicit `chromosomes=` list is honoured verbatim. Record the value in `md.uns["dmr_params"]`.
2. `feat(convert,io): canonical_only at ingestion, default off`
   - `convert_sample`, `ensure_converted_sample`, `_can_reuse_sample` gain `canonical_only=False`; the per-sample manifest records it so a changed value invalidates the conversion cache. `read_bismark`, `read_methyldackel`, `read_combined_strand_bed` forward it.
3. `refactor(pl): compute_manhattan_data uses the shared canonical list`
   - Behaviour unchanged; the hardcoded copy goes.
4. `test: canonical filter, opt-in`
   - Carry `tests/test_canonical_chrom_filter.py` from the branch with the default-on assertions inverted; keep the override, explicit-list, audit-line and ingestion tests. Drop the DMC-default tests; R6 carries those.
5. `docs: canonical_only`
   - The read-bismark, dmr and architecture pages; CHANGELOG under Added.

## Contract

Default off everywhere, so no existing number moves. The synth fixtures use chr1 to chr5, so the suite sees no change even when on. Regen hashes unchanged.

## Deliver

PR title: `Opt-in canonical chromosome filter for ingestion and DMR calling`. Then `worker_done`.
