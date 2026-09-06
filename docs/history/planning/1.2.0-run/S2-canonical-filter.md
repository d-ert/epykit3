# S2: add opt-in canonical filtering at ingestion and tile DMR calling

Start from `salvage-1-region-fdr` after its implementation and gates pass.
Use `salvage-2-canonical-filter` and target S1.
Read the [run rules](README.md) and issue [27](https://github.com/d-ert/epykit3/issues/27).

Own `src/epykit/dmr.py`, `tl.dmr`, `convert.py`, `io.py`, `pl/_compute.py`, and their tests and docs.
Use L's helper. R6 owns DMC filtering and the CLI.

## Implement

1. Add keyword-only `canonical_only: bool = False` to `call_dmr_tile_based` and `tl.dmr`. Preserve positional arguments. Filter only auto-detected chromosomes before analysis and multiple-testing correction. An explicit list, including an empty list, takes precedence.
2. Define this DMR option as tile-only. DMC-derived methods use the chromosome universe chosen by the upstream DMC call. For those methods, reject `canonical_only=True` with instructions to run DMC with the filter. Do not silently ignore the option or filter a finished q-value table.
3. Use the same resolved chromosome list for observed tiles and S1 permutations. Record the option in tile `dmr_params`. Keep the default-off path unchanged.
4. Add and forward `canonical_only=False` through `convert_sample`, `ensure_converted_sample`, `_can_reuse_sample`, `read_bismark`, `read_methyldackel`, and `read_combined_strand_bed`.
5. Include the option in the per-sample conversion manifest and reuse check. Treat an old manifest without the key as false. When the option changes, regenerate the sample and remove stale partitions so excluded contigs cannot leak through a glob.
6. Make `pl/_compute.py` import the existing UCSC ordering from the shared helper. Preserve the plot's chromosome order.

## Accept when

Use the relevant cases from commit `cd9f89b3fd7a56bd3f9c7fbd2bbdf78aac5f1676` as input, then adapt them to current fixtures.
Own the ingestion and tile cases in `tests/test_canonical_chrom_filter.py`. R6 adds the DMC cases after S3 merges.

The fixture must contain a canonical chromosome and a noncanonical contig. A chr1-only fixture cannot prove the filter.

- Default false preserves both contigs. True removes only the noncanonical contig during auto-detection.
- Explicit chromosome lists take precedence in the tile caller. Unsupported non-tile use fails clearly.
- False-to-true-to-false conversion in the same cache gives the expected partition set at each step.
- A manifest without the new key remains reusable for false and is invalid for true.
- Tile permutation calls use exactly the observed chromosome universe.
- The helper emits one summary per selection operation, with no per-row log output.

Update the ingestion, DMR, and architecture docs. State the fixed human-style chromosome set and tile-only DMR scope.
The CLI documentation waits for R6. Run all code-layer gates.

PR title: `Add opt-in canonical filtering for ingestion and tile DMR calling`.
