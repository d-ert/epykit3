# S3: correct the permissive preset documentation

Start from `salvage-2-canonical-filter` after its implementation and gates pass.
Use `salvage-3-bench-tune-additive` and target S2.
Read the [run rules](README.md) and issue [28](https://github.com/d-ert/epykit3/issues/28).

Own the `tl.dmr` docstring and corresponding DMR documentation only.

## Correct the donor plan

The smoothing metadata additions in `bench-tune` belong to DMC, not DMR.
Current `DMCConfig.to_uns` already records `smoothing` and `smoothing_span_bp`.
Current `tl.dmr` and `call_dmr_tile_based` have no smoothing parameters.
Do not add DMR parameters or record a smoothing setting that no DMR engine consumed.

Read donor commits `47d787823f7c2d4a509493a621d44823f38591de` and `66ee2ef67fb3e6a7c16a38fa68065110236a6ac0` for context.
Port only the doc correction that matches current `DMR_PRESETS`.
At the reviewed baseline, permissive uses `dis_merge_bp=1000` and `pct_sig=0.5`.
The `tl.dmr` docstring incorrectly says `dis_merge_bp=200`.
Keep the layer-level bare chain_merge minimum at 5.

## Accept when

The DMR docstring and docs match `DMR_PRESETS` and `resolve_layer_min_cpgs`.
No runtime code, defaults, metadata, or lockfile changes.
Run the documentation-only gates in the run rules.

DMC CLI smoothing is R6. DMR smoothing and the donor's changed defaults are outside this run.
If the integrated base already corrected the docstring, report this layer as already satisfied and omit an empty PR.

PR title: `Correct the documented permissive DMR preset`.
