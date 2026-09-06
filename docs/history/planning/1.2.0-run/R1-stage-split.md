# R1: split DMC orchestration into stages

Start from main after L merges. Use `refactor-1-stage-split` and target main.
Read the [run rules](README.md) and issues [13](https://github.com/d-ert/epykit3/issues/13) and [26](https://github.com/d-ert/epykit3/issues/26).
Read `src/epykit/_dmc_stages.py` at prototype commit `3836a6d88709d88bd51c552310e1b37c76cb82b5`.
Use its stage outline as a guide. Its line ranges and warning stacklevels are not executable specifications.

Own `tl.dmc`, `_run_dmc_contrast`, DMC-only helpers, `_dmc_config.py`, the new `_dmc_stages.py`, `.github/workflows/test.yml`, and their tests.
Preserve `tl.dmr` and shared helpers in `tl.py`.

## Implement

1. Add `--extra bam` to every Ubuntu entry in both dynamic matrix JSON arrays and to the slow job. Keep Windows entries without the BAM extra. Preserve the existing PR and main matrix sizes, event conditions, and Python versions.
2. Move DMC bodies into the nine stages: `plan_run`, `run_contrast`, `lookup_resume`, `open_input_store`, `run_engine`, `post_process`, `publish`, `persist_resume`, and `finish`.
3. Keep the prototype's five records, `TsvPlan`, `DMCPlan`, `ContrastDesign`, `ResumeTicket`, and `DMCOutcome`, with concrete annotations. Keep the public `tl.dmc` signature and defaults.
4. Make the contrast, resume-hit, and ordinary binary paths explicit in `tl.dmc`. `publish` alone writes the DMC result record through `DMCConfig.to_uns`.
5. Preserve TSV resolution before config validation, then contrast dispatch. Preserve materialization restrictions, return behavior, result keys, best-effort exports, and resume persistence after publication.
6. Keep contrast resume behavior unchanged. Run `post_process` while the temporary input store is still alive. Clean up that temporary store on both normal exit and exceptions.
7. Emit the transitional column warning through `finish` on a resume hit too, as recorded in issue 13. Preserve its existing behavior on other paths.
8. Check warning filenames against a caller outside epykit. Do not change every warning to stacklevel 3: shared helpers and contextmanager entry add different stack depths.
9. Keep shared TSV, warning, and sample-count helpers callable by DMR. Stages can import these helpers locally after `tl.py` initializes. Do not add a top-level circular import or duplicate their implementations.

## Accept when

- Existing metadata, resume, multigroup, CLI parity, and streamed DMC tests pass. Update mock targets only where the import location moved.
- One smoothed-store test checks coverage preservation, rounded and clipped pseudo-counts, raw-count fallback for missing smooth values, and directory cleanup. Include exceptional exit in the same fixture.
- Binary, contrast, and resume-hit metadata keep the same key set and values except the recorded resume-hit warning.
- `materialize=False` does not call `DMCStore.to_dataframe()` on the ordinary binary path.
- Both matrix branches parse as JSON with the intended extras. CI executes `test_asm.py`, `test_bam_io.py`, and `test_entropy.py` on Ubuntu.
- The code-layer gates pass with the BAM extra installed locally.
- Report the measured complexity of `tl.dmc` before and after. Do not change the global ceiling in this layer.

PR title: `Split DMC orchestration into stages and run BAM tests in CI`.
