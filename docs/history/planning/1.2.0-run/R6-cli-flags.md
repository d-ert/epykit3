# R6: add DMC filtering and supported CLI options

Start from R5 on the integrated base that already contains S1 through S3.
Use `refactor-6-cli-flags` and target `refactor-5-compatibility`.
Read the [run rules](README.md) and issues [27](https://github.com/d-ert/epykit3/issues/27) and [28](https://github.com/d-ert/epykit3/issues/28).

Own DMC config, stages, engine, related cache and metadata tests, the CLI package, README, and DMC and CLI docs.
S2's ingestion and tile DMR APIs must exist before this layer starts.

## Add the DMC filter

1. Add keyword-only `canonical_only: bool = False` to `tl.dmc` and `process_chromosomes_dmc`, including its overloads. Preserve existing positional argument meaning.
2. Carry the field through `DMCConfig` and both stage paths, `run_engine` and `run_contrast`.
3. Filter auto-detected chromosomes before engine execution and global multiple-testing correction. Preserve explicit chromosome lists, including an empty list. Retain valid empty-store behavior if no canonical contigs remain.
4. Record `canonical_only` through `to_uns` on binary, resume-hit, and contrast paths. Update `CANONICAL_UNS_KEYS` in `tests/test_dmc_metadata.py` and its path assertions.
5. Include `canonical_only` in `resume_signature_params`. Also key the now-exposed count-smoothing setting and effective smoothing span. Changing a setting must not return a prior full-table resume result.
6. Ensure the resolved chromosome list remains part of the low-level cache identity. Forward the same list and settings into DMC permutations and preserve S1's observed/null chromosome contract.

## Wire the actual CLI call paths

Current binary CLI DMC calls `process_chromosomes_dmc` directly. Only the formula path calls `tl.dmc`.
Current CLI DMR calls the low-level region engines. Do not assume either handler is a universal high-level wrapper.

- Add `--canonical-only` with default false to convert, dmc, and dmr.
- Convert forwards to `convert_sample`.
- Binary DMC forwards to `process_chromosomes_dmc`. Formula DMC forwards to `tl.dmc`; do not silently discard it during contrast dispatch.
- Tile DMR forwards to the tile caller and its empirical permutations. Reject the flag for DMC-derived DMR methods, as S2 specifies. Tell users to filter upstream DMC.
- Add `--smoothing` and `--smoothing-span-bp` to dmc only, with defaults false and 500.
- Forward both smoothing options into the binary DMC engine. Require a positive span when smoothing is enabled.
- Document smoothing as an LR option. If a CLI user enables it with a non-LR engine or a formula/contrast call that does not consume smoothing, reject the combination instead of silently ignoring it.
- Do not add DMR smoothing flags, an `--all-contigs` flag, or lr+ CLI options.

## Accept when

Use a fixture with chr1 and a noncanonical contig.

- DMC default false retains both. True auto-detection retains chr1. Explicit lists override filtering.
- Both high-level and CLI binary and contrast paths honor canonical filtering.
- False-to-true-to-false resumable calls invalidate when required. A repeated identical call resumes.
- Smoothing enabled with two different spans cannot reuse the same high-level result cache.
- Metadata has the same key set on all three DMC paths.
- CLI parity checks compare results for binary DMC smoothing and canonical filtering, conversion filtering, and tile DMR filtering.
- Unsupported DMR and smoothing combinations fail with a useful message.
- Default outputs and LR hashes are unchanged. All code-layer gates pass.

PR title: `Add canonical filtering and supported DMC smoothing options to the CLI`.
