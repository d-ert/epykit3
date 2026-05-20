# Changelog

All notable changes to **epykit** are tracked here. Format roughly follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project uses
SemVer (`MAJOR.MINOR.PATCH`).

## [Unreleased] — 0.7.2

Eight targeted fixes identified by a full benchmark comparison against
methylKit, DSS, RADMeth, BiSeq, methylSig, and Fisher (Piao et al.
2021). Two bugs fixed, three default changes, three guardrail
improvements. Existing user code that passes explicit parameters is
unaffected; only implicit defaults change.

### Fixed

- **Fisher calibration bug (`dmc.py`).** `fisher_exact_vectorized`
  computed only the upper tail via `hypergeom.sf(meth_a - 1, ...)`.
  When group A had less methylation (hypo direction), `meth_a - 1`
  could be -1, making `sf(-1, ...)` return 1.0. Fixed to compute both
  tails and pick the correct one based on whether observed >= expected.
  Fisher TPR went from 0.000 to 0.668-0.998 across all coverage levels.

- **chain_merge DMR fragmentation (`dmr.py`).** Default `dis_merge_bp=100`
  was too narrow for real genomic CpG spacing (median 100-200 bp in
  dense regions, 500-1000 bp in intergenic/intronic). DMRs spanning
  500-2000 bp with CpGs >100 bp apart were fragmented into 4-5 pieces.
  Widened defaults: `strict` 100->250, `default` 100->500, `permissive`
  200->1000, function default 100->500, `_SIG_DEFAULTS` 100->500.
  chain_merge TPR went from 0.086 to 0.971-1.000.

### Added

- **`power_stack` parameter on `ep.tl.dmc` (default `False`).**
  Pass `power_stack=True` to enable `neighbour_combine=True` and
  `sep_fallback=True` (the lr+ stack) in one switch. Pass
  `power_stack="auto"` to auto-enable the stack at n <= 2 replicates
  per group. Default is `False` — no silent auto-engagement.

- **Adjacent tile merging in `ep.tl.dmr(method="tile")`.**
  `merge_adjacent` parameter (default `True`) on `call_dmr_tile_based`
  merges adjacent significant tiles on the same chromosome with the same
  direction via Stouffer Z-combination. P-values are re-corrected with
  BH after merging. Pass `merge_adjacent=False` for exact old behaviour.

- **bb_lr auto-shrink guardrail (`dmc.py`).** When `n_samples < 6` and
  `dispersion="site"`, bb_lr auto-promotes to `dispersion="shrink"` with
  a warning, since per-site dispersion from Pearson residuals with
  `df_resid <= 4` is extremely noisy.

- **`dispersion="eb"` support in bb_lr (`_glm.py`).** The bb_lr path
  now accepts `dispersion="eb"` (empirical-Bayes shrinkage), matching
  the lr engine. Previously bb_lr only accepted `site/chrom/shrink`
  and raised `ValueError` on `"eb"`.

### Changed

- **`dispersion` default on `ep.tl.dmc` changed from `"site"` to
  `"eb"`.** Empirical-Bayes shrinkage of per-site quasi-binomial
  dispersion toward a chromosome-wide inverse-Gamma prior. At high n
  (large per-site df), the weight on the per-site estimate dominates
  and `"eb"` reduces to `"site"`. At low n, `"eb"` shrinks toward the
  chromosome mean, stabilising noisy dispersion estimates without losing
  site-level resolution. The change improves power at n <= 3 per group
  on real WGBS data with heterogeneous overdispersion; on underdispersed
  simulation data (Piao et al. 2021) it is a no-op (per-site phi is
  clamped at 1.0 either way). Pass `dispersion="site"` to restore
  exact pre-0.7.2 behaviour.

- **`method` default on `ep.tl.dmr` changed from `"tile"` to
  `"chain_merge"`.** On the Piao et al. 2021 DMR simulation,
  chain_merge recovers 97-100% of truth DMRs at all coverage levels
  while tile recovers 54-100%. chain_merge is also more adaptive to
  irregular CpG spacing. The tile engine remains available via
  `method="tile"`. Pass `method="tile"` to restore pre-0.7.2 behaviour.

- **welch_t warning tiers (`dmc.py`).** Split `_validate_sample_size_and_warn`
  into severity tiers: CRITICAL at `min_n <= 2` (degenerate
  Welch-Satterthwaite DOF, near-zero power), softer warning at `min_n < 6`.

- **bb_lr low-n warning (`dmc.py`).** Added specific warning when
  `test="bb_lr"` and `min_n < 3`: recommends `test="lr"` instead.

- **`dis_merge_bp` default in `ep.tl.dmr`** changed from 100 to 500 to
  match the updated chain_merge defaults.

- **Docstring tuning guidance in `dmr.py`** updated: "loosen dis_merge_bp
  100 -> 200" changed to "500 -> 1000" to reflect the new defaults.

### Tests

- `test_primitives.py`: `test_fisher_reverse_separation_significant` and
  `test_fisher_symmetry` covering the dual-tail Fisher fix.
- `test_dmr_chain_merge.py`: `test_chain_merge_default_merges_300bp_gap`
  verifying sig CpGs 300 bp apart chain at the new 500 bp default.
- New `test_dmr_tile_merge.py` (6 tests): adjacent tile merging,
  direction-aware non-merging, gap handling, multi-chromosome, three-way
  merge, empty input.

---

## [Unreleased] — 0.7.1

Targeted improvements to the `lr` DMC engine that close its
asymptotic-quasi-binomial gap to methylKit / RADMeth / DSS at low
coverage and small cohorts. All four are opt-in keyword arguments on
`ep.tl.dmc(test="lr", ...)`; defaults preserve 0.7.0 behaviour
exactly.

### Added

- **`dispersion="eb"` on `ep.tl.dmc`.** Empirical-Bayes shrinkage of
  per-site quasi-binomial dispersion toward a chromosome-wide
  inverse-Gamma prior whose pseudo-df is estimated from the
  per-site phi distribution via method-of-moments. Generalises the
  existing `dispersion="shrink"` mode (which uses a fixed
  pseudo-df = 4). No-op on data without genuine overdispersion.

- **`neighbour_combine=True`, `neighbour_bp=200` on `ep.tl.dmc`.**
  Signed-Stouffer Z combiner over neighbouring CpGs. Gated by
  `min_sign_agreement=0.6` (focal site's neighbours must agree on
  direction by ≥ 60 %) and `require_focal_signal=True` (focal raw
  p must be < `focal_p_thresh=0.5`) so spatially isolated false
  positives are not amplified. The output `pvalue` becomes the
  combined p; the raw is preserved as `pvalue_raw`. Two helper
  columns `pvalue_combined` and `pvalue_combined_n_neighbours`
  are added for audit. Exposed standalone as
  `ep.dmc.combine_neighbour_pvalues(dmc_df, ...)`.

- **`sep_fallback=True`, `sep_threshold=0.9` on `ep.tl.dmc`.**
  Separation-aware Fisher fallback inside `_score_finalize`: for
  sites where `|meth_diff| >= sep_threshold` AND the LR p-value
  failed to reject (p > 0.05), re-test with
  `scipy.stats.fisher_exact` on pooled counts and take the more
  powerful of the two. Never inflates p. Affects only sites the
  LR missed, so the overall FPR is unchanged.

- **`fdr_method="fdr_tsbh"` / `"fdr_storey"` on `ep.tl.dmc`.**
  Selects the FDR procedure run after the per-CpG test. New
  options: `"fdr_tsbh"` (Benjamini-Krieger-Yekutieli two-stage BH;
  statsmodels' adaptive variant), `"fdr_storey"` (Storey-Tibshirani
  q-values with π₀ estimated at lam = 0.5). Defaults to
  `"fdr_bh"` for back-compat. Threaded through
  `ep.apply_multiple_testing_correction` and the streaming
  `DMCStore` path; the manifest now records the method used.

- **`ep.dmc._storey_pi0`, `ep.dmc._apply_storey_qvalues`** helpers
  surfaced for re-use from custom pipelines.

- **Benchmark and write-up.** End-to-end reproduction of the Piao
  et al. 2021 (IJERPH 18:7975) simulated benchmark under
  `comparison_test/benchmark/`, including ground-truth
  reconstruction (35 reference DMRs / 19,999 true DMCs, both
  matching the paper's design exactly), baseline-table
  transcription, 7 figures, and a paper-style manuscript.

### Changed

- `_dmc_input_signature` now includes `sep_fallback` and
  `sep_threshold` in the SHA-256 fingerprint, so toggling them
  forces a recompute rather than serving a stale cache.

- `DMCStore.mark_bh_applied(method=...)` records the FDR method
  used so back-to-back calls with different `fdr_method` do not
  silently reuse the cached q-values.

### Tests

- New `tests/test_lr_improvements.py` (10 tests) covering Storey
  π₀ behaviour at known mixtures, fdr_method validation, the
  neighbour combiner's sign-agreement guard, the
  never-inflates-p property, and NaN preservation.

### Bug surfaced (fixed in 0.7.2)

- The pooled `fisher` backend in v0.7.0 returns `pvalue ≈ 1.0` on
  near-perfect-separation 2 × 2 tables. Fixed in 0.7.2 with a
  dual-tail hypergeometric computation.

## [Unreleased] — 0.7.0

Adds a DSS-compatible DMR caller, DSS-style raw-count smoothing for DMC,
annotatr-style multi-annotation (nearest TSS + all overlapping genes /
features), and a UCSC `refGene.txt` parser. Existing 0.6.x defaults are
unchanged — every new capability is opt-in via an explicit kwarg.

### Added

#### DMR
- **`ep.call_dmr_chain_merge` / `ep.tl.dmr(method="chain_merge")`.** DSS
  `callDMR` semantics on top of an epykit DMC table: mark sig CpGs by
  `(pvalue < alpha) AND (|meth_diff| >= min_abs_meth_diff)`, chain
  contiguous sig CpGs within `dis_merge_bp`, then filter by `min_cpgs`,
  `pct_sig`, and `minlen_bp`. `use_q_for_sig=True` switches the gate to
  the q-value column when one is present. Cached per chrom by the same
  DMC `input_sig` fingerprint as the sliding-window caller.
- **`ep.DMR_PRESETS` parameter bundles** for chain-merge: `"strict"`
  (alpha=1e-6, min_cpgs=5, |Δβ|≥0.20 — validation-ready), `"default"`
  (alpha=1e-4, |Δβ|≥0.10, min_cpgs=3 — balanced; one order looser than
  DSS to capture real-but-moderate signal without crashing PPV), and
  `"permissive"` (alpha=1e-4, dis_merge_bp=200, |Δβ|≥0.05 — recall-
  oriented). Surfaced through `ep.tl.dmr(..., preset="...")`; any
  explicit kwarg overrides the bundle value.
- **`ep.tl.diagnose_dmr_calling(md, reference_dmrs, ...)`.** Bucket
  reference DMRs into actionable categories — `SUCCESS_OVERLAP`,
  `H1_NO_CPGS` (lost in coverage / unite), `H2_NO_SIG_CPGS` (test
  too conservative), `H3a_WEAK_ALPHA` (loosen `alpha`), `H3b_STRUCTURE`
  (loosen `dis_merge_bp` / `min_cpgs`) — so a low-recall number maps to
  a specific knob instead of guesswork.

#### DMC
- **DSS-style raw-count smoothing** in `ep.tl.dmc(..., smoothing=True,
  smoothing_span_bp=500)`. Applies DSS's uniform-box ±`span/2`
  moving-average to each sample's per-CpG `(meth, cov)` before the test
  hits them — matches `DMLfit.multiFactor(smoothing=TRUE)` semantics.
  Cached separately from the un-smoothed run (`dmc/<test>_smooth/` vs
  `dmc/<test>/`) and signature-versioned so `False → True` correctly
  invalidates a stale cache.

#### Annotation
- **annotatr-style multi-annotation** in `ep.annotate_features` (default
  on, also `ep.tl.annotate(..., multi_annotation=True)`). Adds
  `nearest_tss_gene` / `nearest_tss_distance` (HOMER-style signed
  distance, flipped on `-` strand so positive is downstream of the TSS)
  and `all_overlapping_genes` / `all_overlapping_features` (one-to-many,
  so a site that's intronic for one gene AND promoter for another is
  faithfully represented instead of collapsed by the best-pick rule).
- **UCSC `refGene.txt(.gz)` annotation source.** Pass
  `ep.tl.annotate(md, refgene=...)` (or `annotate_features(..., source=
  "refgene")`) to use HOMER's default catalog — curated, protein-coding-
  biased, and gives the highest paper-gene recall on methylation work.
  Schema-compatible with the GTF path; same downstream consumers.
- **`gene_type_filter` kwarg.** Restrict the gene catalog before
  building overlap intervals (`"protein_coding"` to drop lincRNAs /
  pseudogenes). Works on both GTF (`gene_type` / `gene_biotype`) and
  refGene (derived from accession prefix).
- **GENCODE / Ensembl GTF gene-type parity.** `_parse_gtf_streaming`
  now accepts both `gene_type` (GENCODE) and `gene_biotype` (Ensembl);
  files that omit it stay annotatable.

### Deprecated
- **`ep.tl.dmc(..., use_smoothed=True)`** (the BSmooth pseudo-count
  transform) now emits a `DeprecationWarning` pointing at the new
  `smoothing=True` path. The pseudo-count path is too aggressive
  (replaces the count signal entirely with the smoothed version, washing
  out per-CpG resolution at default BSmooth parameters) and will be
  removed in a future minor release. Existing callers keep working
  for now.

### Changed
- **Internal cleanup: ASCII docstrings.** Library docstrings, log
  messages, and inline comments have been rewritten in plain ASCII
  (`β → beta`, `μ → mu`, `Σ → Sigma`, `² → ^2`, `→ → ->`, `— → --`,
  `× → x`, `≥ → >=`, `≤ → <=`, `≈ → ~=`). No semantic changes — this
  fixes `argparse --help` crashing with `UnicodeEncodeError` on Windows
  consoles running the default `cp1252` codec and removes the need for
  the CLI's stdout/stderr `utf-8` reconfigure to be load-bearing on
  every path. The CLI still reconfigures to UTF-8 defensively.

### Tests
- New: `test_dmr_chain_merge.py`, `test_dmr_presets_and_diagnose.py`,
  `test_dmc_smooth_dispersion.py`, `test_annotate_multi.py`. Existing
  tests realigned with the ASCII docstring rewrite (no behavioural
  changes).

### Deferred to 0.8+

CLI surface for `chain_merge` (`epykit dmr --method chain_merge ...`),
RefSeq / UCSC functional-element annotation tracks beyond gene models,
DMR caller for mixed designs (random-effects), Zarr storage backing,
single-cell sparse stores.

---

## [0.6.0] — HMM segmentation

One shared HMM engine, three callers on top. All three operate on the
existing methylstore (no new I/O path) and route through the 0.4
chrom-streaming dispatcher, so the distributed backend works for free.
No new dependencies — the HMM is hand-rolled in ~200 LoC of numpy.

### Added

#### Shared HMM engine
- **`epykit._hmm.segment(observations, n_states, ...)`**. Hand-rolled
  forward-backward + Viterbi for either Bernoulli or Gaussian
  emissions, with a sticky-chain transition prior (configurable
  ``self_loop`` and full ``transition_priors`` overrides).
- **`epykit._hmm.runs_of_state(viterbi, target_state, positions)`**
  extracts contiguous runs of a target state from the Viterbi path,
  optionally translated to genomic positions.

#### Callers
- **`ep.tl.pmd(md)`** — partially methylated domains. Per-sample,
  megabase-scale 2-state HMM on coverage-weighted smoothed β.
  Output in ``md.uns["pmd"]``.
- **`ep.tl.hmr(md)`** — hypo- and low-methylated regions
  (MethylSeekR-style). Per-sample 2-state HMM on raw per-CpG β.
  Tagging splits the runs into ``md.uns["hmr"]`` (dense, CpG-island-
  like) and ``md.uns["lmr"]`` (sparse, distal regulatory).
- **`ep.tl.dmr(md, method="hmm")`** — HMM-based DMR caller. Three-state
  Gaussian HMM on the per-CpG ``meth_diff`` signal from any DMC table.
  Schema-compatible with ``method="tile"``, so existing plot / export
  paths work unchanged.

### Deferred to 0.7+

Zarr storage backing, single-cell sparse stores, eQTM, motif/TFBS
enrichment, matrix-completion imputation beyond kNN, pyGenomeTracks
track plot, differential entropy CpG-window test.

---

## [0.5.0] — BAM-based read-level analyses

Adds the first analyses that need read-level methylation information.
Builds on the 0.4 engine lifts but doesn't change any existing default
behaviour. New optional `bam` extra (pysam, Linux/macOS only) gates
the new analyses; everything else stays installable on Windows.

### Added

#### BAM ingestion
- **`epykit.bam_io.read_methylation_calls(bam, ...)`**. Returns a
  long-form polars DataFrame with one row per (read, covered CpG):
  `(read_id, chrom, pos, methylation_status, base_qual, mapq,
  mate_pair_id, strand, allele_base)`. Two BAM dialects: Bismark `XM`
  tags and SAM-standard `MM`/`ML` tags (MethylDackel).
- **`bam` optional extra** in `pyproject.toml` (`pysam>=0.22`).
  Linux/macOS only — pysam has no Windows wheel.

#### Allele-specific methylation (ASM)
- **`ep.tl.asm(md, bam=..., vcf=...)`**. Per-CpG Fisher exact test of
  H1 vs H2 read methylation, with heterozygous SNVs from a VCF as
  phasing anchors. Result lands in `md.varm["asm"]` with columns
  matching the `dmc_*` family (`pvalue`, `qvalue`, `meth_diff`) so
  `pl.volcano(md, key="asm")` works without modification.
- Reuses `fisher_exact_vectorized` and `apply_multiple_testing_correction`
  from the DMC stack.

#### Methylation entropy
- **`ep.tl.entropy(md, bam=..., window_cpgs=4)`**. Per-CpG-window
  Shannon entropy over the observed read methylation patterns. Reads
  that cover all CpGs in the window contribute one binary pattern; the
  full distribution's Shannon entropy is normalised to `[0, 1]`.
  Result lands in `md.varm["entropy"]`.

### Deferred to 0.6

PMD / HMR / LMR callers and HMM-DMR all share an HMM segmentation
engine; they land together in the 0.6 release.

---

## [0.4.0] — engine lifts

Infrastructure release. Default behaviour unchanged (every 0.3 test
passes verbatim); every new capability is reached by an explicit kwarg
or optional extras install. No new analyses — those land in 0.5.

### Added

#### Compute backends
- **Distributed compute via Dask** (`tl.dmc(..., backend="dask", n_workers=4)`,
  `tl.dmr`, `tl.dvc`). Per-chromosome work submitted to a local cluster
  or an existing `dask.distributed.Client`. Requires the new
  `pip install 'epykit[distributed]'` extra. Results are bit-identical
  to the sequential path.
- **Distributed compute via Ray** (`backend="ray"`). Same surface as
  the Dask backend; uses `ray.remote` actors. Requires
  `pip install 'epykit[ray]'`.
- **GPU IRLS via CuPy** (`tl.dmc(test="glm", glm_backend="gpu")`). The
  batched binomial IRLS hot path in `_glm.py` now has a CuPy mirror
  (`_glm_gpu.py`). Requires `pip install 'epykit[gpu]'` (CUDA 12). The
  closed-form `lr` / `score` tests stay CPU-only by design.

#### Pipeline manifest + resume
- **Formal checkpoint / resume API.** Each `MethylData` analysis root
  now hosts a top-level `.epykit_manifest.json` recording completed
  pipeline stages with their input signatures and sidecar parquet
  paths. Call `ep.tl.dmc(md, ..., resumable=True)` twice with the same
  inputs and the second call loads the cached result instead of
  recomputing. `md.completed_stages` reports the recorded list;
  `md.resume_from("dmc_lr")` re-hydrates a fresh `MethylData` from the
  on-disk manifest.

#### Tabix-on-Parquet random access
- **`ep.query` module.** Three entry points —
  `query_region(store, chrom, start, end)`,
  `query_regions(store, regions_df)`, and
  `query_sites(store, sites_df)` — return long-form
  `(sample_id, chrom, pos, strand, N_meth, coverage, beta)` frames for
  arbitrary genomic loci. No new dependency: built on
  `pl.scan_parquet` predicate pushdown over the existing
  hive-partitioned store.

### Internal

- **Per-chromosome compute dispatcher** at
  `src/epykit/_compute.py:run_chrom_pipeline`. The chrom loops in
  `dmc.py`, `dmr.py` (tile), and `dvc.py` were refactored to route
  through this shared dispatcher. `backend="sequential"` (default) is
  bit-identical to the pre-0.4 in-line loop.
- **`irls_dispatch`** in `_glm.py` routes `irls_binomial_batch`
  between CPU (numpy) and GPU (CuPy via `_glm_gpu.py`).

### Deferred to 0.5+

Out of scope for 0.4: ASM, methylation entropy, PMD, HMR/LMR, HMM-DMR,
Zarr backing, single-cell sparse stores, eQTM, motif/TFBS enrichment,
matrix-completion imputation, pyGenomeTracks-style track plot.

---

## [Unreleased] — 0.3.0

### Added

#### Visualization
- **Karyogram / chromosome painter** (`ep.pl.karyogram`). One row per
  chromosome, megabase-binned mean of any per-CpG metric (`meth_diff`,
  `-log10_qvalue`, raw β). RdBu_r by default; symmetric colour limits
  for signed metrics.
- **DMR UpSet / Venn overlap** (`ep.pl.dmr_overlap`). 2-set inputs fall
  back to a 2-circle Venn; 3-6 set inputs render an UpSet plot
  (bar chart + dot matrix + per-set totals). Matplotlib-only; no
  upsetplot or matplotlib_venn dependency.
- **Gene-body metaplot** (`ep.pl.gene_body_metaplot`). Three-zone TSS /
  body / TES plot with flanking windows; body length is normalised
  across genes so a 50 kb and a 500 kb gene contribute equally.

#### Statistical features
- **DVR — Differentially Variable Regions** (`ep.tl.dvr`,
  `ep.call_dvr_density`). Region-level aggregation of `tl.dvc` output
  via per-tile DVC-density enrichment (one-sided binomial vs the
  genome-wide rate, BH-corrected). Avoids the variance-statistic
  combining problem that defeats Fisher / Stouffer's at the
  per-CpG level.
- **Effect-size shrinkage** (`ep.shrink_meth_diff`). Empirical-Bayes
  Normal-prior James-Stein-style shrinkage of `meth_diff` toward zero.
  Adds `meth_diff_shrunk`, `meth_diff_se`, `shrinkage_factor` columns
  to a DMC table. Same spirit as `ashr` / `apeglm`; pulls
  low-coverage sites harder than well-powered ones.
- **kNN β imputation** (`ep.impute_knn_beta`, `ep.impute_knn_anndata`).
  Per-chromosome inverse-distance-weighted kNN over genomic position;
  optional `max_distance_bp` cap so cross-CGI gaps don't pull
  long-distance neighbours.

#### Clocks / deconvolution scaffolding
- **Generic linear age-clock runner** (`ep.tl.age_clock`,
  `ep.age_clock`). Takes a user-supplied `(cpg_id, coefficient)` table
  and a probe → `(chrom, pos)` manifest, computes per-sample age, and
  writes the result into `md.obs`. Supports an optional `transform`
  argument (`"horvath"` for the standard ≥20-year piecewise
  anti-transform); other published clocks (Hannum, PhenoAge,
  DunedinPACE) plug in via their own coefficient CSVs. Coefficient
  tables and manifests are **not bundled** — licensing and probe
  vendor specifics differ per clock.
- **Reference-based cell-type deconvolution** (`ep.tl.deconvolve`,
  `ep.deconvolve`). Non-negative least squares solve against a
  user-supplied reference β matrix (EpiDISH / CIBERSORT / Houseman
  style). Long-format result on `md.uns['deconvolution']`; wide
  per-cell-type columns (`frac_<celltype>`) joined onto `md.obs`.

### Tests
- New: `test_viz_new.py` (8 tests), `test_dvr.py` (6 tests),
  `test_stats_new.py` (12 tests) covering shrinkage, kNN imputation
  end-to-end on AnnData, age-clock recovery on synthetic data, and
  deconvolution NNLS round-trip.

### Notes on what's *still* left on the 0.3+ roadmap

Architectural lifts that didn't fit this round and remain unimplemented:
ASM (per-read haplotype phasing from BAM), PMD / HMR / LMR callers,
methylation entropy, single-cell methylation, distributed compute
(Dask / Ray), GPU IRLS, HMM-based DMR, Zarr backing, Tabix-on-Parquet,
formal checkpoint / resume API. These each warrant their own focused
release; the karyogram + UpSet + gene-body trio plus DVR / shrinkage /
imputation / clocks / deconvolution covers the highest-ROI subset of
the roadmap.

---

## [0.2.0]

### Added
- **MethylDackel input adapter.** `ep.read_methyldackel(samplesheet, ...)` and
  `epykit convert --format methyldackel` ingest MethylDackel `.bedGraph[.gz]`
  output through the same partitioned-Parquet pipeline as Bismark. The
  conversion cache is format-aware so a Bismark store cannot be silently
  reused for MethylDackel input.
- **Permutation empirical FDR for DMC.** `ep.tl.dmc(..., empirical_fdr=True,
  n_perm=100)` shuffles treatment / control labels, re-runs the per-CpG
  DMC engine, and emits `empirical_pvalue` / `empirical_qvalue` columns —
  parity with the existing DMR permutation FDR. Refused with `formula=` /
  `contrast=` designs (label shuffling invalidates stratified models).
- **Bismark M-bias parser and plot.** `nfcore_qc.parse_bismark_mbias(path)`
  parses the Bismark M-bias text format into a long table; `pl.mbias_plot(
  {sample_id: df_or_path})` renders percent methylation per read position
  with context / R1 / R2 lines.
- **CLI `--version` flag.** `epykit --version` prints the installed
  `__version__` from `importlib.metadata`.
- **DVC engine re-export.** `ep.process_chromosomes_dvc` is now part of the
  public surface alongside the high-level `ep.tl.dvc` orchestrator.

### Changed
- **GLM degeneracy is visible.** When the batched `(X'WX)⁻¹` solve in
  `_glm.py` falls back to a per-site solve, the helper emits a
  `logger.warning` summarising affected sites instead of failing silently.
  Per-site GLM separation is logged at `info` (or `warning` when ≥5 % of
  sites separate) rather than `debug`.
- **Bisulfite conversion rate is reported, not applied.** Doc clarification
  in `qc.bisulfite_conversion_rate` and the README: epykit follows
  `bsseq` / `methylKit` defaults and does not rescale per-CpG counts by the
  conversion rate. The rate is surfaced through QC, MultiQC export, and
  the HTML report so users can gate on it.
- **Multi-group DMC accuracy test.** `tests/test_dmc_multigroup.py` now
  asserts power (≥30 %) and FDR (≤15 %) on the 3-group joint F-test
  fixture instead of only checking column presence. Continuous-covariate
  test is now an explicit structural check (engine runs, finite
  p-values, FDR isn't catastrophic) with the fixture limitation
  documented in the test.
- **DMR sensitivity / specificity tests.** `test_accuracy.py` gains
  `test_dmr_tile_sensitivity_and_fdp` and
  `test_dmr_sliding_window_sensitivity_and_fdp`, pinning both recovery
  rate and false-discovery proportion per method. Conditional
  `pytest.skip` calls in the existing DMR direction test are now hard
  assertions to surface regressions instead of masking them.

### Fixed
- `tests/test_dmc_multigroup.py` no longer silently swallows a `ValueError`
  in the contrast-resolution test; on `ValueError` the test now asserts
  the error message is informative and on success it asserts finite
  p-values are produced.

## [0.1.0]

- Initial release. Bismark `.cov` → partitioned Parquet methylstore;
  scanpy-style `pp` / `tl` / `pl` API; 8 DMC test backends; two DMR
  engines plus permutation FDR; DVC (iEVORA-style); covariate-aware
  contrasts; AnnData / MuData / methylKit / MultiQC interop;
  self-contained HTML report.
