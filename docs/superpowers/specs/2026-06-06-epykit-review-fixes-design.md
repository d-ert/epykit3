# epykit review-remediation — design

**Status:** approved for Tier 1 execution
**Date:** 2026-06-06
**Source of findings:** [docs/review/2026-06-06-epykit-peer-review.md](../../review/2026-06-06-epykit-peer-review.md) (8-agent pre-submission audit; all Critical claims hand-verified)
**Decision record:** phased remediation; "fix the science properly" posture; C1 via start/end auto-detection; C4/C5 deferred to Tier 2.

---

## Goal

Remediate the defects found in the pre-submission review so that epykit 1.0 is correct for end users and reproducible for reviewers, without disturbing the benchmark paper's published numbers. Work proceeds in three approve-and-execute tiers; this document specifies **Tier 1 in full** and enumerates Tiers 2–3.

## Non-goals / out of scope

- No change to the `combined_strand_bed` or MethylDackel ingestion coordinate handling (both genuinely 0-based; the paper's GSE263850 path uses `combined_strand_bed`).
- No re-running of the benchmark ablations in Tier 1 (no `lr+` knob-default changes in Tier 1).
- No new features. This is remediation only.
- Putting a bug-fix manifest into the manuscript is explicitly out of scope — fixes live in the repo; the paper presents the current stable version.

## Guiding constraints

- Windows + Linux, py3.9–3.12 must stay green (CI matrix is load-bearing).
- Library code emits via `logging`, never `print` (CLI owns stdout).
- Preserve the streaming contract (O(largest chromosome) peak memory); don't materialize full per-CpG frames.
- Every fix lands with a test that would have caught the bug. Where the existing test encodes the bug (C1 fixtures, C3 transform test), the test is corrected as part of the fix.

---

## Tier 1 — paper-blockers + silent-wrong-science criticals

Tier 1 acceptance: `uv run pytest -m "not slow"` green on both OSes; the new slow-tier job green on Linux; CLI and API produce identical DMC q-values on the parity fixture; a fresh `uv sync` from the committed lock yields a working `import epykit` + `pivot` path; real Bismark `.cov` round-trips to the correct 0-based `pos`.

### 1.1 — Dependency floor correction (M-PKG1)
- **Change:** `pyproject.toml` — `polars>=0.20.0` → `polars>=1.0,<2`. Re-derive every other floor against a clean resolve; in particular reconcile `numba`/`numpy` (numba≥0.60 for numpy-2) and relax the auto-pinned `psutil>=7.2.2` / `matplotlib>=3.9.4` / `seaborn>=0.13.2` / `scikit-learn>=1.6.1` to genuinely-tested minimums.
- **Files:** `pyproject.toml`.
- **Test:** a test that exercises a `DataFrame.pivot(on=...)` path (e.g. via `pl/_compute.py`) so a wrong floor can't pass silently; CI `uv sync` from the committed lock then `python -c "import epykit"`.
- **Acceptance:** fresh resolve installs polars≥1.0; the pivot smoke test passes.

### 1.2 — CLI/API q-value parity (M-PKG2)
- **Change:** `_cmd_dmc` (`cli.py`) currently calls `process_chromosomes_dmc` with no `dispersion`, inheriting that function's `dispersion="site"` default while `tl.dmc` uses `"eb"`. Make the CLI default match `tl.dmc`: pass `dispersion` (default `"eb"`), `reference` (default `"adaptive"`), and `fdr_method` (default `"fdr_bh"`), exposed as `--dispersion`, `--reference`, `--fdr-method`. Update the README CLI table accordingly (also fixes m-PKG6 for the dmc row).
- **Files:** `cli.py`, `README.md`.
- **Test:** new `tests/test_cli_api_parity.py` — run `epykit dmc` (via the CLI entry) and `tl.dmc` with matched args on one synthetic store; assert `qvalue` columns equal within 1e-12.
- **Acceptance:** parity test green; `epykit dmc --help` lists the new flags.

### 1.3 — Commit the environment lock (M-PKG4)
- **Change:** remove the `uv.lock` ignore line in `.gitignore`; commit the existing `uv.lock`. Document in `benchmark/README.md` (or the paper's repro section) the thread env-vars (`POLARS_MAX_THREADS`, `OMP_NUM_THREADS`) used to produce the published numbers.
- **Files:** `.gitignore`, `uv.lock` (tracked), `benchmark/README.md`.
- **Acceptance:** `uv.lock` tracked; a documented `uv sync --frozen` path exists.

### 1.4 — Slow-tier CI + real null-calibration test (M-PKG3)
- **Change A (CI):** add a job to `.github/workflows/test.yml` that runs `uv run pytest -m slow` on `ubuntu-latest` / one Python (3.12). May be `if: github.event_name == 'push'` or nightly to bound cost; must run somewhere automated.
- **Change B (test):** new `tests/test_null_calibration.py` (slow-marked) that builds a true-null store (both groups drawn from one beta-binomial at realistic φ≈2; or label-permuted from a real-ish fixture), runs `ep.tl.dmc(test=...)` for `lr`, `glm`, `welch_t`, and asserts: (i) KS test of raw p-values vs Uniform(0,1) is **not** rejected at α=0.01; (ii) empirical FPR at raw p<0.05 ∈ [0.03, 0.07]. This is the test that anchors the paper's calibration claim and currently does not exist (`test_calibration.py` exercises mock engines only).
- **Files:** `.github/workflows/test.yml`, `tests/test_null_calibration.py`.
- **Acceptance:** the calibration test passes for all three engines on the new CI job. *If an engine fails calibration here, that is a real finding* — surface it and rescope (do not loosen the test to pass).

### 1.5 — C1: Bismark coordinate correction (auto-detect via start/end)
- **Change:**
  1. In `convert.py`, the `else` (bismark/methyldackel) branch: stop discarding `end`. Add a `coordinate_base` resolution step keyed on `format` and an explicit override:
     - new param `coordinate_base: Literal["auto","one_based","zero_based"] = "auto"` on `convert`/`read_bismark`.
     - `"auto"` for `format="bismark"`: sample the first N rows; if `start == end` predominates → 1-based → `pos = start - 1`; if `start == end - 1` → 0-based → `pos = start`. Emit one `logger.info` stating the detected convention.
     - `format="methyldackel"` → always 0-based (`pos = start`); `combined_strand_bed` branch unchanged.
     - explicit `one_based`/`zero_based` override bypasses detection.
  2. Rewrite the `convert.py:3-9` module docstring to state the truth: standard Bismark `.cov` is 1-based (`start==end`); MethylDackel bedGraph and combined-strand BED are 0-based; auto-detection + override.
  3. Bump the raw-store manifest version (`RAW_MANIFEST_NAME` payload) so stores written under the old (shifted) convention are detectable; on load of an old-version raw store, warn that positions may be +1 shifted.
- **Files:** `convert.py`, `io.py` (thread the param through `read_bismark`), `_cache.py`/manifest writer.
- **Tests:**
  - rewrite `tests/fixtures/synth.py` Bismark writer to emit real 1-based `.cov` (`start == end`); keep/clearly-name a 0-based variant for the MethylDackel-equivalence test.
  - rewrite `tests/test_methyldackel.py` so the equivalence test feeds a Bismark row at `start=N` and a MethylDackel row at `start=N-1` for the *same* cytosine and asserts identical `pos`.
  - new end-to-end regression: a `.cov` row at a known 1-based position → `read_bismark` → assert `pos` is the 0-based coordinate → `annotate_features` against a known feature → assert the expected `feature_type`/distance (pins the 0/1-based + strand chain).
- **Acceptance:** real 1-based `.cov` ingests to correct 0-based `pos`; pre-shifted 0-based input is detected and not double-shifted; cross-format unite matches; the coord→annotation regression pins the convention.
- **Migration note:** existing methylstores built from real Bismark `.cov` under the old behavior were +1 shifted; the manifest-version bump + load warning flags them. Document "re-`read_bismark` to rebuild" in CHANGELOG.

### 1.6 — C2: Annotation chromosome-name mismatch guard
- **Change:** in `annotate.py`, after building `features_by_chrom` / `tss_by_chrom` and after loading the CpG-island BED, compute `set(site_chroms) ∩ set(feature_chroms)`. If empty → raise `ValueError` showing sampled names from each side (`sites: ['chr1','chr2'] vs features: ['1','2']`). If non-empty but below a fraction (e.g. <50% of site chroms covered) → `logger.warning`. Applies to both the gene-feature and CpG-island joins. Optionally accept `harmonize_chr=True` to add/strip the `chr` prefix on one side.
- **Files:** `annotate.py`.
- **Test:** annotate a `chr1`-named store against a `1`-named GTF; assert it raises (and that `harmonize_chr=True`, if implemented, succeeds).
- **Acceptance:** mismatch no longer silently yields 100% intergenic.

### 1.7 — C3: Horvath clock anti-transform
- **Change:** `clocks.py:238` negative branch → `21.0*np.exp(linear) - 1.0`; parameterize `adult_age: float = 20.0` and compute both branches from it (`(1+adult_age)*exp(x)-1` and `(1+adult_age)*x+adult_age`).
- **Files:** `clocks.py`.
- **Test:** rewrite `tests/test_stats_new.py:230-240` to call `age_clock(..., transform="horvath")` and assert against a hand-computed reference (e.g. `linear=-0.5 → 11.74`, `linear=+0.1 → 22.1`).
- **Acceptance:** young/under-20 ages are correct; the test exercises the real function.

---

## Tier 2 — remaining correctness (enumerated; detailed at tier start)

"Fix the science properly" posture; C4/C5 first.

- **C4 — ASM bisulfite confound** (`asm.py:167-176`): restrict phasing anchors to bisulfite-safe SNV classes (exclude C/T, A/G, and strand-context-confounded types) or document + gate behind `epykit.experimental`; apply the same dup/QC/base-quality filters as `read_methylation_calls`; anchor each CpG to a single het SNV rather than a global ref/alt label.
- **C5 — `qc.power`** (`qc.py:721-793`): non-central t (or t critical values, df=2(n−1)); add overdispersion φ multiplying the binomial variance; accept `n_tests`/target-FDR to derive a multiple-testing-aware α; document single-locus default.
- **Statistical engines:** M-STAT1 (GLM CI via `delta_method_meth_diff_ci` + φ-scaled `coef_se`), M-STAT2 (emit the F denominator df actually used), M-STAT3 (share the `DF_PHI_FLOOR` convention across `lr`/`glm`), M-STAT4 (BH/Storey over finite p-values only), M-STAT5 (implement documented pooled-null empirical FDR or rename to FWER).
- **DMR:** M-DMR1 (combine raw p, not q — verified), M-DMR2 (correlation-aware combine / honest docstring), M-DMR3 (within-stratum permutation), M-DMR4 ("DSS-style" relabel + document missing smoothing), M-DMR5 (FDR over all candidate regions; no double-BH), M-DMR6 (signed combine + genome-wide BH in `dmr_segment`), M-DMR7 (standardize DMR width convention + document 0-based half-open), M-DMR8 (`_MAX_DMR_BP` configurable), M-DMR9 (PMD NaN-safe smoother).
- **DVC:** M-DVC1 (warn/NaN at n<3 instead of silent zero), M-DVC2 (rename to `brown_forsythe`, warn on the `bartlett` alias, fix docstring — verified), M-DVC3 (`min_coverage` floor).
- **QC:** M-QC1 (sex_check guards + warn when diptest absent), M-QC2 (rename `contamination_estimate`→`intermediate_beta_fraction` or gate experimental + document the heterogeneity confound; a validated estimator needs spike-in data and is out of scope), M-QC3 (pairwise-complete `sample_correlation` + report site counts).
- **Secondary modules:** M-SEC1 (clocks `impute_missing=False` true drop, not zero-fill), M-SEC2 (gold-standard-mean imputation), M-SEC3 (impute returns a `was_imputed` mask + docs), M-SEC4 (entropy Miller-Madow bias correction + matched-coverage doc).
- **Annotation:** M-ANN1 (intron via cumulative-max end), M-ANN2 (stable deterministic gene tie-break), M-ANN3 (vectorize nearest-TSS), M-ANN4 (per-transcript TSS or document gene-level limitation).
- **API/packaging:** M-PKG5 (move clocks/asm/entropy/impute/pmd/hmr behind `epykit.experimental.*` or warn-on-use; stop exporting `_glm.build_design`).

## Tier 3 — minor/nitpick polish (enumerated)

All Minor + Nitpick IDs from the review (m-STAT6–11, m-IO1–6, m-DVC4–5, m-QC4–6, m-ANN5, m-DMR10–11, m-PKG6–9, and the nitpick list): CI/lint ratchet (`select=["F"]`→`F,E,W,I,B`), `python -m epykit` entry, CHANGELOG ordering, FutureWarning noise + version reconciliation, docs-URL fix, save/load symmetry + self-containment, samplesheet covariate dtype, region_beta half-open + multi-part glob, and the rest.

---

## Risks & rollback

- **C1 changes outputs for real-Bismark users.** Mitigated by: detection (pre-shifters unaffected), explicit override, manifest-version bump + load warning, CHANGELOG migration note. Rollback = revert the offset; stores remain readable.
- **Null-calibration test may fail a real engine.** That is a discovery, not a blocker to loosen — it would reshape Tier 2. Surface immediately.
- **polars≥1.0 floor** could surface latent API usage elsewhere; the pivot smoke test + full `not slow` run on both OSes covers it.

## Verification per tier

Each fix = atomic commit with its test. Tier gate = full `not slow` suite green on Windows+Linux, plus the new slow job green, plus the tier's specific acceptance criteria above. No tier advances until its gate passes.
