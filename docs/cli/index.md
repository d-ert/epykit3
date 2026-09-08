# CLI Reference

The `epykit` command-line interface provides access to the core analysis pipeline
without writing Python code. It is installed as a console script with the package.

```bash
epykit --version
epykit <subcommand> [options]
```

---

## Subcommands

| Subcommand | Purpose |
|------------|---------|
| `convert` | Convert Bismark or MethylDackel output to Parquet store |
| `filter` | Coverage and blacklist filtering |
| `summary` | Per-sample summary statistics |
| `dmc` | Per-CpG differential methylation calling |
| `dmr` | Differentially methylated region calling |
| `annotate` | Gene-feature and CpG-island annotation |
| `qc-report` | QC and coverage report |
| `smooth` | Gaussian-kernel smoothing |
| `report` | Interactive HTML report |
| `aggregate-regions` | Aggregate methylation to BED regions |
| `export` | Export to external formats (sub-commands below) |

---

## Common Flags

These flags are shared across most subcommands:

| Flag | Description |
|------|-------------|
| `--methylstore PATH` | Path to the Parquet methylation store |
| `--samplesheet PATH` | Path to the samplesheet CSV (`sample_id`, `group`, `path`) |
| `--treatment-group GROUP` | Treatment group name |
| `--control-group GROUP` | Control group name |
| `--canonical-only` | Restrict `convert`, `dmc` and `dmr --method tile` to the fixed human-style chromosome set; see each subcommand |

## Sibling TSV auto-emit

The `dmc`, `dmr`, `annotate`, and `qc-report` subcommands write a sibling
TSV next to their primary output by default. For example,
`epykit dmc ... --output dmc.parquet` also writes
`dmc.significant.tsv` (and, with `--tsv-full`, `dmc.tsv`).

| Flag / env var | Effect |
|----------------|--------|
| `--tsv PATH` | Override the auto-derived path. A `.csv` suffix selects comma delimiter; anything else uses tab. Implies the file is written. |
| `--no-tsv` | Suppress the sibling TSV auto-emit entirely. |
| `--tsv-alpha FLOAT` *(`dmc`)* | q-value threshold for the significant-only TSV. Default `0.05`. |
| `--tsv-full` *(`dmc`)* | Also write the full (unfiltered) table alongside the significant one. |
| `EPYKIT_NO_AUTO_TSV=1` *(env var)* | Suppress the auto-emit globally for the current shell session, handy for batch scripts that only want parquet output. |

The older `--csv`, `--no-csv`, `--csv-alpha` and `--csv-full` flags and the
`EPYKIT_NO_AUTO_CSV` variable are deprecated aliases of the `tsv` names. They
still work through the same code path, are hidden from `--help`, and log a
deprecation warning when used. The name never selects the delimiter; the
path suffix does. See [Deprecations](../reference/deprecations.md).

---

## convert

Convert one Bismark `.cov[.gz]` or MethylDackel `.bedGraph[.gz]` file into
the epykit Parquet store format. Run it once per sample; `ep.read_bismark()`
in the Python API does the same for a whole samplesheet.

```bash
epykit convert \
    --input tumor_1.cov.gz \
    --sample-id tumor_1 \
    --output-dir methylstore/ \
    --format bismark
```

| Option | Description |
|--------|-------------|
| `--input PATH` | Required. The `.cov[.gz]` or `.bedGraph[.gz]` file. |
| `--sample-id ID` | Required. Sample identifier written into the `sample` partition. |
| `--output-dir PATH` | Required. Output Parquet store directory. |
| `--format {bismark,methyldackel}` | Input format. Default `bismark`. |
| `--context {CpG,CHG,CHH}` | Cytosine context label. Default `CpG`. |
| `--reference-fasta PATH` | Indexed reference FASTA used to infer strand. |
| `--merge-cpg` / `--no-merge-cpg` | Merge symmetric CpG dyads into one record (the default) or keep per-strand records. |
| `--canonical-only` | Keep only the fixed human-style chromosome set (`1`-`22`, `X`, `Y`, `M`/`MT`, with or without `chr`) and drop every other contig before the partition write. The setting is part of the per-sample conversion cache. Default off. Same as `convert_sample(canonical_only=True)`; see [Canonical chromosomes only](../io/read-bismark.md#canonical-chromosomes-only). |

---

## filter

Apply coverage and blacklist filters to the store.

```bash
epykit filter \
    --methylstore methylstore/ \
    --output-dir methylstore_filtered/ \
    --min-coverage 10 \
    --blacklist-bed encode_blacklist.bed
```

| Option | Description |
|--------|-------------|
| `--output-dir PATH` | Required. Directory for the filtered store. |
| `--min-coverage INT` | Minimum read coverage per CpG (default: 10) |
| `--max-coverage-quantile FLOAT` | Drop CpGs above this coverage quantile (default: 0.999) |
| `--blacklist-bed PATH` | BED file of regions to exclude |
| `--sample ID` | Filter one sample only |

---

## summary

Print per-sample summary statistics.

```bash
epykit summary --methylstore methylstore/ --sample tumor_1
```

---

## dmc

Run per-CpG differential methylation calling. The binary treatment /
control path streams the engine over the methylstore; `--formula` or
`--contrast` switches to the GLM-contrast path, which reads every sample in
the samplesheet.

```bash
# Binary treatment / control path, the lr engine.
epykit dmc \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --treatment-group tumor --control-group normal \
    --test lr \
    --output dmc.parquet

# Covariate-adjusted GLM contrast on the samplesheet columns.
epykit dmc \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --treatment-group tumor --control-group normal \
    --formula "~ group + age" --contrast group \
    --output dmc_glm.parquet
```

| Option | Description |
|--------|-------------|
| `--output PATH` | Required. DMC parquet output path. |
| `--test {lr,glm,welch_t,fisher}` | Statistical test to use. Default `lr` (quasi-binomial LR). `glm` needs a design via `--formula`. `fisher` is the n=1 fallback. The pre-1.0 engines `logit_t`, `bb_lr`, `score`, `cmh` were removed in 0.7.5. |
| `--formula TEXT` | patsy formula on samplesheet columns (e.g. `"~ group"`, `"~ group + age"`). Selects the GLM-contrast path; pair with `--contrast`. |
| `--contrast SPEC` | A column name (continuous covariate effect), a factor name for a joint F-test (`group`), or a patsy linear combination (`'group[T.KO] - group[T.WT]'`). |
| `--covariates a,b` | Comma-separated nuisance covariate columns from the samplesheet. |
| `--unite` / `--no-unite` | Intersect (sites covered in every sample) or union (the default, sites covered in at least one sample). |
| `--min-samples-treatment N` / `--min-samples-control N` | Per-site minimum number of samples with coverage in each group. Default 0. |
| `--allow-n1` | Permit n=1 per group (falls back to Fisher exact on pooled reads). Off by default; p-values become anti-conservative. |
| `--dispersion {site,eb,shrink,chrom}` | Dispersion estimator for `lr`. Default `eb`, matching `ep.tl.dmc`. |
| `--reference {adaptive,F,chi2}` | Reference distribution for the `lr` statistic. Default `adaptive`. |
| `--fdr-method NAME` | Multiple-testing correction. Default `fdr_bh`. |
| `--smoothing` | DSS-style per-sample count smoothing for `--test lr`: each sample's raw counts are replaced by a uniform-box average over the CpGs within half the span on each side before the test, as in `DMLfit.multiFactor(smoothing=TRUE)`. Default off. Rejected with the other engines (including the `--allow-n1` Fisher fallback) and on the `--formula` / `--contrast` path, which do not read it. Same as `ep.tl.dmc(smoothing=True)`. |
| `--smoothing-span-bp INT` | Full smoothing window in bp for `--smoothing`. Default `500`, the DSS default. Must be positive while `--smoothing` is set. |
| `--canonical-only` | Test only the fixed human-style chromosome set (`1`-`22`, `X`, `Y`, `M`/`MT`, with or without `chr`) of the store's partitions; other contigs are dropped before the test and the FDR correction, on the binary and the `--formula` / `--contrast` path. Default off. Same as `ep.tl.dmc(canonical_only=True)`; see [DMC calling](../analysis/dmc.md#canonical-chromosomes-only). |
| `--no-tsv` / `--tsv PATH` / `--tsv-alpha 0.05` / `--tsv-full` | TSV auto-emit controls (see [Sibling TSV auto-emit](#sibling-tsv-auto-emit)). |

---

## dmr

Call differentially methylated regions. `chain_merge` (the default since 1.0)
chains contiguous significant CpGs from a precomputed DMC parquet into
DSS-style regions. `tile` pools reads across CpGs within each fixed-size
tile and tests directly from the methylstore (no DMC parquet required).
`sliding_window` and `segment` both read a precomputed DMC parquet.

```bash
# Default: DSS-style chain_merge from a precomputed DMC parquet.
epykit dmr \
    --method chain_merge \
    --dmc-results md_out/varm/dmc_lr.parquet \
    --preset strict \
    --output dmrs.parquet

# Tile-based DMR (reads the methylstore directly), canonical chromosomes only.
epykit dmr \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --treatment-group tumor --control-group normal \
    --method tile \
    --canonical-only \
    --empirical-fdr --n-perm 1000 \
    --output dmrs.parquet
```

| Option | Description |
|--------|-------------|
| `--method {chain_merge,tile,sliding_window,segment}` | DMR detection method. Default `chain_merge`. `hmm` was renamed `segment` in 1.0 (same engine). |
| `--output PATH` | Required. DMR parquet output path. |
| `--dmc-results PATH` | (`chain_merge`, `sliding_window`, `segment`) Parquet from `epykit dmc`. |
| `--preset {strict,default,permissive}` | (`chain_merge`) Parameter bundle. Explicit knob flags override the bundled value. |
| `--dis-merge-bp INT` | (`chain_merge`) Max bp gap between consecutive significant CpGs. Default 500. |
| `--pct-sig FLOAT` | (`chain_merge`) Min fraction of CpGs in a span that must be significant. Default 0.5. |
| `--minlen-bp INT` | (`chain_merge`) Min DMR span in bp. Default 50. |
| `--use-q-for-sig` | (`chain_merge`) Gate significance on q-value rather than p-value. |
| `--tile-size-bp INT` | (`tile`) Tile width. Default 1000. |
| `--min-cpgs-per-tile INT` | (`tile`) Minimum CpGs per tile per sample. Default 5. |
| `--test {lr,glm,welch_t,fisher}` | (`tile`) Statistical test applied to tile-level counts. Default `lr`. |
| `--empirical-fdr` | (`tile`) Permutation-based empirical FDR. |
| `--n-perm INT` | (`tile`) Number of permutations. Default 100. |
| `--perm-seed INT` | (`tile`) Seed for the label shuffles. Default 42. |
| `--canonical-only` | (`tile`) Test only the fixed human-style chromosome set (`1`-`22`, `X`, `Y`, `M`/`MT`, with or without `chr`) of the store's partitions; the same set is used by every `--empirical-fdr` permutation. Default off. The other methods inherit the chromosomes of the DMC parquet and exit with an error: run `epykit dmc --canonical-only` instead. Same as `ep.tl.dmr(method="tile", canonical_only=True)`. |
| `--min-mean-qvalue FLOAT` | (`chain_merge`, `sliding_window`, `tile`) Region-level q-value post-filter, matching `ep.tl.dmr`. Default 0.05; a value above 1.0 disables it. |
| `--no-tsv` / `--tsv PATH` | TSV auto-emit controls. |

---

## annotate

Annotate CpGs with gene features and CpG island context.

```bash
epykit annotate \
    --methylstore methylstore/ \
    --gtf gencode.v44.gtf.gz \
    --cpg-islands cpg_islands.bed
```

| Option | Description |
|--------|-------------|
| `--gtf PATH` | GTF file for gene-feature annotation (promoter, exon, intron, intergenic) |
| `--cpg-islands PATH` | BED file for CpG island/shore/shelf annotation |
| `--no-tsv` / `--tsv PATH` | TSV auto-emit controls. |

---

## qc-report

Generate a QC and coverage report.

```bash
epykit qc-report \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    -o qc_report.html
```

| Option | Description |
|--------|-------------|
| `--no-tsv` | Suppress the sibling TSV auto-emit. |

---

## smooth

Apply Gaussian-kernel smoothing to methylation values.

```bash
epykit smooth --methylstore methylstore/ --bandwidth 300
```

| Option | Description |
|--------|-------------|
| `--bandwidth INT` | Smoothing bandwidth in base pairs (default: 300) |

---

## report

Generate an interactive HTML report.

```bash
epykit report \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    -o report.html
```

---

## aggregate-regions

Aggregate per-CpG methylation values to user-defined regions.

```bash
epykit aggregate-regions \
    --methylstore methylstore/ \
    --regions promoters.bed \
    -o region_means.tsv
```

| Option | Description |
|--------|-------------|
| `--regions PATH` | BED file defining regions to aggregate over |
| `-o PATH` | Output TSV file |

---

## export

Export data to external formats. Each export type is a sub-command:

```bash
epykit export bedgraph --methylstore methylstore/ --sample tumor_1 -o tumor_1.bedgraph
epykit export bigwig --methylstore methylstore/ --sample tumor_1 -o tumor_1.bw
epykit export dmcs-bed --methylstore methylstore/ -o dmcs.bed
epykit export dmrs-bed --methylstore methylstore/ -o dmrs.bed
epykit export mudata --methylstore methylstore/ --samplesheet samplesheet.csv -o data.h5mu
epykit export methylkit-tabix --methylstore methylstore/ --samplesheet samplesheet.csv -o methylkit_out/
epykit export multiqc --methylstore methylstore/ --samplesheet samplesheet.csv -o multiqc_input/
```

---

## Example Pipeline

A complete CLI pipeline from raw Bismark output to an annotated HTML report:

```bash
# 1. Convert each sample's Bismark output to the Parquet store,
#    keeping the canonical chromosomes only.
for sample in tumor_1 tumor_2 normal_1 normal_2; do
    epykit convert \
        --input "$sample.cov.gz" \
        --sample-id "$sample" \
        --output-dir methylstore/ \
        --format bismark \
        --canonical-only
done

# 2. Filter low-coverage CpGs
epykit filter \
    --methylstore methylstore/ \
    --output-dir methylstore_filtered/ \
    --min-coverage 10

# 3. Differential methylation calling
#    Auto-emits dmc.significant.tsv next to dmc.parquet.
epykit dmc \
    --methylstore methylstore_filtered/ \
    --samplesheet samplesheet.csv \
    --treatment-group tumor --control-group normal \
    --test lr \
    --output dmc.parquet

# 4. DMR detection (DSS-style chain_merge over the DMC parquet)
epykit dmr \
    --method chain_merge \
    --dmc-results dmc.parquet \
    --preset default \
    --output dmrs.parquet

# 5. Annotate with gene features and CpG islands
epykit annotate \
    --methylstore methylstore/ \
    --gtf gencode.v44.gtf.gz \
    --cpg-islands cpg_islands.bed

# 6. Generate HTML report
epykit report \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    -o methylation_report.html
```

---

## Python-API-Only Features

The following features are currently available only through the Python API
and do not have CLI equivalents:

- `tl.dmc(power_stack=...)` and the lr+ knobs `neighbour_combine`, `neighbour_bp`, `sep_fallback` and `sep_threshold`. The CLI runs bare `lr` (with `--dispersion`, `--reference`, `--fdr-method` and `--smoothing`); there is no schedule for lr+ flags.
- `tl.dmc(resumable=True)`, `materialize=False` and the `backend` / `n_workers` execution knobs.
- AnnData / MuData export with custom modalities
- `ep.query` random-access queries

`chain_merge` DMRs (previously API-only) ship in the CLI since 1.0 via
`epykit dmr --method chain_merge --dmc-results <dmc.parquet>`.

See the [Python API documentation](../export/index.md) for these features.
