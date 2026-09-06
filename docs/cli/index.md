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
| `--samplesheet PATH` | Path to the samplesheet CSV |
| `--treatment GROUP` | Treatment group name |
| `--control GROUP` | Control group name |
| `--assembly NAME` | Genome assembly (e.g., `hg38`, `mm10`) |

## Sibling TSV auto-emit

The `dmc`, `dmr`, `annotate`, and `qc-report` subcommands write a sibling
TSV next to their primary output by default. For example,
`epykit dmc ... --output dmc.parquet` also writes
`dmc.significant.tsv` (and, with `--csv-full`, `dmc.tsv`).

| Flag / env var | Effect |
|----------------|--------|
| `--csv PATH` | Override the auto-derived path. A `.csv` suffix selects comma delimiter; anything else uses tab. Implies the file is written. |
| `--no-csv` | Suppress the sibling TSV auto-emit entirely. |
| `--csv-alpha FLOAT` *(`dmc`, `dvc`)* | q-value threshold for the significant-only TSV. Default `0.05`. |
| `--csv-full` *(`dmc`, `dvc`)* | Also write the full (unfiltered) table alongside the significant one. |
| `EPYKIT_NO_AUTO_CSV=1` *(env var)* | Suppress the auto-emit globally for the current shell session -- handy for batch scripts that only want parquet output. |

---

## convert

Convert alignment tool output into the epykit Parquet store format.

```bash
epykit convert \
    --samplesheet samplesheet.csv \
    --methylstore methylstore/ \
    --format bismark
```

| Option | Description |
|--------|-------------|
| `--format {bismark,methyldackel}` | Input format |
| `--samplesheet PATH` | CSV with columns: `sample_id`, `file`, `treatment` |
| `--methylstore PATH` | Output Parquet store directory |

---

## filter

Apply coverage and blacklist filters to the store.

```bash
epykit filter \
    --methylstore methylstore/ \
    --min-coverage 10 \
    --blacklist encode_blacklist.bed
```

| Option | Description |
|--------|-------------|
| `--min-coverage INT` | Minimum read coverage per CpG (default: 10) |
| `--max-coverage INT` | Maximum read coverage per CpG |
| `--blacklist PATH` | BED file of regions to exclude |

---

## summary

Print per-sample summary statistics.

```bash
epykit summary --methylstore methylstore/ --samplesheet samplesheet.csv
```

---

## dmc

Run per-CpG differential methylation calling.

```bash
epykit dmc \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --test lr \
    --formula "~ treatment" \
    --contrast treatment tumor normal
```

| Option | Description |
|--------|-------------|
| `--test {lr,glm,welch_t,fisher}` | Statistical test to use. Default `lr` (quasi-binomial LR). `glm` is auto-selected when `--formula` is set. `fisher` is the n=1 fallback. The pre-1.0 engines `logit_t`, `bb_lr`, `score`, `cmh` were removed in 0.7.5. |
| `--formula TEXT` | Model formula (e.g., `"~ treatment"`, `"~ treatment + age"`) |
| `--contrast COL TREAT CTRL` | Contrast specification: column name, treatment level, control level |
| `--allow-n1` | Permit n=1 per group (falls back to Fisher exact on pooled reads). Off by default -- p-values become anti-conservative. |
| `--no-csv` / `--csv PATH` / `--csv-alpha 0.05` / `--csv-full` | TSV auto-emit controls (see [Sibling TSV auto-emit](#sibling-tsv-auto-emit)). |

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

# Tile-based DMR (reads the methylstore directly).
epykit dmr \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --treatment-group tumor --control-group normal \
    --method tile \
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
| `--no-csv` / `--csv PATH` | TSV auto-emit controls. |

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
| `--no-csv` / `--csv PATH` | TSV auto-emit controls. |

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
| `--no-csv` / `--csv PATH` | Suppress or override the sibling TSV auto-emit. |

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
# 1. Convert Bismark output to Parquet store
epykit convert \
    --samplesheet samplesheet.csv \
    --methylstore methylstore/ \
    --format bismark

# 2. Filter low-coverage CpGs
epykit filter \
    --methylstore methylstore/ \
    --min-coverage 10

# 3. Differential methylation calling
#    Auto-emits dmc.significant.tsv next to dmc.parquet.
epykit dmc \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --test lr \
    --formula "~ treatment" \
    --contrast treatment tumor normal \
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

- `tl.dmc(power_stack=...)` and the `neighbour_combine` / `sep_fallback` knobs of the lr+ stack. The CLI runs bare `lr` (`--fdr-method` and `--dispersion` are exposed); the lr+ stack is Python API only and no CLI flags are planned.
- AnnData / MuData export with custom modalities
- `ep.query` random-access queries

`chain_merge` DMRs (previously API-only) ship in the CLI since 1.0 via
`epykit dmr --method chain_merge --dmc-results <dmc.parquet>`.

See the [Python API documentation](../export/index.md) for these features.
