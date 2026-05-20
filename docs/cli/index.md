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
| `--test {lr,score,glm,logit_t,welch_t,bb_lr,cmh,fisher}` | Statistical test to use |
| `--formula TEXT` | Model formula (e.g., `"~ treatment"`, `"~ treatment + age"`) |
| `--contrast COL TREAT CTRL` | Contrast specification: column name, treatment level, control level |

---

## dmr

Call differentially methylated regions from DMC results.

```bash
epykit dmr \
    --methylstore methylstore/ \
    --method tile \
    --empirical-fdr \
    --n-perm 1000
```

| Option | Description |
|--------|-------------|
| `--method {tile,sliding_window}` | DMR detection method |
| `--empirical-fdr` | Use permutation-based empirical FDR |
| `--n-perm INT` | Number of permutations for empirical FDR (default: 1000) |

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

---

## qc-report

Generate a QC and coverage report.

```bash
epykit qc-report \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    -o qc_report.html
```

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
epykit dmc \
    --methylstore methylstore/ \
    --samplesheet samplesheet.csv \
    --test lr \
    --formula "~ treatment" \
    --contrast treatment tumor normal

# 4. DMR detection
epykit dmr \
    --methylstore methylstore/ \
    --method tile

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

The following features are currently available only through the Python API and do not
have CLI equivalents:

- `lr` test with `power_stack` options
- `chain_merge` DMR method
- AnnData / MuData export with custom modalities
- `ep.query` random-access queries

See the [Python API documentation](../export/index.md) for these features.
