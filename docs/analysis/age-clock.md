# Epigenetic Age Clocks

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.age_clock(md)` runs linear epigenetic-age clock models and stores a
per-sample age estimate in `md.obs["epigenetic_age"]`. The user supplies the
clock coefficient table and an Illumina CpG manifest that maps CpG probe IDs
to genomic coordinates.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

# Run the Horvath multi-tissue clock
ep.tl.age_clock(
    md,
    coefficients="horvath_coefficients.csv",
    manifest="HM450.hg38.manifest.tsv.gz",
    clock_name="horvath",
)

# Per-sample predicted age is now on md.obs
print(md.obs.select("sample_id", "epigenetic_age"))
```

## How It Works

1. The coefficient table is loaded. Each row maps a CpG probe ID to a weight.
2. The manifest translates probe IDs (e.g. `cg00000029`) into genomic
   coordinates (`chr1:15,000,000`).
3. `ep.query.query_sites()` fetches the beta values at the required positions
   from the methylstore.
4. For each sample, the predicted age is computed as a weighted linear
   combination of beta values plus the intercept term.
5. The result is joined onto `md.obs` as the column `epigenetic_age`.

## Supported Clocks

| Clock | CpGs | Tissue | Reference |
|-------|-------|--------|-----------|
| Horvath multi-tissue | 353 | Pan-tissue | Horvath 2013 |
| Hannum blood | 71 | Blood | Hannum et al. 2013 |
| PhenoAge | 513 | Blood | Levine et al. 2018 |
| DunedinPACE | ~20,000 | Blood (pace-of-aging) | Belsky et al. 2022 |

The `clock_name` parameter selects the intercept handling and any
clock-specific normalization. The actual weights always come from the
user-supplied `coefficients` file.

## Coefficient Table Format

The coefficient CSV must contain at least two columns:

| Column | Description |
|--------|-------------|
| `probe_id` | Illumina CpG probe ID (e.g. `cg00000029`) |
| `coefficient` | Clock weight for this probe |

An optional `intercept` row (with `probe_id = "intercept"`) supplies the
model intercept. If absent, the intercept defaults to zero.

```
probe_id,coefficient
intercept,0.6955
cg00000029,0.0123
cg00000165,-0.0042
...
```

## Manifest

The manifest maps probe IDs to genomic positions so that epykit can locate
the corresponding CpG sites in the methylstore. Standard Illumina manifest
files work directly:

- **HM450** -- `HM450.hg38.manifest.tsv.gz`
- **EPIC v1** -- `EPIC.hg38.manifest.tsv.gz`
- **EPIC v2** -- `EPICv2.hg38.manifest.tsv.gz`

The manifest must contain `probe_id`, `chrom`, and `pos` columns (or the
Illumina-standard equivalents `IlmnID`, `CHR`, `MAPINFO`).

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `coefficients` | str / Path / DataFrame | required | Clock coefficient table (CSV, TSV, or Polars DataFrame) |
| `manifest` | str / Path / DataFrame | required | Illumina CpG manifest mapping probe IDs to positions |
| `clock_name` | str | `"horvath"` | Clock model name (controls intercept handling) |
| `missing` | str | `"drop"` | How to handle missing CpGs: `"drop"` skips them, `"impute_mean"` uses the training-set mean |

## Comparing Predicted vs. Chronological Age

```python
import polars as pl

obs = md.obs.select("sample_id", "age", "epigenetic_age")

# Age acceleration = residual from regressing epigenetic age on chronological age
obs = obs.with_columns(
    (pl.col("epigenetic_age") - pl.col("age")).alias("age_acceleration")
)
print(obs)
```

## Notes

- Clock accuracy depends on CpG coverage. Low-coverage WGBS data may lack
  many of the required probe positions, leading to unreliable estimates.
  Check `md.obs["age_clock_cpgs_found"]` for the number of matched sites.
- DunedinPACE measures the *pace* of aging (years of biological aging per
  calendar year), not absolute age. Its output scale differs from the other
  clocks.
- The function does not ship coefficient files. Obtain them from the
  original publications or from public repositories such as the
  `methylclock` R package.
