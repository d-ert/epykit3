# Samplesheet Format

All `ep.read_*` functions take a samplesheet CSV as their first argument. This
page documents the expected format.

## Required columns

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | Unique identifier for each sample. Used as the partition key in the methylstore (`sample=<sample_id>/`). |
| `group` | string | Group label (e.g. `"tumor"`, `"normal"`, `"KO"`, `"WT"`). Used for treatment/control assignment and multi-group contrasts. |
| `path` | string | Filesystem path to the input coverage file (Bismark `.cov[.gz]` or MethylDackel `.bedGraph[.gz]`). Absolute or relative to the working directory. |

## Minimal example

```csv
sample_id,group,path
tumor_1,tumor,/data/bismark/tumor_1.cov.gz
tumor_2,tumor,/data/bismark/tumor_2.cov.gz
normal_1,normal,/data/bismark/normal_1.cov.gz
normal_2,normal,/data/bismark/normal_2.cov.gz
```

## Extra columns (covariates)

Any column beyond the three required ones is preserved on `md.obs` and can be
used as a GLM covariate in `ep.tl.dmc()`.

```csv
sample_id,group,path,sex,batch,age,donor
tumor_1,tumor,/data/tumor_1.cov.gz,M,batch_1,62,D001
tumor_2,tumor,/data/tumor_2.cov.gz,F,batch_2,55,D002
tumor_3,tumor,/data/tumor_3.cov.gz,M,batch_1,71,D003
normal_1,normal,/data/normal_1.cov.gz,M,batch_1,62,D001
normal_2,normal,/data/normal_2.cov.gz,F,batch_2,55,D002
normal_3,normal,/data/normal_3.cov.gz,M,batch_1,71,D003
```

Use these covariates in the DMC call via the `formula` parameter:

```python
ep.tl.dmc(md, formula="~ group + age + sex", contrast="group")
```

Or via the `covariates` convenience parameter:

```python
ep.tl.dmc(md, covariates=["age", "sex"])
```

Covariates are also available for the tile-based DMR method:

```python
ep.tl.dmr(md, method="tile", design="~ group + batch", covariates=["batch"])
```

## Treatment and control groups

### Binary mode (two groups)

Pass `treatment_group` and `control_group` to designate which group is
case (treatment=1) and which is control (treatment=0):

```python
md = ep.read_bismark(
    "samples.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
)
```

Only rows whose `group` column matches one of these two values are included.
A `treatment` column (integer 0 or 1) is added to `md.obs` automatically.

### Multi-group mode (3+ groups)

Pass `groups=[...]` instead of the treatment/control pair:

```python
md = ep.read_bismark(
    "samples.csv",
    groups=["WT", "KO", "rescue"],
    assembly="hg38",
)
```

In multi-group mode, no binary `treatment` column is created. Use the
`formula` and `contrast` parameters in `ep.tl.dmc()` for pairwise or
joint tests:

```python
# Joint F-test across all three groups
ep.tl.dmc(md, formula="~ group", contrast="group")

# Pairwise: KO vs WT
ep.tl.dmc(md, formula="~ group", contrast="group[T.KO] - group[T.WT]")
```

## Validation

`ep.read_bismark()` validates the samplesheet on load:

- Raises `ValueError` if any of the three required columns are missing.
- Raises `ValueError` if the samplesheet contains no rows.
- Raises `ValueError` if no rows match the requested `treatment_group` /
  `control_group` or `groups`.
- Rows whose `group` value does not match the requested groups are silently
  skipped.

## File format notes

- The CSV must have a header row.
- Standard Python CSV parsing is used (comma-separated, optional quoting).
- Paths may be absolute or relative. Relative paths are resolved from the
  current working directory at the time `ep.read_bismark()` is called.
- Bismark `.cov` and `.cov.gz` files are both supported. MethylDackel
  `.bedGraph` and `.bedGraph.gz` files are supported by
  `ep.read_methyldackel()`.
