# methylKit Tabix

Export methylation data as tabix-indexed TSV files compatible with the methylKit R
package.

## ep.to_methylkit_tabix

```python
ep.to_methylkit_tabix(md, dir)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with loaded store |
| `dir` | `str` or `Path` | Output directory. One TSV + tabix index per sample |

**Output**

For each sample, the function writes:

- `<sample_id>.txt.gz` -- bgzipped TSV with methylKit columns
- `<sample_id>.txt.gz.tbi` -- tabix index

The TSV columns follow the methylKit format:

| Column | Description |
|--------|-------------|
| `chrBase` | Chromosome and position (e.g., `chr1.1234`) |
| `chr` | Chromosome |
| `base` | Genomic position |
| `strand` | `F` or `R` |
| `coverage` | Total read coverage |
| `freqC` | Percent methylation (0--100) |
| `freqT` | Percent unmethylated (0--100) |

## Platform Requirements

This function depends on `pysam` for bgzip compression and tabix indexing. `pysam` is
available on Linux and macOS only -- there is no Windows build.

```bash
pip install pysam
```

## Example

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.pp.filter_coverage(md, min_cov=10)

ep.to_methylkit_tabix(md, dir="methylkit_output/")
```

This produces one file pair per sample:

```
methylkit_output/
  tumor_1.txt.gz
  tumor_1.txt.gz.tbi
  tumor_2.txt.gz
  tumor_2.txt.gz.tbi
  normal_1.txt.gz
  normal_1.txt.gz.tbi
  ...
```

## Cross-Validation with methylKit in R

The exported files can be read directly in R using `methylKit::methRead()`:

```r
library(methylKit)

file_list <- list(
  "methylkit_output/tumor_1.txt.gz",
  "methylkit_output/tumor_2.txt.gz",
  "methylkit_output/normal_1.txt.gz",
  "methylkit_output/normal_2.txt.gz"
)

sample_ids <- list("tumor_1", "tumor_2", "normal_1", "normal_2")
treatment <- c(1, 1, 0, 0)

obj <- methRead(
  file_list,
  sample.id = sample_ids,
  treatment = treatment,
  assembly = "hg38",
  pipeline = "bismarkCytosineReport",
  header = TRUE
)

# Proceed with methylKit analysis
meth <- unite(obj)
diff <- calculateDiffMeth(meth)
```

This workflow is useful for cross-validating epykit results against methylKit, or for
handing off data to collaborators who prefer the R ecosystem.
