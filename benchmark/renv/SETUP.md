# Benchmark R environment

The R-side benchmarks (methylKit, DSS, dmrseq, BSmooth) need a frozen
environment so a third reader can reproduce the exact numbers in the
paper. We use [renv](https://rstudio.github.io/renv/) for that.

This directory holds two files (one of which is generated):

  * [install_packages.R](install_packages.R) -- one-shot bootstrap
    that installs every package the benchmarks need at the versions
    pinned in PROTOCOL.md, then snapshots the result into
    `benchmark/renv.lock`. **Run this once on the Linux machine
    before the first benchmark run.**
  * `benchmark/renv.lock` (generated) -- the snapshot that
    [Dockerfile.r](../../Dockerfile.r) restores into the container.
    Committed to git so the lockfile is part of the submission
    bundle.

## Running install_packages.R

On the Linux machine (or inside the rocker/r-ver:4.5.0 container):

```bash
Rscript benchmark/renv/install_packages.R
```

The script:

1. Installs `renv` from CRAN if absent.
2. Initialises an `renv` project rooted at the repo top-level.
3. Installs Bioconductor 3.21 (matching R 4.5.0) and the package
   set the benchmarks use:
   - methylKit (>= 1.36.0; PROTOCOL.md §1)
   - DSS (>= 2.52.0)
   - dmrseq (>= 1.26.0)
   - bsseq (>= 1.42.0)
   - data.table, optparse, GenomicRanges, DelayedMatrixStats (transitives)
4. Runs `renv::snapshot()` to write `benchmark/renv.lock`.

After this, [Dockerfile.r](../../Dockerfile.r) will pick up the
lockfile during `docker build` and produce a container that
restores the exact frozen environment.

## Why this is bootstrapped, not committed pre-built

A renv lockfile pins every transitive dependency by exact version
(plus checksum). The right way to generate one is on the target OS
(Linux x86_64) so the system-specific binaries are recorded
correctly. Shipping a hand-written stub lockfile from Windows
would be misleading and would not actually reproduce the paper
numbers.
