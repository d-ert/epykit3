#!/usr/bin/env Rscript
# install_packages.R -- one-shot bootstrap for the benchmark R stack.
#
# Run once on the Linux machine (or inside rocker/r-ver:4.5.0) to
# install every package the benchmarks need at the versions pinned in
# PROTOCOL.md, then snapshot the result into benchmark/renv.lock.
# Subsequent runs of Dockerfile.r will pick up the lockfile and
# restore the exact frozen environment.
#
# Usage:
#   Rscript benchmark/renv/install_packages.R
#
# This is intentionally a once-per-environment bootstrap, not part
# of regen_all.py -- the goal is to produce the lockfile, after
# which everything else flows through renv::restore() inside the
# container.

cat("[install_packages.R] starting at", as.character(Sys.time()), "\n")

# ---- renv ----------------------------------------------------------------
if (!"renv" %in% rownames(installed.packages())) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialise the renv project at the repo root so the lockfile lands at
# benchmark/renv.lock (we set restore-anchor below).
repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."))
setwd(repo_root)
cat("[install_packages.R] repo root:", repo_root, "\n")

if (!dir.exists("benchmark/renv")) {
  dir.create("benchmark/renv", recursive = TRUE)
}

# renv::init writes .Rprofile + an `renv/` skeleton in the project
# root. We avoid that here -- we only want the LOCKFILE under
# benchmark/, not a full project takeover -- by going straight to
# install + snapshot with an explicit project hint.
options(renv.config.startup.quiet = TRUE)
renv::activate(project = repo_root)

# ---- Bioconductor 3.21 (matches R 4.5.0) ---------------------------------
if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(version = "3.21", update = FALSE, ask = FALSE)

# ---- Package versions pinned in PROTOCOL.md sec 1 -----------------------
# We install the names; the eventual lockfile records the exact built
# versions for renv::restore(). Bioconductor 3.21 ships with the
# versions referenced in PROTOCOL.md (methylKit 1.36.0+, DSS 2.52.0+,
# dmrseq 1.26.0+, bsseq 1.42.0+).
bioc_pkgs <- c(
  "methylKit",
  "DSS",
  "dmrseq",
  "bsseq",
  "GenomicRanges",
  "DelayedMatrixStats"
)
cran_pkgs <- c(
  "data.table",
  "optparse"
)

cat("[install_packages.R] installing CRAN deps:", paste(cran_pkgs, collapse = ", "), "\n")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

cat("[install_packages.R] installing Bioconductor deps:", paste(bioc_pkgs, collapse = ", "), "\n")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# ---- Snapshot to benchmark/renv.lock -------------------------------------
# renv::snapshot() writes the lockfile at the project root by default;
# we override the path to keep it under benchmark/ for parity with the
# rest of the benchmark assets.
cat("[install_packages.R] snapshotting environment...\n")
renv::snapshot(
  project = repo_root,
  lockfile = "benchmark/renv.lock",
  prompt = FALSE,
  type = "all"
)

cat("[install_packages.R] DONE. Lockfile at benchmark/renv.lock\n")
cat("[install_packages.R] To restore in a fresh container:\n")
cat("    R -q -e \"renv::restore(lockfile = 'benchmark/renv.lock', prompt = FALSE)\"\n")
