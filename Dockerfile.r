# Dockerfile.r --- epykit benchmark R stack for reproducible competitor runs.
#
# Phase 1.4 of the GB resubmission plan: ship the R-side counterpart to
# Dockerfile.python, with methylKit, DSS, dmrseq, BSmooth (bsseq), and
# all their Bioconductor + CRAN transitive deps locked via a renv.lock
# file. M7 (no renv.lock / Docker for the R stack) was a real
# reproducibility critique in the Nature-tier review; this container
# closes it.
#
# Usage:
#   docker build -t epykit-r -f Dockerfile.r .
#   docker run --rm -v "$PWD":/work -w /work epykit-r \
#       Rscript benchmark/scripts/run_dss_simulator.R --in-dir ... --out ...
#   docker run --rm -v "$PWD":/work -w /work epykit-r \
#       Rscript benchmark/scripts/run_dmrseq.R --in-dir ... --out ...
#
# The rocker/r-ver:4.5.0 base is the canonical "frozen R + system libs"
# anchor; on top of it we restore the project renv lockfile which pins
# every R package (CRAN + Bioconductor + dev) to the exact version used
# by the paper.
FROM rocker/r-ver:4.5.0

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    RENV_PATHS_LIBRARY=/work/benchmark/renv/library

# System deps:
#   - libcurl/openssl/libxml2: every Bioconductor HTTP-using package
#   - libgit2/libssh2: renv pulls some packages from GitHub
#   - libz/libbz2/liblzma: Rsamtools / GenomicAlignments
#   - libpng/libtiff: bsseq plotting (occasionally invoked from scripts)
#   - libgsl/libfftw3/libxt: dmrseq's optional plotting deps
#   - tini: PID-1 reaping
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates curl git tini \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libgit2-dev libssh2-1-dev \
        libz-dev libbz2-dev liblzma-dev \
        libpng-dev libtiff5-dev \
        libgsl-dev libfftw3-dev libxt-dev \
        libcairo2-dev libpango1.0-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Install renv before copying the project so the renv binary layer is
# stable across iterations. renv pins ITSELF in the lockfile, so the
# bootstrap version here just needs to be recent enough to restore.
RUN R -q -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# Restore the locked R environment. We copy ONLY the lockfile here so
# Docker's layer cache only busts on lockfile changes; the project
# source is volume-mounted at run time.
#
# If benchmark/renv.lock is absent (early-stage development), this
# section degrades to "no R packages installed beyond rocker base"
# rather than failing the build -- so the container is usable for
# bootstrapping a fresh lockfile.
COPY benchmark/renv.lock* benchmark/renv.lock
RUN if [ -f benchmark/renv.lock ]; then \
        R -q -e "renv::restore(lockfile = 'benchmark/renv.lock', prompt = FALSE)" ; \
    else \
        echo "WARNING: benchmark/renv.lock not present; R deps NOT installed." ; \
    fi

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
