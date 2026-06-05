#!/usr/bin/env bash
# run_tabix_vs_epykit_gse.sh
# Staged + RESUMABLE full-genome GSE263850 head-to-head: epykit lr vs methylKit
# {ram,tabix} at {1,8} cores. Each config is an isolated driver call appending to
# one CSV. The CSV is NEVER wiped: a config already present (returncode 0) is
# skipped, so an accidental kill only ever costs the single config in flight.
# Ordered by information value + safety: hero + fair-fast first, OOM-risky
# RAM-8core isolated in the middle, the long single-core runs last (overnight).
#
# Launch DETACHED so it survives a client disconnect:
#   setsid nohup bash run_tabix_vs_epykit_gse.sh > gse_run.log 2>&1 < /dev/null &
#
# mincov=1 matches epykit (no coverage filter; unite intersect keeps sites
# covered in all 6 samples) -> fair + max-data memory stress test.
set -u
cd /home/mlegrand/Desktop/Deniz/epykit_bechmarking/epykit3

IN=/home/mlegrand/Desktop/Deniz/raw_data/GSE263850/covs
SCRATCH=/home/mlegrand/Desktop/Deniz/tve_scratch_gse
OUT=/home/mlegrand/Desktop/Deniz/phi_sweep_export_2026-06-05/summaries/tabix_vs_epykit_gse_per_run.csv
DRV=benchmark/scripts/run_tabix_vs_epykit.py
COMMON="--in-dir $IN --glob *.cov.gz --n-per-group 3 --mincov 1 --assembly hg38 \
        --dataset gse263850 --cell-id whole_genome --scratch $SCRATCH --out-csv $OUT"

mkdir -p "$SCRATCH"
echo "### GSE full-genome tabix-vs-epykit (resumable) started: $(date -Is)"

have () {  # $1=tool $2=cores -> 0 if a successful row already exists
  [ -f "$OUT" ] || return 1
  awk -F, -v t="$1" -v c="$2" 'NR>1 && $3==t && $5==c && $14==0 {f=1} END{exit !f}' "$OUT"
}

run () {  # $1=label $2=tool $3=cores, rest = extra driver args
  local label="$1" tool="$2" cores="$3"; shift 3
  if have "$tool" "$cores"; then
    echo "### [$label] SKIP (already done) $(date -Is)"; return 0
  fi
  echo "### [$label] start $(date -Is)"
  uv run --python 3.12 python $DRV $COMMON "$@"
  echo "### [$label] end   $(date -Is) rc=$?"
}

# 1) epykit lr (single-thread, out-of-core parquet) -- hero number (~80 s)
run epykit   epykit_lr       1 --cores 1 --skip-ram --skip-tabix
# 2) methylKit tabix, 8 cores -- fair fast out-of-core baseline (~1-1.5 h)
run tabix_c8 methylkit_tabix 8 --cores 8 --skip-ram --skip-epykit
# 3) methylKit RAM, 8 cores -- OOM-risk (~48 GB vs ~50 GB free); isolated (~1 h)
run ram_c8   methylkit_ram   8 --cores 8 --skip-tabix --skip-epykit
# 4) methylKit RAM, 1 core -- honest single-process peak RSS (~5 h)
run ram_c1   methylkit_ram   1 --cores 1 --skip-tabix --skip-epykit
# 5) methylKit tabix, 1 core -- completeness (~6 h)
run tabix_c1 methylkit_tabix 1 --cores 1 --skip-ram --skip-epykit

echo "### GSE full-genome tabix-vs-epykit finished: $(date -Is)"
echo "### results: $OUT"
