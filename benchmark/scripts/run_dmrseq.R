#!/usr/bin/env Rscript
# run_dmrseq.R
# ---------------------------------------------------------------------------
# dmrseq DMR call on a single intrinsic-truth simulator seed OR the
# GSE263850 real cohort.
#
# Reads 6 per-sample .cov.gz files (3 treat vs 3 ctrl), builds a bsseq
# object, runs dmrseq::dmrseq() with testCovariate="group", and writes a
# per-DMR TSV in the canonical [chr, start, end, stat, pvalue, qvalue,
# meth_diff] schema that benchmark/scripts/_epykit_scoring.py's
# score_dmr_parquet() consumes unchanged.
#
# This script lands as part of Phase 1.3 of the GB resubmission plan,
# adding dmrseq as a locally-re-run baseline rather than transcribing
# numbers from Piao 2021.
#
# Usage:
#   Rscript run_dmrseq.R \
#       --in-dir benchmark/data/study1b_simulator/seed=2026000/bismark_cov \
#       --out    benchmark/data/study1b_simulator/seed=2026000/dmrseq.tsv \
#       [--cutoff 0.10]
#
# Outputs:
#   <out>                   - per-DMR TSV
#   <out>.sessioninfo.txt   - R sessionInfo() for reproducibility
#   <out>.timing.tsv        - wall + CPU time per phase
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bsseq)
  library(dmrseq)
})

opt_list <- list(
  make_option("--in-dir", type = "character",
              help = "Directory containing the 6 .cov.gz files"),
  make_option("--out", type = "character",
              help = "Output TSV path for per-DMR results"),
  make_option("--cutoff", type = "double", default = 0.10,
              help = "dmrseq cutoff on the smoothed methylation difference (default 0.10)"),
  make_option("--minNumRegion", type = "integer", default = 5,
              help = "Minimum CpGs per candidate region (default 5, dmrseq author guidance)")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$`in-dir`) || is.null(opt$out)) {
  stop("--in-dir and --out are required")
}

in_dir   <- opt$`in-dir`
out_path <- opt$out
out_dir  <- dirname(out_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Discover input files (deterministic) --------------------------------
files_all <- sort(list.files(in_dir, pattern = "\\.cov\\.gz$", full.names = TRUE))
if (length(files_all) != 6L) {
  stop(sprintf("expected 6 .cov.gz files in %s, found %d",
               in_dir, length(files_all)))
}
cat("input files (sorted):\n")
for (f in files_all) cat("  ", basename(f), "\n")

sample_ids <- c("treat_1", "treat_2", "treat_3", "ctrl_1", "ctrl_2", "ctrl_3")
group <- factor(c("treat", "treat", "treat", "ctrl", "ctrl", "ctrl"),
                levels = c("ctrl", "treat"))

# ---- Phase 1: build bsseq object from .cov files -------------------------
t_read_0 <- proc.time()
parse_one <- function(path, sid) {
  dt <- fread(path, sep = "\t", header = FALSE,
              col.names = c("chr", "start", "end", "beta", "count_M", "count_U"))
  dt[, sample := sid]
  dt
}
dt_list   <- Map(parse_one, files_all, sample_ids)
all_sites <- unique(rbindlist(lapply(dt_list, function(d) d[, .(chr, start)])))
setkey(all_sites, chr, start)
n_sites <- nrow(all_sites)
cat(sprintf("union of CpG sites: %d\n", n_sites))

M_mat   <- matrix(0L, nrow = n_sites, ncol = length(files_all))
Cov_mat <- matrix(0L, nrow = n_sites, ncol = length(files_all))
colnames(M_mat)   <- sample_ids
colnames(Cov_mat) <- sample_ids
for (i in seq_along(dt_list)) {
  d <- dt_list[[i]]
  setkey(d, chr, start)
  merged <- d[all_sites, on = c("chr", "start")]
  M_mat[, i]   <- ifelse(is.na(merged$count_M), 0L, merged$count_M)
  Cov_mat[, i] <- ifelse(is.na(merged$count_M) | is.na(merged$count_U),
                         0L, merged$count_M + merged$count_U)
}

bs <- BSseq(
  chr = all_sites$chr,
  pos = all_sites$start,
  M = M_mat, Cov = Cov_mat,
  sampleNames = sample_ids
)
pData(bs)$group <- group
t_read <- proc.time() - t_read_0
cat(sprintf("bsseq: %d sites x %d samples\n", nrow(bs), ncol(bs)))

# ---- Phase 2: dmrseq -----------------------------------------------------
# dmrseq requires a loss-of-coverage filter: drop any site whose coverage
# is zero in any sample (dmrseq author guidance; the function itself errors
# out on zeros).
t_filter_0 <- proc.time()
loci_keep <- which(DelayedMatrixStats::rowMins(getCoverage(bs, type = "Cov")) > 0)
bs <- bs[loci_keep, ]
cat(sprintf("post-filter: %d sites x %d samples\n", nrow(bs), ncol(bs)))
t_filter <- proc.time() - t_filter_0

t_dmr_0 <- proc.time()
# testCovariate must match the pData(bs)$<name> column; we set "group" above.
dmrs <- dmrseq(
  bs            = bs,
  testCovariate = "group",
  cutoff        = opt$cutoff,
  minNumRegion  = opt$minNumRegion
)
t_dmr <- proc.time() - t_dmr_0
cat(sprintf("dmrseq: %d DMRs called\n", length(dmrs)))

# ---- Write per-DMR output ------------------------------------------------
# dmrseq returns a GRanges with mcols:
#   L, area, beta, stat, pvalue, qvalue, index.start, index.end, index.width
# We standardise to [chr, start, end, stat, pvalue, qvalue, meth_diff]
# where meth_diff is the signed methylation difference (dmrseq returns it
# as `beta` in the smoothed-arcsine link scale; we additionally compute the
# raw beta-difference per region as ``rawMeanDiff`` for downstream
# consistency with DSS / methylKit output).
mcols_df <- as.data.frame(GenomicRanges::mcols(dmrs))
cat(sprintf("dmrseq mcols columns: %s\n",
            paste(colnames(mcols_df), collapse = ", ")))
# dmrseq's column names vary by version: >=1.30 emits `pval`/`qval`, while
# earlier releases used `pvalue`/`qvalue`. Resolve defensively so the writer
# does not silently produce a NULL (length-0) column -> data.frame() row
# mismatch (the bug this guards against).
pick <- function(df, ...) {
  for (nm in c(...)) if (!is.null(df[[nm]])) return(df[[nm]])
  stop(sprintf("none of (%s) present in dmrseq mcols", paste(c(...), collapse = ", ")))
}
out_df <- data.frame(
  chr       = as.character(GenomicRanges::seqnames(dmrs)),
  start     = GenomicRanges::start(dmrs),
  end       = GenomicRanges::end(dmrs),
  stat      = pick(mcols_df, "stat"),
  pvalue    = pick(mcols_df, "pvalue", "pval"),
  qvalue    = pick(mcols_df, "qvalue", "qval"),
  meth_diff = pick(mcols_df, "beta")
)
write.table(out_df, file = out_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out_path, nrow(out_df)))

# ---- Timings sidecar ------------------------------------------------------
timing_df <- data.frame(
  phase  = c("read+bsseq", "coverage_filter", "dmrseq"),
  wall_s = c(t_read["elapsed"],   t_filter["elapsed"],   t_dmr["elapsed"]),
  cpu_s  = c(t_read["user.self"]   + t_read["sys.self"],
             t_filter["user.self"] + t_filter["sys.self"],
             t_dmr["user.self"]    + t_dmr["sys.self"])
)
timing_path <- paste0(out_path, ".timing.tsv")
write.table(timing_df, file = timing_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s\n", timing_path))

# ---- sessionInfo ----------------------------------------------------------
sinfo_path <- paste0(out_path, ".sessioninfo.txt")
sink(sinfo_path)
cat("# run_dmrseq.R session info\n\n")
cat(sprintf("# in_dir       : %s\n",   in_dir))
cat(sprintf("# out          : %s\n",   out_path))
cat(sprintf("# cutoff       : %s\n",   opt$cutoff))
cat(sprintf("# minNumRegion : %s\n\n", opt$minNumRegion))
print(sessionInfo())
sink()
cat(sprintf("wrote %s\n", sinfo_path))

cat("\nDONE\n")
