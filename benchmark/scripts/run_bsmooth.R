#!/usr/bin/env Rscript
# run_bsmooth.R
# ---------------------------------------------------------------------------
# BSmooth DMR call on a single intrinsic-truth simulator seed OR the
# GSE263850 real cohort.
#
# Reads 6 per-sample .cov.gz files (3 treat vs 3 ctrl), builds a bsseq
# object, runs bsseq::BSmooth() + BSmooth.tstat() + dmrFinder(), and
# writes a per-DMR TSV in the canonical scoring schema.
#
# Note on the scoring contract:
#   BSmooth's dmrFinder() is a cutoff-based caller -- it does NOT emit
#   p-values or q-values. To be consumable by
#   _epykit_scoring.score_dmr_parquet() (which thresholds on qvalue), we
#   emit a synthesised ``qvalue`` derived from the rank of the absolute
#   t-statistic (smaller q == more confident). This preserves DMR
#   ranking but is NOT calibrated FDR; the paper text must call this
#   out alongside the BSmooth numbers. dmrseq is the calibrated
#   p-value/q-value reference among DMR callers; BSmooth is the
#   smoothing baseline.
#
# This script lands as part of Phase 1.3 of the GB resubmission plan,
# adding BSmooth as a locally-re-run baseline rather than transcribing
# numbers from Piao 2021.
#
# Usage:
#   Rscript run_bsmooth.R \
#       --in-dir benchmark/data/study1b_simulator/seed=2026000/bismark_cov \
#       --out    benchmark/data/study1b_simulator/seed=2026000/bsmooth.tsv \
#       [--tstat-cutoff "-4.6,4.6"] [--min-cpgs 3]
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
})

opt_list <- list(
  make_option("--in-dir", type = "character",
              help = "Directory containing the 6 .cov.gz files"),
  make_option("--out", type = "character",
              help = "Output TSV path for per-DMR results"),
  make_option("--tstat-cutoff", type = "character", default = "-4.6,4.6",
              help = "Comma-separated low,high t-stat thresholds (default: BSmooth defaults -4.6, 4.6)"),
  make_option("--min-cpgs", type = "integer", default = 3,
              help = "Minimum CpGs per DMR (default 3)"),
  make_option("--bsmooth-h", type = "integer", default = 1000,
              help = "BSmooth window half-width in bp (default 1000)"),
  make_option("--bsmooth-ns", type = "integer", default = 70,
              help = "BSmooth minimum number of CpGs per window (default 70)")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$`in-dir`) || is.null(opt$out)) {
  stop("--in-dir and --out are required")
}

in_dir   <- opt$`in-dir`
out_path <- opt$out
out_dir  <- dirname(out_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cutoff_parts <- as.numeric(strsplit(opt$`tstat-cutoff`, ",", fixed = TRUE)[[1]])
if (length(cutoff_parts) != 2L || any(is.na(cutoff_parts))) {
  stop(sprintf("--tstat-cutoff must be 'low,high' (got %s)",
               opt$`tstat-cutoff`))
}
cutoff_low  <- cutoff_parts[1]
cutoff_high <- cutoff_parts[2]

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

# ---- Phase 2: BSmooth -----------------------------------------------------
t_smooth_0 <- proc.time()
bs_smoothed <- BSmooth(bs, h = opt$`bsmooth-h`, ns = opt$`bsmooth-ns`,
                       verbose = FALSE)
t_smooth <- proc.time() - t_smooth_0
cat("BSmooth: done\n")

# ---- Phase 3: BSmooth.tstat ----------------------------------------------
t_tstat_0 <- proc.time()
group1_ids <- sample_ids[group == "treat"]
group2_ids <- sample_ids[group == "ctrl"]
# BSmooth.tstat requires both groups have >= 2 samples with non-zero
# coverage at the tested CpG. Drop low-coverage sites first.
keep_idx <- which(
  DelayedMatrixStats::rowSums2(getCoverage(bs_smoothed) >= 1L,
                               cols = which(group == "treat")) >= 2L &
  DelayedMatrixStats::rowSums2(getCoverage(bs_smoothed) >= 1L,
                               cols = which(group == "ctrl"))  >= 2L
)
bs_smoothed <- bs_smoothed[keep_idx, ]
tstat <- BSmooth.tstat(
  bs_smoothed,
  group1   = group1_ids,
  group2   = group2_ids,
  estimate.var = "same",
  local.correct = TRUE,
  verbose       = FALSE
)
t_tstat <- proc.time() - t_tstat_0
cat("BSmooth.tstat: done\n")

# ---- Phase 4: dmrFinder ---------------------------------------------------
t_dmr_0 <- proc.time()
dmrs <- dmrFinder(
  tstat,
  cutoff = c(cutoff_low, cutoff_high),
  qcutoff = NULL,
  maxGap = 300L,
  verbose = FALSE
)
t_dmr <- proc.time() - t_dmr_0
if (is.null(dmrs)) dmrs <- data.frame()
cat(sprintf("dmrFinder: %d candidate DMRs\n", nrow(dmrs)))

# Filter by minimum number of CpGs
if (nrow(dmrs) > 0L) {
  dmrs <- dmrs[dmrs$n >= opt$`min-cpgs`, , drop = FALSE]
  cat(sprintf("dmrFinder + min-cpgs >= %d: %d DMRs\n",
              opt$`min-cpgs`, nrow(dmrs)))
}

# ---- Write per-DMR output ------------------------------------------------
# Standardise to [chr, start, end, stat, pvalue, qvalue, meth_diff].
# BSmooth has no native pvalue/qvalue; we synthesise a rank-based qvalue
# from |tstat| so downstream qvalue-threshold scoring is well-defined.
# This is documented above and must be called out in the paper.
if (nrow(dmrs) == 0L) {
  out_df <- data.frame(
    chr = character(), start = integer(), end = integer(),
    stat = numeric(), pvalue = numeric(), qvalue = numeric(),
    meth_diff = numeric()
  )
} else {
  ranks <- rank(-abs(dmrs$areaStat), ties.method = "average")
  synth_qvalue <- ranks / nrow(dmrs)  # monotone in -|tstat|; in range 0..1
  out_df <- data.frame(
    chr       = as.character(dmrs$chr),
    start     = as.integer(dmrs$start),
    end       = as.integer(dmrs$end),
    stat      = dmrs$areaStat,
    pvalue    = synth_qvalue,  # synthesised; not a calibrated p-value
    qvalue    = synth_qvalue,
    meth_diff = dmrs$meanDiff
  )
}
write.table(out_df, file = out_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out_path, nrow(out_df)))

# ---- Timings sidecar ------------------------------------------------------
timing_df <- data.frame(
  phase  = c("read+bsseq", "BSmooth", "BSmooth.tstat", "dmrFinder"),
  wall_s = c(t_read["elapsed"],   t_smooth["elapsed"],
             t_tstat["elapsed"],  t_dmr["elapsed"]),
  cpu_s  = c(t_read["user.self"]   + t_read["sys.self"],
             t_smooth["user.self"] + t_smooth["sys.self"],
             t_tstat["user.self"]  + t_tstat["sys.self"],
             t_dmr["user.self"]    + t_dmr["sys.self"])
)
timing_path <- paste0(out_path, ".timing.tsv")
write.table(timing_df, file = timing_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s\n", timing_path))

# ---- sessionInfo ----------------------------------------------------------
sinfo_path <- paste0(out_path, ".sessioninfo.txt")
sink(sinfo_path)
cat("# run_bsmooth.R session info\n\n")
cat(sprintf("# in_dir       : %s\n",   in_dir))
cat(sprintf("# out          : %s\n",   out_path))
cat(sprintf("# tstat_cutoff : (%s, %s)\n", cutoff_low, cutoff_high))
cat(sprintf("# min_cpgs     : %s\n",   opt$`min-cpgs`))
cat(sprintf("# bsmooth_h    : %s\n",   opt$`bsmooth-h`))
cat(sprintf("# bsmooth_ns   : %s\n\n", opt$`bsmooth-ns`))
print(sessionInfo())
sink()
cat(sprintf("wrote %s\n", sinfo_path))

cat("\nDONE\n")
