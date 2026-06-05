#!/usr/bin/env Rscript
# run_methylkit_backend.R
# ---------------------------------------------------------------------------
# methylKit DMC run with a SELECTABLE storage backend:
#   --dbtype memory  : in-memory (data.frame/data.table) -- the default backend
#   --dbtype tabix   : on-disk, bgzipped + tabix-indexed flat files under --dbdir
#                      (methylRawListDB / methylBaseDB / methylDiffDB), streamed
#                      through instead of held in RAM.
#
# Both paths call the IDENTICAL methRead -> unite -> calculateDiffMeth pipeline
# with the same mincov / treatment / overdispersion, so the ONLY difference is
# where the data lives. tabix changes storage/IO, NOT the statistics: the
# per-CpG p/q/meth.diff are identical to memory mode. The point of comparison is
# therefore peak RSS and wall time (and n_sites, as a correctness check).
#
# Reads sorted bismarkCoverage .cov.gz in --in-dir (first --n-per-group =
# treatment, next = control), matching run_epykit_cell.py so inputs are
# byte-identical across tools.
#
# Usage:
#   Rscript run_methylkit_backend.R --dbtype tabix \
#       --in-dir <dir> --out <per-CpG.tsv> --dbdir <scratch tabix dir> \
#       [--cores 1] [--mincov 1] [--n-per-group 3] [--assembly sim] \
#       [--overdispersion none|MN] [--test default|Chisq|F]
#
# Outputs: <out> (chr,start,end,strand,pvalue,qvalue,meth.diff),
#          <out>.timing.tsv, <out>.sessioninfo.txt
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(methylKit)
})

opt_list <- list(
  make_option("--dbtype",  type = "character", default = "memory",
              help = "storage backend: 'memory' (default) or 'tabix'"),
  make_option("--in-dir",  type = "character", help = "Directory containing the .cov.gz files"),
  make_option("--out",     type = "character", help = "Output TSV path for per-CpG results"),
  make_option("--dbdir",   type = "character", default = NA_character_,
              help = "Scratch dir for tabix DB files (required when --dbtype tabix; wiped)"),
  make_option("--cores",   type = "integer",   default = 1L, help = "mc.cores for calculateDiffMeth"),
  make_option("--mincov",  type = "integer",   default = 1L, help = "Min per-CpG coverage in methRead"),
  make_option("--n-per-group", type = "integer", default = 3L,
              help = "Samples per group; first N = treatment, next N = control"),
  make_option("--overdispersion", type = "character", default = "none",
              help = "calculateDiffMeth overdispersion: 'none' (default) or 'MN'"),
  make_option("--test",    type = "character", default = "default",
              help = "test: 'default' (F / Chisq-when-MN) or explicit 'Chisq'/'F'"),
  make_option("--assembly", type = "character", default = "sim", help = "assembly label")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$`in-dir`) || is.null(opt$out)) stop("--in-dir and --out are required")
use_tabix <- identical(opt$dbtype, "tabix")
if (use_tabix && is.na(opt$dbdir)) stop("--dbdir is required when --dbtype tabix")

in_dir   <- opt$`in-dir`
out_path <- opt$out
npg      <- opt$`n-per-group`
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

if (use_tabix) {
  # Fresh tabix scratch every run (methylKit errors if a target db exists).
  if (dir.exists(opt$dbdir)) unlink(opt$dbdir, recursive = TRUE, force = TRUE)
  dir.create(opt$dbdir, recursive = TRUE, showWarnings = FALSE)
}

files_all <- sort(list.files(in_dir, pattern = "\\.cov\\.gz$", full.names = TRUE))
if (length(files_all) != 2L * npg)
  stop(sprintf("expected %d .cov.gz in %s, found %d", 2L * npg, in_dir, length(files_all)))
cat(sprintf("backend=%s  cores=%d  mincov=%d\n", opt$dbtype, opt$cores, opt$mincov))
for (f in files_all) cat("  ", basename(f), "\n")

sample_ids <- as.list(c(sprintf("treat_%d", seq_len(npg)), sprintf("ctrl_%d", seq_len(npg))))
treatment  <- c(rep(1L, npg), rep(0L, npg))

# ---- read + unite (backend-specific only in the dbtype/dbdir args) --------
read_args <- list(
  as.list(files_all), sample.id = sample_ids, treatment = treatment,
  assembly = opt$assembly, pipeline = "bismarkCoverage", context = "CpG",
  mincov = opt$mincov
)
if (use_tabix) { read_args$dbtype <- "tabix"; read_args$dbdir <- opt$dbdir }

t_read_0 <- proc.time()
obj <- do.call(methRead, read_args)
t_read <- proc.time() - t_read_0
cat(sprintf("methRead: %d samples\n", length(obj)))

t_unite_0 <- proc.time()
united <- methylKit::unite(obj, destrand = FALSE)  # DB input -> methylBaseDB (save.db default TRUE)
t_unite <- proc.time() - t_unite_0
cat(sprintf("unite: %d sites covered in all samples\n", nrow(united)))

# ---- calculateDiffMeth ----------------------------------------------------
t_diff_0 <- proc.time()
if (identical(opt$overdispersion, "MN")) {
  mk_test <- if (identical(opt$test, "default")) "Chisq" else opt$test
  diff <- calculateDiffMeth(united, overdispersion = "MN", test = mk_test, mc.cores = opt$cores)
} else {
  diff <- calculateDiffMeth(united, mc.cores = opt$cores)
}
t_diff <- proc.time() - t_diff_0
cat(sprintf("calculateDiffMeth: %d sites tested\n", nrow(diff)))

out_df <- getData(diff)
write.table(out_df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out_path, nrow(out_df)))

timing_df <- data.frame(
  phase  = c("read", "unite", "diffmeth"),
  wall_s = c(t_read["elapsed"], t_unite["elapsed"], t_diff["elapsed"]),
  cpu_s  = c(t_read["user.self"] + t_read["sys.self"],
             t_unite["user.self"] + t_unite["sys.self"],
             t_diff["user.self"] + t_diff["sys.self"])
)
write.table(timing_df, file = paste0(out_path, ".timing.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(paste0(out_path, ".sessioninfo.txt"))
cat(sprintf("# run_methylkit_backend.R  dbtype=%s  in_dir=%s  out=%s\n\n",
            opt$dbtype, in_dir, out_path))
print(sessionInfo()); sink()
cat("\nDONE\n")
