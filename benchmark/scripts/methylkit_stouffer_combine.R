#!/usr/bin/env Rscript
# methylkit_stouffer_combine.R
#
# Apply adjacent-CpG Stouffer combination to a methylKit per-scenario
# TSV. Mirrors epykit's neighbour_combine so Phase 4 head-to-head is
# tuned-vs-tuned per PROTOCOL R1.
#
# Usage:
#   Rscript methylkit_stouffer_combine.R \
#     --in <tsv> --out <tsv> [--max-gap-bp 1000] [--window 3]

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--in",         type = "character", help = "Input TSV path"),
  make_option("--out",        type = "character", help = "Output TSV path"),
  make_option("--max-gap-bp", type = "integer",   default = 1000L),
  make_option("--window",     type = "integer",   default = 3L)
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$`in`) || is.null(opt$out)) stop("--in and --out are required")

df <- fread(opt$`in`, sep = "\t", header = TRUE)
stopifnot(all(c("chr", "start", "pvalue") %in% names(df)))

stouffer <- function(pv) {
  pv <- pmin(pmax(pv, 1e-300), 1 - 1e-15)
  z  <- qnorm(1 - pv / 2)
  z_comb <- sum(z) / sqrt(length(z))
  2 * (1 - pnorm(abs(z_comb)))
}

df[, pvalue_combined := NA_real_]
setkey(df, chr, start)

for (chrom in unique(df$chr)) {
  rows  <- df[chr == chrom]
  pos   <- rows$start
  pv    <- rows$pvalue
  md    <- if ("meth.diff" %in% names(rows)) rows$`meth.diff` else rep(0, nrow(rows))
  half  <- (opt$window - 1L) %/% 2L
  comb  <- numeric(nrow(rows))

  for (i in seq_along(pv)) {
    lo <- max(1L, i - half)
    hi <- min(length(pv), i + half)
    # Keep only CpGs within max-gap-bp of the focal.
    in_gap <- which(abs(pos[lo:hi] - pos[i]) <= opt$`max-gap-bp`) + lo - 1L
    if (length(in_gap) < 2L) { comb[i] <- pv[i]; next }
    # Direction check: skip mixed-direction windows (non-zero signs must agree).
    signs <- sign(md[in_gap])
    nonzero_signs <- signs[signs != 0]
    if (length(unique(nonzero_signs)) > 1L) { comb[i] <- pv[i]; next }
    comb[i] <- stouffer(pv[in_gap])
  }
  df[chr == chrom, pvalue_combined := comb]
}

# BH per chromosome.
df[, qvalue_combined := NA_real_]
for (chrom in unique(df$chr)) {
  m  <- df$chr == chrom & !is.na(df$pvalue_combined)
  if (any(m)) df[m, qvalue_combined := p.adjust(pvalue_combined, method = "BH")]
}

fwrite(df, opt$out, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d rows)\n", opt$out, nrow(df)))
