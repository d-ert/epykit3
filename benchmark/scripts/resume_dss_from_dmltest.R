#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Resume the DSS replication from the cached dmltest_per_cpg.tsv.gz produced
# by the prior run. Skips the 34-minute DMLfit + DMLtest step.
#
# Steps:
#   1. Load cached DMLtest table.
#   2. Call DSS::callDMR with paper-matched parameters.
#   3. Read the six per-CpG BED files; compute per-DMR per-group mean
#      methylation -> derive diff.Methy.
#   4. HOMER-equivalent annotation via UCSC refGene.
#   5. Build the 100 kb DMR-gene linkage table.
#   6. Write all DSS outputs (parquet not needed for callDMR output; CSV/TSV).
#
# Note: the resume_dss_from_dmltest.py wrapper handles the psutil-based
# resource sampling around this script.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
})

opt_list <- list(
  make_option("--samplesheet",   type = "character"),
  make_option("--out-dir",       type = "character"),
  make_option("--refgene",       type = "character"),
  make_option("--coverage-min",  type = "integer", default = 5L),
  make_option("--p-threshold",   type = "double",  default = 1e-5),
  make_option("--delta",         type = "double",  default = 0),
  make_option("--minlen",        type = "integer", default = 50L),
  make_option("--minCG",         type = "integer", default = 3L),
  make_option("--dis-merge",     type = "integer", default = 100L),
  make_option("--pct-sig",       type = "double",  default = 0.5),
  make_option("--gene-link-bp",  type = "integer", default = 100000L)
)
opts <- parse_args(OptionParser(option_list = opt_list))

OUT <- normalizePath(opts[["out-dir"]], winslash = "/", mustWork = TRUE)
LOG <- file.path(OUT, "resume_log.txt")
sink(LOG, split = TRUE)

cat("Resuming DSS replication from cached DMLtest\n")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Out dir: ", OUT, "\n\n", sep = "")

timings <- list()
record_step <- function(label, expr) {
  cat(sprintf("\n[step] %s ...\n", label))
  invisible(gc(reset = TRUE, full = TRUE))
  t0 <- Sys.time(); c0 <- proc.time()
  val <- eval(expr, envir = parent.frame())
  c1 <- proc.time(); t1 <- Sys.time()
  g <- gc()
  dt_w <- as.numeric(difftime(t1, t0, units = "secs"))
  dt_c <- (c1 - c0)[["user.self"]] + (c1 - c0)[["sys.self"]]
  peak <- sum(g[, "max used"])
  timings[[length(timings) + 1L]] <<- data.table(
    step = label, wall_seconds = round(dt_w, 3),
    total_cpu_seconds = round(dt_c, 3),
    r_mem_peak_mb = round(peak, 1)
  )
  cat(sprintf("       %.2fs wall / %.2fs CPU / R-mem peak %.0f MB\n",
              dt_w, dt_c, peak))
  val
}

# ---- 1. Load cached DMLtest ------------------------------------------------

test_path <- file.path(OUT, "dmltest_per_cpg.tsv.gz")
stopifnot(file.exists(test_path))
test_res <- record_step("load_dmltest_cache", quote({
  dt <- fread(test_path, sep = "\t")
  cat(sprintf("       loaded %s rows, %d cols: %s\n",
              format(nrow(dt), big.mark = ","), ncol(dt),
              paste(colnames(dt), collapse = ", ")))
  # callDMR wants a data.frame with the DSS schema. Keep it.
  as.data.frame(dt)
}))

# ---- 2. callDMR (multi-factor result returns chr,start,end,length,nCG,areaStat) ---

dmrs <- record_step("callDMR", quote({
  callDMR(test_res,
          delta       = opts$delta,
          p.threshold = opts[["p-threshold"]],
          minlen      = opts$minlen,
          minCG       = opts$minCG,
          dis.merge   = opts[["dis-merge"]],
          pct.sig     = opts[["pct-sig"]])
}))
stopifnot(!is.null(dmrs), nrow(dmrs) > 0L)
dmrs <- as.data.table(dmrs)
setnames(dmrs, "chr", "chrom")
setnames(dmrs, "nCG", "n_cpgs", skip_absent = TRUE)
cat(sprintf("\ncallDMR -> %d DMRs; columns: %s\n",
            nrow(dmrs), paste(colnames(dmrs), collapse = ", ")))

fwrite(dmrs, file.path(OUT, "dmr_dss_raw.tsv"), sep = "\t")

# callDMR's diff.Methy on the multifactor path is derived from the
# smoothed/fit model — call it diff_Methy_DSSfit for clarity. We still
# compute the from-raw-counts version below to match the paper's Table 5
# convention (per-sample mean within each DMR, then averaged within group).
stopifnot("diff.Methy" %in% colnames(dmrs))
setnames(dmrs, "diff.Methy", "diff_Methy_DSSfit")
if ("meanMethy1" %in% colnames(dmrs)) setnames(dmrs, "meanMethy1", "meanMethy1_DSSfit")
if ("meanMethy2" %in% colnames(dmrs)) setnames(dmrs, "meanMethy2", "meanMethy2_DSSfit")
dmrs[, dmr_id := seq_len(.N)]

# ---- 3. Compute per-DMR per-group mean methylation -------------------------
#
# Multi-factor callDMR doesn't return meanMethy1/2 or diff.Methy because the
# concept is only well-defined for two-group tests. The paper's Table 5 has
# diff.meth_mean — they derived it themselves from the BSseq object. We do
# the same: for each DMR, average (sum M / sum T) per sample across the CpGs
# inside the DMR, then average across samples within each group.
#
# Sample <-> group mapping from samplesheet.

samplesheet <- fread(opts$samplesheet)
samplesheet[, group := factor(group)]
stopifnot("group" %in% colnames(samplesheet),
          "path"  %in% colnames(samplesheet),
          "sample_id" %in% colnames(samplesheet))

# Read each sample's per-CpG counts (cols 10/11 are combined-strand M, T).
# We index by (chrom, pos) and only keep CpGs that fall inside any DMR.

# Build a set of (chrom, pos) we care about by expanding each DMR to the
# CpG positions it contains. We don't know the CpG positions a priori, so
# instead we'll do the per-sample work region-by-region after loading.

# Build a DMR interval table for foverlaps (interval join). DSS callDMR
# coords are 1-based inclusive; BED is 0-based start, 1-based end. We
# convert the DMR interval to BED's coordinate convention so a BED row
# hits DMR i when bed.start in [DMR.start - 1, DMR.end - 1].
dmr_intervals <- dmrs[, .(chrom = as.character(chrom),
                           start = as.integer(start) - 1L,
                           end   = as.integer(end),
                           dmr_id = dmr_id)]
setkey(dmr_intervals, chrom, start, end)

compute_means_for_sample <- function(path, sample_id) {
  cat(sprintf("       reading %s\n", sample_id))
  bed <- fread(path, header = FALSE, select = c(1L, 2L, 10L, 11L),
               col.names = c("chrom", "start", "M", "T"))
  bed[, chrom := as.character(chrom)]
  bed[, start := as.integer(start)]
  bed <- bed[T >= opts[["coverage-min"]]]
  bed[, end := start + 1L]
  setkey(bed, chrom, start, end)
  # foverlaps semantics: rows in x that overlap intervals in y (keyed).
  ov <- foverlaps(bed, dmr_intervals, nomatch = NULL,
                  type = "any", which = FALSE)
  if (nrow(ov) == 0L) {
    return(data.table(dmr_id = integer(), sum_M = numeric(),
                      sum_T = numeric(), n_cpg = integer()))
  }
  ov[, .(sum_M = sum(as.numeric(M)),
         sum_T = sum(as.numeric(T)),
         n_cpg = .N),
     by = dmr_id]
}

per_sample_means <- record_step("per_sample_dmr_means", quote({
  lst <- vector("list", nrow(samplesheet))
  for (i in seq_len(nrow(samplesheet))) {
    row <- samplesheet[i]
    res <- compute_means_for_sample(row$path, row$sample_id)
    res[, sample_id := row$sample_id]
    res[, group := as.character(row$group)]
    lst[[i]] <- res
  }
  rbindlist(lst)
}))

# Per-DMR per-sample beta = sum_M / sum_T
per_sample_means[, beta := sum_M / pmax(sum_T, 1)]

# Wide per-sample table (one column per sample) — matches paper Supp Table 5
sample_wide <- dcast(per_sample_means, dmr_id ~ sample_id, value.var = "beta")
sample_cols <- setdiff(colnames(sample_wide), "dmr_id")
setnames(sample_wide, sample_cols, paste0("beta_", sample_cols))

# Per-group mean across samples
group_means <- per_sample_means[, .(beta_group = mean(beta, na.rm = TRUE)),
                                 by = .(dmr_id, group)]
group_wide <- dcast(group_means, dmr_id ~ group, value.var = "beta_group")
groups <- levels(samplesheet$group)
cat("Group levels (in order):", paste(groups, collapse = ", "), "\n")
treatment_level <- groups[!grepl("^WT$", groups)][1L]
control_level   <- "WT"
stopifnot(treatment_level %in% colnames(group_wide),
          control_level   %in% colnames(group_wide))
setnames(group_wide,
         old = c(treatment_level, control_level),
         new = c("meanMethy_treatment", "meanMethy_control"))
group_wide[, diff_Methy_fromCounts := meanMethy_treatment - meanMethy_control]

# Merge per-sample + per-group methylation back onto the DMR table by dmr_id.
dmrs <- merge(dmrs, group_wide, by = "dmr_id", all.x = TRUE, sort = FALSE)
dmrs <- merge(dmrs, sample_wide, by = "dmr_id", all.x = TRUE, sort = FALSE)

# Primary direction: use the from-counts diff (paper-equivalent). Sign of
# diff_Methy_DSSfit (model-based) is cross-checked below.
dmrs[, dmr_type := fifelse(diff_Methy_fromCounts > 0, "hyper",
                   fifelse(diff_Methy_fromCounts < 0, "hypo", "none"))]

n_h <- sum(dmrs$dmr_type == "hyper"); n_o <- sum(dmrs$dmr_type == "hypo")
n_disagree <- sum(sign(dmrs$diff_Methy_fromCounts) != sign(dmrs$diff_Methy_DSSfit),
                  na.rm = TRUE)
cat(sprintf("\nDMRs after per-sample mean derivation: %d (%d hyper, %d hypo)\n",
            nrow(dmrs), n_h, n_o))
cat(sprintf("Sign agreement DSSfit vs fromCounts: %d / %d disagree\n",
            n_disagree, nrow(dmrs)))

# ---- 4. HOMER-equivalent annotation (UCSC refGene) -------------------------

annotate_homer_refseq <- function(dmrs, refgene_path) {
  cat(sprintf("\nLoading refGene catalog: %s\n", refgene_path))
  rg <- fread(refgene_path, header = FALSE)
  setnames(rg, paste0("V", 1:ncol(rg)),
           c("bin","acc","chrom","strand","txStart","txEnd",
             "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
             "score","gene","cdsStartStat","cdsEndStat","exonFrames")[1:ncol(rg)])
  rg <- rg[grepl("^chr", chrom) & !grepl("_", chrom)]
  rg[, tss := ifelse(strand == "+", txStart, txEnd)]
  rg[, tts := ifelse(strand == "+", txEnd,   txStart)]
  rg[, coding := grepl("^NM_", acc)]

  PROMOTER_UP   <- 1000L; PROMOTER_DOWN <- 100L
  TTS_UP        <- 100L;  TTS_DOWN      <- 1000L
  FEAT_PRIO <- c("promoter-TSS" = 0, "TTS" = 1, "5UTR" = 2, "3UTR" = 3,
                 "exon" = 4, "intron" = 5, "non-coding" = 6, "intergenic" = 7)

  classify_one <- function(center, chrom_) {
    cands <- rg[chrom == chrom_]
    if (nrow(cands) == 0L) return(c("intergenic", ""))
    best_prio <- 8L; best_feat <- "intergenic"; best_gene <- ""
    for (i in seq_len(nrow(cands))) {
      g <- cands[i]
      if (g$strand == "+") {
        prom_lo <- g$tss - PROMOTER_UP;   prom_hi <- g$tss + PROMOTER_DOWN
        tts_lo  <- g$tts - TTS_UP;        tts_hi  <- g$tts + TTS_DOWN
      } else {
        prom_lo <- g$tss - PROMOTER_DOWN; prom_hi <- g$tss + PROMOTER_UP
        tts_lo  <- g$tts - TTS_DOWN;      tts_hi  <- g$tts + TTS_UP
      }
      this_feat <- NA_character_
      if (center >= prom_lo && center <= prom_hi) {
        this_feat <- "promoter-TSS"
      } else if (center >= tts_lo && center <= tts_hi) {
        this_feat <- "TTS"
      } else if (center >= g$txStart && center <= g$txEnd) {
        if (!g$coding) {
          this_feat <- "non-coding"
        } else {
          exonStarts <- as.integer(strsplit(sub(",$", "", g$exonStarts), ",")[[1]])
          exonEnds   <- as.integer(strsplit(sub(",$", "", g$exonEnds),   ",")[[1]])
          in_exon <- any(center >= exonStarts & center <= exonEnds)
          if (!in_exon) {
            this_feat <- "intron"
          } else if (g$cdsStart == g$cdsEnd) {
            this_feat <- "non-coding"
          } else if (g$strand == "+") {
            if      (center <  g$cdsStart) this_feat <- "5UTR"
            else if (center >= g$cdsEnd)   this_feat <- "3UTR"
            else                            this_feat <- "exon"
          } else {
            if      (center >= g$cdsEnd)   this_feat <- "5UTR"
            else if (center <  g$cdsStart) this_feat <- "3UTR"
            else                            this_feat <- "exon"
          }
        }
      }
      if (!is.na(this_feat)) {
        p <- FEAT_PRIO[[this_feat]]
        if (p < best_prio) {
          best_prio <- p; best_feat <- this_feat; best_gene <- g$gene
        }
      }
    }
    c(best_feat, best_gene)
  }

  rg_by_chr <- split(rg, by = "chrom")
  for (ch in names(rg_by_chr)) setorder(rg_by_chr[[ch]], tss)
  nearest_tss <- function(center, chrom_) {
    cands <- rg_by_chr[[chrom_]]
    if (is.null(cands) || nrow(cands) == 0L) return(list(gene = "", dist = NA_integer_))
    tsspos <- cands$tss
    i <- findInterval(center, tsspos)
    pick <- c()
    if (i >= 1L)                  pick <- c(pick, i)
    if (i + 1L <= length(tsspos)) pick <- c(pick, i + 1L)
    if (length(pick) == 0L)       return(list(gene = "", dist = NA_integer_))
    pick <- pick[which.min(abs(tsspos[pick] - center))]
    g    <- cands[pick]
    d    <- center - g$tss
    if (g$strand == "-") d <- -d
    list(gene = g$gene, dist = as.integer(d))
  }

  cat(sprintf("Annotating %d DMRs...\n", nrow(dmrs)))
  feat <- character(nrow(dmrs)); feat_gene <- character(nrow(dmrs))
  near_gene <- character(nrow(dmrs)); near_dist <- integer(nrow(dmrs))
  for (i in seq_len(nrow(dmrs))) {
    center <- as.integer((dmrs$start[i] + dmrs$end[i]) %/% 2L)
    ff <- classify_one(center, dmrs$chrom[i])
    feat[i] <- ff[1]; feat_gene[i] <- ff[2]
    nn <- nearest_tss(center, dmrs$chrom[i])
    near_gene[i] <- nn$gene
    near_dist[i] <- if (is.null(nn$dist) || is.na(nn$dist)) NA_integer_ else nn$dist
  }
  dmrs[, `:=`(feature_type = feat,
              feature_gene = feat_gene,
              nearest_tss_gene = near_gene,
              nearest_tss_distance = near_dist)]
  list(dmrs = dmrs, rg_by_chr = rg_by_chr)
}

ann <- record_step("annotate_refgene_homer", quote({
  annotate_homer_refseq(dmrs, opts$refgene)
}))
dmrs   <- ann$dmrs
rg_chr <- ann$rg_by_chr

fwrite(dmrs, file.path(OUT, "dmr_dss.csv"))
fwrite(dmrs, file.path(OUT, "dmr_dss.tsv"), sep = "\t")

# ---- 5. 100 kb DMR-gene linkage --------------------------------------------

build_100kb_links <- function(dmrs, rg_by_chr, max_bp) {
  rg_all <- rbindlist(rg_by_chr)
  canon <- rg_all[, .(tss = ifelse(strand[1] == "+", min(tss), max(tss)),
                      strand = strand[1]),
                  by = .(chrom, gene)]
  canon_by_chr <- split(canon, by = "chrom")
  for (ch in names(canon_by_chr)) setorder(canon_by_chr[[ch]], tss)

  out <- vector("list", nrow(dmrs))
  for (i in seq_len(nrow(dmrs))) {
    ch <- dmrs$chrom[i]
    if (is.null(canon_by_chr[[ch]])) next
    midpoint <- as.integer((dmrs$start[i] + dmrs$end[i]) %/% 2L)
    cands <- canon_by_chr[[ch]]
    keep  <- abs(cands$tss - midpoint) <= max_bp
    if (!any(keep)) next
    sub <- cands[keep]
    out[[i]] <- data.table(
      dmr_index                = i - 1L,
      chrom                    = ch,
      dmr_start                = dmrs$start[i],
      dmr_end                  = dmrs$end[i],
      dmr_midpoint             = midpoint,
      dmr_type                 = dmrs$dmr_type[i],
      mean_meth_diff           = dmrs$diff_Methy_fromCounts[i],
      gene                     = sub$gene,
      tss                      = sub$tss,
      distance_tss_to_midpoint = sub$tss - midpoint,
      abs_distance             = abs(sub$tss - midpoint)
    )
  }
  out <- rbindlist(out)
  out <- unique(out, by = c("dmr_index", "gene"))
  setorder(out, dmr_index, abs_distance)
  out
}

links <- record_step("build_100kb_gene_links", quote({
  build_100kb_links(dmrs, rg_chr, opts[["gene-link-bp"]])
}))
fwrite(links, file.path(OUT, "dmr_gene_links_100kb.csv"))
cat(sprintf("Wrote dmr_gene_links_100kb.csv  (%d rows, %d unique genes)\n",
            nrow(links), uniqueN(links$gene)))

# ---- 6. Step timings + session info ---------------------------------------

timing_dt <- rbindlist(timings)
fwrite(timing_dt, file.path(OUT, "step_timings_resume.tsv"), sep = "\t")
cat("\nstep_timings_resume.tsv:\n"); print(timing_dt)

writeLines(capture.output(sessionInfo()),
           file.path(OUT, "dss_session_info.txt"))

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
sink()
