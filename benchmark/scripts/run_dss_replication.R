#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# DSS replication of GSE263850 with paper-matched parameters.
#
# Paper (Farhangdoost et al., 2024), Methods, DMR section:
#
#   "Differential methylation analyses were performed using the R package
#    DSS (DMLfit.multiFactor function with parameter smoothing = TRUE).
#    DMRs were called using the function callDMR at default parameters
#    with a minimum length of 50 base pairs and 3 CpG sites (delta = 0,
#    p.threshold = 1e-5, minlen = 50, minCG = 3, dis.merge = 100,
#    pct.sig = 0.5). Annotation determined using Homer (hg38). [...]
#    For DMR-gene expression correlations, we only included genes
#    associated within a 100 kb distance from the DMR (distance of 100 kb
#    from the TSS to the middle of the DMR)."
#
# This script implements that call set exactly.
#
# Outputs (written into --out-dir):
#   dmr_dss_raw.tsv                   - raw callDMR output (one row per DMR)
#   dmr_dss.parquet / .csv            - DMR table with HOMER-style refGene annotation
#   dmr_gene_links_100kb.csv          - long-form (DMR, gene) pairs within 100 kb
#   step_timings.tsv                  - per-step wall + CPU + R-side memory
#   dss_session_info.txt              - sessionInfo() for reproducibility
#   run_log.txt                       - stdout/stderr of this script
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(bsseq)
  library(DSS)
})

# ---- CLI -------------------------------------------------------------------

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
parser <- OptionParser(option_list = opt_list)
opts <- parse_args(parser)

stopifnot(!is.null(opts$samplesheet),
          !is.null(opts[["out-dir"]]),
          !is.null(opts$refgene))

OUT  <- normalizePath(opts[["out-dir"]], winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(OUT, "run_log.txt")
sink(LOG, split = TRUE)   # also echoes to console

cat("DSS replication of GSE263850 — paper-matched parameters\n")
cat(strrep("=", 72), "\n", sep = "")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Out dir: ", OUT, "\n", sep = "")

# ---- Step timing helper ----------------------------------------------------

timings <- list()
record_step <- function(label, expr) {
  cat(sprintf("\n[step] %s ...\n", label))
  invisible(gc(reset = TRUE, full = TRUE))
  t0_wall <- Sys.time()
  t0_cpu  <- proc.time()
  val <- eval(expr, envir = parent.frame())
  t1_cpu  <- proc.time()
  t1_wall <- Sys.time()
  gcinfo  <- gc()  # captures max used in MB
  dt_wall <- as.numeric(difftime(t1_wall, t0_wall, units = "secs"))
  dt_user <- (t1_cpu - t0_cpu)[["user.self"]]
  dt_sys  <- (t1_cpu - t0_cpu)[["sys.self"]]
  peak_mb <- sum(gcinfo[, "max used"])  # MB across Vcells+Ncells
  used_mb <- sum(gcinfo[, "used"])
  timings[[length(timings) + 1L]] <<- data.table(
    step              = label,
    wall_seconds      = round(dt_wall, 3),
    user_cpu_seconds  = round(dt_user, 3),
    sys_cpu_seconds   = round(dt_sys, 3),
    total_cpu_seconds = round(dt_user + dt_sys, 3),
    cpu_pct_of_wall   = round(100 * (dt_user + dt_sys) / max(dt_wall, 1e-9), 1),
    r_mem_used_mb     = round(used_mb, 1),
    r_mem_peak_mb     = round(peak_mb, 1)
  )
  cat(sprintf("       %.2fs wall / %.2fs CPU / R-mem peak %.0f MB\n",
              dt_wall, dt_user + dt_sys, peak_mb))
  val
}

# ---- 1. Load samplesheet + read BEDs ---------------------------------------

samplesheet <- fread(opts$samplesheet)
cat("\nSamplesheet:\n"); print(samplesheet)
stopifnot(all(c("sample_id", "group", "path") %in% colnames(samplesheet)))

samplesheet[, group := factor(group)]
sample_ids <- samplesheet$sample_id
n_samples  <- length(sample_ids)

# 12-col combined-strand BED: chrom start end fwd_M fwd_T fwd_pct
# rev_M rev_T rev_pct M T pct      (cols 10/11 = combined-strand M/T)
read_one <- function(path) {
  dt <- fread(path, header = FALSE, select = c(1L, 2L, 10L, 11L),
              col.names = c("chrom", "start", "M", "T"))
  dt[T >= opts[["coverage-min"]]]
}

per_sample <- record_step("ingest_beds", quote({
  lst <- vector("list", n_samples)
  names(lst) <- sample_ids
  for (i in seq_len(n_samples)) {
    cat(sprintf("       [%d/%d] %s\n", i, n_samples, sample_ids[i]))
    lst[[i]] <- read_one(samplesheet$path[i])
  }
  lst
}))

# ---- 2. Intersect CpGs across samples + build BSseq ------------------------

bsseq_obj <- record_step("intersect_and_bsseq", quote({
  # Inner-join across samples on (chrom, start). Result: one row per CpG
  # present in all 6 samples with coverage >= coverage_min.
  key_cols <- c("chrom", "start")
  joined <- per_sample[[1L]][, .(chrom, start)]
  for (i in 2:n_samples) {
    joined <- merge(joined, per_sample[[i]][, .(chrom, start)],
                    by = key_cols)
  }
  setkeyv(joined, key_cols)
  cat(sprintf("       intersected CpGs: %s across %d samples\n",
              format(nrow(joined), big.mark = ","), n_samples))

  # Build M and T matrices aligned to `joined` order
  n_cpg <- nrow(joined)
  M_mat <- matrix(0L, nrow = n_cpg, ncol = n_samples,
                  dimnames = list(NULL, sample_ids))
  T_mat <- matrix(0L, nrow = n_cpg, ncol = n_samples,
                  dimnames = list(NULL, sample_ids))
  for (i in seq_len(n_samples)) {
    sub <- per_sample[[i]][joined, on = key_cols, nomatch = NULL]
    M_mat[, i] <- as.integer(sub$M)
    T_mat[, i] <- as.integer(sub$T)
  }
  # Free the per-sample list to reclaim RAM before smoothing.
  rm(per_sample, envir = parent.frame())
  invisible(gc(full = TRUE))

  BSseq(chr  = joined$chrom,
        pos  = joined$start + 1L,   # BSseq is 1-based; BED start is 0-based
        M    = M_mat,
        Cov  = T_mat,
        sampleNames = sample_ids)
}))

# pData with group factor (paper's WT vs Het_AKAP11_KO)
pData(bsseq_obj)$group <- factor(
  samplesheet$group[match(sampleNames(bsseq_obj), samplesheet$sample_id)]
)
cat("\nBSseq object built:\n"); print(bsseq_obj)
cat("Group levels:\n"); print(table(pData(bsseq_obj)$group))

# ---- 3. DMLfit.multiFactor(smoothing = TRUE) -------------------------------

design <- data.frame(group = pData(bsseq_obj)$group)
formula_str <- "~ group"
cat(sprintf("\nDesign (formula = %s):\n", formula_str)); print(design)

fit <- record_step("DMLfit.multiFactor_smooth", quote({
  DMLfit.multiFactor(bsseq_obj,
                     design   = design,
                     formula  = as.formula(formula_str),
                     smoothing = TRUE)
}))

# Which coefficient does "group" map to? Take the last fitted coefficient
# (the non-intercept term, i.e. the contrast for the second factor level).
coef_names <- colnames(fit$X)
cat("Design matrix coefficients:", paste(coef_names, collapse = ", "), "\n")
contrast_coef <- tail(coef_names, 1L)
cat("Testing coefficient:", contrast_coef, "\n")

# ---- 4. DMLtest.multiFactor -----------------------------------------------

test_res <- record_step("DMLtest.multiFactor", quote({
  DMLtest.multiFactor(fit, coef = contrast_coef)
}))
cat(sprintf("\nDMLtest sites: %s\n", format(nrow(test_res), big.mark = ",")))
n_sig_p <- sum(test_res$pvals < opts[["p-threshold"]], na.rm = TRUE)
n_sig_q <- sum(test_res$fdrs  < 0.05, na.rm = TRUE)
cat(sprintf("  sig at p < %g: %s\n", opts[["p-threshold"]],
            format(n_sig_p, big.mark = ",")))
cat(sprintf("  sig at FDR < 0.05: %s\n", format(n_sig_q, big.mark = ",")))

# Save the per-CpG test table for downstream use.
fwrite(as.data.table(test_res),
       file.path(OUT, "dmltest_per_cpg.tsv.gz"),
       sep = "\t", compress = "gzip")
cat("Wrote dmltest_per_cpg.tsv.gz\n")

# ---- 5. callDMR with paper-matched parameters ------------------------------

dmrs <- record_step("callDMR", quote({
  callDMR(test_res,
          delta       = opts$delta,
          p.threshold = opts[["p-threshold"]],
          minlen      = opts$minlen,
          minCG       = opts$minCG,
          dis.merge   = opts[["dis-merge"]],
          pct.sig     = opts[["pct-sig"]])
}))

if (is.null(dmrs) || nrow(dmrs) == 0L) {
  stop("DSS callDMR returned no DMRs — check inputs / parameters.")
}
dmrs <- as.data.table(dmrs)
setnames(dmrs, old = c("nCG"), new = c("n_cpgs"), skip_absent = TRUE)
# Standard DSS callDMR columns: chr, start, end, length, nCG, meanMethy1,
#   meanMethy2, diff.Methy, areaStat. We keep all + derive dmr_type.
dmrs[, dmr_type := ifelse(diff.Methy > 0, "hyper",
                   ifelse(diff.Methy < 0, "hypo", "none"))]
setnames(dmrs, "chr", "chrom")
fwrite(dmrs, file.path(OUT, "dmr_dss_raw.tsv"), sep = "\t")
cat(sprintf("DMRs called: %d  (hyper: %d, hypo: %d)\n",
            nrow(dmrs),
            sum(dmrs$dmr_type == "hyper"),
            sum(dmrs$dmr_type == "hypo")))

# ---- 6. HOMER-equivalent annotation (refGene nearest TSS + feature) --------

annotate_homer_refseq <- function(dmrs, refgene_path) {
  cat(sprintf("\nLoading refGene catalog: %s\n", refgene_path))
  rg <- fread(refgene_path, header = FALSE)
  # UCSC refGene fields: bin name chrom strand txStart txEnd cdsStart
  #   cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat
  #   cdsEndStat exonFrames
  setnames(rg, paste0("V", 1:ncol(rg)),
           c("bin","acc","chrom","strand","txStart","txEnd",
             "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
             "score","gene","cdsStartStat","cdsEndStat","exonFrames")[1:ncol(rg)])
  rg <- rg[grepl("^chr", chrom) & !grepl("_", chrom)]
  rg[, tss := ifelse(strand == "+", txStart, txEnd)]
  rg[, tts := ifelse(strand == "+", txEnd, txStart)]
  rg[, coding := grepl("^NM_", acc)]

  # ---- Feature classification (HOMER priorities) ----
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
        prom_lo <- g$tss - PROMOTER_UP; prom_hi <- g$tss + PROMOTER_DOWN
        tts_lo  <- g$tts - TTS_UP;      tts_hi  <- g$tts + TTS_DOWN
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

  # Nearest TSS per chromosome (signed distance, strand-aware)
  rg_by_chr <- split(rg, by = "chrom")
  for (ch in names(rg_by_chr)) {
    setorder(rg_by_chr[[ch]], tss)
  }
  nearest_tss <- function(center, chrom_) {
    cands <- rg_by_chr[[chrom_]]
    if (is.null(cands) || nrow(cands) == 0L) return(list(gene = "", dist = NA_integer_))
    tsspos <- cands$tss
    i <- findInterval(center, tsspos)
    pick <- c()
    if (i >= 1L)              pick <- c(pick, i)
    if (i + 1L <= length(tsspos)) pick <- c(pick, i + 1L)
    if (length(pick) == 0L)   return(list(gene = "", dist = NA_integer_))
    pick <- pick[which.min(abs(tsspos[pick] - center))]
    g    <- cands[pick]
    d    <- center - g$tss
    if (g$strand == "-") d <- -d
    list(gene = g$gene, dist = as.integer(d))
  }

  cat(sprintf("Annotating %d DMRs (HOMER refGene rules)...\n", nrow(dmrs)))
  feat <- character(nrow(dmrs)); feat_gene <- character(nrow(dmrs))
  near_gene <- character(nrow(dmrs)); near_dist <- integer(nrow(dmrs))
  for (i in seq_len(nrow(dmrs))) {
    center <- as.integer((dmrs$start[i] + dmrs$end[i]) %/% 2L)
    ch     <- dmrs$chrom[i]
    ff <- classify_one(center, ch)
    feat[i] <- ff[1]; feat_gene[i] <- ff[2]
    nn <- nearest_tss(center, ch)
    near_gene[i] <- nn$gene
    near_dist[i] <- if (is.null(nn$dist) || is.na(nn$dist)) NA_integer_ else nn$dist
  }
  dmrs[, `:=`(feature_type = feat,
              feature_gene = feat_gene,
              nearest_tss_gene = near_gene,
              nearest_tss_distance = near_dist)]
  list(dmrs = dmrs, rg = rg, rg_by_chr = rg_by_chr)
}

ann <- record_step("annotate_refgene_homer", quote({
  annotate_homer_refseq(dmrs, opts$refgene)
}))
dmrs    <- ann$dmrs
rg      <- ann$rg
rg_chr  <- ann$rg_by_chr

fwrite(dmrs, file.path(OUT, "dmr_dss.csv"))
fwrite(dmrs, file.path(OUT, "dmr_dss.tsv"), sep = "\t")
cat("Wrote dmr_dss.csv / dmr_dss.tsv\n")

# ---- 7. 100 kb DMR–gene linkage table --------------------------------------

build_100kb_links <- function(dmrs, rg_by_chr, max_bp) {
  # Collapse to canonical TSS per gene per chromosome (most upstream).
  rg_all <- rbindlist(rg_by_chr)
  canon <- rg_all[, .(tss = ifelse(strand[1] == "+", min(tss), max(tss)),
                      strand = strand[1]),
                  by = .(chrom, gene)]
  canon_by_chr <- split(canon, by = "chrom")
  for (ch in names(canon_by_chr)) {
    setorder(canon_by_chr[[ch]], tss)
  }

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
      dmr_index                = i - 1L,           # 0-based for parity w/ chain_merge
      chrom                    = ch,
      dmr_start                = dmrs$start[i],
      dmr_end                  = dmrs$end[i],
      dmr_midpoint             = midpoint,
      dmr_type                 = dmrs$dmr_type[i],
      mean_meth_diff           = dmrs$diff.Methy[i],
      gene                     = sub$gene,
      tss                      = sub$tss,
      distance_tss_to_midpoint = sub$tss - midpoint,
      abs_distance             = abs(sub$tss - midpoint)
    )
  }
  out <- rbindlist(out)
  # Dedup (gene-level)
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

# ---- 8. Write step_timings + session info ---------------------------------

timing_dt <- rbindlist(timings)
fwrite(timing_dt, file.path(OUT, "step_timings.tsv"), sep = "\t")
cat("\nstep_timings.tsv:\n"); print(timing_dt)

writeLines(capture.output(sessionInfo()),
           file.path(OUT, "dss_session_info.txt"))
cat("\nWrote dss_session_info.txt\n")

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
sink()
