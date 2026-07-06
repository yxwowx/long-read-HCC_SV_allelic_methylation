#!/usr/bin/env Rscript
# const_amr_threshold_sensitivity.R
# Priority 1 + 2 analyses (GRADE gap remediation):
#   1. Constitutional AMR overlap threshold sensitivity (0–50% width overlap)
#      → Does C14 somatic aDMR OR stay null regardless of threshold?
#   2. Normal tissue aDMR statistics (per-patient counts, overlap fractions)
#
# Outputs:
#   SV_aDMR/const_amr_threshold_or.csv     — OR table per threshold
#   SV_aDMR/normal_admr_stats.csv          — per-patient normal/tumor/const/somatic aDMR counts

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(stringr)
})

UTILS_DIR <- "/home/kachungk/script/SV-DMR/shared_file/pipeline"
source(file.path(UTILS_DIR, "shared_utils.R"))

DATA_ROOT    <- "/node200data/kachungk/hcc_data"
DMR_DIR      <- file.path(DATA_ROOT, "DMR_minimap2.out_hg38/DSS")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")

REF_ROOT   <- "/node200data/kachungk/reference/GRCh38"
SEGDUP_BED <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD_BED    <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW     <- file.path(REF_ROOT, "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")

MAIN_CHROMS  <- paste0("chr", c(1:22, "X"))
CTRL_MULT    <- 10L
THRESHOLDS   <- c(0, 0.05, 0.10, 0.25, 0.50)  # width overlap fractions

# hg38 main chromosome lengths (GRCh38.p13)
CHROM_LENS <- setNames(
  as.list(c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
            170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
            135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
            90338345L, 83257441L, 80373285L, 58617616L, 64444167L,
            46709983L, 50818468L, 156040895L)),
  c(paste0("chr", 1:22), "chrX")
)

pmap <- fread(PATIENT_MAP_PATH)
patients <- unique(pmap$patient_code)

# ── Reference features ────────────────────────────────────────────────────────
cat("[0] Loading reference features...\n")
segdup_gr <- import.bed(SEGDUP_BED) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
lad_gr    <- import.bed(LAD_BED)    |> keepStandardChromosomes(pruning.mode = "coarse")
lad_gr    <- lad_gr[seqnames(lad_gr) %in% MAIN_CHROMS]
has_pc1   <- file.exists(PC1_BW)
if (has_pc1) pc1_rle <- import(PC1_BW, as = "RleList")

annotate_gr <- function(gr) {
  mcols(gr)$segdup <- gr %over% segdup_gr
  mcols(gr)$lad    <- gr %over% lad_gr
  mcols(gr)$b_compartment <- NA_real_
  if (has_pc1) {
    shared <- intersect(seqlevels(gr), names(pc1_rle))
    if (length(shared)) {
      gs  <- keepSeqlevels(gr, shared, pruning.mode = "coarse")
      ps  <- pc1_rle[seqlevels(gs)]
      sc  <- binnedAverage(gs, ps, varname = "PC1")
      idx <- which(as.character(seqnames(gr)) %in% shared)
      mcols(gr)$b_compartment[idx] <- mcols(sc)$PC1 < 0
    }
  }
  gr
}

make_controls <- function(case_gr, n_mult, seed) {
  set.seed(seed)
  ctrl_list <- lapply(unique(as.character(seqnames(case_gr))), function(chr) {
    cc  <- case_gr[seqnames(case_gr) == chr]
    if (length(cc) == 0) return(NULL)
    chr_len <- CHROM_LENS[[chr]]
    if (is.null(chr_len)) return(NULL)
    ww  <- width(cc)
    starts <- unlist(lapply(ww, function(w) {
      ms <- chr_len - w; if (ms < 1) return(NULL)
      sample.int(ms, n_mult, replace = TRUE)
    }))
    if (!length(starts)) return(NULL)
    GRanges(seqnames = chr,
            ranges   = IRanges(starts, starts + rep(ww, each = n_mult) - 1L),
            is_case  = 0L, obs_weight = 1 / n_mult)
  })
  do.call(c, Filter(Negate(is.null), ctrl_list))
}

run_glm_gr <- function(case_gr, threshold) {
  case_gr <- annotate_gr(case_gr)
  ctrl_gr <- annotate_gr(make_controls(case_gr, CTRL_MULT, 42L))
  mcols(case_gr)$is_case    <- 1L
  mcols(case_gr)$obs_weight <- 1
  mcols(ctrl_gr)$is_case    <- 0L
  df <- rbind(
    as.data.frame(mcols(case_gr)[, c("is_case","segdup","lad","b_compartment","obs_weight")]),
    as.data.frame(mcols(ctrl_gr)[, c("is_case","segdup","lad","b_compartment","obs_weight")])
  )
  df <- df[!is.na(df$b_compartment), ]
  df$segdup        <- as.integer(df$segdup)
  df$lad           <- as.integer(df$lad)
  df$b_compartment <- as.integer(df$b_compartment)
  m <- glm(is_case ~ segdup + lad + b_compartment, data = df,
           family = binomial, weights = df$obs_weight)
  ci  <- confint.default(m)["segdup", ]
  co  <- coef(m)["segdup"]
  pv  <- summary(m)$coefficients["segdup", "Pr(>|z|)"]
  data.frame(threshold = threshold,
             n_somatic  = sum(df$is_case == 1),
             OR  = exp(co), CI_lo = exp(ci[1]), CI_hi = exp(ci[2]),
             p_val = pv)
}

load_admr <- function(files) {
  dt <- lapply(files, function(f) {
    d <- fread(f)
    setnames(d, c("chr","start","end"), c("seqnames","start","end"), skip_absent = TRUE)
    if ("meanMethy1" %in% names(d)) setnames(d, c("meanMethy1","meanMethy2"),
                                               c("HP1.Methy","HP2.Methy"), skip_absent = TRUE)
    d
  }) |> rbindlist(fill = TRUE)
  gr <- makeGRangesFromDataFrame(dt, seqnames.field = "seqnames",
                                 keep.extra.columns = TRUE)
  gr[seqnames(gr) %in% MAIN_CHROMS]
}

# ── Per-patient data collection ───────────────────────────────────────────────
cat("[1] Loading per-patient tumor/normal aDMR...\n")
tumor_files  <- list.files(DMR_DIR, pattern = "\\.tumor_aDMR\\.sorted\\.txt$",  full.names = TRUE)
normal_files <- list.files(DMR_DIR, pattern = "\\.normal_aDMR\\.sorted\\.txt$", full.names = TRUE)
names(tumor_files)  <- str_remove(basename(tumor_files),  "\\.tumor_aDMR\\.sorted\\.txt$")
names(normal_files) <- str_remove(basename(normal_files), "\\.normal_aDMR\\.sorted\\.txt$")

stats_list <- list()
# For each threshold: collect somatic aDMR GRanges
# We store in a named list: key = threshold (as char), value = combined GRanges
somatic_gr_by_thresh <- setNames(vector("list", length(THRESHOLDS)), as.character(THRESHOLDS))
for (thr in THRESHOLDS) somatic_gr_by_thresh[[as.character(thr)]] <- list()

for (pt in patients) {
  sample_ids <- pmap$Samples_ID[pmap$patient_code == pt]
  tf <- tumor_files [names(tumor_files)  %in% sample_ids]
  nf <- normal_files[names(normal_files) %in% sample_ids]
  if (!length(tf) || !length(nf)) {
    cat(sprintf("  %s: SKIP (missing files)\n", pt)); next
  }

  tumor_gr  <- load_admr(tf)
  normal_gr <- load_admr(nf)
  cat(sprintf("  %s: tumor=%d, normal=%d\n", pt, length(tumor_gr), length(normal_gr)))

  # Compute overlap fraction for each tumor aDMR vs normal aDMR
  hits <- findOverlaps(tumor_gr, normal_gr, ignore.strand = TRUE)
  overlap_frac <- numeric(length(tumor_gr))
  if (length(hits) > 0) {
    ol_len  <- width(pintersect(tumor_gr[queryHits(hits)], normal_gr[subjectHits(hits)],
                                 ignore.strand = TRUE))
    hit_frac <- ol_len / width(tumor_gr[queryHits(hits)])
    # For each tumor aDMR, take max overlap fraction across all normal hits
    max_frac <- tapply(hit_frac, queryHits(hits), max)
    overlap_frac[as.integer(names(max_frac))] <- max_frac
  }
  mcols(tumor_gr)$overlap_frac_norm <- overlap_frac

  # Count constitutional vs somatic at each threshold
  pt_stats <- data.frame(
    patient_code   = pt,
    n_tumor_admr   = length(tumor_gr),
    n_normal_admr  = length(normal_gr)
  )
  for (thr in THRESHOLDS) {
    const_mask  <- overlap_frac >= thr
    somatic_mask <- !const_mask
    pt_stats[[paste0("n_const_thr",  sprintf("%02d", thr*100))]] <- sum(const_mask)
    pt_stats[[paste0("n_somatic_thr", sprintf("%02d", thr*100))]] <- sum(somatic_mask)
    somatic_gr_by_thresh[[as.character(thr)]][[pt]] <- tumor_gr[somatic_mask]
  }
  stats_list[[pt]] <- pt_stats
  rm(tumor_gr, normal_gr); gc(verbose = FALSE)
}

stats_dt <- rbindlist(stats_list, fill = TRUE)
fwrite(stats_dt, file.path(OUT_ROOT, "normal_admr_stats.csv"))
cat(sprintf("\n[Priority 2] Normal aDMR stats saved: %d patients\n", nrow(stats_dt)))
print(stats_dt[, .(patient_code, n_tumor_admr, n_normal_admr,
                   n_const_thr10, n_somatic_thr10)])

# ── GLM per threshold ──────────────────────────────────────────────────────────
cat("\n[2] Running GLM for each threshold...\n")
or_results <- list()

for (thr in THRESHOLDS) {
  cat(sprintf("  Threshold %.0f%%...", thr * 100))
  combined <- do.call(c, unname(Filter(Negate(is.null),
                                        somatic_gr_by_thresh[[as.character(thr)]])))
  if (length(combined) < 100) {
    cat(sprintf(" SKIP (n=%d too small)\n", length(combined))); next
  }
  cat(sprintf(" n=%d\n", length(combined)))
  res <- tryCatch(run_glm_gr(combined, thr),
                  error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
  if (!is.null(res)) or_results[[as.character(thr)]] <- res
}

or_dt <- rbindlist(or_results)
or_dt[, sig := fcase(p_val < 0.001, "***", p_val < 0.01, "**",
                     p_val < 0.05, "*", default = "ns")]
or_dt[, threshold_label := paste0(threshold * 100, "% overlap")]
fwrite(or_dt, file.path(OUT_ROOT, "const_amr_threshold_or.csv"))
cat("\n[Priority 1] Threshold sensitivity results:\n")
print(or_dt[, .(threshold_label, n_somatic, OR, CI_lo, CI_hi, p_val, sig)])

log_decision(sprintf(
  "const_amr_threshold_sensitivity: C14 OR at 5%%=%.3f, 10%%=%.3f, 25%%=%.3f, 50%%=%.3f; all %s",
  or_dt[threshold == 0.05, OR], or_dt[threshold == 0.10, OR],
  or_dt[threshold == 0.25, OR], or_dt[threshold == 0.50, OR],
  if (all(or_dt$p_val > 0.05)) "ns" else "some significant"
))
cat("Done.\n")
