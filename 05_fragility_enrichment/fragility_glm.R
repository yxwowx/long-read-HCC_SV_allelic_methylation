#!/usr/bin/env Rscript
# 02_fragility_glm.R
# Weighted GLM for SegDup enrichment in somatic aDMR and SV breakpoints
# Claims: C13 (SV → SegDup), C14 (somatic aDMR → SegDup), C22 (+log10_nCG), C23 (SV-far)
# Pseudoreplication fix (§4.1): weight controls by 1/N_CTRL_MULT
# Output: SV_aDMR/segdup_or_table.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

DATA_ROOT    <- Sys.getenv("HCC_DATA_DIR")
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")

REF_ROOT  <- Sys.getenv("REFERENCE_DIR")
SEGDUP_BED <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD_BED    <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW     <- file.path(REF_ROOT, "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")

CTRL_MULT_ADMR <- 10L   # controls per aDMR locus
CTRL_MULT_SV   <- 5L    # controls per SV breakpoint
SEEDS          <- c(42L, 7L, 13L, 100L, 999L)   # 5-seed sensitivity
DIST_FAR_KB    <- 10000L  # bp; aDMR >10kb from any SV for C23

MAIN_CHROMS <- paste0("chr", c(1:22, "X"))

# Load reference features ======================================================
cat("[0] Loading reference features...\n")
segdup_gr <- import.bed(SEGDUP_BED) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
lad_gr    <- import.bed(LAD_BED)    |> keepStandardChromosomes(pruning.mode = "coarse")
lad_gr    <- lad_gr[seqnames(lad_gr) %in% MAIN_CHROMS]

has_pc1 <- file.exists(PC1_BW)
if (has_pc1) pc1_rle <- import(PC1_BW, as = "RleList")

annotate_gr <- function(gr) {
  mcols(gr)$segdup <- gr %over% segdup_gr
  mcols(gr)$lad    <- gr %over% lad_gr
  mcols(gr)$b_compartment <- NA
  if (has_pc1) {
    shared_chroms <- intersect(seqlevels(gr), names(pc1_rle))
    if (length(shared_chroms) > 0) {
      gr_sub  <- keepSeqlevels(gr, shared_chroms, pruning.mode = "coarse")
      # Use seqlevels(gr_sub) order to ensure names match for binnedAverage
      pc1_sub <- pc1_rle[seqlevels(gr_sub)]
      scores  <- binnedAverage(gr_sub, pc1_sub, varname = "PC1")
      idx_in  <- which(as.character(seqnames(gr)) %in% shared_chroms)
      mcols(gr)$b_compartment[idx_in] <- mcols(scores)$PC1 < 0
    }
  }
  gr
}

# Generate width+chr-matched controls ==========================================
make_controls <- function(case_gr, n_mult, seed) {
  set.seed(seed)
  ctrl_list <- lapply(seqlevels(case_gr)[seqlevels(case_gr) %in% MAIN_CHROMS], function(chr) {
    cases_chr <- case_gr[seqnames(case_gr) == chr]
    if (length(cases_chr) == 0) return(NULL)
    chr_len <- CHROM_LENS[[chr]]
    if (is.null(chr_len) || is.na(chr_len)) return(NULL)
    widths <- width(cases_chr)
    starts <- unlist(lapply(widths, function(w) {
      max_start <- chr_len - w
      if (max_start < 1L) return(NULL)
      sample.int(max_start, size = n_mult, replace = TRUE)
    }))
    if (length(starts) == 0) return(NULL)
    rep_widths <- rep(widths, each = n_mult)
    GRanges(seqnames = chr,
            ranges   = IRanges(starts, starts + rep_widths - 1L),
            is_case  = 0L,
            obs_weight = 1 / n_mult)
  })
  do.call(c, Filter(Negate(is.null), ctrl_list))
}

# GLM helper ===================================================================
run_glm <- function(case_gr, ctrl_gr, formula_str, claim, seed, granularity = "window") {
  case_gr <- annotate_gr(case_gr)
  ctrl_gr <- annotate_gr(ctrl_gr)

  mcols(case_gr)$is_case     <- 1L
  mcols(case_gr)$obs_weight  <- 1
  mcols(case_gr)$log10_nCG   <- if (!is.null(mcols(case_gr)$log10_nCG)) mcols(case_gr)$log10_nCG else NA_real_
  # Width-expected genome background CpG density (hg38: ~28M CpGs / 2.9 Gb = 0.00966 CpG/bp)
  mcols(ctrl_gr)$log10_nCG   <- log10(pmax(width(ctrl_gr) * 0.00966, 1))

  df <- rbind(
    as.data.frame(mcols(case_gr))[, c("is_case", "segdup", "lad", "b_compartment",
                                       "log10_nCG", "obs_weight")],
    as.data.frame(mcols(ctrl_gr))[, c("is_case", "segdup", "lad", "b_compartment",
                                       "log10_nCG", "obs_weight")]
  )
  df$segdup        <- as.integer(df$segdup)
  df$lad           <- as.integer(df$lad)
  df$b_compartment <- as.integer(df$b_compartment)
  df <- df[!is.na(df$b_compartment), ]  # drop NA PC1 rows

  fml <- as.formula(formula_str)
  m   <- glm(fml, data = df, family = binomial, weights = df$obs_weight)

  coefs  <- coef(m)
  ci     <- tryCatch(confint.default(m), error = function(e) matrix(NA, nrow = length(coefs), ncol = 2))
  pvals  <- summary(m)$coefficients[, "Pr(>|z|)"]

  data.frame(
    claim       = claim,
    seed        = seed,
    granularity = granularity,
    term        = names(coefs),
    beta        = as.numeric(coefs),
    se          = as.numeric(sqrt(diag(vcov(m)))),
    OR          = exp(coefs),
    CI_lo       = exp(ci[, 1]),
    CI_hi       = exp(ci[, 2]),
    p_val       = pvals,
    n_case      = sum(df$is_case == 1L),
    n_ctrl      = sum(df$is_case == 0L)
  )
}

# Load data ====================================================================
cat("[1] Loading somatic aDMR...\n")
admr_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
admr_dt <- admr_dt[seqnames %in% MAIN_CHROMS]
admr_gr_all <- makeGRangesFromDataFrame(admr_dt, keep.extra.columns = TRUE)
cat(sprintf("  %d somatic aDMR\n", length(admr_gr_all)))

cat("[2] Loading SV breakpoints...\n")
sv_dt <- fread(file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz"))
sv_dt <- sv_dt[seqnames %in% MAIN_CHROMS]
# Each row is a breakpoint (use start position as bp)
sv_gr_all <- makeGRangesFromDataFrame(
  sv_dt[, .(seqnames, start, end = start + 1L, sample, geom_type,
             cnv_class, stratification)],
  keep.extra.columns = TRUE
)
cat(sprintf("  %d SV breakpoints\n", length(sv_gr_all)))

# C23: aDMR far from SV (>10kb) ================================================
sv_flanked_gr <- suppressWarnings(trim(resize(sv_gr_all, width = 2 * DIST_FAR_KB + 1L, fix = "center")))
admr_far_mask <- !IRanges::overlapsAny(admr_gr_all, sv_flanked_gr)
cat(sprintf("  C23: aDMR >%dkb from SV: %d / %d (%.1f%%)\n",
            DIST_FAR_KB %/% 1000L, sum(admr_far_mask),
            length(admr_gr_all), mean(admr_far_mask)*100))
admr_gr_far <- admr_gr_all[admr_far_mask]

# Run GLMs across seeds ========================================================
results <- list()

for (seed in SEEDS) {
  cat(sprintf("\n--- Seed %d ---\n", seed))

  # C13: SV → SegDup (5× controls, point breakpoints) — always "region" granularity
  ctrl_sv <- make_controls(sv_gr_all, CTRL_MULT_SV, seed)
  res <- try(run_glm(sv_gr_all, ctrl_sv,
                     "is_case ~ segdup + lad + b_compartment",
                     "C13_sv_segdup", seed, granularity = "region"), silent = FALSE)
  if (inherits(res, "try-error")) cat("  C13 FAILED\n") else {
    results[[length(results)+1]] <- res
    cat(sprintf("  C13 OK (n_case=%d, n_ctrl=%d)\n", length(sv_gr_all), length(ctrl_sv)))
  }

  # C13-marginal: SV vs random, one feature at a time (seed 42 only) — feeds the
  # 3-population x fragility-feature comparison (figS9); avoids the sign-flip
  # that mutual adjustment introduces for LAD (adjusted OR=0.79 vs marginal ~1.36)
  if (seed == 42L) {
    for (feat in c("segdup", "lad", "b_compartment")) {
      res <- try(run_glm(sv_gr_all, ctrl_sv,
                         sprintf("is_case ~ %s", feat),
                         sprintf("C13_sv_%s_marginal", feat), seed, granularity = "region"),
                 silent = FALSE)
      if (inherits(res, "try-error")) cat(sprintf("  C13-marginal (%s) FAILED\n", feat)) else {
        results[[length(results)+1]] <- res
        cat(sprintf("  C13-marginal (%s) OK\n", feat))
      }
    }
  }

  # C14: somatic aDMR → SegDup (10× controls)
  ctrl_admr <- make_controls(admr_gr_all, CTRL_MULT_ADMR, seed)
  mcols(admr_gr_all)$log10_nCG <- log10(pmax(mcols(admr_gr_all)$nCG, 1))
  cat(sprintf("  C14 running (n_case=%d, n_ctrl=%d)...\n",
              length(admr_gr_all), length(ctrl_admr)))
  res <- try(run_glm(admr_gr_all, ctrl_admr,
                     "is_case ~ segdup + lad + b_compartment",
                     "C14_admr_segdup", seed), silent = FALSE)
  if (inherits(res, "try-error")) cat("  C14 FAILED\n") else {
    results[[length(results)+1]] <- res
    cat("  C14 OK\n")
  }

  # C22: somatic aDMR → SegDup + log10_nCG (CpG-density adjusted)
  res <- try(run_glm(admr_gr_all, ctrl_admr,
                     "is_case ~ segdup + lad + b_compartment + log10_nCG",
                     "C22_admr_segdup_cpg", seed), silent = FALSE)
  if (inherits(res, "try-error")) cat("  C22 FAILED\n") else {
    results[[length(results)+1]] <- res
    cat("  C22 OK\n")
  }

  # C23: SV-far aDMR → SegDup
  if (length(admr_gr_far) > 100L) {
    ctrl_far <- make_controls(admr_gr_far, CTRL_MULT_ADMR, seed)
    res <- try(run_glm(admr_gr_far, ctrl_far,
                       "is_case ~ segdup + lad + b_compartment",
                       "C23_admr_sv_far_segdup", seed), silent = FALSE)
    if (inherits(res, "try-error")) cat("  C23 FAILED\n") else {
      results[[length(results)+1]] <- res
      cat("  C23 OK\n")
    }
  }
}

# Merged-region analyses (cross-patient union) =================================
# Purpose: collapse 2.88M per-patient per-window loci into non-redundant genomic
#   regions for an apples-to-apples comparison with the per-SV-breakpoint C13 GLM.
#   Window-level OR (~1.13) is a pseudoreplication artifact of counting CpG windows
#   per region; merged-region OR (~2–4) is biologically interpretable.
cat("\n[Merged-region] Building cross-patient union...\n")
admr_gr_merged_raw <- GenomicRanges::reduce(admr_gr_all)

# Assign max nCG from constituent windows per merged region
hits_m  <- findOverlaps(admr_gr_merged_raw, admr_gr_all)
ncg_dt  <- data.table(qhit = queryHits(hits_m),
                       ncg  = mcols(admr_gr_all)$nCG[subjectHits(hits_m)])
max_ncg_dt <- ncg_dt[, .(max_ncg = max(ncg, na.rm = TRUE)), by = qhit][order(qhit)]
ncg_full    <- rep(1L, length(admr_gr_merged_raw))
ncg_full[max_ncg_dt$qhit] <- max_ncg_dt$max_ncg
mcols(admr_gr_merged_raw)$nCG       <- ncg_full
mcols(admr_gr_merged_raw)$log10_nCG <- log10(pmax(ncg_full, 1L))
admr_gr_merged <- admr_gr_merged_raw
cat(sprintf("  %d merged regions from %d windows (%.1f%% reduction)\n",
            length(admr_gr_merged),
            length(admr_gr_all),
            (1 - length(admr_gr_merged)/length(admr_gr_all)) * 100))

# C23 mask for merged regions (>10kb from any SV)
admr_merged_far_mask <- !IRanges::overlapsAny(admr_gr_merged, sv_flanked_gr)
admr_gr_merged_far   <- admr_gr_merged[admr_merged_far_mask]
cat(sprintf("  C23-region: %d / %d merged regions far from SV (%.1f%%)\n",
            sum(admr_merged_far_mask), length(admr_gr_merged),
            mean(admr_merged_far_mask) * 100))

results_merged <- list()
for (seed in SEEDS) {
  cat(sprintf("\n--- Merged-region seed %d ---\n", seed))

  ctrl_m <- make_controls(admr_gr_merged, CTRL_MULT_ADMR, seed)

  # C14-region: merged aDMR → SegDup
  res <- try(run_glm(admr_gr_merged, ctrl_m,
                     "is_case ~ segdup + lad + b_compartment",
                     "C14_admr_segdup", seed, granularity = "region"), silent = FALSE)
  if (inherits(res, "try-error")) cat("  C14-r FAILED\n") else {
    results_merged[[length(results_merged)+1]] <- res; cat("  C14-r OK\n")
  }

  # C22-region: merged aDMR + CpG-density adj.
  res <- try(run_glm(admr_gr_merged, ctrl_m,
                     "is_case ~ segdup + lad + b_compartment + log10_nCG",
                     "C22_admr_segdup_cpg", seed, granularity = "region"), silent = FALSE)
  if (inherits(res, "try-error")) cat("  C22-r FAILED\n") else {
    results_merged[[length(results_merged)+1]] <- res; cat("  C22-r OK\n")
  }

  # C14-region-marginal: merged aDMR vs random, one feature at a time (seed 42 only)
  if (seed == 42L) {
    for (feat in c("segdup", "lad", "b_compartment")) {
      res <- try(run_glm(admr_gr_merged, ctrl_m,
                         sprintf("is_case ~ %s", feat),
                         sprintf("C14_admr_%s_marginal", feat), seed, granularity = "region"),
                 silent = FALSE)
      if (inherits(res, "try-error")) cat(sprintf("  C14-r-marginal (%s) FAILED\n", feat)) else {
        results_merged[[length(results_merged)+1]] <- res
        cat(sprintf("  C14-r-marginal (%s) OK\n", feat))
      }
    }
  }

  # C23-region: SV-far merged aDMR → SegDup
  if (length(admr_gr_merged_far) > 100L) {
    ctrl_mf <- make_controls(admr_gr_merged_far, CTRL_MULT_ADMR, seed)
    res <- try(run_glm(admr_gr_merged_far, ctrl_mf,
                       "is_case ~ segdup + lad + b_compartment",
                       "C23_admr_sv_far_segdup", seed, granularity = "region"), silent = FALSE)
    if (inherits(res, "try-error")) cat("  C23-r FAILED\n") else {
      results_merged[[length(results_merged)+1]] <- res; cat("  C23-r OK\n")
    }
  }
}

# Summarize ====================================================================
if (length(results) == 0) stop("All window-granularity run_glm calls returned try-error")
results_all <- c(results, results_merged)
results_dt  <- rbindlist(results_all)

# Primary results: region-granularity, seed=42, segdup term
primary <- results_dt[granularity == "region" & seed == 42L & term == "segdup"]
cat("\n=== Primary SegDup OR — region granularity (seed=42) ===\n")
print(primary[, .(claim, granularity, OR, CI_lo, CI_hi, p_val, n_case, n_ctrl)])

# Comparison: window vs region for C14 (seed=42, segdup term)
cat("\n=== Window vs Region OR comparison (C14, seed=42) ===\n")
comp <- results_dt[claim == "C14_admr_segdup" & seed == 42L & term == "segdup",
                    .(granularity, OR, CI_lo, CI_hi, p_val, n_case)]
print(comp)

# Seed-to-seed OR stability — region granularity
stability <- results_dt[granularity == "region" & term == "segdup", .(
  OR_mean  = mean(OR),
  OR_sd    = sd(OR),
  OR_cv    = sd(OR)/mean(OR),
  pct_chg  = (max(OR) - min(OR)) / mean(OR) * 100
), by = .(claim, granularity)]
cat("\n=== Seed stability — region (segdup term) ===\n")
print(stability)

fwrite(results_dt, file.path(OUT_ROOT, "segdup_or_table.csv"))
cat(sprintf("\nSaved: segdup_or_table.csv (%d rows, both window + region granularity)\n",
            nrow(results_dt)))

cat("=== Done: 02_fragility_glm.R ===\n")

