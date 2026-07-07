#!/usr/bin/env Rscript
# 05_somatic_vs_constitutional.R
# Three-way SegDup OR comparison:
#   (A) Somatic aDMR: tumor aDMR minus normal aDMR overlap ≥10%
#   (B) Recurrently remodeled constitutional aDMR: all confident_dmr-derived tumor aDMR
#       (bulk tumor-normal DMR ∩ tumor aDMR, admr_* coords)
#   (C) Recurrently remodeled normal aDMR: subset of (B) with ov_norm_admr=TRUE
# Method: per-set weighted GLM + stacked interaction (A vs B) + bootstrap ΔOR
# Output: SV_aDMR/somatic_vs_constitutional_segdup.csv

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
CONF_DMR_PATH <- file.path(DMR_SVS_ROOT, "01.DMR_recurrence/confident_dmr_per_patient.csv.gz")

REF_ROOT   <- Sys.getenv("REFERENCE_DIR")
SEGDUP_BED <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD_BED    <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW     <- file.path(REF_ROOT, "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")

CTRL_MULT  <- 10L
N_BOOT     <- 1000L
MAIN_CHROMS <- paste0("chr", c(1:22, "X"))
set.seed(42L)

# Load reference features ======================================================
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
  mcols(gr)$b_compartment <- NA
  if (has_pc1) {
    shared_chroms <- intersect(seqlevels(gr), names(pc1_rle))
    if (length(shared_chroms) > 0) {
      gr_sub  <- keepSeqlevels(gr, shared_chroms, pruning.mode = "coarse")
      pc1_sub <- pc1_rle[seqlevels(gr_sub)]
      sc      <- binnedAverage(gr_sub, pc1_sub, varname = "PC1")
      idx_in  <- which(as.character(seqnames(gr)) %in% shared_chroms)
      mcols(gr)$b_compartment[idx_in] <- mcols(sc)$PC1 < 0
    }
  }
  gr
}

make_controls <- function(case_gr, n_mult, seed) {
  set.seed(seed)
  ctrl_list <- lapply(unique(as.character(seqnames(case_gr))), function(chr) {
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
    GRanges(seqnames = chr,
            ranges   = IRanges(starts, starts + rep(widths, each = n_mult) - 1L),
            is_case  = 0L)
  })
  do.call(c, Filter(Negate(is.null), ctrl_list))
}

# Load somatic aDMR (set A) ====================================================
cat("[1] Loading somatic aDMR...\n")
somatic_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
somatic_dt <- somatic_dt[seqnames %in% MAIN_CHROMS]
somatic_gr <- makeGRangesFromDataFrame(somatic_dt, seqnames.field = "seqnames")
cat(sprintf("  Somatic aDMR: %d\n", length(somatic_gr)))

# Load recurrently remodeled constitutional aDMR (set B) =======================
# Source: confident_dmr_per_patient.csv.gz (bulk tumor-normal DMR ∩ tumor aDMR)
# Coordinates: admr_chr/admr_start/admr_end = tumor aDMR locus
# Set B (all): ALL confident_dmr-derived tumor aDMR = recurrently remodeled constitutional aDMR
# Set C (subset): ov_norm_admr=TRUE = recurrently remodeled normal aDMR
cat("[2] Loading recurrently remodeled constitutional aDMR from confident_dmr...\n")
conf_dt <- fread(CONF_DMR_PATH)
conf_dt <- conf_dt[admr_chr %in% MAIN_CHROMS]
# Rename bulk DMR coord cols — GRanges forbids metadata named seqnames/start/end/width/strand
setnames(conf_dt,
  old = c("seqnames","start","end","width","strand"),
  new = c("bulk_chr","bulk_start","bulk_end","bulk_width","bulk_strand"))
# Use admr_* coords (tumor aDMR locus) as GRanges coordinates
conf_gr_all <- filter_dmr_region(
  makeGRangesFromDataFrame(conf_dt,
    seqnames.field     = "admr_chr",
    start.field        = "admr_start",
    end.field          = "admr_end",
    keep.extra.columns = TRUE),
  meth_diff_cutoff = 0.2
)
conf_gr_all <- conf_gr_all[seqnames(conf_gr_all) %in% MAIN_CHROMS]

constitutional_gr <- conf_gr_all
normal_admr_gr    <- conf_gr_all[mcols(conf_gr_all)$ov_norm_admr == TRUE]

cat(sprintf("  Recurrently remodeled constitutional aDMR (all): %d\n",
            length(constitutional_gr)))
cat(sprintf("  Recurrently remodeled normal aDMR (ov_norm_admr=TRUE): %d (%.1f%%)\n",
            length(normal_admr_gr),
            length(normal_admr_gr) / length(constitutional_gr) * 100))

# Annotate and build analysis data frames ======================================
cat("[3] Annotating all three sets...\n")
somatic_gr        <- annotate_gr(somatic_gr)
constitutional_gr <- annotate_gr(constitutional_gr)
normal_admr_gr    <- annotate_gr(normal_admr_gr)

ctrl_som    <- annotate_gr(make_controls(somatic_gr,        CTRL_MULT, 42L))
ctrl_const  <- annotate_gr(make_controls(constitutional_gr, CTRL_MULT, 42L))
ctrl_normal <- annotate_gr(make_controls(normal_admr_gr,    CTRL_MULT, 42L))

build_df <- function(case_gr, ctrl_gr, set_label) {
  mc <- function(gr) {
    df <- as.data.frame(mcols(gr))
    df$is_case <- if ("is_case" %in% names(df)) df$is_case else 1L
    df
  }
  df_case <- mc(case_gr); df_case$is_case <- 1L; df_case$obs_weight <- 1.0
  df_ctrl <- mc(ctrl_gr); df_ctrl$is_case <- 0L; df_ctrl$obs_weight <- 1 / CTRL_MULT
  df <- rbind(df_case[, c("is_case","segdup","lad","b_compartment","obs_weight")],
              df_ctrl[, c("is_case","segdup","lad","b_compartment","obs_weight")])
  df$set           <- set_label
  df$segdup        <- as.integer(df$segdup)
  df$lad           <- as.integer(df$lad)
  df$b_compartment <- as.integer(df$b_compartment)
  df[!is.na(df$b_compartment), ]
}

df_som    <- build_df(somatic_gr,        ctrl_som,    "somatic")
df_const  <- build_df(constitutional_gr, ctrl_const,  "constitutional")
df_normal <- build_df(normal_admr_gr,    ctrl_normal, "normal_admr")

# Per-set GLM ==================================================================
cat("[4] Fitting per-set GLMs...\n")
fit_glm <- function(df) {
  glm(is_case ~ segdup + lad + b_compartment, data = df, family = binomial,
      weights = df$obs_weight)
}

m_som    <- fit_glm(df_som)
m_const  <- fit_glm(df_const)
m_normal <- fit_glm(df_normal)

get_or <- function(m, term = "segdup") {
  coefs <- coef(m)
  ci    <- tryCatch(confint.default(m), error = function(e) matrix(NA, length(coefs), 2))
  pv    <- summary(m)$coefficients[, "Pr(>|z|)"]
  data.frame(term = term, OR = exp(coefs[[term]]),
             CI_lo = exp(ci[term, 1]), CI_hi = exp(ci[term, 2]),
             p_val = pv[[term]])
}

or_somatic <- get_or(m_som)
or_const   <- get_or(m_const)
or_normal  <- get_or(m_normal)

cat(sprintf("  Somatic aDMR                       SegDup OR=%.3f [%.3f–%.3f] p=%.2e\n",
            or_somatic$OR, or_somatic$CI_lo, or_somatic$CI_hi, or_somatic$p_val))
cat(sprintf("  Recurrently remodeled constitutional SegDup OR=%.3f [%.3f–%.3f] p=%.2e\n",
            or_const$OR, or_const$CI_lo, or_const$CI_hi, or_const$p_val))
cat(sprintf("  Recurrently remodeled normal aDMR   SegDup OR=%.3f [%.3f–%.3f] p=%.2e\n",
            or_normal$OR, or_normal$CI_lo, or_normal$CI_hi, or_normal$p_val))

# Multi-feature OR table (SegDup/LAD/B-comp, adjusted + marginal) ==============
# Feeds the 3-population x fragility-feature comparison (viz/v1/figS9): the
# adjusted (mutually-conditioned) GLM above is reused for the "multivariate"
# framework; a separate one-feature-at-a-time GLM gives the "marginal" (vs
# matched-random) framework used as the primary panel encoding, since SegDup/
# LAD/B-comp collinearity can flip adjusted coefficient signs (documented for
# SV in 02_fragility_glm.R: LAD adjusted OR=0.79 vs marginal ~1.36).
cat("[4b] Multi-feature (SegDup/LAD/B-comp) OR table, adjusted + marginal...\n")
fit_glm_uni <- function(df, feature) {
  glm(as.formula(sprintf("is_case ~ %s", feature)), data = df, family = binomial,
      weights = df$obs_weight)
}

FEATURES <- c("segdup", "lad", "b_compartment")
SET_DFS  <- list(somatic = df_som, constitutional = df_const, normal_admr = df_normal)
SET_M    <- list(somatic = m_som,  constitutional = m_const,  normal_admr = m_normal)

multifeature_rows <- rbindlist(lapply(names(SET_DFS), function(set_label) {
  df <- SET_DFS[[set_label]]
  m_adj <- SET_M[[set_label]]
  rbindlist(lapply(FEATURES, function(feat) {
    adj <- get_or(m_adj, term = feat)
    m_uni <- fit_glm_uni(df, feat)
    uni   <- get_or(m_uni, term = feat)
    rbind(
      data.frame(population = set_label, feature = feat, framework = "multivariate",
                 OR = adj$OR, CI_lo = adj$CI_lo, CI_hi = adj$CI_hi, p_val = adj$p_val,
                 n_case = sum(df$is_case == 1L), n_ctrl = sum(df$is_case == 0L)),
      data.frame(population = set_label, feature = feat, framework = "marginal",
                 OR = uni$OR, CI_lo = uni$CI_lo, CI_hi = uni$CI_hi, p_val = uni$p_val,
                 n_case = sum(df$is_case == 1L), n_ctrl = sum(df$is_case == 0L))
    )
  }))
}))
multifeature_rows[, sig := fcase(
  p_val < 0.001, "***", p_val < 0.01, "**", p_val < 0.05, "*", default = "ns"
)]

cat("  Multi-feature OR table:\n")
print(multifeature_rows[, .(population, feature, framework, OR = round(OR,3),
                             CI_lo = round(CI_lo,3), CI_hi = round(CI_hi,3), sig)])

fwrite(multifeature_rows, file.path(OUT_ROOT, "fragility_or_by_population.csv"))
cat(sprintf("  Saved: %s\n", file.path(OUT_ROOT, "fragility_or_by_population.csv")))

# Stacked GLM: interaction test ================================================
cat("[5] Stacked GLM interaction test (segdup × set)...\n")
df_stacked <- rbind(
  cbind(df_som,   set_indicator = 0L),  # somatic = reference
  cbind(df_const, set_indicator = 1L)   # constitutional = 1
)
df_stacked$obs_weight_adj <- df_stacked$obs_weight * 0.5  # equal weight per set

m_interact <- glm(is_case ~ segdup * set_indicator + lad + b_compartment,
                  data    = df_stacked,
                  family  = binomial,
                  weights = df_stacked$obs_weight_adj)

interact_coef <- coef(m_interact)[["segdup:set_indicator"]]
interact_ci   <- tryCatch(confint.default(m_interact)["segdup:set_indicator", ],
                           error = function(e) c(NA, NA))
interact_p    <- summary(m_interact)$coefficients["segdup:set_indicator", "Pr(>|z|)"]

cat(sprintf("  Interaction β=%.3f, OR_ratio=%.2f [%.2f–%.2f], p=%.4f\n",
            interact_coef, exp(interact_coef),
            exp(interact_ci[1]), exp(interact_ci[2]), interact_p))

# Bootstrap CI on delta-OR =====================================================
cat(sprintf("[6] Bootstrapping ΔOR (%d resamples)...\n", N_BOOT))

boot_delta_or <- replicate(N_BOOT, {
  n_s <- nrow(df_som)
  n_c <- nrow(df_const)
  d_s <- df_som[sample.int(n_s, replace = TRUE), ]
  d_c <- df_const[sample.int(n_c, replace = TRUE), ]
  m1  <- tryCatch(fit_glm(d_s), error = function(e) NULL)
  m2  <- tryCatch(fit_glm(d_c), error = function(e) NULL)
  if (is.null(m1) || is.null(m2)) return(NA_real_)
  exp(coef(m1)[["segdup"]]) - exp(coef(m2)[["segdup"]])
})
boot_delta_or <- boot_delta_or[!is.na(boot_delta_or)]

delta_or_obs  <- or_somatic$OR - or_const$OR
boot_ci       <- quantile(boot_delta_or, c(0.025, 0.975))
boot_p        <- 2 * min(mean(boot_delta_or > 0), mean(boot_delta_or < 0))

cat(sprintf("  ΔOR (somatic−constitutional) = %.3f [%.3f–%.3f] boot_p=%.4f\n",
            delta_or_obs, boot_ci[1], boot_ci[2], boot_p))

# Save =========================================================================
result_dt <- data.table(
  set            = c("somatic", "constitutional", "normal_admr"),
  label          = c("Somatic aDMR",
                     "Recurrently remodeled constitutional aDMR",
                     "Recurrently remodeled normal aDMR"),
  OR_segdup      = c(or_somatic$OR,  or_const$OR,  or_normal$OR),
  CI_lo          = c(or_somatic$CI_lo, or_const$CI_lo, or_normal$CI_lo),
  CI_hi          = c(or_somatic$CI_hi, or_const$CI_hi, or_normal$CI_hi),
  p_val          = c(or_somatic$p_val, or_const$p_val, or_normal$p_val),
  n_case         = c(length(somatic_gr), length(constitutional_gr), length(normal_admr_gr)),
  n_ctrl         = c(length(ctrl_som),   length(ctrl_const),        length(ctrl_normal)),
  delta_or_vs_somatic   = c(NA_real_,
                             or_const$OR  - or_somatic$OR,
                             or_normal$OR - or_somatic$OR),
  delta_or_ci_lo        = c(NA_real_, boot_ci[1], NA_real_),
  delta_or_ci_hi        = c(NA_real_, boot_ci[2], NA_real_),
  delta_or_boot_p       = c(NA_real_, boot_p,     NA_real_),
  interaction_OR_ratio  = c(NA_real_, exp(interact_coef), NA_real_),
  interaction_p         = c(NA_real_, interact_p,         NA_real_)
)

fwrite(result_dt, file.path(OUT_ROOT, "somatic_vs_constitutional_segdup.csv"))
cat(sprintf("\nSaved: somatic_vs_constitutional_segdup.csv\n"))
print(result_dt[, .(set, OR_segdup, CI_lo, CI_hi, p_val, n_case)])

cat("=== Done: 05_somatic_vs_constitutional.R ===\n")
