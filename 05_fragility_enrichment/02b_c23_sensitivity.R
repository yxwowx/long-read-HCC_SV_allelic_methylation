#!/usr/bin/env Rscript
# 02b_c23_sensitivity.R
# Sensitivity analysis: C23 re-run with 1kb, 10kb (already done), 50kb, 100kb cutoffs
# Uses already-annotated somatic_admr_annotated.csv.gz (has segdup/lad/b_compartment)
# Output: SV_aDMR/c23_sensitivity.csv

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

UTILS_DIR <- "/home/kachungk/script/SV-DMR/shared_file/pipeline"
source(file.path(UTILS_DIR, "shared_utils.R"))

DATA_ROOT    <- "/node200data/kachungk/hcc_data"
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")

CTRL_MULT  <- 10L
SEEDS      <- c(42L, 7L, 13L, 100L, 999L)
MAIN_CHROMS <- paste0("chr", c(1:22, "X"))

DIST_CUTOFFS_KB <- c(1L, 10L, 50L, 100L)  # distances to test

REF_ROOT   <- "/node200data/kachungk/reference/GRCh38"
SEGDUP_BED <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD_BED    <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")

cat("[0] Loading reference features...\n")
segdup_gr <- import.bed(SEGDUP_BED) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
lad_gr    <- import.bed(LAD_BED)    |> keepStandardChromosomes(pruning.mode = "coarse")
lad_gr    <- lad_gr[seqnames(lad_gr) %in% MAIN_CHROMS]

# ── Load annotated somatic aDMR (already has segdup/lad/b_compartment) ────────
cat("[1] Loading annotated somatic aDMR...\n")
admr_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
admr_dt <- admr_dt[seqnames %in% MAIN_CHROMS]
admr_dt[, log10_nCG := log10(pmax(nCG, 1))]
cat(sprintf("  %d somatic aDMR\n", nrow(admr_dt)))

admr_gr_all <- makeGRangesFromDataFrame(admr_dt, keep.extra.columns = TRUE)

cat("[2] Loading SV breakpoints...\n")
sv_dt <- fread(file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz"))
sv_dt <- sv_dt[seqnames %in% MAIN_CHROMS]
sv_gr_all <- makeGRangesFromDataFrame(
  sv_dt[, .(seqnames, start, end = start + 1L)],
  keep.extra.columns = FALSE
)
cat(sprintf("  %d SV breakpoints\n", length(sv_gr_all)))

# ── Pre-compute aDMR subsets for each distance cutoff ────────────────────────
cat("[3] Computing aDMR subsets per distance cutoff...\n")
admr_far <- lapply(DIST_CUTOFFS_KB, function(d_kb) {
  sv_flanked <- suppressWarnings(
    trim(resize(sv_gr_all, width = 2L * d_kb * 1000L + 1L, fix = "center"))
  )
  mask <- !IRanges::overlapsAny(admr_gr_all, sv_flanked)
  n    <- sum(mask)
  cat(sprintf("  >%dkb from SV: %d / %d (%.1f%%)\n",
              d_kb, n, length(admr_gr_all), n/length(admr_gr_all)*100))
  admr_gr_all[mask]
})
names(admr_far) <- paste0(DIST_CUTOFFS_KB, "kb")

# ── Width+chr-matched control generator ───────────────────────────────────────
make_controls <- function(case_gr, n_mult, seed) {
  set.seed(seed)
  ctrl_list <- lapply(unique(as.character(seqnames(case_gr))), function(chr) {
    gc_chr <- case_gr[seqnames(case_gr) == chr]
    if (length(gc_chr) == 0) return(NULL)
    chr_len <- CHROM_LENS[[chr]]
    if (is.null(chr_len) || is.na(chr_len)) return(NULL)
    widths  <- width(gc_chr)
    starts  <- unlist(lapply(widths, function(w) {
      max_s <- chr_len - w
      if (max_s < 1L) return(NULL)
      sample.int(max_s, size = n_mult, replace = TRUE)
    }))
    if (length(starts) == 0) return(NULL)
    rep_w <- rep(widths, each = n_mult)
    GRanges(seqnames = chr,
            ranges   = IRanges(starts, starts + rep_w - 1L),
            log10_nCG = log10(pmax(rep_w * 0.00966, 1)),
            is_case  = 0L, obs_weight = 1 / n_mult)
  })
  ctrl_gr <- do.call(c, Filter(Negate(is.null), ctrl_list))
  # Annotate segdup/lad overlap (fast %over%, no binnedAverage needed)
  mcols(ctrl_gr)$segdup <- ctrl_gr %over% segdup_gr
  mcols(ctrl_gr)$lad    <- ctrl_gr %over% lad_gr
  ctrl_gr
}

# ── GLM helper (uses pre-computed annotations on case_gr) ────────────────────
run_c23_glm <- function(case_gr, ctrl_gr, seed) {
  df_case <- as.data.frame(mcols(case_gr))
  df_case$is_case    <- 1L
  df_case$obs_weight <- 1

  df_ctrl <- as.data.frame(mcols(ctrl_gr))
  df_ctrl$is_case <- 0L

  needed <- c("is_case", "segdup", "lad", "log10_nCG", "obs_weight")
  df <- rbind(df_case[, needed], df_ctrl[, needed])
  df$segdup <- as.integer(df$segdup)
  df$lad    <- as.integer(df$lad)

  m <- tryCatch(
    glm(is_case ~ segdup + lad, data = df,
        family = binomial, weights = df$obs_weight),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  coefs <- coef(m)
  ci    <- tryCatch(confint.default(m), error = function(e) matrix(NA, nrow=length(coefs), ncol=2))
  pvals <- summary(m)$coefficients[, "Pr(>|z|)"]

  data.frame(
    seed    = seed,
    term    = names(coefs),
    OR      = exp(coefs),
    CI_lo   = exp(ci[, 1]),
    CI_hi   = exp(ci[, 2]),
    p_val   = pvals,
    n_case  = sum(df$is_case == 1L),
    n_ctrl  = sum(df$is_case == 0L)
  )
}

# ── Run across cutoffs and seeds ──────────────────────────────────────────────
cat("[4] Running GLMs...\n")
results <- list()

for (d_label in names(admr_far)) {
  ag <- admr_far[[d_label]]
  cat(sprintf("\n  Cutoff: >%s (n=%d)\n", d_label, length(ag)))
  for (seed in SEEDS) {
    ctrl <- make_controls(ag, CTRL_MULT, seed)
    res  <- run_c23_glm(ag, ctrl, seed)
    if (!is.null(res)) {
      res$cutoff_kb <- as.integer(gsub("kb", "", d_label))
      res$claim     <- paste0("C23_admr_sv_far_", d_label)
      results[[length(results)+1]] <- res
    }
  }
  cat(sprintf("  Done: %d\n", length(ag)))
}

# ── Also pull existing 10kb results from segdup_or_table.csv for comparison ──
existing <- fread(file.path(OUT_ROOT, "segdup_or_table.csv"))
c23_existing <- existing[claim == "C23_admr_sv_far_segdup" & term == "segdup"]

# ── Summarize ─────────────────────────────────────────────────────────────────
results_dt <- rbindlist(results, fill = TRUE)
segdup_dt  <- results_dt[term == "segdup"]

cat("\n=== C23 Sensitivity: SegDup OR by distance cutoff (seed=42) ===\n")
print(segdup_dt[seed == 42L, .(cutoff_kb, n_case, n_ctrl, OR, CI_lo, CI_hi, p_val)])

cat("\n=== Seed stability by cutoff ===\n")
stab <- segdup_dt[, .(
  OR_mean = mean(OR), OR_sd = sd(OR), pct_chg = (max(OR)-min(OR))/mean(OR)*100
), by = .(cutoff_kb, claim)]
print(stab[order(cutoff_kb)])

# Append 10kb from existing for side-by-side
c23_ref <- data.table(
  cutoff_kb = 10L, claim = "C23_admr_sv_far_10kb (02_glm)",
  OR_mean = mean(c23_existing$OR), OR_sd = sd(c23_existing$OR),
  pct_chg = (max(c23_existing$OR)-min(c23_existing$OR))/mean(c23_existing$OR)*100
)

cat("\n=== Reference C23 from 02_fragility_glm (10kb, from segdup_or_table.csv) ===\n")
print(c23_existing[, .(seed, OR, CI_lo, CI_hi, p_val, n_case)])

fwrite(results_dt, file.path(OUT_ROOT, "c23_sensitivity.csv"))
cat(sprintf("\nSaved: c23_sensitivity.csv (%d rows)\n", nrow(results_dt)))

log_decision("02b_c23_sensitivity: C23 distance sensitivity 1/10/50/100kb; all ORs reported")
cat("=== Done: 02b_c23_sensitivity.R ===\n")
