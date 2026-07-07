suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  library(stringr)
  library(IRanges)
  library(scales)
  library(ggridges)
  library(ggplot2)
  library(S4Vectors)
  library(optparse)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character"),
  make_option(c("-o", "--outdir"),    type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Basic Objects & Functions ----------------------------------------------
## Common Functions ------------------------------------------------------

#' Filter direction-consistent hits
#' @description gr1의 diff.Methy 방향이 gr2의 diff.Methy 일치하는 hit만 반환
direction_consistent_hits <- function(hits, gr1, gr2) {
  if (length(hits) == 0) return(hits)
  dir1 <- sign(gr1$diff.Methy[queryHits(hits)])
  dir2 <- sign(gr2$diff.Methy[subjectHits(hits)])
  hits[dir1 == dir2]
}

#' Merged region 기준 recurrence 카운트 (환자별 binary matrix)
count_recurrence <- function(merged, per_patient_list, min_bp = 1L) {
  count_mat <- vapply(per_patient_list, function(gr) {
    as.integer(countOverlaps(merged, gr, minoverlap = min_bp) > 0)
  }, integer(length(merged)))
  rowSums(count_mat)
}

#' mingapwidth 결정을 위한 DMR 크기 및 간격 분포 분석 후 결과 report
diagnose_mingapwidth <- function(dmr_list,
                               candidate_gaps = c(200, 500, 1000, 2000)) {
  all_dmr <- do.call(c, unname(dmr_list))

  # 1. DMR size distribution  --------------------------------------------
  dmr_widths <- width(all_dmr)
  cat("=== DMR size distribution ===\n")
  cat(sprintf("  median : %d bp\n",  as.integer(median(dmr_widths))))
  cat(sprintf("  mean   : %d bp\n",  as.integer(mean(dmr_widths))))
  cat(sprintf("  Q1/Q3  : %d / %d bp\n",
              as.integer(quantile(dmr_widths, 0.25)),
              as.integer(quantile(dmr_widths, 0.75))))
  cat(sprintf("  90th   : %d bp\n",  as.integer(quantile(dmr_widths, 0.90))))
  cat(sprintf("  max    : %d bp\n",  max(dmr_widths)))

  # 2. Intergenic DMR gap distribution  -----------------------------------
  gaps <- unlist(lapply(split(all_dmr, seqnames(all_dmr)), function(chr_dmr) {
    if (length(chr_dmr) < 2) return(NULL)
    chr_sorted <- sort(chr_dmr)
    # Distance between adjacent DMRs
    as.integer(start(chr_sorted)[-1] - end(chr_sorted)[-length(chr_sorted)])
  }))
  gaps <- gaps[gaps > 0]   # exclude overlapping DMRs (gap=0)

  cat("\n=== Intergenic DMR gap distribution ===\n")
  cat(sprintf("  median : %d bp\n",  as.integer(median(gaps))))
  cat(sprintf("  Q1/Q3  : %d / %d bp\n",
              as.integer(quantile(gaps, 0.25)),
              as.integer(quantile(gaps, 0.75))))
  cat(sprintf("  10th   : %d bp\n",  as.integer(quantile(gaps, 0.10))))
  cat(sprintf("  5th    : %d bp\n",  as.integer(quantile(gaps, 0.05))))

  # 3. Simulation of candidate gap effects --------------------------------
  cat("\n=== mingapwidth candidate gap effects ===\n")
  cat(sprintf("%-12s %10s %10s %10s\n",
              "mingapwidth", "merged_n", "comparsion to original(%)", "median_w(bp)"))

  results <- lapply(candidate_gaps, function(g) {
    merged <- GenomicRanges::reduce(all_dmr,
                                    min.gapwidth  = g,
                                    ignore.strand = TRUE)
    data.frame(
      mingapwidth   = g,
      n_merged      = length(merged),
      pct_of_orig   = round(length(merged) / length(all_dmr) * 100, 1),
      median_width  = as.integer(median(width(merged)))
    )
  })

  res_df <- do.call(rbind, results)
  for (i in seq_len(nrow(res_df))) {
    cat(sprintf("%-12d %10d %10.1f %10d\n",
                res_df$mingapwidth[i],
                res_df$n_merged[i],
                res_df$pct_of_orig[i],
                res_df$median_width[i]))
  }

  # 4. Set cutoff based on statistics ----------------------------------------
  # Recommended gap: 10th percentile
  # 10% gap: high posibility of located in same regulatory unit
  recommended <- as.integer(quantile(gaps, 0.10))
  cat(sprintf("\nRecommended mingapwidth: %d bp (gap distribution 10th percentile)\n",
              recommended))

  invisible(list(
    dmr_widths  = dmr_widths,
    gaps        = gaps,
    candidates  = res_df,
    recommended = recommended
  ))
}

## Load datasets for each sample -----------------------------------------
setwd(opt$input_dir)
# Load DMRs for each sample and convert lists of data.frames into a GRangelist
dmr_list <- list.files(
  "DSS",
  pattern = "*.normal_vs_tumor_DMR.sorted.txt$",
  full.names = TRUE
)
per_pt_dmr <- lapply(dmr_list, function(f) {
  prefix <- str_remove(basename(f), "\\.normal_vs_tumor_DMR.sorted.txt")
  fread(f, header = TRUE) %>%
    mutate(sample = prefix) %>%
    dplyr::rename(normal.Methy = meanMethy1,
           tumor.Methy  = meanMethy2)
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  GRanges() %>%
  split(mcols(.)$sample)

per_pt_dmr <- endoapply(per_pt_dmr, filter_dmr_region)

# Load aDMR for each sample
admr_list <- list.files(
  "DSS", full.names = TRUE,
  pattern = "*tumor_aDMR.sorted.txt$",
)

per_pt_admr <- lapply(admr_list, function(f) {
  prefix <- str_remove(basename(f), ".tumor_aDMR.sorted.txt")
  fread(f, header = TRUE) %>%
    mutate(sample = prefix) %>%
    dplyr::rename(HP1.Methy  = meanMethy1,
           HP2.Methy  = meanMethy2)
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  GRanges() %>%
  split(mcols(.)$sample)
per_pt_admr <- endoapply(per_pt_admr, filter_dmr_region)

# Load aDMR for normal tissue (used to flag germline-like allelic imbalance)
normal_admr_files <- list.files(
  "DSS", full.names = TRUE,
  pattern = "*.normal_aDMR.sorted.txt$"
)
per_pt_normal_admr <- lapply(normal_admr_files, function(f) {
  prefix <- str_remove(basename(f), ".normal_aDMR.sorted.txt")
  fread(f, header = TRUE) %>%
    mutate(sample = prefix) %>%
    dplyr::rename(HP1.Methy = meanMethy1,
                  HP2.Methy = meanMethy2)
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  GRanges() %>%
  split(mcols(.)$sample)
per_pt_normal_admr <- endoapply(per_pt_normal_admr, filter_dmr_region)

## Parameters and Common Theme --------------------------------------------
PATIENT_IDS <- names(per_pt_dmr)

COLORS <- list(
  blue   = "#3B8BD4", green  = "#1D9E75",
  amber  = "#BA7517", red    = "#E24B4A",
  purple = "#7F77DD", gray   = "#888780"
)

NORM_ADMR_OVERLAP_PCT <- 0.10   # normal aDMR exclusion standard (based on tumor aDMR)

# Somatic aDMR: Exclude tumor aDMR regions that overlap with normal aDMR by ≥10% -----------
per_pt_somatic_admr <- lapply(PATIENT_IDS, function(pt) {
  tadmr     <- per_pt_admr[[pt]]
  norm_admr <- per_pt_normal_admr[[pt]]
  if (length(norm_admr) == 0) return(tadmr)
  hits <- findOverlaps(tadmr, norm_admr, ignore.strand = TRUE)
  if (length(hits) == 0) return(tadmr)
  ov_len  <- width(pintersect(tadmr[queryHits(hits)], norm_admr[subjectHits(hits)],
                               ignore.strand = TRUE))
  frac    <- ov_len / width(tadmr[queryHits(hits)])
  exclude <- unique(queryHits(hits)[frac >= NORM_ADMR_OVERLAP_PCT])
  tadmr[setdiff(seq_along(tadmr), exclude)]
}) %>% setNames(PATIENT_IDS)

{
  n_t <- sapply(PATIENT_IDS, function(pt) length(per_pt_admr[[pt]]))
  n_s <- sapply(PATIENT_IDS, function(pt) length(per_pt_somatic_admr[[pt]]))
  cat("\n[Somatic aDMR: tumor aDMR - normal aDMR (exclude ≥10% overlap, patient-specific)]\n")
  print(data.frame(
    patient_id     = PATIENT_IDS,
    n_tumor_admr   = n_t,
    n_somatic_admr = n_s,
    pct_removed    = round((1 - n_s / pmax(n_t, 1)) * 100, 1)
  ))
  cat(sprintf("\nTotal: tumor aDMR %d → somatic aDMR %d (%.1f%% removed)\n",
              sum(n_t), sum(n_s), (1 - sum(n_s) / sum(n_t)) * 100))

  somatic_df <- lapply(per_pt_somatic_admr, as.data.frame) %>% bind_rows()
  fwrite(somatic_df, file.path(opt$outdir, "somatic_admr_per_patient.csv.gz"),
         row.names = FALSE, quote = FALSE)
  cat(sprintf("저장: %s\n", file.path(opt$outdir, "somatic_admr_per_patient.csv.gz")))
  rm(n_t, n_s, somatic_df)
}

# STEP 1: tumor-normal DMR ∩ tumor aDMR in individual patients -----
# output: confident_dmr, layer_counts

## Reciprocal overlap of tumor-normal DMR ∩ tumor aDMR ---------------
# Parameters for reciprocal overlap
MIN_OVERLAP_BP  <- 100L     # minimum absolute overlap (bp)
RECIPROCAL_PCT  <- 0.30    # minimum 30% overlap in both regions

# dataframe of statistical raw dmr summary
layer_counts <- data.frame(
  patient_id   = PATIENT_IDS,
  n_tumor_dmr  = sapply(PATIENT_IDS, function(pt) length(per_pt_dmr[[pt]])),
  n_admr       = sapply(PATIENT_IDS, function(pt) length(per_pt_admr[[pt]]))
)

# Reciprocal overlap of tumor-normal DMR ∩ tumor aDMR in each patient
confident_dmr <- lapply(PATIENT_IDS, function(pt) {
  dmr  <- per_pt_dmr[[pt]]
  admr <- per_pt_admr[[pt]]

  # Step 1: reciprocal overlap
  hits <- reciprocal_overlap_hits(dmr, admr,
                                  min_pct = RECIPROCAL_PCT,
                                  min_bp  = MIN_OVERLAP_BP)
  if (length(hits) == 0) return(GRanges())

  out <- dmr[queryHits(hits)]

  # Add aDMR metadata
  mcols(out)$has_admr_support <- TRUE
  mcols(out)$admr_chr   <- as.character(seqnames(admr))[subjectHits(hits)]
  mcols(out)$admr_start <- start(admr)[subjectHits(hits)]
  mcols(out)$admr_end   <- end(admr)[subjectHits(hits)]
  mcols(out)$HP1.Methy  <- mcols(admr)$HP1.Methy[subjectHits(hits)]
  mcols(out)$HP2.Methy  <- mcols(admr)$HP2.Methy[subjectHits(hits)]
  mcols(out)$admr.stat  <- mcols(admr)$areaStat[subjectHits(hits)]

  out
}) %>% setNames(PATIENT_IDS)

# Flag confident DMRs overlapping normal aDMR (≥10% of tumor aDMR width → possible germline)
confident_dmr <- lapply(PATIENT_IDS, function(pt) {
  gr        <- confident_dmr[[pt]]
  norm_admr <- per_pt_normal_admr[[pt]]
  if (length(gr) == 0) {
    mcols(gr)$ov_norm_admr <- logical(0)
    return(gr)
  }
  if (is.null(norm_admr) || length(norm_admr) == 0) {
    mcols(gr)$ov_norm_admr <- FALSE
    return(gr)
  }
  # confident_dmr coordiate based on tumor-normal DMR; aDMR coordinates in admr_* columns
  tumor_admr_gr <- GRanges(
    seqnames = mcols(gr)$admr_chr,
    ranges   = IRanges(mcols(gr)$admr_start, mcols(gr)$admr_end)
  )
  hits    <- findOverlaps(tumor_admr_gr, norm_admr, ignore.strand = TRUE)
  ov_len  <- width(pintersect(tumor_admr_gr[queryHits(hits)], norm_admr[subjectHits(hits)],
                               ignore.strand = TRUE))
  frac    <- ov_len / width(tumor_admr_gr[queryHits(hits)])
  flagged <- unique(queryHits(hits)[frac >= NORM_ADMR_OVERLAP_PCT])
  mcols(gr)$ov_norm_admr <- seq_along(gr) %in% flagged
  gr
}) %>% setNames(PATIENT_IDS)

# Report ov_norm_admr distribution
{
  n_conf <- sapply(PATIENT_IDS, function(pt) length(confident_dmr[[pt]]))
  n_flag <- sapply(PATIENT_IDS, function(pt)
    sum(mcols(confident_dmr[[pt]])$ov_norm_admr, na.rm = TRUE))
  cat("\n[ov_norm_admr distribution (aDMR_normal ≥10% overlap, patient-specific)]\n")
  print(data.frame(
    patient_id     = PATIENT_IDS,
    n_confident    = n_conf,
    n_ov_norm_admr = n_flag,
    pct_ov         = round(n_flag / pmax(n_conf, 1) * 100, 1)
  ))
  all_tmp <- do.call(c, unname(confident_dmr))
  cat(sprintf("\nTotal: ov_norm_admr=TRUE %d / %d (%.1f%%)\n",
              sum(mcols(all_tmp)$ov_norm_admr),
              length(all_tmp),
              sum(mcols(all_tmp)$ov_norm_admr) / length(all_tmp) * 100))
  rm(all_tmp, n_conf, n_flag)
}

# Extract best-match normal aDMR on ov_norm_admr=TRUE pairs
{
  diag_list <- lapply(PATIENT_IDS, function(pt) {
    gr        <- confident_dmr[[pt]]
    norm_admr <- per_pt_normal_admr[[pt]]
    if (length(gr) == 0 || is.null(norm_admr) || length(norm_admr) == 0) return(NULL)

    tumor_admr_gr <- GRanges(
      seqnames = mcols(gr)$admr_chr,
      ranges   = IRanges(mcols(gr)$admr_start, mcols(gr)$admr_end)
    )
    hits   <- findOverlaps(tumor_admr_gr, norm_admr, ignore.strand = TRUE)
    ov_len <- width(pintersect(tumor_admr_gr[queryHits(hits)], norm_admr[subjectHits(hits)],
                               ignore.strand = TRUE))
    frac   <- ov_len / width(tumor_admr_gr[queryHits(hits)])

    keep  <- frac >= NORM_ADMR_OVERLAP_PCT
    q_idx <- queryHits(hits)[keep]
    s_idx <- subjectHits(hits)[keep]
    f_val <- frac[keep]
    if (length(q_idx) == 0) return(NULL)

    # Select one best-match normal aDMR for each tumor aDMR with maximum overlap
    best_i <- tapply(seq_along(q_idx), q_idx, function(i) i[which.max(f_val[i])])
    bq     <- as.integer(names(best_i))
    bs     <- s_idx[unlist(best_i)]

    db_t <- mcols(gr)$HP1.Methy[bq]  - mcols(gr)$HP2.Methy[bq]
    db_n <- mcols(norm_admr)$HP1.Methy[bs] - mcols(norm_admr)$HP2.Methy[bs]

    data.frame(
      patient_id     = pt,
      delta_tumor    = db_t,
      delta_normal   = db_n,
      abs_diff       = abs(db_t) - abs(db_n),   # >0: tumor에서 증폭
      direction_flip = sign(db_t) != sign(db_n)
    )
  })
  diag_df <- do.call(rbind, diag_list)

  cat("\n[P1: HP |Δβ| change (|tumor| - |normal|)]\n")
  cat(sprintf("  median |Δβ|_tumor : %.3f\n",  median(abs(diag_df$delta_tumor),  na.rm = TRUE)))
  cat(sprintf("  median |Δβ|_normal: %.3f\n",  median(abs(diag_df$delta_normal), na.rm = TRUE)))
  cat(sprintf("  median diff (T-N) : %.3f\n",  median(diag_df$abs_diff, na.rm = TRUE)))
  for (thr in c(0.10, 0.15, 0.20)) {
    n <- sum(diag_df$abs_diff > thr, na.rm = TRUE)
    cat(sprintf("  amplified >%.2f : %d / %d (%.1f%%)\n",
                thr, n, nrow(diag_df), n / nrow(diag_df) * 100))
  }

  cat("\n[P2: HP direction flip (tumor vs normal)]\n")
  n_flip <- sum(diag_df$direction_flip, na.rm = TRUE)
  cat(sprintf("  flip=TRUE : %d / %d (%.1f%%)\n",
              n_flip, nrow(diag_df), n_flip / nrow(diag_df) * 100))
  cat(sprintf("  flip=FALSE: %d / %d (%.1f%%)\n",
              nrow(diag_df) - n_flip, nrow(diag_df),
              (nrow(diag_df) - n_flip) / nrow(diag_df) * 100))

  cat("\n[P1 & P2 (flip=TRUE AND abs_diff>0.15)]\n")
  n_both <- sum(diag_df$direction_flip & diag_df$abs_diff > 0.15, na.rm = TRUE)
  cat(sprintf("  %d / %d (%.1f%%)\n", n_both, nrow(diag_df), n_both / nrow(diag_df) * 100))

  rm(diag_list, diag_df)
}

# Summary
layer_counts$n_confident <- sapply(confident_dmr, length)
layer_counts$pct_dmr_with_admr <- round(
  layer_counts$n_confident / pmax(layer_counts$n_tumor_dmr, 1) * 100, 1
)
layer_counts$pct_admr_with_dmr <- round(
  layer_counts$n_confident / pmax(layer_counts$n_admr, 1) * 100, 1
)

cat("\n[Stratify DMR # (patient-specific)]\n")
print(layer_counts)
cat(sprintf(
  "\n[Cohort Summary]\n  tumor-normal DMR: total %d (Average %.0f per patient)\n",
  sum(layer_counts$n_tumor_dmr), mean(layer_counts$n_tumor_dmr)
))
cat(sprintf(
  "  aDMR:             total %d (Average %.0f per patient)\n",
  sum(layer_counts$n_admr), mean(layer_counts$n_admr)
))
cat(sprintf(
  "  allelic imbalance DMR:    total %d (tumor-normal DMR's average %.1f%% support aDMR)\n",
  sum(layer_counts$n_confident), mean(layer_counts$pct_dmr_with_admr)
))

# STEP 2: Cohort-level recurrence analysis of confident DMRs --------
## Merge confident DMRs to make standardized regions ----------------

per_pt_confident_distinct <- confident_dmr %>%
  endoapply(function(gr) {
    GRanges(
      seqnames = seqnames(gr),
      ranges   = IRanges(start(gr), end(gr)),
      strand   = strand(gr),
      sample   = mcols(gr)$sample
    )
 })

# Diagnose mingapwidth for merging confident DMRs
diag <- diagnose_mingapwidth(
  dmr_list       = per_pt_confident_distinct,
  candidate_gaps = c(200, 500, 1000, 2000)
)

# Standardize merged regions across all patients
all_confident <- do.call(c, unname(per_pt_confident_distinct))

merged_regions <- GenomicRanges::reduce(
  all_confident,
  ignore.strand = TRUE,
  min.gapwidth = 500
)

## Count Recurrence between merged regions 
merged_regions$n_patients <- count_recurrence(
  merged_regions, per_pt_confident_distinct, min_bp = MIN_OVERLAP_BP
)
merged_regions$pct_patients <- merged_regions$n_patients / length(PATIENT_IDS) * 100

# Step 2: Add recurrence column on confident_dmr --------------------

confident_dmr <- lapply(PATIENT_IDS, function(pt) {
  gr <- confident_dmr[[pt]]
  if (length(gr) == 0) return(gr)

  # Count overlap on other patients' DMR
  other_pts <- setdiff(PATIENT_IDS, pt)

  recurrence_count <- rowSums(vapply(other_pts, function(other) {
    as.integer(countOverlaps(gr, confident_dmr[[other]],
                             minoverlap = MIN_OVERLAP_BP) > 0)
  }, integer(length(gr))))

  # Count recurrence including self (n_patients)
  mcols(gr)$n_patients <- recurrence_count + 1L
  mcols(gr)$pct_patients <- mcols(gr)$n_patients / length(PATIENT_IDS) * 100

  gr
})
names(confident_dmr) <- PATIENT_IDS

# 확인
cat("recurrence distribution (Total patients pooled):\n")
all_annotated <- do.call(c, unname(confident_dmr))
print(table(mcols(all_annotated)$n_patients))

confident_dmr_df <- lapply(confident_dmr, function(gr) {
  as.data.frame(gr)
}) %>% bind_rows()

fwrite(
  confident_dmr_df,
  file.path(opt$outdir, "confident_dmr_per_patient.csv.gz")
)

rm(per_pt_admr, per_pt_dmr, per_pt_normal_admr, per_pt_somatic_admr,
   diag, all_confident, confident_dmr_df)
gc()

# STEP 3: Check independency of pooled DMR & state consensus dmr =========
pooled_dmr <- fread("DSS/DMR/total.hap.DMR.csv.gz", nThread = 4) %>%
  GRanges() %>%
  filter_dmr_region()


## Explore recurrence cutoff --------------------------------------
RECOMMENDED_N   <- max(3L, ceiling(length(PATIENT_IDS) * 0.25))  # 25%
DIR_CONSISTENCY_CUTOFF <- 0.80   # Direction consistency cutoff: 80%

# 3-B. Final Consensus DMR based on DSS coordinates ------------
#   Filter 1: n_patients ≥ RECOMMENDED_N
#   Filter 2:  reciprocal overlap with pooled DMR
#   Filter 3: direction consistency (Optional, pooled DMR & mean_diff direction concordance into new column)

consensus_dmr <- lapply(PATIENT_IDS, function(pt) {
  gr <- confident_dmr[[pt]]
  if (length(gr) == 0) return(gr)

  # Filter 1
  gr <- gr[mcols(gr)$n_patients >= RECOMMENDED_N]
  if (length(gr) == 0) return(gr)

  # Filter 2
  hits <- reciprocal_overlap_hits(gr, pooled_dmr,
                                  min_pct = RECIPROCAL_PCT,
                                  min_bp  = MIN_OVERLAP_BP)
  if (length(hits) == 0) return(GRanges())

  # Filter 3 (Optional)
  dir_hits <- direction_consistent_hits(hits, gr, pooled_dmr)
  if (length(dir_hits) == 0) return(GRanges())
  mcols(gr)$pooled_consistency <- FALSE
  mcols(gr)$pooled_consistency[queryHits(dir_hits)] <- TRUE
  gr %>% as.data.frame()
}) %>%
  bind_rows() %>%
  GRanges()  %>%
  split(mcols(.)$sample)

all_consensus <- do.call(c, unname(consensus_dmr))
cat(sprintf("Direction consistency distribution (pooled DMR and direction consistency):
  %i consistent, %i inconsistent, NA %i\n",
            sum(mcols(all_consensus)$pooled_consistency, na.rm = TRUE),
            sum(!mcols(all_consensus)$pooled_consistency, na.rm = TRUE),
            sum(is.na(mcols(all_consensus)$pooled_consistency))))

## 3-C. Funnel Statistics summary ---------------------------------------
funnel_stats <- data.frame(
  step = c(
    "1. tumor-normal DMR",
    "2. Confident DMR (aDMR reciprocal overlap)",
    "3. Recurrence n ≥ threshold",
    "4. Pooled DMR Validation (Consensus DMR)"
  ),
  n_dmr = c(
    sum(layer_counts$n_tumor_dmr),
    sum(layer_counts$n_confident),
    sum(mcols(all_annotated)$n_patients >= RECOMMENDED_N),
    length(all_consensus)
  ),
  note = c(
    "Patient-specific sum (with duplicates)",
    sprintf("reciprocal ≥%.0f%%, ≥%d bp | no direction filter",
            RECIPROCAL_PCT * 100, MIN_OVERLAP_BP),
    sprintf("n ≥ %d (%.0f%%)", RECOMMENDED_N,
            RECOMMENDED_N / length(PATIENT_IDS) * 100),
    sprintf("reciprocal ≥%.0f%% vs pooled DMR | DSS original coordinates",
            RECIPROCAL_PCT * 100)
  )
)

cat("\n[Funnel Statistics Summary]\n")
print(funnel_stats, row.names = FALSE)
cat(sprintf("\nFinal Consensus DMR: %d\n", length(all_consensus)))

final_df <- as.data.frame(all_consensus)
fwrite(final_df, file.path(opt$outdir, "consensus_dmrs_per_patient.csv.gz"),
       row.names = FALSE, quote = FALSE)
fwrite(funnel_stats, file.path(opt$outdir, "dmr_funnel_summary.csv.gz"),
       row.names = FALSE, quote = FALSE)
fwrite(layer_counts, file.path(opt$outdir, "per_patient_layer_counts.csv.gz"),
       row.names = FALSE, quote = FALSE)

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Final Consensus DMR: %d\n", nrow(final_df)))
cat(sprintf("  Conditions: recurrence n≥%d + direction consistency ≥%.0f%% + pooled DMR reciprocal overlap\n",
            RECOMMENDED_N, DIR_CONSISTENCY_CUTOFF * 100))

cat("\n=== Analysis Parameters Summary ===\n")
cat(sprintf("Number of patients: %d\n", length(PATIENT_IDS)))
cat(sprintf("Minimum reciprocal overlap: %.0f%%\n", RECIPROCAL_PCT * 100))
cat(sprintf("Minimum absolute overlap: %d bp\n", MIN_OVERLAP_BP))
cat(sprintf("Recurrence threshold: n ≥ %d (%.0f%%)\n",
            RECOMMENDED_N,
            RECOMMENDED_N / length(PATIENT_IDS) * 100))
cat(sprintf("Merged region min gap: 500 bp\n"))

cat("\n[Funnel Statistics Summary]\n")
print(funnel_stats[, c("step", "n_dmr")], row.names = FALSE)

# Output final results
# Confident DMR: confident_dmr_per_patient.csv.gz
# Consensus DMR: consensus_dmrs_per_patient.csv.gz
# Funnel summary: dmr_funnel_summary.csv.gz
# Layer counts: per_patient_layer_counts.csv.gz
