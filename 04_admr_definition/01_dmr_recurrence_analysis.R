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
  make_option(c("-i", "--input_dir"), type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_minimap2.out_hg38/",
              help = "각 환자별 후보 DMR 정보가 저장된 디렉토리 [default: %default]"),
  make_option(c("-o", "--outdir"),      type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/", #nolint
              help = "출력 디렉토리 [default: %default]")
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

  # ── 1. DMR 크기 분포 ────────────────────────────────────────────────────────
  dmr_widths <- width(all_dmr)
  cat("=== DMR 크기 분포 ===\n")
  cat(sprintf("  median : %d bp\n",  as.integer(median(dmr_widths))))
  cat(sprintf("  mean   : %d bp\n",  as.integer(mean(dmr_widths))))
  cat(sprintf("  Q1/Q3  : %d / %d bp\n",
              as.integer(quantile(dmr_widths, 0.25)),
              as.integer(quantile(dmr_widths, 0.75))))
  cat(sprintf("  90th   : %d bp\n",  as.integer(quantile(dmr_widths, 0.90))))
  cat(sprintf("  max    : %d bp\n",  max(dmr_widths)))

  # ── 2. 염색체 내 DMR 간 gap 분포 ────────────────────────────────────────────
  gaps <- unlist(lapply(split(all_dmr, seqnames(all_dmr)), function(chr_dmr) {
    if (length(chr_dmr) < 2) return(NULL)
    chr_sorted <- sort(chr_dmr)
    # 인접 DMR 간 거리
    as.integer(start(chr_sorted)[-1] - end(chr_sorted)[-length(chr_sorted)])
  }))
  gaps <- gaps[gaps > 0]   # 겹치는 경우 제외

  cat("\n=== 인접 DMR 간 gap 분포 ===\n")
  cat(sprintf("  median : %d bp\n",  as.integer(median(gaps))))
  cat(sprintf("  Q1/Q3  : %d / %d bp\n",
              as.integer(quantile(gaps, 0.25)),
              as.integer(quantile(gaps, 0.75))))
  cat(sprintf("  10th   : %d bp\n",  as.integer(quantile(gaps, 0.10))))
  cat(sprintf("  5th    : %d bp\n",  as.integer(quantile(gaps, 0.05))))

  # ── 3. candidate gap별 합쳐지는 DMR 비율 시뮬레이션 ─────────────────────────
  cat("\n=== mingapwidth 후보별 영향 ===\n")
  cat(sprintf("%-12s %10s %10s %10s\n",
              "mingapwidth", "merged_n", "원본대비(%)", "median_w(bp)"))

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

  # ── 4. 권장값 자동 탐지 ──────────────────────────────────────────────────────
  # gap 분포의 10th percentile을 권장 mingapwidth로 사용
  # 근거: 하위 10% gap은 같은 regulatory unit 내 조각일 가능성이 높음
  recommended <- as.integer(quantile(gaps, 0.10))
  cat(sprintf("\n권장 mingapwidth: %d bp (gap 분포 10th percentile)\n",
              recommended))
  cat("근거: 인접 DMR의 하위 10% gap은 같은 regulatory unit 내 단편화로 판단\n")

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

NORM_ADMR_OVERLAP_PCT <- 0.10   # normal aDMR 제거 기준 (tumor aDMR 너비 대비)

# Somatic aDMR: tumor aDMR에서 normal aDMR과 ≥10% 겹치는 영역 제거 -----------
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
  cat("\n[Somatic aDMR: tumor aDMR - normal aDMR (≥10% overlap 제거, 환자별)]\n")
  print(data.frame(
    patient_id     = PATIENT_IDS,
    n_tumor_admr   = n_t,
    n_somatic_admr = n_s,
    pct_removed    = round((1 - n_s / pmax(n_t, 1)) * 100, 1)
  ))
  cat(sprintf("\n전체: tumor aDMR %d → somatic aDMR %d (%.1f%% 제거)\n",
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
MIN_OVERLAP_BP  <- 100L     # 최소 절대 overlap (bp)
RECIPROCAL_PCT  <- 0.30    # 양쪽 영역의 최소 30% 겹침

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
  # confident_dmr 좌표는 tumor-normal DMR; aDMR 좌표는 admr_* 컬럼에 저장됨
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
  cat("\n[ov_norm_admr 분포 (aDMR_normal ≥10% overlap, 환자별)]\n")
  print(data.frame(
    patient_id     = PATIENT_IDS,
    n_confident    = n_conf,
    n_ov_norm_admr = n_flag,
    pct_ov         = round(n_flag / pmax(n_conf, 1) * 100, 1)
  ))
  all_tmp <- do.call(c, unname(confident_dmr))
  cat(sprintf("\n전체: ov_norm_admr=TRUE %d / %d (%.1f%%)\n",
              sum(mcols(all_tmp)$ov_norm_admr),
              length(all_tmp),
              sum(mcols(all_tmp)$ov_norm_admr) / length(all_tmp) * 100))
  rm(all_tmp, n_conf, n_flag)
}

# 두 관점 진단: ov_norm_admr=TRUE 쌍에서 best-match normal aDMR 추출
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

    # 각 tumor aDMR에서 overlap 가장 큰 normal aDMR 1개 선택
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

  cat("\n[관점 1: HP |Δβ| 변화량 (|tumor| - |normal|)]\n")
  cat(sprintf("  median |Δβ|_tumor : %.3f\n",  median(abs(diag_df$delta_tumor),  na.rm = TRUE)))
  cat(sprintf("  median |Δβ|_normal: %.3f\n",  median(abs(diag_df$delta_normal), na.rm = TRUE)))
  cat(sprintf("  median diff (T-N) : %.3f\n",  median(diag_df$abs_diff, na.rm = TRUE)))
  for (thr in c(0.10, 0.15, 0.20)) {
    n <- sum(diag_df$abs_diff > thr, na.rm = TRUE)
    cat(sprintf("  amplified >%.2f : %d / %d (%.1f%%)\n",
                thr, n, nrow(diag_df), n / nrow(diag_df) * 100))
  }

  cat("\n[관점 2: HP 방향 전환 (direction flip, tumor vs normal)]\n")
  n_flip <- sum(diag_df$direction_flip, na.rm = TRUE)
  cat(sprintf("  flip=TRUE : %d / %d (%.1f%%)\n",
              n_flip, nrow(diag_df), n_flip / nrow(diag_df) * 100))
  cat(sprintf("  flip=FALSE: %d / %d (%.1f%%)\n",
              nrow(diag_df) - n_flip, nrow(diag_df),
              (nrow(diag_df) - n_flip) / nrow(diag_df) * 100))

  cat("\n[두 조건 동시 (flip=TRUE AND abs_diff>0.15)]\n")
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

cat("\n[층위별 DMR 수 (환자별)]\n")
print(layer_counts)
cat(sprintf(
  "\n[코호트 전체 요약]\n  tumor-normal DMR: 합계 %d개 (환자당 평균 %.0f개)\n",
  sum(layer_counts$n_tumor_dmr), mean(layer_counts$n_tumor_dmr)
))
cat(sprintf(
  "  aDMR:             합계 %d개 (환자당 평균 %.0f개)\n",
  sum(layer_counts$n_admr), mean(layer_counts$n_admr)
))
cat(sprintf(
  "  allelic imbalance DMR:    합계 %d개 (tumor-normal DMR의 평균 %.1f%%가 aDMR 지지)\n",
  sum(layer_counts$n_confident), mean(layer_counts$pct_dmr_with_admr)
))

# STEP 2: Cohort-level recurrence analysis of confident DMRs --------
## Merge confident DMRs to make standardized regions ----------------
# 환자별로 메타데이터 일부 제거 후 중복 제거
per_pt_confident_distinct <- confident_dmr %>%
  endoapply(function(gr) {
    GRanges(
      seqnames = seqnames(gr),
      ranges   = IRanges(start(gr), end(gr)),
      strand   = strand(gr),
      sample   = mcols(gr)$sample
    )
 })

# 실행
diag <- diagnose_mingapwidth(
  dmr_list       = per_pt_confident_distinct,           # Step 1 산출물
  candidate_gaps = c(200, 500, 1000, 2000)
)

# 전체 환자들에서 500bp 이내로 병합하여 표준화된 영역 생성
all_confident <- do.call(c, unname(per_pt_confident_distinct))

merged_regions <- GenomicRanges::reduce(
  all_confident,
  ignore.strand = TRUE,
  min.gapwidth = 500
)

## Count Recurrence between merged regions --------------------------
merged_regions$n_patients <- count_recurrence(
  merged_regions, per_pt_confident_distinct, min_bp = MIN_OVERLAP_BP
)
merged_regions$pct_patients <- merged_regions$n_patients / length(PATIENT_IDS) * 100

# Step 2: confident_dmr에 recurrence metadata 추가
# 각 환자의 DMR이 다른 환자들과 얼마나 겹치는지 카운트

confident_dmr <- lapply(PATIENT_IDS, function(pt) {
  gr <- confident_dmr[[pt]]
  if (length(gr) == 0) return(gr)

  # 다른 환자들의 DMR과 overlap 카운트
  other_pts <- setdiff(PATIENT_IDS, pt)

  recurrence_count <- rowSums(vapply(other_pts, function(other) {
    as.integer(countOverlaps(gr, confident_dmr[[other]],
                             minoverlap = MIN_OVERLAP_BP) > 0)
  }, integer(length(gr))))

  # 자기 자신 포함한 recurrence (n_patients)
  mcols(gr)$n_patients <- recurrence_count + 1L
  mcols(gr)$pct_patients <- mcols(gr)$n_patients / length(PATIENT_IDS) * 100

  gr
})
names(confident_dmr) <- PATIENT_IDS

# 확인
cat("recurrence 분포 (전체 환자 pooled):\n")
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

# STEP 3: Pooled DMR 독립 검증, consensus dmr 정의 ======================
pooled_dmr <- fread("DSS/DMR/total.hap.DMR.csv.gz", nThread = 4) %>%
  GRanges() %>%
  filter_dmr_region()


## Recurrence cutoff 탐색 --------------------------------------
RECOMMENDED_N   <- max(3L, ceiling(length(PATIENT_IDS) * 0.25))  # 25%
DIR_CONSISTENCY_CUTOFF <- 0.80   # recurrent 환자 중 80% 이상 동일 방향

# ── 3-B. 최종 Consensus DMR 확정 (DSS 원본 좌표 기반) ───────────────────────
#   조건 1: n_patients ≥ RECOMMENDED_N
#   조건 2: pooled DMR과 reciprocal overlap
#   조건 3: 방향 일관성 (선택, pooled DMR과 mean_diff 방향 일치 여부 metadata로 삽입)

consensus_dmr <- lapply(PATIENT_IDS, function(pt) {
  gr <- confident_dmr[[pt]]
  if (length(gr) == 0) return(gr)

  # 조건 1
  gr <- gr[mcols(gr)$n_patients >= RECOMMENDED_N]
  if (length(gr) == 0) return(gr)

  # 조건 2
  hits <- reciprocal_overlap_hits(gr, pooled_dmr,
                                  min_pct = RECIPROCAL_PCT,
                                  min_bp  = MIN_OVERLAP_BP)
  if (length(hits) == 0) return(GRanges())

  # 조건 3 (옵션)
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
cat(sprintf("Direction consistency 분포 (pooled DMR과 방향 일치 여부):
  %i개 일치, %i개 불일치, NA %i개\n",
            sum(mcols(all_consensus)$pooled_consistency, na.rm = TRUE),
            sum(!mcols(all_consensus)$pooled_consistency, na.rm = TRUE),
            sum(is.na(mcols(all_consensus)$pooled_consistency))))

## 3-C. Funnel 통계 요약 ---------------------------------------
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
    "환자별 합산 (중복 포함)",
    sprintf("reciprocal ≥%.0f%%, ≥%d bp | direction 필터 없음",
            RECIPROCAL_PCT * 100, MIN_OVERLAP_BP),
    sprintf("n ≥ %d (%.0f%%)", RECOMMENDED_N,
            RECOMMENDED_N / length(PATIENT_IDS) * 100),
    sprintf("reciprocal ≥%.0f%% vs pooled DMR | DSS 원본 좌표",
            RECIPROCAL_PCT * 100)
  )
)

cat("\n[Funnel 요약]\n")
print(funnel_stats, row.names = FALSE)
cat(sprintf("\n최종 Consensus DMR: %d개\n", length(all_consensus)))

final_df <- as.data.frame(all_consensus)
fwrite(final_df, file.path(opt$outdir, "consensus_dmrs_per_patient.csv.gz"),
       row.names = FALSE, quote = FALSE)
fwrite(funnel_stats, file.path(opt$outdir, "dmr_funnel_summary.csv.gz"),
       row.names = FALSE, quote = FALSE)
fwrite(layer_counts, file.path(opt$outdir, "per_patient_layer_counts.csv.gz"),
       row.names = FALSE, quote = FALSE)

cat("\n=== 분석 완료 ===\n")
cat(sprintf("최종 Consensus DMR: %d개\n", nrow(final_df)))
cat(sprintf("  조건: recurrence n≥%d + 방향 일관성 ≥%.0f%% + pooled DMR reciprocal overlap\n",
            RECOMMENDED_N, DIR_CONSISTENCY_CUTOFF * 100))

cat("\n=== 분석 파라미터 요약 ===\n")
cat(sprintf("환자 수: %d명\n", length(PATIENT_IDS)))
cat(sprintf("Reciprocal overlap 최소 비율: %.0f%%\n", RECIPROCAL_PCT * 100))
cat(sprintf("최소 절대 overlap: %d bp\n", MIN_OVERLAP_BP))
cat(sprintf("Recurrence threshold: n ≥ %d (%.0f%%)\n",
            RECOMMENDED_N,
            RECOMMENDED_N / length(PATIENT_IDS) * 100))
cat(sprintf("Merged region min gap: 500 bp\n"))

cat("\n[Funnel 최종 요약]\n")
print(funnel_stats[, c("step", "n_dmr")], row.names = FALSE)

# Output final results
# Confident DMR: confident_dmr_per_patient.csv.gz
# Consensus DMR: consensus_dmrs_per_patient.csv.gz
# Funnel summary: dmr_funnel_summary.csv.gz
# Layer counts: per_patient_layer_counts.csv.gz
# 최종 결과물은 opt$outdir에 저장됨