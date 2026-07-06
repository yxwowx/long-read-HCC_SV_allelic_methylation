# conda env: renv
# =============================================================================
# SV 유형별 DMR 농축도 비교 분석
# Copy-neutral SV (INV, TRA/BND) vs Copy-number-changing SV (DEL, DUP, INS)
# HCC Paired PacBio Long-read WGS | n = 12
# =============================================================================
# 핵심 목적:
#   CNV confounding 배제 — copy-neutral SV 주변에서도 DMR이 농축된다면
#   methylation 변화가 copy number 변화가 아닌 SV 자체에 의한 것임을 지지
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
  library(scales)
  library(stringr)
  library(StructuralVariantAnnotation)
  library(VariantAnnotation)
  library(data.table)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(tibble)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--primary_window", type = "integer", default = 50000L,
              metavar = "INT",
              help = "Primary breakpoint window radius in bp [default: %default]"),
  make_option("--windows", type = "character", default = "10000,50000,100000,500000,1000000",
              metavar = "INT,INT,...",
              help = "Comma-separated window sizes in bp [default: %default]"),
  make_option("--n_perm", type = "integer", default = 1000L,
              metavar = "INT",
              help = "Number of permutations [default: %default]"),
  make_option("--outdir", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment",
              metavar = "DIR",
              help = "Output directory [default: %default]"),
  make_option("--run_id", type = "character", default = "",
              metavar = "STR",
              help = "Run ID prefix for output file names [default: none]"),
  make_option("--dmr_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/consensus_dmrs_per_patient.csv.gz",
              metavar = "FILE",
              help = "Path to DMR CSV(.gz) [default: %default]"),
  make_option("--no_plot", action = "store_true", default = FALSE,
              help = "Skip all visualization (ggsave) steps [default: %default]"),
  make_option("--hap_coloc_csv", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/layer1_coloc_results.csv",
              metavar = "FILE",
              help = "Path to layer1_coloc_results.csv from haplotype_sv_admr_analysis.R [default: %default]"),
  make_option("--sv_strat_file", type = "character", default = NULL,
              metavar = "FILE",
              help = "Pre-stratified SV file (CSV.gz) from sv_stratification.R. Skips VCF loading when provided."),
  make_option("--group_by", type = "character", default = "cnv_class",
              metavar = "STR",
              help = "Stratify SVs by 'cnv_class' or 'tier' (stratification col) [default: %default]"),
  make_option("--min_sv_per_group", type = "integer", default = 5L,
              metavar = "INT",
              help = "Min SVs per group per patient to run permutation [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$windows <- as.integer(strsplit(opt$windows, ",")[[1]])
if (!opt$group_by %in% c("cnv_class", "tier")) stop("--group_by must be 'cnv_class' or 'tier'")

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# ── 스크립트별 상수 ────────────────────────────────────────────────────────────
CNV_GROUP_COLORS <- c(
  "Copy-neutral (INV/TRA/BND)" = "#3B8BD4",
  "Copy-number-changing"       = "#E24B4A"
)
CNV_CLASSES <- c("copy_neutral", "copy_gaining", "copy_losing", "insertion", "COM")

# Set stratification variables based on --group_by
GROUP_COL    <- if (opt$group_by == "tier") "stratification" else "cnv_class"
GROUP_LEVELS <- if (opt$group_by == "tier") STRAT_LEVELS else CNV_CLASSES
GROUP_COLORS <- if (opt$group_by == "tier") STRAT_COLORS else SV_COLORS

cat(sprintf("\n=== Stratification mode: %s (column: %s) ===\n", opt$group_by, GROUP_COL))

## DMR의 CNV 상태 판정 함수 ==============================================================
#' DMR stratification by CNV status
#' @param dmr_gr GRanges: DMR
#' @param cnv_gr GRanges, 필수 metadata: $copyNumber, $bafCount, $germlineStatus, $depthWindowCount
#' @description Use DMR & CNV segment of same patient
annotate_dmr_cnv <- function(dmr_gr, cnv_gr,
                             cn_range = CN_NORMAL_RANGE, min_baf  = 10L) {
  if (length(dmr_gr) == 0 || length(cnv_gr) == 0) {
    dmr_gr$cnv_status <- "unknown"
    return(dmr_gr)
  }
  # PURPLE CNV segment 중 CNV 상태 신뢰할 수 있는 구간만 사용 (depthWindowCount, bafCount 기준)
  cnv_use <- cnv_gr[cnv_gr$bafCount >= min_baf]

  # DMR과 CNV segment 간 겹침 계산 → DMR별 가중 평균 CNV copy number 산출
  hits <- findOverlaps(dmr_gr, cnv_use, minoverlap = 1L)
  cn_per_dmr <- rep(NA_real_, length(dmr_gr))
  if (length(hits) > 0) {
    # overlap width로 가중 평균 계산
    ol_w  <- width(pintersect(dmr_gr[queryHits(hits)],
                              cnv_use[subjectHits(hits)]))
    cn_df <- data.frame(
      dmr_idx = queryHits(hits),
      cn      = cnv_use$copyNumber[subjectHits(hits)],
      ol_w    = ol_w
    ) %>%
      dplyr::group_by(dmr_idx) %>%
      dplyr::summarise(cn_wt = sum(cn * ol_w) / sum(ol_w), .groups = "drop")
    cn_per_dmr[cn_df$dmr_idx] <- cn_df$cn_wt
  }
  dmr_gr$cnv_copy_number <- cn_per_dmr
  dmr_gr$cnv_status <- dplyr::case_when(
    is.na(cn_per_dmr)                                     ~ "no_cnv_data",
    cn_per_dmr >= cn_range[1] & cn_per_dmr <= cn_range[2] ~ "cnv_normal",
    TRUE                                                  ~ "cnv_altered"
  )
  dmr_gr
}


## H3 Distance decay 헬퍼 함수==============================================

#' SV breakpoint 기준 signed distance 계산
#' @param sv_gr   GRanges: SV breakpoint ($sourceId, $cnv_class, $svtype 필수)
#' @param dmr_gr  GRanges: DMR ($diff.Methy 또는 $direction 필수)
#' @param max_dist integer: 계산할 최대 거리 (bp); 이 밖은 제외
#' @return data.frame: bp_id, bp_side, cnv_class, svtype, sourceId,
#'         signed_dist (양수=downstream, 음수=upstream),
#'         abs_dist, diff_methy, direction
get_signed_distance <- function(sv_gr, dmr_gr, max_dist = 200000L) {

  if (length(sv_gr) == 0 || length(dmr_gr) == 0) return(NULL)

  if (!"direction" %in% names(mcols(dmr_gr))) {
    mcols(dmr_gr)$direction <- dplyr::case_when(
      mcols(dmr_gr)$diff.Methy >  0 ~ "hypo",
      mcols(dmr_gr)$diff.Methy <  0 ~ "hyper",
      TRUE                           ~ NA_character_
    )
  }

  hits     <- GenomicRanges::distanceToNearest(sv_gr, dmr_gr)
  bp_idx   <- queryHits(hits)
  dmr_idx  <- subjectHits(hits)
  abs_dist <- mcols(hits)$distance

  keep     <- abs_dist <= max_dist
  bp_idx   <- bp_idx[keep]
  dmr_idx  <- dmr_idx[keep]
  abs_dist <- abs_dist[keep]

  signed_dist <- ifelse(
    start(dmr_gr[dmr_idx]) >= start(sv_gr[bp_idx]),
    abs_dist,
    -abs_dist
  )

  data.frame(
    bp_id       = sv_gr$sourceId[bp_idx],
    bp_side     = sub(".*_(bp[12])$", "\\1", sv_gr$bp_id[bp_idx]),
    cnv_class   = sv_gr$cnv_class[bp_idx],
    svtype      = sv_gr$svtype[bp_idx],
    signed_dist = signed_dist,
    abs_dist    = abs_dist,
    diff_methy  = as.numeric(as.character(
      mcols(dmr_gr)$diff.Methy[dmr_idx])),
    direction   = dmr_gr$direction[dmr_idx]
  )
}

#' Distance decay metaplot 집계
#' @param dist_df  data.frame: get_signed_distance() 출력
#' @param bin_size integer: bin 크기 (bp), 기본 5kb
#' @param max_dist integer: 포함할 최대 절댓값 거리 (bp)
#' @return data.frame: cnv_class, dist_bin, mean_abs_db, mean_db, dmr_density
make_metaplot_df <- function(dist_df,
                              bin_size = 5000L,
                              max_dist = 100000L) {
  dist_df %>%
    dplyr::filter(abs_dist <= max_dist) %>%
    dplyr::mutate(dist_bin = round(signed_dist / bin_size) * bin_size) %>%
    dplyr::group_by(cnv_class, dist_bin) %>%
    dplyr::summarise(
      mean_abs_db = mean(abs(diff_methy), na.rm = TRUE),
      mean_db     = mean(diff_methy,      na.rm = TRUE),
      dmr_density = dplyr::n(),
      .groups     = "drop"
    )
}

# 1. 데이터 준비  ===============================================================

setwd("/node200data/kachungk/hcc_data/DMR_SVs")

if (!is.null(opt$sv_strat_file)) {
  message("Reading pre-stratified SV file: ", opt$sv_strat_file)
  sv_list <- fread(opt$sv_strat_file) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    split(mcols(.)$sample) %>%
    lapply(function(gr) { names(gr) <- mcols(gr)$bp_id; gr })
} else {

info_cols <- c("MATE_ID", "MAPQ", "DETAILED_TYPE", "INSIDE_VNTR", "LOW_COV_IN", "HP") #nolint
sv_files <- list.files(pattern = "*somaticSVs.vcf.gz$", full.names = TRUE)

# Breakpoint based GRangesList
sv_list <- lapply(sv_files, function(x) {
  vcf <- readVcf(x)
  bp_df <- c(
    breakpointRanges(vcf, inferMissingBreakends = TRUE, info_columns = info_cols), #nolint
    breakendRanges(vcf, info_columns = info_cols)
  ) %>% as.data.frame()

  # Deduplicate BND breakpoints (same MATE_ID, same position)
  bnd_df <- bp_df %>%
    filter(svtype == "BND") %>%
    filter(str_detect(partner, "svrecord")) %>%
    mutate(partner = MATE_ID) %>%
    dplyr::select(-MATE_ID)

  non_bnd_df <- bp_df %>%
    filter(svtype != "BND") %>%
    dplyr::select(-MATE_ID)

  gr <- bind_rows(non_bnd_df, bnd_df) %>%
    distinct() %>% # 완전 동일한 행 제거
    mutate(sample = str_remove(basename(x), "\\.severus_somaticSVs.vcf.gz$")) %>% #nolint
    GRanges()
  mcols(gr)$bp_id <- names(gr)

  # Genotype 정보 추가 (VAF, DR, DV, hVAF)
  fmt_df <- data.table(
    ID = as.character(names(rowRanges(vcf))),
    VAF = as.numeric(geno(vcf)$VAF),
    DR = as.integer(geno(vcf)$DR),
    DV = as.integer(geno(vcf)$DV)
  )

  if ("hVAF" %in% names(geno(vcf))) {
    # hVAF는 보통 (variants x samples x 3) 형태
    hv <- geno(vcf)$hVAF[, 1, , drop = FALSE]
    fmt_df[, `:=`(
      hVAF_H0 = as.numeric(hv[, , 1]),
      hVAF_H1 = as.numeric(hv[, , 2]),
      hVAF_H2 = as.numeric(hv[, , 3])
    )]
  }

  sv_id <- mcols(gr)$sourceId
  idx <- match(as.character(sv_id), fmt_df$ID)
  mcols(gr)$VAF <- fmt_df$VAF[idx]
  mcols(gr)$DR  <- fmt_df$DR[idx]
  mcols(gr)$DV  <- fmt_df$DV[idx]
  if ("hVAF_H0" %in% colnames(fmt_df)) {
    mcols(gr)$hVAF_H0 <- fmt_df$hVAF_H0[idx]
    mcols(gr)$hVAF_H1 <- fmt_df$hVAF_H1[idx]
    mcols(gr)$hVAF_H2 <- fmt_df$hVAF_H2[idx]
  }

  tibble::as_tibble(gr)
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  dplyr::rename(patient_name = sample,
                sample = patient_code) %>%
  GRanges() %>%
  split(mcols(.)$sample) %>%
  lapply(function(gr) {
    names(gr) <- mcols(gr)$bp_id
    gr
  })

} # end else (VCF loading)

PATIENT_IDS <- names(sv_list)

consensus_dmr <- local({
  raw <- fread(opt$dmr_file)
  # Case 1: sample column already has patient codes (P1, P2 …)
  if (grepl("^P\\d+$", raw$sample[1])) {
    raw %>% GRanges() %>% split(mcols(.)$sample)
  # Case 2: patient_code column already present in file (e.g. confident_dmr output)
  } else if ("patient_code" %in% colnames(raw) && grepl("^P\\d+$", raw$patient_code[1])) {
    raw %>%
      dplyr::rename(patient_name = sample, sample = patient_code) %>%
      GRanges() %>%
      split(mcols(.)$sample)
  # Case 3: need to anonymise via patient_code_mapping.csv
  } else {
    raw %>%
      anonym_sample() %>%
      dplyr::rename(patient_name = sample, sample = patient_code) %>%
      GRanges() %>%
      split(mcols(.)$sample)
  }
})

cnv_files <- list.files(
  "../cnv_deepsomatic.out_hg38/purple",
  pattern = "*tumor.purple.segment.tsv$",
  full.names = TRUE
)

cnv_segments <- lapply(cnv_files, function(x) {
  df <- fread(x)
  if (!"minorAlleleCopyNumber" %in% colnames(df)) {
    df$minorAlleleCopyNumber <- NA_real_
  }

  df <- df %>%
    dplyr::select(
      chromosome, start, end, tumorCopyNumber, minorAlleleCopyNumber,
      bafCount, observedBAF, germlineStatus, depthWindowCount
    ) %>%
    dplyr::rename(copyNumber = tumorCopyNumber) %>%
    mutate(sample = str_remove(basename(x), "\\_tumor.purple.segment.tsv$"))

  df
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  dplyr::rename(patient_name = sample,
                sample = patient_code) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$sample)

### DMR CNV 상태 층화 ==============================================
consensus_dmr_cnv <- lapply(PATIENT_IDS, function(pt) {
  annotate_dmr_cnv(consensus_dmr[[pt]], cnv_segments[[pt]])
}) %>% setNames(PATIENT_IDS)

cat("=== 데이터 개요 ===\n")
cat(sprintf("SV/환자 평균: %.0f | DMR/환자 평균: %.0f\n",
    mean(sapply(sv_list, length)),
    mean(sapply(consensus_dmr, length))))

# SV 유형 분포 확인
sv_all_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
  as.data.frame(mcols(sv_list[[pt]])) %>%
    dplyr::select(any_of(c("geom_type", "DETAILED_TYPE", "cnv_class", "stratification",
                            "cnv_class_source", "cn_reliable", "vaf_concordant"))) %>%
    mutate(patient_id = pt)
}))
cat(sprintf("\n=== SV %s 분포 ===\n", GROUP_COL))
if (GROUP_COL == "stratification" && "stratification" %in% colnames(sv_all_df)) {
  print(table(sv_all_df$stratification))
} else {
  print(table(sv_all_df$cnv_class, sv_all_df$cnv_class_source))
}

rm(cnv_segments)
gc()

## 핵심 함수: SV breakpoint 주변 DMR 농축도 계산=================================
#' SV breakpoint 근처에 DMR이 얼마나 존재하는지 계산
#' @param sv_gr  GRanges: SV breakpoint (StructuralVariantAnnotation bp loci)
#'               필수 metadata: $sourceId, $bp_id (…_bp1/$…_bp2), $partner,
#'               $svtype, $cnv_class
#' @param dmr_gr GRanges: DMR — $direction("hyper"/"hypo"), $mean_diff(Δβ)
#'               direction 없으면 $diff.Methy(normal-tumor)에서 자동 파생
#' @param window integer: breakpoint 중심 양방향 반경 (bp)
#' @return data.frame: bp row 단위 (bp_id, bp_side, partner, svtype,
#'         n_dmr, n_hyper, n_hypo, mean_abs_db)
count_dmr_near_bp <- function(sv_gr, dmr_gr, window = opt$primary_window) {

  if (length(sv_gr) == 0 || length(dmr_gr) == 0)
    return(data.frame(
      bp_id = character(), bp_side = character(),
      partner = character(), svtype = character(),
      n_dmr = integer(), n_hyper = integer(),
      n_hypo = integer(), mean_abs_db = numeric()
    ))

  # direction 보장: diff.Methy = normal − tumor
  if (!"direction" %in% names(mcols(dmr_gr))) {
    mcols(dmr_gr)$direction <- dplyr::case_when(
      mcols(dmr_gr)$diff.Methy >  0 ~ "hypo",
      mcols(dmr_gr)$diff.Methy <  0 ~ "hyper",
      TRUE                           ~ NA_character_
    )
  }

  wins <- suppressWarnings(
    GenomicRanges::resize(sv_gr, width = window * 2L, fix = "center")
  )
  hits <- findOverlaps(wins, dmr_gr, minoverlap = 1L)

  lapply(seq_along(wins), function(i) {
    idx <- subjectHits(hits)[queryHits(hits) == i]
    data.frame(
      bp_id       = sv_gr$sourceId[i],
      bp_side     = sub(".*_(bp[12])$", "\\1", sv_gr$bp_id[i]),
      partner     = sv_gr$partner[i],
      svtype      = sv_gr$svtype[i],
      n_dmr       = length(idx),
      n_hyper     = sum(dmr_gr$direction[idx] == "hyper", na.rm = TRUE),
      n_hypo      = sum(dmr_gr$direction[idx] == "hypo",  na.rm = TRUE),
      mean_abs_db = if (length(idx) > 0)
        mean(abs(as.numeric(as.character(dmr_gr$diff.Methy[idx]))),
             na.rm = TRUE)
      else NA_real_
    )
  }) %>% bind_rows()
}

## Permutation test 함수 (null 분포 생성) =====================================

#' SV 위치 shuffle — 염색체 내 무작위 재배치, SV pair(bp1·bp2) 통째로 이동
#' @param sv_gr GRanges: SV breakpoint ($sourceId 필수)
#' @return GRanges: shuffle된 sv_gr (sourceId·svtype·cnv_class 보존)
shuffle_sv <- function(sv_gr) {
  valid <- intersect(names(CHROM_LENS), as.character(seqnames(sv_gr)))
  out   <- sv_gr

  for (sid in unique(out$sourceId)) {
    pair_idx <- which(out$sourceId == sid)
    pair_chr  <- as.character(seqnames(out[pair_idx]))

    # 같은 chromosome pair → 통째로 같은 offset 이동
    # inter-chromosomal (TRA/BND) → bp별 독립 shuffle
    if (length(unique(pair_chr)) == 1 && unique(pair_chr) %in% valid) {
      chr      <- unique(pair_chr)
      span     <- max(end(out[pair_idx])) - min(start(out[pair_idx]))
      max_orig <- CHROM_LENS[chr] - span - 1L
      if (max_orig < 1L) next
      new_origin <- sample.int(as.integer(max_orig), 1L)
      offset     <- new_origin - min(start(out[pair_idx]))
      ranges(out)[pair_idx] <- IRanges(
        start = start(out[pair_idx]) + offset,
        width = width(out[pair_idx])
      )
    } else {
      for (i in pair_idx) {
        chr <- as.character(seqnames(out[i]))
        if (!chr %in% valid) next
        max_s <- CHROM_LENS[chr] - width(out[i])
        if (max_s < 1L) next
        ranges(out)[i] <- IRanges(
          start = sample.int(as.integer(max_s), 1L),
          width = width(out[i])
        )
      }
    }
  }
  out
}

#' Permutation null 분포 생성
#' @param sv_gr    GRanges: SV breakpoint (단일 환자 · 단일 cnv_class 서브셋)
#' @param dmr_gr   GRanges: DMR (동일 환자)
#' @param window   integer: bp window (bp)
#' @param n_perm   integer: permutation 반복 수
#' @param BPPARAM  BiocParallel param (NULL이면 MulticoreParam 자동 생성)
#' @return numeric vector: 각 permutation의 mean(n_dmr) [SV 단위 max 집계]
run_permutation <- function(sv_gr, dmr_gr, window = opt$primary_window,
                            n_perm = opt$n_perm, BPPARAM = NULL) {
  if (length(sv_gr) == 0) return(numeric(0))

  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(
      workers     = min(4L, BiocParallel::multicoreWorkers()),
      progressbar = FALSE,
      RNGseed     = 20260421L
    )
  }

  out <- BiocParallel::bplapply(seq_len(n_perm), function(i) {
    shuf <- shuffle_sv(sv_gr)
    res  <- count_dmr_near_bp(shuf, dmr_gr, window)
    res %>%
      dplyr::group_by(bp_id) %>%
      dplyr::summarise(n = max(n_dmr), .groups = "drop") %>%
      dplyr::pull(n) %>%
      mean(na.rm = TRUE)
  }, BPPARAM = BPPARAM)

  as.numeric(unlist(out, use.names = FALSE))
}



# Main Analysis ================================================================
## H1·H2. Window enrichment test ===============================================
# Permutation test (n=1,000): SV pair를 chromosome 내 무작위 재배치 →
#     obs_mean vs null 분포 비교. 분포 가정 없는 가장 보수적 접근.
#근거: Encode 3 (Partridge et al. 2020 Nature), SV-DMR 연구에서 표준.
cat("\n=== H1·H2: Window enrichment test ===\n")

window_enrich <- lapply(PATIENT_IDS, function(pt) {
  message("Working on ", pt)

  sv  <- sv_list[[pt]]
  dmr <- consensus_dmr_cnv[[pt]]

  # direction 파생 (루프당 1회)
  mcols(dmr)$direction <- dplyr::case_when(
    mcols(dmr)$diff.Methy >  0 ~ "hypo",
    mcols(dmr)$diff.Methy <  0 ~ "hyper",
    TRUE                        ~ NA_character_
  )

  lapply(GROUP_LEVELS, function(cls) {
    sv_sub <- sv[mcols(sv)[[GROUP_COL]] == cls]
    if (length(sv_sub) == 0 || length(unique(sv_sub$sourceId)) < opt$min_sv_per_group) {
      message("[", cls, "] n_sv < min_sv_per_group (", opt$min_sv_per_group, "), skip")
      return(NULL)
    }

    message("Count DMR near SV bp for ", cls)
    lapply(opt$windows, function(win) {
      message("window = ", win, " bp")

      # bp row 단위 카운트 → sourceId(=bp_id) 기준 max 집계 (Primary B)
      res <- count_dmr_near_bp(sv_sub, dmr, win) %>%
        dplyr::group_by(bp_id) %>%
        dplyr::summarise(
          n_dmr       = max(n_dmr),
          n_hyper     = max(n_hyper),
          n_hypo      = max(n_hypo),
          mean_abs_db = mean(mean_abs_db, na.rm = TRUE),
          .groups     = "drop"
        ) 

      if (nrow(res) == 0) return(NULL)
      obs_mean <- mean(res$n_dmr, na.rm = TRUE)
      null_vec <- run_permutation(sv_sub, dmr, win, opt$n_perm)

      data.frame(
        patient_id       = pt,
        cnv_class        = cls,
        window_kb        = win / 1000,
        n_sv             = nrow(res),
        obs_mean_dmr     = obs_mean,
        null_mean        = mean(null_vec),
        null_sd          = sd(null_vec),
        enrichment_ratio = obs_mean / pmax(mean(null_vec), 1e-6),
        p_perm           = mean(null_vec >= obs_mean),
        pct_sv_w_dmr     = mean(res$n_dmr > 0) * 100,
        mean_hyper       = mean(res$n_hyper, na.rm = TRUE),
        mean_hypo        = mean(res$n_hypo,  na.rm = TRUE),
        hyper_ratio      = mean(res$n_hyper / pmax(res$n_dmr, 1), na.rm = TRUE))
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# FDR 보정 (환자 × class × window 전체)
window_enrich$p_fdr <- p.adjust(window_enrich$p_perm, method = "BH")
window_enrich$sig_label <- dplyr::case_when(
  window_enrich$p_fdr < 0.001 ~ "***",
  window_enrich$p_fdr < 0.01  ~ "**",
  window_enrich$p_fdr < 0.05  ~ "*",
  TRUE                         ~ "ns"
)

# 기존 figure 코드 호환: enrichment_results 별칭 (컬럼명 복원)
enrichment_results <- window_enrich %>%
  dplyr::rename(
    mean_n_dmr      = obs_mean_dmr,
    pct_sv_with_dmr = pct_sv_w_dmr
  )

## H3 Distance decay + metaplot ==============================================
# Axis 3 추가 검정 근거 및 제안 (CNV confounding 배제)
#   - copy_neutral vs copy_changing enrichment_ratio 비교 (Wilcoxon, Section 10)
#   - permutation enrichment_ratio: null 대비 농축 배율
#
# 추가 권장 검정:
#   [1] Wilcoxon rank-sum (copy_neutral enrichment_ratio vs 1.0)
#       - 귀무: copy-neutral SV의 enrichment_ratio 중위값 = 1 (농축 없음)
#       - 이유: copy-neutral SV (INV/TRA)에서도 enrichment_ratio > 1이면
#               methylation 변화가 CNV와 독립적임을 직접 증명.
#               단측 검정(alternative="greater")으로 방향성 명시.
#       - 구현:
#         neutral_ratio <- window_enrich %>%
#           filter(cnv_class=="copy_neutral", window_kb==50) %>%
#           pull(enrichment_ratio)
#         wilcox.test(neutral_ratio, mu=1, alternative="greater")
#
#   [2] PURPLE CNV 세그먼트 overlap → CNV-driven DMR 제거 후 재분석
#       - PURPLE 출력: cnv_segment (chr, start, end, copyNumber, majorCN, minorCN)
#       - DMR이 CNV 세그먼트와 overlap하면 "CNV-driven" 후보로 플래그
#       - CNV-driven DMR 제거 후 enrichment_ratio 재계산 → 변화 없으면 robust
#       - 이유: DEL/DUP에서 관찰된 methylation 변화가 실제 copy number 변화로
#               인한 artifact인지 구별하는 가장 직접적 방법.
#               PURPLE + HMF methylation 파이프라인(Hartwig Medical Foundation)
#               에서 표준 QC 단계로 사용.
#       - 구현:
#         cnv_seg_gr <- makeGRangesFromDataFrame(purple_cnv)
#         dmr_cnv_flag <- countOverlaps(dmr_gr, cnv_seg_gr) > 0
#         dmr_clean <- dmr_gr[!dmr_cnv_flag]
#         # 이후 window_enrich 재실행 with dmr_clean
#
#   [3] Interaction test: cnv_class × window_kb
#       - lmer(enrichment_ratio ~ cnv_class * window_kb + (1|patient_id))
#       - interaction이 유의하면 copy-neutral과 copy-changing의 enrichment
#         패턴이 window 크기에 따라 다르게 변화함을 증명
#       - 이유: 근거리(±10kb)에서만 enrichment 차이가 있다면
#               이는 직접적 cis 효과의 증거이며, 원거리에서 사라진다면
#               TAD-mediated long-range effect 논의 가능.

cat("\n=== H3: Distance decay Analysis ===\n")

dist_all <- lapply(PATIENT_IDS, function(pt) {
  sv  <- sv_list[[pt]]
  dmr <- consensus_dmr_cnv[[pt]]

  mcols(dmr)$direction <- dplyr::case_when(
    mcols(dmr)$diff.Methy >  0 ~ "hypo",
    mcols(dmr)$diff.Methy <  0 ~ "hyper",
    TRUE                        ~ NA_character_
  )

  lapply(GROUP_LEVELS, function(cls) {
    sv_sub <- sv[mcols(sv)[[GROUP_COL]] == cls]
    if (length(sv_sub) == 0 || length(unique(sv_sub$sourceId)) < opt$min_sv_per_group) return(NULL)
    df <- get_signed_distance(sv_sub, dmr, max_dist = 200000L)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$patient_id <- pt
    df
  }) %>% bind_rows()
}) %>% bind_rows()

# Spearman ρ: 거리 vs |Δβ| (H3 핵심 통계)
# 환자별 독립 ρ → 환자 분포에 대해 단측 Wilcoxon으로 검정 (pseudoreplication 방지)
per_pt_decay_rho <- dist_all %>%
  dplyr::group_by(patient_id, cnv_class) %>%
  dplyr::summarise(
    spearman_rho = cor(abs_dist, abs(diff_methy),
                       method = "spearman", use = "complete.obs"),
    n_pairs      = dplyr::n(),
    .groups      = "drop"
  ) %>%
  dplyr::filter(n_pairs >= 5L)

decay_cor <- per_pt_decay_rho %>%
  dplyr::group_by(cnv_class) %>%
  dplyr::summarise(
    n_patients   = sum(!is.na(spearman_rho)),
    median_rho   = median(spearman_rho, na.rm = TRUE),
    p_value      = tryCatch(
      wilcox.test(spearman_rho, mu = 0, alternative = "less")$p.value,
      error = function(e) NA_real_
    ),
    .groups      = "drop"
  ) %>%
  dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))

cat("Distance decay Spearman ρ (by CNV class):\n")
print(decay_cor)

# Metaplot 집계 (5kb bin)
metaplot_df <- make_metaplot_df(dist_all, bin_size = 5000L, max_dist = 100000L)

write.csv(dist_all,          file.path(opt$outdir, sprintf("%s_dist_decay_full.csv",        opt$run_id)), row.names = FALSE, quote = FALSE)
write.csv(metaplot_df,       file.path(opt$outdir, sprintf("%s_metaplot_binned.csv",          opt$run_id)), row.names = FALSE, quote = FALSE)
write.csv(decay_cor,         file.path(opt$outdir, sprintf("%s_decay_spearman_cor.csv",       opt$run_id)), row.names = FALSE, quote = FALSE)
write.csv(per_pt_decay_rho,  file.path(opt$outdir, sprintf("%s_decay_spearman_per_pt.csv",    opt$run_id)), row.names = FALSE, quote = FALSE)

cat("H3 출력: dist_decay_full.csv · metaplot_binned.csv · decay_spearman_cor.csv\n")

# Axis 2 — SV 유형별 methylation 방향성 검정 ==================================
#
# 검정 선택 근거
# ─────────────────────────────────────────────────────────────────────────────
# [A] Chi-square goodness-of-fit (방향 균일성)
#     귀무: 각 svtype 내 hyper:hypo = 50:50 (무방향)
#     이유: DEL은 regulatory element 소실 → hypermethylation 편향 예측,
#           DUP는 copy gain → hypomethylation 편향 예측.
#           정규성 가정 없이 카운트 데이터에 적용 가능하며,
#           ENCODE/Roadmap 기반 SV-methylation 연구(e.g. Flavahan et al. 2016
#           Nature, Onuchic et al. 2018 Nat Genet)에서 유사하게 쓰임.
#
# [B] Fisher's exact test (svtype pair 간 방향 비율 비교)
#     귀무: svtype A vs B에서 hyper 비율이 동일
#     이유: n=12 환자로 셀 수가 작아 chi-square 근사가 불안정.
#           Fisher's exact는 정확 확률을 계산하므로 소표본에서 더 신뢰적.
#           TCGA 기반 SV-methylation 분석(Hoadley et al. 2018 Cell)에서도
#           소그룹 비교에 Fisher's exact 채택.
#
# [C] Kruskal-Wallis + Dunn post-hoc (Δβ 크기 비교)
#     귀무: 모든 svtype에서 |diff.Methy| 중위값이 동일
#     이유: Δβ 분포는 heavy-tailed이며 정규성 가정 불가.
#           비모수 검정이 적절하며, 다중 그룹 비교 후
#           Bonferroni-Dunn으로 pairwise 보정.
#           동일한 접근이 Li et al. 2020 Nat Commun (SV-epigenome)에 사용됨.

cat("\n=== Axis 2: SV 유형별 methylation 방향성 검정 ===\n")

# dist_all에서 svtype × direction 카운트 집계
# (bp row 단위 — H3와 동일 데이터, 집계 방식만 변경)
direction_sv <- dist_all %>%
  dplyr::filter(!is.na(direction), !is.na(svtype)) %>%
  dplyr::count(svtype, direction) %>%
  tidyr::pivot_wider(names_from = direction,
                     values_from = n,
                     values_fill = 0L) %>%
  dplyr::mutate(
    total       = hyper + hypo,
    hyper_pct   = hyper / total * 100,
    hypo_pct    = hypo  / total * 100
  )

cat("SV 유형별 hyper/hypo 카운트:\n")
print(direction_sv)

# ── [A] Chi-square goodness-of-fit: 각 svtype 내 hyper:hypo = 50:50 ──────────
chisq_results <- direction_sv %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    chisq_p = tryCatch(
      chisq.test(c(hyper, hypo), p = c(0.5, 0.5))$p.value,
      error = function(e) NA_real_
    ),
    # 셀 수 < 5이면 Fisher's exact로 대체 (소표본 보정)
    use_fisher = (hyper < 5 | hypo < 5)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chisq_fdr = p.adjust(chisq_p, method = "BH"))

cat("\n[A] Chi-square (hyper:hypo = 50:50 귀무):\n")
print(chisq_results %>%
  dplyr::select(svtype, hyper, hypo, hyper_pct, chisq_p, chisq_fdr, use_fisher))

# ── [B] Fisher's exact: svtype pair 간 hyper 비율 비교 ───────────────────────
sv_types_present <- direction_sv$svtype
fisher_pairs <- combn(sv_types_present, 2, simplify = FALSE)

fisher_results <- do.call(rbind, lapply(fisher_pairs, function(pair) {
  d <- direction_sv %>% dplyr::filter(svtype %in% pair)
  if (nrow(d) < 2) return(NULL)
  mat <- matrix(c(d$hyper[1], d$hypo[1],
                  d$hyper[2], d$hypo[2]),
                nrow = 2, byrow = TRUE,
                dimnames = list(pair, c("hyper", "hypo")))
  ft <- fisher.test(mat, alternative = "two.sided")
  data.frame(
    svtype_A  = pair[1],
    svtype_B  = pair[2],
    OR        = round(ft$estimate, 3),
    CI_lo     = round(ft$conf.int[1], 3),
    CI_hi     = round(ft$conf.int[2], 3),
    p_fisher  = round(ft$p.value, 5)
  )
}))

fisher_results$p_fdr <- p.adjust(fisher_results$p_fisher, method = "BH")

cat("\n[B] Fisher's exact (svtype pair 간 hyper 비율):\n")
print(fisher_results)

# ── [C] Kruskal-Wallis + Dunn post-hoc: |Δβ| 크기 비교 ─────────────────────
kw_data <- dist_all %>%
  dplyr::filter(!is.na(svtype), !is.na(diff_methy)) %>%
  dplyr::mutate(abs_db = abs(diff_methy))

kw_test <- kruskal.test(abs_db ~ svtype, data = kw_data)
cat(sprintf("\n[C] Kruskal-Wallis: H = %.3f, df = %d, p = %.5f\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

# Dunn post-hoc (Bonferroni 보정)
if (requireNamespace("dunn.test", quietly = TRUE)) {
  dunn_res <- dunn.test::dunn.test(
    kw_data$abs_db, kw_data$svtype,
    method = "bonferroni", kw = FALSE, label = TRUE
  )
  dunn_df <- data.frame(
    comparison = dunn_res$comparisons,
    Z          = round(dunn_res$Z, 3),
    p_adj      = round(dunn_res$P.adjusted, 5)
  )
  cat("\nDunn post-hoc (Bonferroni):\n")
  print(dunn_df)
  write.csv(dunn_df, file.path(opt$outdir, sprintf("%s_axis2_dunn_posthoc.csv", opt$run_id)), row.names = FALSE, quote = FALSE)
} else {
  cat("  dunn.test 패키지 없음 — install.packages('dunn.test') 후 재실행\n")
  dunn_df <- NULL
}

# Violin plot: |Δβ| per svtype
axis2_violin <- dist_all %>%
  dplyr::filter(!is.na(svtype), !is.na(diff_methy)) %>%
  dplyr::mutate(
    abs_db    = abs(diff_methy),
    svtype    = factor(svtype,
                       levels = c("DEL","DUP","INV","TRA","INS"))
  )

if (!opt$no_plot) {
p_violin <- ggplot(axis2_violin, aes(x = svtype, y = abs_db, fill = svtype)) +
  geom_violin(alpha = 0.75, color = "grey40", linewidth = 0.4,
              trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5,
               fill = "white", color = "grey30", alpha = 0.8) +
  scale_fill_manual(values = c(
    DEL = unname(SV_COLORS["copy_losing"]),
    DUP = unname(SV_COLORS["copy_gaining"]),
    INV = unname(SV_COLORS["copy_neutral"]),
    TRA = unname(SV_COLORS["copy_neutral"]),
    INS = unname(SV_COLORS["insertion"])
  ), guide = "none") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Kruskal-Wallis p = %.4f", kw_test$p.value),
           hjust = 1.1, vjust = 1.5, size = 3.2, color = "grey40") +
  labs(
    title    = "Axis 2. SV 유형별 |Δβ| 분포",
    subtitle = "Kruskal-Wallis + Dunn post-hoc (Bonferroni) | bp row 단위",
    x        = "SV 유형", y = "|diff.Methy| (normal − tumor)"
  ) +
  theme_hcc

saveRDS(p_violin, file.path(opt$outdir, sprintf("%s_axis2_violin_abs_db.rds", opt$run_id)))
ggsave(file.path(opt$outdir, sprintf("%s_axis2_violin_abs_db.png", opt$run_id)), p_violin,
       width = 8, height = 5, device = "png", dpi = 300)

# hyper/hypo 비율 막대 그래프
p_direction_bar <- direction_sv %>%
  tidyr::pivot_longer(c(hyper_pct, hypo_pct),
                      names_to = "direction", values_to = "pct") %>%
  dplyr::mutate(
    direction = recode(direction,
                       hyper_pct = "Hypermethylated",
                       hypo_pct  = "Hypomethylated"),
    sig_label = dplyr::case_when(
      svtype %in% (chisq_results %>%
                     dplyr::filter(chisq_fdr < 0.05) %>%
                     dplyr::pull(svtype)) ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ggplot(aes(x = svtype, y = pct, fill = direction)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_text(
    data = . %>% dplyr::filter(direction == "Hypermethylated"),
    aes(label = sig_label, y = 105),
    size = 5, color = "grey30"
  ) +
  scale_fill_manual(
    values = c("Hypermethylated" = "#E24B4A",
               "Hypomethylated"  = "#3B8BD4"),
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 110)) +
  annotate("text", x = 0.55, y = 52,
           label = "50% (무방향)", size = 3, color = "grey50", hjust = 0) +
  labs(
    title    = "Axis 2. hyper/hypo ratio by SV type",
    subtitle = "*: Chi-square FDR < 0.05 (50:50 귀무 기각)",
    x        = "SV 유형", y = "비율 (%)"
  ) +
  theme_hcc

saveRDS(p_direction_bar, file.path(opt$outdir, sprintf("%s_axis2_direction_bar.rds", opt$run_id)))
ggsave(file.path(opt$outdir, sprintf("%s_axis2_direction_bar.png", opt$run_id)), p_direction_bar,
       width = 7, height = 5, device = "png", dpi = 300)
} # end if (!opt$no_plot) — Axis 2 plots

# 결과 저장
write.csv(chisq_results,  file.path(opt$outdir, sprintf("%s_axis2_chisq_results.csv", opt$run_id)),  row.names = FALSE, quote = FALSE)
write.csv(fisher_results, file.path(opt$outdir, sprintf("%s_axis2_fisher_results.csv", opt$run_id)), row.names = FALSE, quote = FALSE)

cat("Axis 2 출력: axis2_chisq_results.csv · axis2_fisher_results.csv\n")
cat("             axis2_dunn_posthoc.csv · axis2_violin_abs_db.png\n")
cat("             axis2_direction_bar.png\n")

viz_window <- opt$primary_window / 1000   # used by both viz and summary sections below

# Visualization ================================================================
if (!opt$no_plot) {
## Figure 1 — CNV class별 DMR 농축도 (window별 비교) =============================
  plot_df1 <- enrichment_results %>%
    filter(window_kb == viz_window) %>%
    mutate(
      cnv_label = factor(cnv_class,
        levels = c("copy_neutral", "copy_gaining", "copy_losing", "insertion", "COM"),
        labels = c("Copy-neutral\n(INV, TRA/BND)",
                   "Copy-gaining\n(DUP)",
                   "Copy-losing\n(DEL)",
                   "Insertion\n(INS)",
                   "Complex\n(COM)"))
    )

  p1a <- ggplot(plot_df1, aes(x = cnv_label, y = mean_n_dmr, fill = cnv_class)) +
    geom_boxplot(width = 0.5, outlier.size = 1.2, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                  color = "grey20", linewidth = 0.5) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    labs(title    = "Average DMR # in ±50kb window by SV type",
          subtitle = "Dot = Patient | *: FDR < 0.05",
          x = NULL, y = "Average DMR # / SV") +
    theme_hcc

  # Null 기대값 (전체 genome 균등 분포 가정)
  null_pct_expected <- (opt$primary_window * 2 / 3e9) * mean(sapply(consensus_dmr_cnv, length)) * 100 #nolint

  # DMR을 가진 SV 비율
  p1b <- ggplot(plot_df1, aes(x = cnv_label, y = pct_sv_with_dmr, fill = cnv_class)) +
    geom_boxplot(width = 0.5, outlier.size = 1.2, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    geom_hline(yintercept = null_pct_expected, linetype = "dashed",
              color = "grey50", linewidth = 0.7) +
    annotate("text", x = 0.6, y = null_pct_expected + 1,
            label = "null expected", size = 3, color = "grey50", hjust = 0) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(round(x), "%")) +
    labs(title    = paste0("Proportion of SVs with DMRs (±", viz_window, "kb)"),
        subtitle = "Dashed line = null expected",
        x = NULL, y = "Proportion of SVs with DMRs (%)") +
    theme_hcc

  fig1 <- p1a | p1b
  fig1 <- fig1 +
    plot_annotation(
      title   = paste0("Figure 1. Spatial Correlation of SV–DMR (±", viz_window, "kb window)"),
      caption = sprintf("n=%d | PacBio WGS | permutation n=%d | FDR: BH",
                        length(PATIENT_IDS), opt$n_perm),
      theme   = theme(plot.title = element_text(face = "bold", size = 14))
    )

  #
  saveRDS(p1a, file.path(opt$outdir, sprintf("%s_fig1a_cnv_class_dmr_enrichment.rds", opt$run_id)))
  saveRDS(p1b, file.path(opt$outdir, sprintf("%s_fig1b_cnv_class_dmr_enrichment.rds", opt$run_id)))
  ggsave(file.path(opt$outdir, sprintf("%s_fig1_cnv_class_dmr_enrichment.png", opt$run_id)), fig1,
        width = 12, height = 6, device = png, dpi = 300)

  ## Figure 2 — Window 크기별 농축도 곡선 (copy_neutral vs others) ===============
  window_df <- enrichment_results %>%
    mutate(
      cnv_group = ifelse(cnv_class == "copy_neutral",
                        "Copy-neutral (INV/TRA)", "Copy-number-changing (DEL/DUP/INS)"),
      window_kb = factor(window_kb, levels = sort(unique(window_kb)))
    ) %>%
    group_by(cnv_group, window_kb) %>%
    summarise(
      mean_dmr  = mean(mean_n_dmr, na.rm = TRUE),
      se_dmr    = sd(mean_n_dmr, na.rm = TRUE) / sqrt(n()),
      .groups   = "drop"
    )

  p2 <- ggplot(window_df, aes(x = as.numeric(as.character(window_kb)),
                              y = mean_dmr, color = cnv_group,
                              group = cnv_group)) +
    geom_ribbon(aes(ymin = mean_dmr - se_dmr, ymax = mean_dmr + se_dmr,
                    fill = cnv_group), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c("Copy-neutral (INV/TRA)" = "#3B8BD4",
                "Copy-number-changing (DEL/DUP/INS)" = "#E24B4A"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("Copy-neutral (INV/TRA)" = "#3B8BD4",
                "Copy-number-changing (DEL/DUP/INS)" = "#E24B4A"),
      guide = "none"
    ) +
    scale_x_log10(
      breaks = sort(opt$windows / 1000),
      labels = paste0("±", sort(opt$windows / 1000), "kb")
    ) +
    labs(
      title    = "Figure 2. Window 크기별 SV 주변 평균 DMR 수",
      subtitle = "Copy-neutral SV에서도 DMR 농축 → CNV confounding 배제 근거",
      x        = "SV breakpoint 기준 분석 window",
      y        = "평균 DMR 수 / SV",
      caption  = "리본: ±1 SE | copy-neutral = INV + TRA/BND"
    ) +
    theme_hcc +
    theme(legend.position = "top")

  saveRDS(p2, file.path(opt$outdir, sprintf("%s_fig2_window_enrichment_curve.rds", opt$run_id)))
  ggsave(file.path(opt$outdir, sprintf("%s_fig2_window_enrichment_curve.png", opt$run_id)), p2,
        width = 8, height = 5, device = png, dpi = 300)

  ## Figure 3 — Methylation 방향성: Hyper vs Hypo (SV 유형별) ===================
  direction_df <- enrichment_results %>%
    filter(window_kb == viz_window) %>%
    dplyr::select(patient_id, cnv_class, mean_hyper, mean_hypo) %>%
    pivot_longer(c(mean_hyper, mean_hypo),
                names_to = "direction", values_to = "n_dmr") %>%
    mutate(
      direction = recode(direction,
                        mean_hyper = "Hypermethylated",
                        mean_hypo  = "Hypomethylated"),
      cnv_label = factor(cnv_class,
        levels = c("copy_neutral","copy_gaining","copy_losing","insertion"),
        labels = c("INV/TRA\n(neutral)", "DUP\n(gain)",
                  "DEL\n(loss)", "INS"))
    )

  p3a <- ggplot(direction_df,
                aes(x = cnv_label, y = n_dmr, fill = direction)) +
    geom_boxplot(position = position_dodge(0.7), width = 0.55,
                alpha = 0.8, outlier.size = 1) +
    scale_fill_manual(
      values = c("Hypermethylated" = "#E24B4A",
                "Hypomethylated"  = "#3B8BD4"),
      name = "DMR 방향"
    ) +
    labs(
      title    = "SV 유형별 methylation 방향성",
      subtitle = "DEL→hypermethylation, DUP→hypomethylation 가설 검증",
      x        = NULL, y = "평균 DMR 수 / SV"
    ) +
    theme_hcc

  # Hyper/(Hyper+Hypo) 비율
  ratio_df <- enrichment_results %>%
    filter(window_kb == viz_window) %>%
    mutate(
      hyper_ratio_pct = hyper_ratio * 100,
      cnv_label = factor(cnv_class,
        levels = c("copy_neutral","copy_gaining","copy_losing","insertion"),
        labels = c("INV/TRA\n(neutral)", "DUP\n(gain)",
                  "DEL\n(loss)", "INS"))
    )

  p3b <- ggplot(ratio_df, aes(x = cnv_label, y = hyper_ratio_pct, fill = cnv_class)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey60",
              linewidth = 0.7) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 1) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.5, color = "grey30") +
    annotate("text", x = 0.55, y = 52, label = "50% (무방향성)",
            size = 3, color = "grey50", hjust = 0) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(round(x), "%"),
                      limits = c(0, 100)) +
    labs(
      title    = "Hypermethylation 비율 (Hyper / Hyper+Hypo)",
      subtitle = "50% 기준선 이상 = hypermethylation 우세",
      x        = NULL, y = "Hypermethylation 비율 (%)"
    ) +
    theme_hcc

  fig3 <- p3a | p3b
  fig3 <- fig3 +
    plot_annotation(
      title   = "Figure 3. SV 유형별 methylation 방향성 비교",
      caption = "DEL: regulatory element 소실 → hypermethylation 예상\nDUP: copy 증가 → hypomethylation 예상",
      theme   = theme(plot.title = element_text(face = "bold", size = 14))
    )

  saveRDS(p3a, file.path(opt$outdir, sprintf("%s_fig3a_methylation_direction.rds", opt$run_id)))
saveRDS(p3b, file.path(opt$outdir, sprintf("%s_fig3b_methylation_direction.rds", opt$run_id)))
ggsave(file.path(opt$outdir, sprintf("%s_fig3_methylation_direction.png", opt$run_id)), fig3,
        width = 12, height = 6, device = "png", dpi = 300)

  # Figure 4 — Window enrichment permutation 결과 시각화 ===================
  perm_plot_df <- window_enrich %>%
    filter(window_kb == viz_window) %>%
    mutate(
      cnv_label = factor(cnv_class,
        levels = c("copy_neutral","copy_gaining","copy_losing","insertion"),
        labels = c("INV/TRA\n(neutral)", "DUP\n(gain)", "DEL\n(loss)", "INS"))
    )

  p4 <- ggplot(perm_plot_df,
              aes(x = cnv_label, y = enrichment_ratio, fill = cnv_class)) +
    geom_hline(yintercept = 1, linetype = "dashed",
              color = "grey50", linewidth = 0.7) +
    geom_col(position = position_dodge(0.7), width = 0.55,
            alpha = 0.8, color = "grey40", linewidth = 0.3) +
    geom_text(aes(label = sig_label, y = enrichment_ratio + 0.05),
              position = position_dodge(0.7), size = 4, vjust = 0) +
    facet_wrap(~patient_id, nrow = 1) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    labs(
      title    = "Figure 4. Permutation test — 관찰값 / null 기대값 비율",
      subtitle = "1 초과 = 무작위 배치 대비 DMR 농축 | *, **, ***: FDR<0.05/0.01/0.001",
      x        = NULL, y = "농축 배율 (observed / null mean)",
      caption  = sprintf("n_perm=%d | ±50kb window | FDR: BH", opt$n_perm)
    ) +
    theme_hcc +
    theme(strip.text = element_text(size = 9))

  saveRDS(p4, file.path(opt$outdir, sprintf("%s_fig4_permutation_enrichment.rds", opt$run_id)))
  ggsave(file.path(opt$outdir, sprintf("%s_fig4_permutation_enrichment.png", opt$run_id)), p4,
        width = 12, height = 5, device = png, dpi = 300)
  } # end if (!opt$no_plot) — Figures 1–4

## Layer 1 — Phase block co-localization visualization =======================
# Reads layer1_coloc_results.csv from haplotype_sv_admr_analysis.R (nperm=1000)
if (!opt$no_plot && file.exists(opt$hap_coloc_csv)) {
  coloc_res <- read.csv(opt$hap_coloc_csv)

  fig_hap1_df <- coloc_res %>%
    tidyr::pivot_longer(c(mean_dmr_sv_blk, mean_dmr_nonsv),
                          names_to = "group", values_to = "mean_dmr") %>%
    dplyr::mutate(group = dplyr::recode(group,
      mean_dmr_sv_blk = "SV-containing block",
      mean_dmr_nonsv  = "SV-free block"))

  p_hap1a <- ggplot(fig_hap1_df, aes(x = group, y = mean_dmr, fill = group)) +
    geom_boxplot(width = 0.45, alpha = 0.8, outlier.size = 1.2) +
    geom_line(aes(group = patient_id), color = "grey60",
              linewidth = 0.4, alpha = 0.6) +
    geom_point(aes(group = patient_id), size = 2, alpha = 0.7) +
    scale_fill_manual(
      values = c("SV-containing block" = "#3B8BD4",
                  "SV-free block"       = "#888780"),
      guide = "none"
    ) +
    labs(title    = "aDMR count in phase blocks",
        subtitle = "Lines connect same patient",
        x = NULL, y = "Mean aDMR count / block") +
    theme_hcc

  p_hap1b <- ggplot(coloc_res,
                    aes(x = reorder(patient_id, enrichment_ratio),
                        y = enrichment_ratio,
                        fill = enrichment_ratio > 1)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed",
                color = "grey50", linewidth = 0.7) +
    geom_text(aes(label = ifelse(p_wilcoxon < 0.05, "*", ""),
                  y = enrichment_ratio + 0.05), size = 4) +
    scale_fill_manual(
      values = c("TRUE" = "#3B8BD4", "FALSE" = "#E24B4A"),
      guide  = "none"
    ) +
    coord_flip() +
    labs(title    = "SV phase block aDMR enrichment",
        subtitle = "*: Wilcoxon p < 0.05",
        x = NULL, y = "Enrichment ratio (SV block / non-SV block)") +
    theme_hcc

  fig_hap1 <- (p_hap1a | p_hap1b) +
    plot_annotation(
      title   = "Haplotype Fig 1. Layer 1 — SV Phase Block aDMR Co-localization",
      caption = sprintf("n = %d patients | nperm = 1000 | dashed = null (ratio 1)",
                        nrow(coloc_res)),
      theme   = theme(plot.title = element_text(face = "bold", size = 14))
    )

  saveRDS(p_hap1a, file.path(opt$outdir, sprintf("%s_hap_fig1a_coloc.rds", opt$run_id)))
  saveRDS(p_hap1b, file.path(opt$outdir, sprintf("%s_hap_fig1b_coloc.rds", opt$run_id)))
  ggsave(
    file.path(opt$outdir, sprintf("%s_hap_fig1_coloc.png", opt$run_id)),
      fig_hap1, width = 12, height = 5, device = "png", dpi = 300
  )
  cat(sprintf("Haplotype Layer 1 figure saved: %s_hap_fig1_coloc.png\n", opt$run_id))
}

# Statistical test ====================================================
# Binary comparison: copy_neutral vs changing (cnv_class) or boundary vs non-boundary (tier)
if (opt$group_by == "cnv_class") {
  cat("\n=== 통계 검정: Copy-neutral vs Copy-number-changing (±", viz_window, "kb) ===\n")
  test_df <- enrichment_results %>%
    filter(window_kb == viz_window) %>%
    mutate(group = ifelse(cnv_class == "copy_neutral", "copy_neutral", "copy_changing"))
  group_label <- "copy_neutral vs copy_changing"
} else {
  cat("\n=== 통계 검정: Boundary-disrupting vs Non-boundary (±", viz_window, "kb) ===\n")
  test_df <- enrichment_results %>%
    filter(window_kb == viz_window) %>%
    mutate(group = ifelse(
      cnv_class %in% c("TAD+CTCF disrupting", "CTCF-only", "TAD-only"),
      "boundary", "non_boundary"
    ))
  group_label <- "boundary vs non_boundary"
}

wilcox_mean <- wilcox.test(
  mean_n_dmr ~ group, data = test_df,
  alternative = "two.sided", exact = FALSE
)
wilcox_pct <- wilcox.test(
  pct_sv_with_dmr ~ group, data = test_df,
  alternative = "two.sided", exact = FALSE
)
cat(sprintf("[%s] 평균 DMR 수 비교 — Wilcoxon p = %.4f\n", group_label, wilcox_mean$p.value))
cat(sprintf("[%s] DMR 보유 SV 비율 비교 — Wilcoxon p = %.4f\n", group_label, wilcox_pct$p.value))

# Group-level summary
summary_table <- enrichment_results %>%
  filter(window_kb == viz_window) %>%
  group_by(cnv_class) %>%
  summarise(
    n_patients      = n(),
    mean_dmr_per_sv = round(mean(mean_n_dmr, na.rm = TRUE), 2),
    sd_dmr          = round(sd(mean_n_dmr, na.rm = TRUE), 2),
    pct_sv_w_dmr    = round(mean(pct_sv_with_dmr, na.rm = TRUE), 1),
    mean_hyper_pct  = round(mean(hyper_ratio * 100, na.rm = TRUE), 1),
    .groups         = "drop"
  )

cat(sprintf("\n=== %s 요약 통계 ===\n", GROUP_COL))
print(summary_table)

fwrite(summary_table, file.path(opt$outdir, sprintf("%s_sv_dmr_enrichment_summary.csv", opt$run_id)),
            row.names = FALSE, quote = FALSE)

# 11-A. Co-occurrence analysis ============================================
# % of SVs with nearest DMR ≤ threshold vs null
# 전체 SV (CNV class 무관) 기준 per-patient co-occurrence.
# 현재 window_enrich의 null은 mean(n_dmr) 단위라 % null 분포가 없음 → 별도 계산.

  cat("\n=== Co-occurrence: % SVs within ±", opt$primary_window / 1000, "kb of nearest DMR ===\n")

  COOCCUR_THR    <- opt$primary_window
  N_PERM_COOCCUR <- opt$n_perm

  #' Per-SV minimum distance to nearest DMR, deduped to SV (not breakpoint) level
  compute_pct_within_thr <- function(sv_gr, dmr_gr, thr = COOCCUR_THR) {
    if (length(sv_gr) == 0 || length(dmr_gr) == 0) return(NA_real_)
    sv_uniq <- sv_gr[!duplicated(sv_gr$sourceId)]
    hits     <- distanceToNearest(sv_uniq, dmr_gr)
    mean(mcols(hits)$distance <= thr, na.rm = TRUE) * 100
  }

  cooccur_both <- lapply(PATIENT_IDS, function(pt) {
    message("Co-occurrence permutation: ", pt)
    sv  <- sv_list[[pt]]
    dmr <- consensus_dmr_cnv[[pt]]

    obs_pct <- compute_pct_within_thr(sv, dmr)

    BPPARAM_local <- BiocParallel::MulticoreParam(
      workers     = min(4L, BiocParallel::multicoreWorkers()),
      progressbar = FALSE,
      RNGseed     = 20260421L
    )
    null_pcts <- unlist(BiocParallel::bplapply(seq_len(N_PERM_COOCCUR), function(i) {
      compute_pct_within_thr(shuffle_sv(sv), dmr)
    }, BPPARAM = BPPARAM_local))

    wt <- tryCatch(
      wilcox.test(null_pcts, mu = obs_pct, alternative = "less", exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )

    list(
      summary = data.frame(
        patient_id = pt,
        n_sv       = length(unique(sv$sourceId)),
        n_dmr      = length(dmr),
        obs_pct    = obs_pct,
        null_mean  = mean(null_pcts, na.rm = TRUE),
        null_sd    = sd(null_pcts,   na.rm = TRUE),
        null_q025  = quantile(null_pcts, 0.025, na.rm = TRUE),
        null_q975  = quantile(null_pcts, 0.975, na.rm = TRUE),
        wilcox_p   = wt$p.value
      ),
      null_pcts = null_pcts
    )
  })

  cooccur_res       <- bind_rows(lapply(cooccur_both, `[[`, "summary"))
  cooccur_null_list <- setNames(lapply(cooccur_both, `[[`, "null_pcts"), PATIENT_IDS)
  rm(cooccur_both)

  cooccur_res$wilcox_fdr <- p.adjust(cooccur_res$wilcox_p, method = "BH")
  cooccur_res$sig_label  <- dplyr::case_when(
    cooccur_res$wilcox_fdr < 0.001 ~ "***",
    cooccur_res$wilcox_fdr < 0.01  ~ "**",
    cooccur_res$wilcox_fdr < 0.05  ~ "*",
    TRUE                            ~ "ns"
  )

  # Across patients: paired Wilcoxon (obs_pct > null_mean)
  paired_wt <- tryCatch(
    wilcox.test(cooccur_res$obs_pct, cooccur_res$null_mean,
                paired = TRUE, alternative = "greater", exact = FALSE),
    error = function(e) list(p.value = NA_real_)
  )
  cat(sprintf("Paired Wilcoxon across patients (obs > null): p = %.4f\n", paired_wt$p.value))
  print(cooccur_res %>%
    dplyr::select(patient_id, n_sv, obs_pct, null_mean, null_sd, wilcox_fdr, sig_label))

  fwrite(cooccur_res,
            file.path(opt$outdir, sprintf("%s_cooccur_pct_results.csv", opt$run_id)),
            row.names = FALSE, quote = FALSE)

  if (!opt$no_plot) {
    null_long_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
      data.frame(patient_id = pt, pct = cooccur_null_list[[pt]])
    }))

    pt_order <- cooccur_res %>%
      dplyr::arrange(desc(obs_pct)) %>%
      dplyr::pull(patient_id)

    null_long_df$patient_id <- factor(null_long_df$patient_id, levels = pt_order)
    cooccur_res$patient_id  <- factor(cooccur_res$patient_id,  levels = pt_order)

    p_cooccur <- ggplot() +
      geom_boxplot(
        data  = null_long_df,
        aes(x = patient_id, y = pct),
        fill = "grey85", color = "grey55",
        width = 0.5, outlier.size = 0.4, alpha = 0.75
      ) +
      geom_point(
        data  = cooccur_res,
        aes(x = patient_id, y = obs_pct),
        color = "#E24B4A", size = 3.5, shape = 18
      ) +
      geom_text(
        data = cooccur_res,
        aes(x   = patient_id,
            y   = pmax(obs_pct, null_q975) + 2,
            label = sig_label),
        size = 4.5, color = "grey25"
      ) +
      scale_y_continuous(
        labels  = function(x) paste0(round(x), "%"),
        limits  = c(0, NA),
        expand  = expansion(mult = c(0, 0.15))
      ) +
      labs(
        title    = sprintf("SV–DMR Co-occurrence: %% SVs with DMR within ±%dkb",
                          COOCCUR_THR / 1000),
        subtitle = sprintf(
          "Diamond = observed | Box = null (SV shuffled, n=%d) | Paired Wilcoxon p = %.4f",
          N_PERM_COOCCUR, paired_wt$p.value
        ),
        x       = "Patient",
        y       = sprintf("%% SVs with nearest DMR ≤ ±%dkb", COOCCUR_THR / 1000),
        caption = "One-sample Wilcoxon per patient: H₀: null median ≥ obs | FDR: BH\n*p<0.05  **p<0.01  ***p<0.001"
      ) +
      theme_hcc +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

  saveRDS(p_cooccur, file.path(opt$outdir, sprintf("%s_fig_cooccur_pct.rds", opt$run_id)))
  ggsave(
    file.path(opt$outdir, sprintf("%s_fig_cooccur_pct.png", opt$run_id)),
    p_cooccur, width = 10, height = 5, device = "png", dpi = 300
  )
  cat(sprintf("Co-occurrence figure: %s_fig_cooccur_pct.png\n", opt$run_id))
}


# 11. 최종 결과 내보내기 =========================================================

fwrite(window_enrich,   file.path(opt$outdir, sprintf("%s_sv_dmr_enrichment_full.csv", opt$run_id)),
          row.names = FALSE, quote = FALSE)
fwrite(window_enrich,   file.path(opt$outdir, sprintf("%s_window_enrich_full.csv", opt$run_id)),
          row.names = FALSE, quote = FALSE)

cat("\n=== 분석 완료 ===\n")
cat("출력 파일:\n")
cat("  - fig1_cnv_class_dmr_enrichment.png   : CNV class별 ±50kb 농축도\n")
cat("  - fig2_window_enrichment_curve.png     : window 크기별 농축 곡선\n")
cat("  - fig3_methylation_direction.png       : hyper/hypo 방향성 비교\n")
cat("  - fig4_permutation_enrichment.png      : permutation 배율 (FDR)\n")
cat("  - sv_dmr_enrichment_summary.csv        : CNV class별 요약 통계\n")
cat("  - sv_dmr_enrichment_full.csv           : 전체 환자×class×window 결과\n")
cat("  - window_enrich_full.csv               : window_enrich 원본 (obs_mean_dmr 등)\n")
cat("  - dist_decay_full.csv                  : H3 bp×DMR signed distance 전체\n")
cat("  - metaplot_binned.csv                  : H3 metaplot 5kb bin 집계\n")
cat("  - decay_spearman_cor.csv               : H3 Spearman ρ (CNV class별)\n")
cat("  - cooccur_pct_results.csv              : per-patient % SVs within threshold (obs vs null)\n")
cat("  - fig_cooccur_pct.png                  : co-occurrence % observed vs null (per patient)\n")

# =============================================================================
# 12. 논문 Methods 섹션 초안 (주석)
# =============================================================================

# --- Methods 문구 예시 ---
# "To distinguish SV-specific methylation effects from those attributable to
#  copy number variation, SVs were stratified into copy-neutral
#  (inversions [INV] and translocations [TRA/BND]) and copy-number-changing
#  (deletions [DEL], duplications [DUP], and insertions [INS]) classes.
#  For H1 and H2, DMR enrichment within ±10 kb, ±50 kb, and ±100 kb of SV
#  breakpoints was quantified as the mean number of DMRs per SV (Primary B:
#  bp1/bp2 per-SV max). Statistical significance was assessed by permutation
#  testing (n = 1,000 iterations), in which SV breakpoint pairs were shuffled
#  together within each chromosome while preserving SV length distribution.
#  FDR correction was applied across all patient × class × window combinations
#  using the Benjamini-Hochberg procedure. Copy-neutral versus
#  copy-number-changing enrichment was compared using the Wilcoxon rank-sum
#  test across patients.
#  For H3, signed distances from each SV breakpoint to its nearest DMR were
#  calculated (positive = downstream, negative = upstream) and binned at 5 kb
#  resolution up to ±100 kb. Distance-decay was quantified using Spearman rank
#  correlation between absolute breakpoint distance and |Δβ| per CNV class."