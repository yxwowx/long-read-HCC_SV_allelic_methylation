# conda env: renv
# =============================================================================
# TAD/CTCF mechanistic validation (Analyses A, B, C)
# =============================================================================
# 목적: ±50kb window 분석에서 TAD+CTCF disrupting SV가 Non-boundary SV보다
#   낮은 DMR 농축을 보인 결과를 보완 검증.
#
#   A. CTCF anchor 정밀 disruption — breakpoint-CTCF 거리(≤2kb/2-10kb/>10kb)별
#      DMR enrichment 재계산. tier 내 직접 CTCF disrupting SV를 분리하여
#      "TAD+CTCF disrupting"의 heterogeneity 평가.
#
#   B. TAD-level DMR enrichment — ±50kb 대신 SV가 속한 TAD body 전체 내
#      DMR 수를 계산. boundary-disrupting SV가 long-range TAD 재편성 효과로
#      TAD 전체에서 DMR을 보유하는지 검정.
#
#   C. Distance-to-boundary stratification — SV bp ~ TAD boundary 거리 bin별
#      DMR enrichment (Spearman trend test). boundary 근접 SV가 더 큰
#      enrichment를 보이는 monotonic 관계 검정.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
  library(data.table)
  library(BiocParallel)
  library(tibble)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--sv_strat_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz",
              metavar = "FILE",
              help = "Pre-stratified SV file with dist_to_TAD/dist_to_CTCF/tad_condition cols"),
  make_option("--dmr_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/consensus_dmrs_per_patient.csv.gz",
              metavar = "FILE", help = "DMR CSV(.gz) per patient"),
  make_option("--tad_bed", type = "character",
              default = "/node200data/kachungk/reference/GRCh38/3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_tad.bed",
              metavar = "FILE", help = "HepG2 Micro-C TAD body BED"),
  make_option("--ctcf_bed", type = "character",
              default = "/node200data/kachungk/reference/GRCh38/ensembl/HepG2_ChIP_optpeaks_ENCFF543WTP.bed.gz",
              metavar = "FILE", help = "ENCODE HepG2 CTCF peak BED"),
  make_option("--outdir", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/tad_ctcf_validation",
              metavar = "DIR", help = "Output directory"),
  make_option("--run_id", type = "character", default = "tier_v2",
              metavar = "STR", help = "Run ID prefix"),
  make_option("--window", type = "integer", default = 50000L,
              metavar = "INT", help = "Window size for A & C [bp]"),
  make_option("--n_perm", type = "integer", default = 1000L,
              metavar = "INT", help = "Permutation iterations"),
  make_option("--min_sv_per_group", type = "integer", default = 5L,
              metavar = "INT", help = "Min SVs per group per patient"),
  make_option("--no_plot", action = "store_true", default = FALSE,
              help = "Skip ggsave steps")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# ── Visualization constants ─────────────────────────────────────────────────
CTCF_DIST_LEVELS <- c("Direct (≤2kb)", "Proximal (2-10kb)",
                      "Distal (10-50kb)", "Far (>50kb)")
CTCF_DIST_COLORS <- c("Direct (≤2kb)"     = "#C0392B",
                      "Proximal (2-10kb)" = "#E67E22",
                      "Distal (10-50kb)"  = "#F1C40F",
                      "Far (>50kb)"       = "#95A5A6")

BOUNDARY_DIST_LEVELS <- c("In boundary", "≤25kb", "25-100kb",
                          "100-500kb", ">500kb")
BOUNDARY_DIST_COLORS <- c("In boundary" = "#C0392B",
                          "≤25kb"       = "#E67E22",
                          "25-100kb"    = "#F1C40F",
                          "100-500kb"   = "#3B8BD4",
                          ">500kb"      = "#95A5A6")

# ── Helpers ─────────────────────────────────────────────────────────────────

#' Stratify CTCF distance into 4 bins
classify_ctcf_dist <- function(d) {
  factor(dplyr::case_when(
    d <= 2000L   ~ "Direct (≤2kb)",
    d <= 10000L  ~ "Proximal (2-10kb)",
    d <= 50000L  ~ "Distal (10-50kb)",
    TRUE         ~ "Far (>50kb)"
  ), levels = CTCF_DIST_LEVELS)
}

#' Stratify dist-to-boundary into 5 bins
classify_boundary_dist <- function(d) {
  factor(dplyr::case_when(
    d == 0L         ~ "In boundary",
    d <= 25000L     ~ "≤25kb",
    d <= 100000L    ~ "25-100kb",
    d <= 500000L    ~ "100-500kb",
    TRUE            ~ ">500kb"
  ), levels = BOUNDARY_DIST_LEVELS)
}

#' Count DMRs within ±window of breakpoint, per-bp row
count_dmr_near_bp <- function(sv_gr, dmr_gr, window = opt$window) {
  if (length(sv_gr) == 0 || length(dmr_gr) == 0)
    return(data.frame(bp_id = character(), n_dmr = integer(),
                      n_hyper = integer(), n_hypo = integer()))
  if (!"direction" %in% names(mcols(dmr_gr))) {
    mcols(dmr_gr)$direction <- dplyr::case_when(
      mcols(dmr_gr)$diff.Methy >  0 ~ "hypo",
      mcols(dmr_gr)$diff.Methy <  0 ~ "hyper",
      TRUE                          ~ NA_character_)
  }
  wins <- suppressWarnings(GenomicRanges::resize(sv_gr, width = window * 2L, fix = "center"))
  hits <- findOverlaps(wins, dmr_gr, minoverlap = 1L)
  lapply(seq_along(wins), function(i) {
    idx <- subjectHits(hits)[queryHits(hits) == i]
    data.frame(
      bp_id   = sv_gr$bp_id[i],
      sourceId = sv_gr$sourceId[i],
      n_dmr   = length(idx),
      n_hyper = sum(dmr_gr$direction[idx] == "hyper", na.rm = TRUE),
      n_hypo  = sum(dmr_gr$direction[idx] == "hypo",  na.rm = TRUE)
    )
  }) %>% bind_rows()
}

#' Count DMRs within the TAD body that contains the bp (or nearest TAD body if none)
count_dmr_in_tad <- function(sv_gr, dmr_gr, tad_body_gr) {
  if (length(sv_gr) == 0 || length(tad_body_gr) == 0)
    return(data.frame(bp_id = character(), tad_idx = integer(),
                      tad_width = integer(), n_dmr_tad = integer()))

  # find primary TAD: overlapping first, else nearest
  hit_idx <- findOverlaps(sv_gr, tad_body_gr, select = "first")
  miss    <- is.na(hit_idx)
  if (any(miss)) {
    near <- distanceToNearest(sv_gr[miss], tad_body_gr)
    hit_idx[miss][queryHits(near)] <- subjectHits(near)
  }

  out <- data.frame(
    bp_id     = sv_gr$bp_id,
    sourceId  = sv_gr$sourceId,
    tad_idx   = hit_idx,
    tad_width = ifelse(is.na(hit_idx), NA_integer_,
                       width(tad_body_gr)[hit_idx])
  )

  # group bp by TAD → count DMRs once per TAD, then map back
  unique_tads <- unique(na.omit(out$tad_idx))
  tad_n <- integer(length(tad_body_gr))
  if (length(unique_tads) > 0 && length(dmr_gr) > 0) {
    ov <- findOverlaps(tad_body_gr[unique_tads], dmr_gr, minoverlap = 1L)
    counts <- tabulate(queryHits(ov), nbins = length(unique_tads))
    tad_n[unique_tads] <- counts
  }
  out$n_dmr_tad <- ifelse(is.na(out$tad_idx), NA_integer_, tad_n[out$tad_idx])
  out
}

#' Shuffle SV breakpoints (same logic as 03 script): pair-aware intra-chrom shuffle
shuffle_sv <- function(sv_gr) {
  valid <- intersect(names(CHROM_LENS), as.character(seqnames(sv_gr)))
  out   <- sv_gr
  for (sid in unique(out$sourceId)) {
    pair_idx <- which(out$sourceId == sid)
    pair_chr <- as.character(seqnames(out[pair_idx]))
    if (length(unique(pair_chr)) == 1 && unique(pair_chr) %in% valid) {
      chr      <- unique(pair_chr)
      span     <- max(end(out[pair_idx])) - min(start(out[pair_idx]))
      max_orig <- CHROM_LENS[chr] - span - 1L
      if (max_orig < 1L) next
      new_origin <- sample.int(as.integer(max_orig), 1L)
      offset     <- new_origin - min(start(out[pair_idx]))
      ranges(out)[pair_idx] <- IRanges(start = start(out[pair_idx]) + offset,
                                       width = width(out[pair_idx]))
    } else {
      for (i in pair_idx) {
        chr <- as.character(seqnames(out[i]))
        if (!chr %in% valid) next
        max_s <- CHROM_LENS[chr] - width(out[i])
        if (max_s < 1L) next
        ranges(out)[i] <- IRanges(start = sample.int(as.integer(max_s), 1L),
                                  width = width(out[i]))
      }
    }
  }
  out
}

#' Permutation null for window enrichment (mean n_dmr per SV)
perm_window <- function(sv_gr, dmr_gr, window = opt$window,
                        n_perm = opt$n_perm, BPPARAM = NULL) {
  if (length(sv_gr) == 0) return(numeric(0))
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(
      workers = min(4L, BiocParallel::multicoreWorkers()),
      progressbar = FALSE, RNGseed = 20260518L)
  }
  out <- BiocParallel::bplapply(seq_len(n_perm), function(i) {
    shuf <- shuffle_sv(sv_gr)
    res  <- count_dmr_near_bp(shuf, dmr_gr, window)
    res %>% dplyr::group_by(sourceId) %>%
      dplyr::summarise(n = max(n_dmr), .groups = "drop") %>%
      dplyr::pull(n) %>% mean(na.rm = TRUE)
  }, BPPARAM = BPPARAM)
  as.numeric(unlist(out, use.names = FALSE))
}

#' Permutation null for TAD-level enrichment (mean n_dmr_tad per SV)
perm_tad <- function(sv_gr, dmr_gr, tad_body_gr,
                     n_perm = opt$n_perm, BPPARAM = NULL) {
  if (length(sv_gr) == 0) return(numeric(0))
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(
      workers = min(4L, BiocParallel::multicoreWorkers()),
      progressbar = FALSE, RNGseed = 20260518L)
  }
  out <- BiocParallel::bplapply(seq_len(n_perm), function(i) {
    shuf <- shuffle_sv(sv_gr)
    res  <- count_dmr_in_tad(shuf, dmr_gr, tad_body_gr)
    res %>% dplyr::group_by(sourceId) %>%
      dplyr::summarise(n = max(n_dmr_tad, na.rm = TRUE), .groups = "drop") %>%
      dplyr::pull(n) %>% mean(na.rm = TRUE)
  }, BPPARAM = BPPARAM)
  as.numeric(unlist(out, use.names = FALSE))
}

# ── 1. Load data ────────────────────────────────────────────────────────────
message("Reading SV stratification: ", opt$sv_strat_file)
sv_df_all <- fread(opt$sv_strat_file)
sv_df_all[, ctcf_dist_class     := classify_ctcf_dist(dist_to_CTCF)]
sv_df_all[, boundary_dist_class := classify_boundary_dist(dist_to_TAD)]

sv_list <- sv_df_all %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$sample)

PATIENT_IDS <- names(sv_list)
message(sprintf("Patients: %d  |  Total SV bp rows: %d",
                length(PATIENT_IDS), nrow(sv_df_all)))

message("Reading DMRs: ", opt$dmr_file)
consensus_dmr <- local({
  raw <- fread(opt$dmr_file)
  if (grepl("^P\\d+$", raw$sample[1])) {
    raw %>% GRanges() %>% split(mcols(.)$sample)
  } else if ("patient_code" %in% colnames(raw) && grepl("^P\\d+$", raw$patient_code[1])) {
    raw %>% dplyr::rename(patient_name = sample, sample = patient_code) %>%
      GRanges() %>% split(mcols(.)$sample)
  } else {
    raw %>% anonym_sample() %>%
      dplyr::rename(patient_name = sample, sample = patient_code) %>%
      GRanges() %>% split(mcols(.)$sample)
  }
})
# add direction once
consensus_dmr <- lapply(consensus_dmr, function(d) {
  mcols(d)$direction <- dplyr::case_when(
    mcols(d)$diff.Methy >  0 ~ "hypo",
    mcols(d)$diff.Methy <  0 ~ "hyper",
    TRUE                     ~ NA_character_)
  d
})

message("Reading TAD bodies: ", opt$tad_bed)
tad_body_gr <- fread(opt$tad_bed) %>%
  dplyr::rename(seqnames = V1, start = V2, end = V3) %>%
  makeGRangesFromDataFrame() %>% sort()
message(sprintf("TAD bodies: %d  |  median width: %.1f kb",
                length(tad_body_gr),
                median(width(tad_body_gr)) / 1000))

# ── 2. Cross-distribution diagnostics ───────────────────────────────────────
cat("\n=== CTCF dist × tier cross-distribution ===\n")
cross_ctcf <- sv_df_all %>%
  dplyr::count(stratification, ctcf_dist_class) %>%
  tidyr::pivot_wider(names_from = ctcf_dist_class, values_from = n, values_fill = 0L)
print(cross_ctcf)
fwrite(cross_ctcf, file.path(opt$outdir, sprintf("%s_cross_tier_ctcf_dist.csv", opt$run_id)))

cat("\n=== Boundary dist × tier cross-distribution ===\n")
cross_bound <- sv_df_all %>%
  dplyr::count(stratification, boundary_dist_class) %>%
  tidyr::pivot_wider(names_from = boundary_dist_class, values_from = n, values_fill = 0L)
print(cross_bound)
fwrite(cross_bound, file.path(opt$outdir, sprintf("%s_cross_tier_boundary_dist.csv", opt$run_id)))

# =============================================================================
# Analysis A: CTCF anchor precise disruption
# =============================================================================
cat("\n\n=== Analysis A: CTCF anchor precise disruption ===\n")

A_res <- lapply(PATIENT_IDS, function(pt) {
  message("[A] ", pt)
  sv  <- sv_list[[pt]]
  dmr <- consensus_dmr[[pt]]
  if (is.null(dmr) || length(dmr) == 0) return(NULL)

  lapply(CTCF_DIST_LEVELS, function(cls) {
    sv_sub <- sv[as.character(mcols(sv)$ctcf_dist_class) == cls]
    if (length(unique(sv_sub$sourceId)) < opt$min_sv_per_group) return(NULL)

    res <- count_dmr_near_bp(sv_sub, dmr, opt$window) %>%
      dplyr::group_by(sourceId) %>%
      dplyr::summarise(n_dmr   = max(n_dmr),
                       n_hyper = max(n_hyper),
                       n_hypo  = max(n_hypo),
                       .groups = "drop")
    if (nrow(res) == 0) return(NULL)

    obs_mean <- mean(res$n_dmr, na.rm = TRUE)
    null_vec <- perm_window(sv_sub, dmr, opt$window, opt$n_perm)

    data.frame(
      patient_id       = pt,
      ctcf_dist_class  = cls,
      n_sv             = nrow(res),
      obs_mean_dmr     = obs_mean,
      null_mean        = mean(null_vec),
      null_sd          = sd(null_vec),
      enrichment_ratio = obs_mean / pmax(mean(null_vec), 1e-6),
      p_perm           = mean(null_vec >= obs_mean),
      pct_sv_w_dmr     = mean(res$n_dmr > 0) * 100,
      mean_hyper       = mean(res$n_hyper, na.rm = TRUE),
      mean_hypo        = mean(res$n_hypo,  na.rm = TRUE)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

A_res$p_fdr <- p.adjust(A_res$p_perm, method = "BH")
A_res$sig_label <- dplyr::case_when(
  A_res$p_fdr < 0.001 ~ "***",
  A_res$p_fdr < 0.01  ~ "**",
  A_res$p_fdr < 0.05  ~ "*",
  TRUE                ~ "ns")
A_res$ctcf_dist_class <- factor(A_res$ctcf_dist_class, levels = CTCF_DIST_LEVELS)

fwrite(A_res, file.path(opt$outdir, sprintf("%s_A_ctcf_disruption.csv", opt$run_id)))

# Trend test: enrichment_ratio across CTCF dist bins (Spearman)
A_trend <- A_res %>%
  dplyr::mutate(rank_class = as.integer(ctcf_dist_class)) %>%
  dplyr::summarise(
    spearman_rho = cor(rank_class, enrichment_ratio,
                       method = "spearman", use = "complete.obs"),
    p_value      = tryCatch(
      cor.test(rank_class, enrichment_ratio,
               method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_),
    n_obs        = dplyr::n()
  )
cat("[A] Spearman ρ (enrichment_ratio ~ CTCF dist bin):\n")
print(A_trend)
fwrite(A_trend, file.path(opt$outdir, sprintf("%s_A_trend_spearman.csv", opt$run_id)))

# Wilcoxon: Direct vs Far
A_wilcox <- tryCatch({
  direct_v <- A_res$enrichment_ratio[A_res$ctcf_dist_class == "Direct (≤2kb)"]
  far_v    <- A_res$enrichment_ratio[A_res$ctcf_dist_class == "Far (>50kb)"]
  if (length(direct_v) >= 3 && length(far_v) >= 3)
    wilcox.test(direct_v, far_v, alternative = "greater", exact = FALSE)
  else list(p.value = NA_real_)
}, error = function(e) list(p.value = NA_real_))
cat(sprintf("[A] Wilcoxon (Direct > Far) p = %.4f\n", A_wilcox$p.value))

if (!opt$no_plot) {
  pA <- ggplot(A_res, aes(x = ctcf_dist_class, y = enrichment_ratio,
                          fill = ctcf_dist_class)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.7) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 1.2) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    scale_fill_manual(values = CTCF_DIST_COLORS, guide = "none") +
    labs(
      title    = sprintf("A. DMR enrichment by breakpoint-CTCF distance (±%dkb window)",
                        opt$window / 1000),
      subtitle = sprintf("Spearman ρ = %.3f (p = %.3g) | Direct vs Far Wilcoxon p = %.3g",
                         A_trend$spearman_rho, A_trend$p_value, A_wilcox$p.value),
      x = NULL, y = "Enrichment ratio (obs / null)",
      caption = sprintf("n_perm = %d | FDR: BH", opt$n_perm)
    ) + theme_hcc

  ggsave(file.path(opt$outdir, sprintf("%s_A_ctcf_disruption.png", opt$run_id)),
         pA, width = 8, height = 5, device = "png", dpi = 300)
  saveRDS(pA, file.path(opt$outdir, sprintf("%s_A_ctcf_disruption.rds", opt$run_id)))
}

# =============================================================================
# Analysis B: TAD-level DMR enrichment
# =============================================================================
cat("\n\n=== Analysis B: TAD-level DMR enrichment ===\n")

B_res <- lapply(PATIENT_IDS, function(pt) {
  message("[B] ", pt)
  sv  <- sv_list[[pt]]
  dmr <- consensus_dmr[[pt]]
  if (is.null(dmr) || length(dmr) == 0) return(NULL)

  lapply(STRAT_LEVELS, function(tier) {
    sv_sub <- sv[as.character(mcols(sv)$stratification) == tier]
    if (length(unique(sv_sub$sourceId)) < opt$min_sv_per_group) return(NULL)

    res <- count_dmr_in_tad(sv_sub, dmr, tad_body_gr) %>%
      dplyr::group_by(sourceId) %>%
      dplyr::summarise(n_dmr_tad = max(n_dmr_tad, na.rm = TRUE),
                       tad_width = max(tad_width, na.rm = TRUE),
                       .groups = "drop")
    if (nrow(res) == 0) return(NULL)

    obs_mean <- mean(res$n_dmr_tad, na.rm = TRUE)
    null_vec <- perm_tad(sv_sub, dmr, tad_body_gr, opt$n_perm)

    data.frame(
      patient_id       = pt,
      tier             = tier,
      n_sv             = nrow(res),
      obs_mean_dmr_tad = obs_mean,
      null_mean        = mean(null_vec),
      null_sd          = sd(null_vec),
      enrichment_ratio = obs_mean / pmax(mean(null_vec), 1e-6),
      p_perm           = mean(null_vec >= obs_mean),
      median_tad_width = median(res$tad_width, na.rm = TRUE),
      dmr_density_tad  = mean(res$n_dmr_tad / pmax(res$tad_width, 1L) * 1e6, na.rm = TRUE)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

B_res$p_fdr <- p.adjust(B_res$p_perm, method = "BH")
B_res$sig_label <- dplyr::case_when(
  B_res$p_fdr < 0.001 ~ "***",
  B_res$p_fdr < 0.01  ~ "**",
  B_res$p_fdr < 0.05  ~ "*",
  TRUE                ~ "ns")
B_res$tier <- factor(B_res$tier, levels = STRAT_LEVELS)

fwrite(B_res, file.path(opt$outdir, sprintf("%s_B_tad_level.csv", opt$run_id)))

# Wilcoxon: TAD+CTCF disrupting (or any boundary) vs Non-boundary at TAD scale
B_wilcox <- tryCatch({
  bd <- B_res$enrichment_ratio[B_res$tier %in%
    c("TAD+CTCF disrupting","CTCF-only","TAD-only")]
  nb <- B_res$enrichment_ratio[B_res$tier == "Non-boundary"]
  wilcox.test(bd, nb, alternative = "greater", exact = FALSE)
}, error = function(e) list(p.value = NA_real_))
cat(sprintf("[B] Wilcoxon (boundary-disrupting > Non-boundary) p = %.4f\n", B_wilcox$p.value))

if (!opt$no_plot) {
  pB <- ggplot(B_res, aes(x = tier, y = enrichment_ratio, fill = tier)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.7) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 1.2) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    scale_fill_manual(values = STRAT_COLORS, guide = "none") +
    labs(
      title    = "B. TAD-level DMR enrichment (DMRs within containing TAD body)",
      subtitle = sprintf("Boundary-disrupting vs Non-boundary Wilcoxon (one-sided) p = %.3g",
                         B_wilcox$p.value),
      x = NULL, y = "Enrichment ratio (obs / null)",
      caption = sprintf("Null = SV shuffled within chromosome | n_perm = %d", opt$n_perm)
    ) + theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(opt$outdir, sprintf("%s_B_tad_level.png", opt$run_id)),
         pB, width = 9, height = 5, device = "png", dpi = 300)
  saveRDS(pB, file.path(opt$outdir, sprintf("%s_B_tad_level.rds", opt$run_id)))

  # Side plot: TAD-level vs window-level enrichment scatter requires window data;
  # instead show DMR density (per Mb of TAD)
  pB_dens <- ggplot(B_res, aes(x = tier, y = dmr_density_tad, fill = tier)) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 1.2) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    scale_fill_manual(values = STRAT_COLORS, guide = "none") +
    labs(title    = "B. DMR density within TAD body (per Mb)",
         subtitle = "Density normalises TAD-size differences",
         x = NULL, y = "DMR / Mb of TAD") +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(opt$outdir, sprintf("%s_B_tad_density.png", opt$run_id)),
         pB_dens, width = 9, height = 5, device = "png", dpi = 300)
  saveRDS(pB_dens, file.path(opt$outdir, sprintf("%s_B_tad_density.rds", opt$run_id)))
}

# =============================================================================
# Analysis C: Distance-to-boundary stratification
# =============================================================================
cat("\n\n=== Analysis C: Distance-to-boundary stratification ===\n")

C_res <- lapply(PATIENT_IDS, function(pt) {
  message("[C] ", pt)
  sv  <- sv_list[[pt]]
  dmr <- consensus_dmr[[pt]]
  if (is.null(dmr) || length(dmr) == 0) return(NULL)

  lapply(BOUNDARY_DIST_LEVELS, function(cls) {
    sv_sub <- sv[as.character(mcols(sv)$boundary_dist_class) == cls]
    if (length(unique(sv_sub$sourceId)) < opt$min_sv_per_group) return(NULL)

    res <- count_dmr_near_bp(sv_sub, dmr, opt$window) %>%
      dplyr::group_by(sourceId) %>%
      dplyr::summarise(n_dmr   = max(n_dmr),
                       n_hyper = max(n_hyper),
                       n_hypo  = max(n_hypo),
                       .groups = "drop")
    if (nrow(res) == 0) return(NULL)

    obs_mean <- mean(res$n_dmr, na.rm = TRUE)
    null_vec <- perm_window(sv_sub, dmr, opt$window, opt$n_perm)

    data.frame(
      patient_id           = pt,
      boundary_dist_class  = cls,
      n_sv                 = nrow(res),
      obs_mean_dmr         = obs_mean,
      null_mean            = mean(null_vec),
      null_sd              = sd(null_vec),
      enrichment_ratio     = obs_mean / pmax(mean(null_vec), 1e-6),
      p_perm               = mean(null_vec >= obs_mean),
      pct_sv_w_dmr         = mean(res$n_dmr > 0) * 100
    )
  }) %>% bind_rows()
}) %>% bind_rows()

C_res$p_fdr <- p.adjust(C_res$p_perm, method = "BH")
C_res$sig_label <- dplyr::case_when(
  C_res$p_fdr < 0.001 ~ "***",
  C_res$p_fdr < 0.01  ~ "**",
  C_res$p_fdr < 0.05  ~ "*",
  TRUE                ~ "ns")
C_res$boundary_dist_class <- factor(C_res$boundary_dist_class,
                                    levels = BOUNDARY_DIST_LEVELS)

fwrite(C_res, file.path(opt$outdir, sprintf("%s_C_boundary_dist.csv", opt$run_id)))

# Spearman trend test on enrichment_ratio across bins
C_trend <- C_res %>%
  dplyr::mutate(rank_class = as.integer(boundary_dist_class)) %>%
  dplyr::summarise(
    spearman_rho = cor(rank_class, enrichment_ratio,
                       method = "spearman", use = "complete.obs"),
    p_value      = tryCatch(
      cor.test(rank_class, enrichment_ratio,
               method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_),
    n_obs        = dplyr::n()
  )
cat("[C] Spearman ρ (enrichment_ratio ~ boundary dist bin):\n")
print(C_trend)
fwrite(C_trend, file.path(opt$outdir, sprintf("%s_C_trend_spearman.csv", opt$run_id)))

# Wilcoxon: In-boundary vs >500kb (predicted: in-boundary > far)
C_wilcox <- tryCatch({
  in_v  <- C_res$enrichment_ratio[C_res$boundary_dist_class == "In boundary"]
  far_v <- C_res$enrichment_ratio[C_res$boundary_dist_class == ">500kb"]
  if (length(in_v) >= 3 && length(far_v) >= 3)
    wilcox.test(in_v, far_v, alternative = "greater", exact = FALSE)
  else list(p.value = NA_real_)
}, error = function(e) list(p.value = NA_real_))
cat(sprintf("[C] Wilcoxon (In boundary > >500kb) p = %.4f\n", C_wilcox$p.value))

if (!opt$no_plot) {
  pC <- ggplot(C_res, aes(x = boundary_dist_class, y = enrichment_ratio,
                          fill = boundary_dist_class)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.7) +
    geom_boxplot(width = 0.5, alpha = 0.8, outlier.size = 1.2) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "grey30") +
    scale_fill_manual(values = BOUNDARY_DIST_COLORS, guide = "none") +
    labs(
      title    = sprintf("C. DMR enrichment by breakpoint-boundary distance (±%dkb window)",
                        opt$window / 1000),
      subtitle = sprintf("Spearman ρ = %.3f (p = %.3g) | In-boundary vs >500kb p = %.3g",
                         C_trend$spearman_rho, C_trend$p_value, C_wilcox$p.value),
      x = NULL, y = "Enrichment ratio (obs / null)",
      caption = sprintf("n_perm = %d | FDR: BH", opt$n_perm)
    ) + theme_hcc

  ggsave(file.path(opt$outdir, sprintf("%s_C_boundary_dist.png", opt$run_id)),
         pC, width = 9, height = 5, device = "png", dpi = 300)
  saveRDS(pC, file.path(opt$outdir, sprintf("%s_C_boundary_dist.rds", opt$run_id)))

  # Composite Figure
  fig_combined <- (pA | pB | pC) +
    plot_annotation(
      title   = "TAD/CTCF mechanistic validation",
      caption = sprintf("n_perm = %d | window = ±%dkb | TAD body from HepG2 Micro-C",
                        opt$n_perm, opt$window / 1000),
      theme   = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(opt$outdir, sprintf("%s_ABC_combined.png", opt$run_id)),
         fig_combined, width = 18, height = 5.5, device = "png", dpi = 300)
}

# ── Summary log ─────────────────────────────────────────────────────────────
cat("\n\n=== Analysis A summary (median enrichment_ratio per CTCF dist bin) ===\n")
print(A_res %>% dplyr::group_by(ctcf_dist_class) %>%
  dplyr::summarise(n_patients = dplyr::n(),
                   median_ratio = round(median(enrichment_ratio, na.rm = TRUE), 3),
                   pct_significant = round(mean(p_fdr < 0.05, na.rm = TRUE) * 100, 1),
                   .groups = "drop"))

cat("\n=== Analysis B summary (median enrichment_ratio per tier) ===\n")
print(B_res %>% dplyr::group_by(tier) %>%
  dplyr::summarise(n_patients = dplyr::n(),
                   median_ratio   = round(median(enrichment_ratio, na.rm = TRUE), 3),
                   median_density = round(median(dmr_density_tad,  na.rm = TRUE), 2),
                   pct_significant = round(mean(p_fdr < 0.05, na.rm = TRUE) * 100, 1),
                   .groups = "drop"))

cat("\n=== Analysis C summary (median enrichment_ratio per boundary dist bin) ===\n")
print(C_res %>% dplyr::group_by(boundary_dist_class) %>%
  dplyr::summarise(n_patients = dplyr::n(),
                   median_ratio = round(median(enrichment_ratio, na.rm = TRUE), 3),
                   pct_significant = round(mean(p_fdr < 0.05, na.rm = TRUE) * 100, 1),
                   .groups = "drop"))

cat("\n=== Analyses A/B/C complete ===\n")
cat(sprintf("Outputs in: %s\n", opt$outdir))
cat("  - <run>_A_ctcf_disruption.{csv,png,rds}\n")
cat("  - <run>_B_tad_level.{csv,png,rds}  +  <run>_B_tad_density.{png,rds}\n")
cat("  - <run>_C_boundary_dist.{csv,png,rds}\n")
cat("  - <run>_ABC_combined.png\n")
cat("  - <run>_cross_tier_{ctcf_dist,boundary_dist}.csv\n")
cat("  - <run>_{A,C}_trend_spearman.csv\n")
