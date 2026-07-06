#!/usr/bin/env Rscript
# Haplotype-SV–aDMR analysis: sv_tier-stratified phase block co-localisation + HP-specific tests
# Primary stratification: sv_tier (TAD_CTCF > TAD_only > CTCF_only > copy_neutral > non_boundary)
# Run: mamba run -n renv Rscript haplotype_sv_admr_analysis.R --sv_strat_file <file.csv.gz>

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(annotatr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(data.table)
  library(GenomicRanges)
  library(stringr)
  library(rtracklayer)
  library(StructuralVariantAnnotation)
  library(VariantAnnotation)
  library(SummarizedExperiment)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(BiocParallel)
  library(optparse)
  library(clinfun)      # jonckheere.test
  library(dunn.test)    # dunn.test post-hoc
  library(lme4)         # lmer
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--sv_strat_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.csv.gz",
              metavar = "FILE",
              help = "Pre-stratified SV file from sv_stratification.R. Required for tier-stratified analyses."),
  make_option("--input_dmr", type = "character",
               default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz", #nolint
              metavar = "FILE",
              help = "Input DMR file with phased beta values (output from dmr_recurrence_analysis.R)."),
  make_option("--ctcf_bed", type = "character",
              default = "/node200data/kachungk/reference/GRCh38/ensembl/HepG2_ChIP_optpeaks_ENCFF543WTP.bed.gz",
              metavar = "FILE",
              help = "CTCF peak BED (narrowPeak; no header)."),
  make_option("--enhancer_bed", type = "character",
              default = "/node200data/kachungk/reference/GRCh38/genomic_element/hg38_genehancer_enhancer.bed",
              metavar = "FILE",
              help = "Enhancer BED (header row: CHROMOSOME,START,END,ELEMENT,SEGMENT)."),
  make_option("--promoter_bed", type = "character",
              default = "/node200data/kachungk/reference/GRCh38/genomic_element/hg38_genes_promoters.bed",
              metavar = "FILE",
              help = "Promoter BED (header row: CHROMOSOME,START,END,SEGMENT,ELEMENT)."),
  make_option("--outdir", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/03.haplotype_sv_admr_analysis",
              metavar = "DIR",
              help = "Output directory for results and figures."),
  make_option("--stratify_by", type = "character", default = "tier",
              metavar = "STRING",
              help = "Stratification variable: 'tier' (default, SV tier hierarchy) or 'svtype' (geom_type column).")
)
opt <- parse_args(OptionParser(option_list = option_list))

BPPARAM <- BiocParallel::MulticoreParam(
  workers     = min(4L, BiocParallel::multicoreWorkers()),
  progressbar = FALSE,
  RNGseed     = 20260421L
)

N_PERM_BLOCK <- 500L   # block-label permutation (Layer 1)
N_PERM_JT    <- 2000L  # Jonckheere-Terpstra permutation


# 1. Data loading =============================================================
setwd("/node200data/kachungk/hcc_data/DMR_SVs")
phase_gtfs <- list.files(
  "../hg38+HBV/clairS/phased_vcf",
  pattern = "*.gtf$", full.names = TRUE
)
phase_blocks <- lapply(phase_gtfs, function(x) {
  gr <- import(x)
  mcols(gr)$sample <- str_remove(basename(x), ".gtf$")
  as.data.frame(gr)
}) %>%
  bind_rows() %>%
  dplyr::select(seqnames, start, end, phase, gene_id, sample) %>%
  anonym_sample() %>%
  dplyr::rename(
    patient_name = sample,
    sample       = patient_code,
    block_id     = gene_id # gene_id is phase block ID
  ) %>%
  GRanges() %>%
  split(mcols(.)$sample)

if (!is.null(opt$sv_strat_file)) {
  message("Reading pre-stratified SV file: ", opt$sv_strat_file)
  sv_phased <- fread(opt$sv_strat_file) %>%
    dplyr::rename(
      hp       = HP,
      svlen    = svLen,
      qual     = QUAL,
      block_id = PHASESETID
    ) %>%
    dplyr::mutate(
      sv_tier = factor(TIER_RECODE[stratification],
                       levels = SV_TIER_LEVELS),
      sv_type = geom_type
    ) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    split(mcols(.)$sample) %>%
    lapply(function(gr) { names(gr) <- mcols(gr)$bp_id; gr })
} else {
  stop("Note: sv_tier not available from VCF. Run with --sv_strat_file for tier-stratified analyses.")
} # end else (VCF loading)
PATIENT_IDS <- names(sv_phased)

# Override stratification to SV type if --stratify_by svtype
if (opt$stratify_by == "svtype") {
  SV_TIER_LEVELS <- c("DEL", "DUP", "INV", "TRA", "INS", "COM")
  SV_TIER_COLORS <- c(
    DEL = "#BA7517", DUP = "#E24B4A",
    INV = "#3B8BD4", TRA = "#7F77DD", INS = "#1D9E75", COM = "#888780"
  )
  sv_phased <- lapply(sv_phased, function(gr) {
    mcols(gr)$sv_tier <- factor(mcols(gr)$sv_type, levels = SV_TIER_LEVELS)
    gr
  })
  message("Stratifying by SV type (geom_type). JT tests replaced with KW (no natural order).")
}
STRAT_LABEL <- if (opt$stratify_by == "svtype") "SV type (geom_type)" else "SV tier"

# Revised 2D tier variables (added alongside sv_tier, not replacing it)
#   sv_arch : boundary (TAD_CTCF|TAD_only|CTCF_only) vs non_boundary (copy_neutral|non_boundary)
#   sv_bio  : CN_changing (DEL|DUP) vs balanced (INV|TRA|INS|COM)
#   arch_bio: 2D cross product — primary axis for revised analyses
ARCH_BIO_LEVELS <- c("boundary_CN_changing", "boundary_balanced",
                      "non_boundary_CN_changing", "non_boundary_balanced")
ARCH_BIO_COLORS <- c(
  boundary_CN_changing     = "#7F77DD",
  boundary_balanced        = "#3B8BD4",
  non_boundary_CN_changing = "#BA7517",
  non_boundary_balanced    = "#888780"
)

sv_phased <- lapply(sv_phased, function(gr) {
  tier_chr <- as.character(mcols(gr)$sv_tier)
  type_chr <- as.character(mcols(gr)$sv_type)

  mcols(gr)$sv_arch <- factor(
    ifelse(tier_chr %in% c("TAD_CTCF", "TAD_only", "CTCF_only"),
           "boundary", "non_boundary"),
    levels = c("boundary", "non_boundary")
  )
  mcols(gr)$sv_bio <- factor(
    ifelse(type_chr %in% c("DEL", "DUP"), "CN_changing", "balanced"),
    levels = c("CN_changing", "balanced")
  )
  mcols(gr)$arch_bio <- factor(
    paste(as.character(mcols(gr)$sv_arch), as.character(mcols(gr)$sv_bio), sep = "_"),
    levels = ARCH_BIO_LEVELS
  )
  gr
})

# Use confident DMRs with phased beta values from dmr_recurrence_analysis.R output
message("Reading aDMR file: ", opt$input_dmr)
admr_phased <- fread(opt$input_dmr) %>%
  dplyr::rename(
    patient_name    = sample,
    sample          = patient_code,
    dmr_chr         = seqnames,
    dmr_start       = start,
    dmr_end         = end,
    normal_meth     = normal.Methy,
    tumor_meth      = tumor.Methy,
    hp1_beta        = HP1.Methy,
    hp2_beta        = HP2.Methy,
    delta_beta_bulk = diff.Methy,
    n_cpg           = nCG
  ) %>%
  dplyr::mutate(hp_delta = hp1_beta - hp2_beta) %>%
  dplyr::select(
    starts_with("admr_"), starts_with("hp"), starts_with("dmr"),
    delta_beta_bulk, n_cpg,
    sample, patient_name, n_patients, pct_patients
  ) %>%
  makeGRangesFromDataFrame(
    seqnames.field     = "admr_chr",
    start.field        = "admr_start",
    end.field          = "admr_end",
    keep.extra.columns = TRUE
  ) %>%
  split(mcols(.)$sample)

## 1-1. Add block_id to admr_phased via overlap with phase_blocks ===============
message("Preprocess aDMR dataset")
admr_phased <- lapply(names(admr_phased), function(pt) {
  dmr_gr <- admr_phased[[pt]]
  blk_gr <- phase_blocks[[pt]]

  if (is.null(blk_gr) || length(blk_gr) == 0) {
    mcols(dmr_gr)$block_id <- NA
    return(dmr_gr)
  }

  hits <- findOverlaps(dmr_gr, blk_gr, select = "first")
  mcols(dmr_gr)$block_id <- blk_gr$block_id[hits]
  dmr_gr[!is.na(mcols(dmr_gr)$block_id)]
}) %>% setNames(PATIENT_IDS)

cat(sprintf("Patients loaded: %d — %s\n", length(PATIENT_IDS),
            paste(PATIENT_IDS, collapse = ", ")))

has_tier <- !is.null(opt$sv_strat_file)
if (!has_tier) {
  stop(
    "ERROR: --sv_strat_file not provided.
    Tier-stratified Layers 1-4 require this argument."
  )
}

# Helper functions ============================================================

#' Infer SV-bearing haplotype per phase block
#'
#' Note: HP1/HP2 labels are arbitrary (not global WT/SV). We therefore define
#' "SV haplotype" per *phase block* using SV evidence, and define WT haplotype
#' as the opposite HP. Blocks with SVs on both haplotypes are treated as
#' ambiguous and excluded from HP-specific WT-vs-SV comparisons.
#'
#' @param hp_vec Integer vector of SV haplotype calls (expected 1/2; may contain 0/NA)
#' @param block_id_vec Character/integer vector of phase block IDs
#' @return Named integer vector: block_id → SV haplotype (1 or 2)
infer_sv_hp_map <- function(hp_vec, block_id_vec) {
  ok <- !is.na(block_id_vec)
  hp_vec       <- as.integer(hp_vec[ok])
  block_id_vec <- as.character(block_id_vec[ok])

  spl <- split(hp_vec, block_id_vec)
  res <- vapply(spl, function(x) {
    x <- x[!is.na(x) & x %in% c(1L, 2L)]
    if (length(x) == 0L) return(NA_integer_)
    ux <- unique(x)
    if (length(ux) != 1L) return(NA_integer_)  # SVs on both HPs → ambiguous
    ux[[1]]
  }, integer(1))

  res <- res[!is.na(res)]
  res
}

#' Extract SV-haplotype vs WT-haplotype beta values for aDMRs in same phase blocks
#' @param admr_gr GRanges of aDMRs with block_id metadata column
#' @param sv_block_ids Character vector of block_ids with SVs (kept for API compatibility)
#' @param sv_hp_map Named integer vector (block_id → SV haplotype 1 or 2)
get_sv_hp_beta <- function(admr_gr, sv_block_ids, sv_hp_map) {

  sv_hp_map <- sv_hp_map[!is.na(sv_hp_map) & sv_hp_map %in% c(1L, 2L)]
  if (length(sv_hp_map) == 0L) return(NULL)

  in_sv <- admr_gr$block_id %in% names(sv_hp_map)
  if (!any(in_sv)) return(NULL)

  admr_sub  <- admr_gr[in_sv]
  sv_hp_vec <- as.integer(sv_hp_map[as.character(admr_sub$block_id)])
  wt_hp_vec <- ifelse(sv_hp_vec == 1L, 2L, 1L)

  sv_beta <- ifelse(sv_hp_vec == 1L, admr_sub$hp1_beta, admr_sub$hp2_beta)
  wt_beta <- ifelse(sv_hp_vec == 1L, admr_sub$hp2_beta, admr_sub$hp1_beta)

  data.frame(
    block_id     = admr_sub$block_id,
    sv_hp        = sv_hp_vec,
    wt_hp        = wt_hp_vec,
    sv_hp_beta   = sv_beta,
    wt_hp_beta   = wt_beta,
    sv_minus_wt  = sv_beta - wt_beta,
    bulk_delta   = admr_sub$delta_beta_bulk,
    n_cpg        = admr_sub$n_cpg,
    hp_delta_abs = abs(admr_sub$hp_delta)
  )
}

#' Annotate aDMR GRanges with logical CRE overlap flags
#' @param dmr_gr GRanges of aDMRs
#' @param promoter_gr GRanges of promoters
#' @param ctcf_gr GRanges of CTCF sites
#' @param enhancer_gr GRanges of enhancers
annotate_dmr_cre <- function(dmr_gr, promoter_gr, ctcf_gr, enhancer_gr) {
  mcols(dmr_gr)$overlaps_promoter <- IRanges::overlapsAny(dmr_gr, promoter_gr)
  mcols(dmr_gr)$overlaps_ctcf     <- IRanges::overlapsAny(dmr_gr, ctcf_gr)
  mcols(dmr_gr)$overlaps_enhancer <- IRanges::overlapsAny(dmr_gr, enhancer_gr)
  mcols(dmr_gr)$overlaps_any_cre  <- mcols(dmr_gr)$overlaps_promoter |
                                      mcols(dmr_gr)$overlaps_ctcf     |
                                      mcols(dmr_gr)$overlaps_enhancer
  dmr_gr
}

# Initialise all result objects (updated inside has_tier blocks)
layer1_results          <- data.frame()
layer2_results          <- data.frame()
all_hp_df               <- data.frame()
layer3_results          <- data.frame()
layer4_gradient         <- data.frame()
layer4_model_comparison <- data.frame()
jt_l1   <- list(p.value = NA_real_)
kw_l1   <- list(p.value = NA_real_)
kw_l2   <- list(p.value = NA_real_)
jt_l2   <- list(p.value = NA_real_)
wt_l3   <- list(p.value = NA_real_)
jt_l3   <- list(p.value = NA_real_)
kw_l3   <- list(p.value = NA_real_)
aic_tier <- NA_real_
aic_dist <- NA_real_


# 2. Layer 1 — Tier-stratified phase block co-localization ===================
#    Width-corrected aDMR density; block-label permutation (n=500) + Wilcoxon;
#    Jonckheere-Terpstra across tier hierarchy

cat("\n=== Layer 1: Tier-stratified phase block co-localization ===\n")

if (!has_tier) {
  message("Skipping Layer 1: --sv_strat_file required.")
} else {

layer1_results <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_phased[[pt]]
  blks   <- phase_blocks[[pt]]
  dmr_gr <- admr_phased[[pt]]
  if (is.null(sv_gr)  || length(sv_gr)  == 0 ||
      is.null(blks)   || length(blks)   == 0 ||
      is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  # Width-corrected aDMR density per block (aDMRs per Mb)
  hits        <- findOverlaps(blks, dmr_gr)
  blk_dmr_n   <- tabulate(queryHits(hits), nbins = length(blks))
  blk_density <- blk_dmr_n / pmax(width(blks) / 1e6, 0.01)
  n_all       <- length(blks)

  do.call(rbind, lapply(SV_TIER_LEVELS, function(tier) {
    sv_tier_gr <- sv_gr[
      !is.na(mcols(sv_gr)$sv_tier) &
      as.character(mcols(sv_gr)$sv_tier) == tier
    ]
    if (length(sv_tier_gr) == 0) return(NULL)

    sv_block_ids <- unique(na.omit(mcols(sv_tier_gr)$block_id))
    has_sv <- blks$block_id %in% sv_block_ids
    n_sv   <- sum(has_sv)
    if (n_sv == 0L || n_sv >= n_all) return(NULL)

    obs_sv  <- mean(blk_density[has_sv],  na.rm = TRUE)
    obs_wt  <- mean(blk_density[!has_sv], na.rm = TRUE)
    obs_ratio <- obs_sv / pmax(obs_wt, 1e-6)

    # Block-label permutation (shuffle which n_sv blocks are "SV blocks")
    perm_ratios <- replicate(N_PERM_BLOCK, {
      idx   <- sample.int(n_all, n_sv)
      p_sv  <- mean(blk_density[idx],  na.rm = TRUE)
      p_wt  <- mean(blk_density[-idx], na.rm = TRUE)
      p_sv  / pmax(p_wt, 1e-6)
    })
    p_perm <- mean(perm_ratios >= obs_ratio)

    wt <- tryCatch(
      wilcox.test(blk_density[has_sv], blk_density[!has_sv],
                  alternative = "greater", exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )

    data.frame(
      patient_id       = pt,
      sv_tier          = tier,
      n_sv_blocks      = n_sv,
      n_nonsv_blocks   = sum(!has_sv),
      obs_sv_density   = round(obs_sv, 4),
      obs_wt_density   = round(obs_wt, 4),
      enrichment_ratio = round(obs_ratio, 3),
      p_perm           = round(p_perm, 4),
      p_wilcoxon       = round(wt$p.value, 4)
    )
  }))
}))

if (!is.null(layer1_results) && nrow(layer1_results) > 0) {
  layer1_results$p_perm_fdr   <- p.adjust(layer1_results$p_perm,     method = "BH")
  layer1_results$p_wilcox_fdr <- p.adjust(layer1_results$p_wilcoxon, method = "BH")
  layer1_results$sv_tier       <- factor(layer1_results$sv_tier, levels = SV_TIER_LEVELS)

  if (opt$stratify_by == "tier") {
    jt_l1 <- tryCatch({
      jt_df <- layer1_results %>% filter(!is.na(enrichment_ratio))
      clinfun::jonckheere.test(
        jt_df$enrichment_ratio,
        as.integer(jt_df$sv_tier),
        alternative = "decreasing",
        nperm       = N_PERM_JT
      )
    }, error = function(e) { message("JT L1 failed: ", e$message); list(p.value = NA_real_) })
    cat(sprintf("Layer 1 JT (enrichment_ratio across tiers): p = %.4f\n", jt_l1$p.value))
  } else {
    kw_l1 <- tryCatch(
      kruskal.test(enrichment_ratio ~ sv_tier, data = layer1_results),
      error = function(e) { message("KW L1 failed: ", e$message); list(p.value = NA_real_) }
    )
    cat(sprintf("Layer 1 KW (enrichment_ratio across SV types): p = %.4f\n", kw_l1$p.value))
  }

  cat(sprintf("Layer 1 summary (median enrichment_ratio by %s):\n", STRAT_LABEL))
  print(layer1_results %>%
    group_by(sv_tier) %>%
    summarise(
      n_patients   = n(),
      median_ratio = median(enrichment_ratio, na.rm = TRUE),
      pct_sig_perm = mean(p_perm < 0.05, na.rm = TRUE) * 100,
      .groups      = "drop"
    ))
}

} # end if has_tier Layer 1


# 3. Layer 2 — Tier-stratified HP-specific Δβ ==================================
#    pooled KW + JT + Dunn post-hoc (Bonferroni)

cat("\n=== Layer 2: Tier-stratified HP-specific Δβ ===\n")

if (!has_tier) {
  message("Skipping Layer 2: --sv_strat_file required.")
} else {
  layer2_results <- lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr) || length(sv_gr) == 0 ||
          is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    lapply(SV_TIER_LEVELS, function(tier) {
      sv_tier_gr <- sv_gr[
        !is.na(mcols(sv_gr)$sv_tier) &
          as.character(mcols(sv_gr)$sv_tier) == tier
      ]
      if (length(sv_tier_gr) == 0) return(NULL)

      sv_hp_map <- infer_sv_hp_map(mcols(sv_tier_gr)$hp, mcols(sv_tier_gr)$block_id) #nolint
      hp_df <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
      if (is.null(hp_df) || nrow(hp_df) < 3L) return(NULL)

      wt <- tryCatch(
        wilcox.test(hp_df$sv_hp_beta, hp_df$wt_hp_beta,
                    paired = TRUE, alternative = "two.sided", exact = FALSE),
        error = function(e) list(p.value = NA_real_, statistic = NA_real_)
      )
      n  <- nrow(hp_df)
      W  <- as.numeric(wt$statistic)
      rb <- if (!is.na(W)) (2 * W) / (n * (n + 1) / 2) - 1 else NA_real_

      data.frame(
        patient_id       = pt,
        sv_tier          = tier,
        n_admr           = n,
        mean_sv_beta     = round(mean(hp_df$sv_hp_beta), 3),
        mean_wt_beta     = round(mean(hp_df$wt_hp_beta), 3),
        median_abs_delta = round(median(hp_df$hp_delta_abs), 3),
        mean_sv_minus_wt = round(mean(hp_df$sv_minus_wt), 3),
        wilcox_p         = round(wt$p.value, 5),
        rank_biserial    = round(rb, 3)
      )
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()

  if (!is.null(layer2_results) && nrow(layer2_results) > 0) {
    layer2_results$wilcox_fdr <- p.adjust(layer2_results$wilcox_p, method = "BH")
    layer2_results$sv_tier    <- factor(layer2_results$sv_tier, levels = SV_TIER_LEVELS)
  }

  # Pool all aDMRs for KW + JT + Dunn
  all_hp_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr) || length(sv_gr) == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    do.call(rbind, lapply(SV_TIER_LEVELS, function(tier) {
      sv_tier_gr <- sv_gr[
        !is.na(mcols(sv_gr)$sv_tier) &
        as.character(mcols(sv_gr)$sv_tier) == tier
      ]
      if (length(sv_tier_gr) == 0) return(NULL)
      sv_hp_map <- infer_sv_hp_map(mcols(sv_tier_gr)$hp, mcols(sv_tier_gr)$block_id) #nolint
      hp_df <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
      if (is.null(hp_df)) return(NULL)
      hp_df$patient_id  <- pt
      hp_df$sv_tier     <- tier
      hp_df$hp_abs_diff <- hp_df$hp_delta_abs
      hp_df
    }))
  }))

  if (!is.null(all_hp_df) && nrow(all_hp_df) > 0) {
    all_hp_df$sv_tier <- factor(all_hp_df$sv_tier, levels = SV_TIER_LEVELS)
    fwrite(all_hp_df, file.path(opt$outdir, "all_hp_admr_tier.csv.gz"))

# Aggregate to patient × tier medians before KW/JT to avoid pseudoreplication
# (aDMRs within a patient are correlated; raw pooling inflates effective N)
    pt_tier_medians <- all_hp_df %>%
      group_by(patient_id, sv_tier) %>%
      summarise(hp_abs_diff = median(hp_abs_diff, na.rm = TRUE), .groups = "drop")

    kw_l2 <- kruskal.test(hp_abs_diff ~ sv_tier, data = pt_tier_medians)
    cat(sprintf("Layer 2 KW (patient-median hp_abs_diff across tiers): p = %.4f\n", kw_l2$p.value))

    if (opt$stratify_by == "tier") {
      jt_l2 <- tryCatch(
        clinfun::jonckheere.test(
          pt_tier_medians$hp_abs_diff,
          as.integer(pt_tier_medians$sv_tier),
          alternative = "decreasing",
          nperm       = N_PERM_JT
        ),
        error = function(e) { message("JT L2 failed: ", e$message); list(p.value = NA_real_) }
      )
      cat(sprintf("Layer 2 JT (patient-median hp_abs_diff across tiers): p = %.4f\n", jt_l2$p.value))
    }

    tryCatch(
      dunn.test::dunn.test(pt_tier_medians$hp_abs_diff, pt_tier_medians$sv_tier,
                          method = "bonferroni", kw = FALSE),
      error = function(e) message("Dunn L2 failed: ", e$message)
    )

    # Direction concordance test: only meaningful for SV types
    if (opt$stratify_by == "svtype") {
      dir_df <- all_hp_df %>%
        mutate(sv_positive = sv_minus_wt > 0) %>%
        group_by(sv_tier) %>%
        summarise(
          n_total      = n(),
          n_sv_hyper   = sum(sv_positive, na.rm = TRUE),
          pct_sv_hyper = round(mean(sv_positive, na.rm = TRUE) * 100, 1),
          binom_p      = tryCatch(
            binom.test(sum(sv_positive, na.rm = TRUE), n(), p = 0.5)$p.value,
            error = function(e) NA_real_
          ),
          .groups = "drop"
        )
      cat("\nLayer 2 direction concordance (% aDMRs where SV-HP β > WT-HP β = hypermeth on SV allele):\n")
      cat("Expected by mechanism: DEL→>50% hyper, DUP→<50% hyper\n")
      print(dir_df)
    }
  }

  cat(sprintf("Layer 2 summary (median_abs_delta by %s):\n", STRAT_LABEL))
  if (nrow(layer2_results) > 0) {
    print(layer2_results %>%
      group_by(sv_tier) %>%
      summarise(
        n_patients       = n(),
        median_abs_delta = median(median_abs_delta, na.rm = TRUE),
        pct_sig_wilcox   = mean(wilcox_fdr < 0.05, na.rm = TRUE) * 100,
        .groups          = "drop"
      ))
  }
} # end if has_tier Layer 2


# 4. Layer 3 — CRE overlap fraction by sv_tier ================================
#    % aDMRs in SV blocks overlapping promoter / CTCF / enhancer
#    Wilcoxon (TAD_CTCF vs non_boundary) + JT across all tiers

cat("\n=== Layer 3: CRE overlap by sv_tier ===\n")

if (!has_tier) {
  message("Skipping Layer 3: --sv_strat_file required.")
} else {

message("Loading CRE references: ", appendLF= FALSE)
message("build promoters/enhancers via annotatr (hg38) when available
         load CTCF from BED.")

# Build promoter + enhancer annotations via annotatr (preferred)
annotations <- tryCatch({
  annotatr::build_annotations(
    genome = "hg38",
    annotations = c("hg38_genes_promoters", "hg38_enhancers_fantom")
  )
}, error = function(e) { message("annotatr build_annotations failed: ", e$message); NULL })

prom_gr <- GRanges()
enh_gr  <- GRanges()
if (!is.null(annotations) && length(annotations) > 0) {
  ann_md <- mcols(annotations)
  ann_col <- NULL
  if ("annotation" %in% colnames(ann_md)) ann_col <- ann_md$annotation
  else if ("annot.type" %in% colnames(ann_md)) ann_col <- ann_md$`annot.type`
  else if ("annot_type" %in% colnames(ann_md)) ann_col <- ann_md$annot_type
  else ann_col <- rep(NA_character_, length(annotations))

  prom_gr <- annotations[grepl("promoter", tolower(as.character(ann_col)), fixed = FALSE)]
  enh_gr  <- annotations[grepl("enhancer", tolower(as.character(ann_col)), fixed = FALSE)]
}

# Attempt canonical-transcript-based promoter definitions via GENCODE v49 gtf file (cached)
canonical_prom_gr <- readRDS(file.path(getwd(), "canonical_promoters_hg38.gencode_v49.rds"))

# If canonical promoters found, prefer them
if (length(canonical_prom_gr) > 0) {
  prom_gr <- canonical_prom_gr
}

# Fallback: if user supplied explicit BEDs, prefer those (keeps compatibility)
ctcf_gr <- tryCatch({
  dt <- fread(opt$ctcf_bed,
              col.names = c("chr","start","end","name","score","strand",
                            "signalValue","pValue","qValue","peak"))
  GRanges(seqnames = dt$chr, ranges = IRanges(dt$start + 1L, dt$end))
}, error = function(e) { message("Failed to load CTCF BED: ", e$message); GRanges() })

if ((!is.null(opt$enhancer_bed) && file.exists(opt$enhancer_bed)) && length(enh_gr) == 0) {
  enh_gr <- tryCatch({
    dt <- fread(opt$enhancer_bed, skip = 1L,
                col.names = c("chr","start","end","element","segment"))
    GRanges(seqnames = dt$chr, ranges = IRanges(dt$start + 1L, dt$end))
  }, error = function(e) { message("Failed to load enhancer BED: ", e$message); GRanges() })
}

if ((!is.null(opt$promoter_bed) && file.exists(opt$promoter_bed)) && length(prom_gr) == 0) {
  prom_gr <- tryCatch({
    dt <- fread(opt$promoter_bed, skip = 1L,
                col.names = c("chr","start","end","segment","element"))
    GRanges(seqnames = dt$chr, ranges = IRanges(dt$start + 1L, dt$end))
  }, error = function(e) { message("Failed to load promoter BED: ", e$message); GRanges() })
}

message(sprintf("CRE loaded: CTCF=%d, Enhancer=%d, Promoter=%d",
                length(ctcf_gr), length(enh_gr), length(prom_gr)))

layer3_results <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_phased[[pt]]
  dmr_gr <- admr_phased[[pt]]
  if (is.null(sv_gr) || length(sv_gr) == 0 ||
      is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  dmr_ann <- annotate_dmr_cre(dmr_gr, prom_gr, ctcf_gr, enh_gr)

  do.call(rbind, lapply(SV_TIER_LEVELS, function(tier) {
    sv_tier_gr <- sv_gr[
      !is.na(mcols(sv_gr)$sv_tier) &
      as.character(mcols(sv_gr)$sv_tier) == tier
    ]
    if (length(sv_tier_gr) == 0) return(NULL)

    sv_block_ids  <- unique(na.omit(mcols(sv_tier_gr)$block_id))
    dmr_in_sv_blk <- dmr_ann[mcols(dmr_ann)$block_id %in% sv_block_ids]
    if (length(dmr_in_sv_blk) == 0) return(NULL)

    data.frame(
      patient_id   = pt,
      sv_tier      = tier,
      n_admr       = length(dmr_in_sv_blk),
      pct_promoter = round(mean(mcols(dmr_in_sv_blk)$overlaps_promoter) * 100, 1),
      pct_ctcf     = round(mean(mcols(dmr_in_sv_blk)$overlaps_ctcf)     * 100, 1),
      pct_enhancer = round(mean(mcols(dmr_in_sv_blk)$overlaps_enhancer) * 100, 1),
      pct_any_cre  = round(mean(mcols(dmr_in_sv_blk)$overlaps_any_cre)  * 100, 1)
    )
  }))
}))

if (!is.null(layer3_results) && nrow(layer3_results) > 0) {
  layer3_results$sv_tier <- factor(layer3_results$sv_tier, levels = SV_TIER_LEVELS)

  if (opt$stratify_by == "tier") {
    wt_l3 <- tryCatch({
      tctcf <- layer3_results$pct_any_cre[layer3_results$sv_tier == "TAD_CTCF"]
      tnb   <- layer3_results$pct_any_cre[layer3_results$sv_tier == "non_boundary"]
      wilcox.test(tctcf, tnb, alternative = "greater", exact = FALSE)
    }, error = function(e) list(p.value = NA_real_))
    jt_l3 <- tryCatch(
      clinfun::jonckheere.test(
        layer3_results$pct_any_cre,
        as.integer(layer3_results$sv_tier),
        alternative = "decreasing",
        nperm       = N_PERM_JT
      ),
      error = function(e) { message("JT L3 failed: ", e$message); list(p.value = NA_real_) }
    )
    cat(sprintf("Layer 3: TAD_CTCF vs non_boundary pct_any_cre Wilcoxon p = %.4f\n", wt_l3$p.value))
    cat(sprintf("Layer 3: JT pct_any_cre across tiers p = %.4f\n", jt_l3$p.value))
  } else {
    kw_l3 <- tryCatch(
      kruskal.test(pct_any_cre ~ sv_tier, data = layer3_results),
      error = function(e) { message("KW L3 failed: ", e$message); list(p.value = NA_real_) }
    )
    cat(sprintf("Layer 3 KW (pct_any_cre across SV types): p = %.4f\n", kw_l3$p.value))
    tryCatch(
      dunn.test::dunn.test(layer3_results$pct_any_cre, layer3_results$sv_tier,
                           method = "bonferroni", kw = FALSE),
      error = function(e) message("Dunn L3 failed: ", e$message)
    )
  }

  cat(sprintf("Layer 3 summary (median %% CRE overlap by %s):\n", STRAT_LABEL))
  print(layer3_results %>%
    group_by(sv_tier) %>%
    summarise(
      n_patients   = n(),
      pct_promoter = median(pct_promoter, na.rm = TRUE),
      pct_ctcf     = median(pct_ctcf,     na.rm = TRUE),
      pct_enhancer = median(pct_enhancer, na.rm = TRUE),
      pct_any_cre  = median(pct_any_cre,  na.rm = TRUE),
      .groups      = "drop"
    ))
}

} # end if has_tier Layer 3

# 5. Layer 4 — Tier effect size gradient + LME model comparison ===============
#    m_tier: hp_abs_diff ~ sv_tier + (1|patient_id)
#    m_dist: hp_abs_diff ~ log10(nearest_bp_dist + 1) + (1|patient_id)
#    Compare AIC; visualise median_abs_delta ± IQR by tier

cat("\n=== Layer 4: Tier gradient + LME model comparison ===\n")

if (!has_tier) {
  message("Skipping Layer 4: --sv_strat_file required.")
} else {

tier_dist_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_phased[[pt]]
  dmr_gr <- admr_phased[[pt]]
  if (is.null(sv_gr) || length(sv_gr) == 0 ||
      is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  do.call(rbind, lapply(SV_TIER_LEVELS, function(tier) {
    sv_tier_gr <- sv_gr[
      !is.na(mcols(sv_gr)$sv_tier) &
      as.character(mcols(sv_gr)$sv_tier) == tier
    ]
    if (length(sv_tier_gr) == 0) return(NULL)

    sv_hp_map <- infer_sv_hp_map(mcols(sv_tier_gr)$hp, mcols(sv_tier_gr)$block_id)
    hp_df <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
    if (is.null(hp_df) || nrow(hp_df) == 0L) return(NULL)

    # Distance from each aDMR to its nearest SV breakpoint (same-block SVs only)
    dmr_blk   <- dmr_gr[mcols(dmr_gr)$block_id %in% names(sv_hp_map)]
    sv_in_blk <- sv_tier_gr[mcols(sv_tier_gr)$block_id %in% mcols(dmr_blk)$block_id]
    dist_vec  <- if (length(sv_in_blk) > 0L) {
      dn <- distanceToNearest(dmr_blk, sv_in_blk)
      d  <- rep(NA_real_, length(dmr_blk))
      d[queryHits(dn)] <- mcols(dn)$distance
      d
    } else {
      rep(NA_real_, length(dmr_blk))
    }

    n <- nrow(hp_df)
    data.frame(
      patient_id      = pt,
      sv_tier         = tier,
      hp_abs_diff     = hp_df$hp_delta_abs,
      nearest_bp_dist = dist_vec[seq_len(n)]
    )
  }))
}))

if (!is.null(tier_dist_df) && nrow(tier_dist_df) > 0L) {
  tier_dist_df <- tier_dist_df %>% filter(!is.na(hp_abs_diff))
  tier_dist_df$sv_tier    <- factor(tier_dist_df$sv_tier, levels = SV_TIER_LEVELS)
  tier_dist_df$log10_dist <- log10(pmax(tier_dist_df$nearest_bp_dist, 1, na.rm = TRUE) + 1)

  lme_df <- tier_dist_df %>% filter(!is.na(log10_dist))

  m_tier <- tryCatch(
    lme4::lmer(hp_abs_diff ~ sv_tier + (1 | patient_id), data = lme_df, REML = FALSE),
    error = function(e) { message("LME tier failed: ", e$message); NULL }
  )
  m_dist <- tryCatch(
    lme4::lmer(hp_abs_diff ~ log10_dist + (1 | patient_id), data = lme_df, REML = FALSE),
    error = function(e) { message("LME dist failed: ", e$message); NULL }
  )

  aic_tier <- if (!is.null(m_tier)) AIC(m_tier) else NA_real_
  aic_dist <- if (!is.null(m_dist)) AIC(m_dist) else NA_real_

  # Full model tests whether tier retains explanatory power after conditioning on distance
  m_full <- tryCatch(
    lme4::lmer(hp_abs_diff ~ sv_tier + log10_dist + (1 | patient_id), data = lme_df, REML = FALSE),
    error = function(e) { message("LME full failed: ", e$message); NULL }
  )
  aic_full <- if (!is.null(m_full)) AIC(m_full) else NA_real_

  cat(sprintf("Layer 4 LME AIC: m_tier = %.2f | m_dist = %.2f | m_full = %.2f\n",
              aic_tier, aic_dist, aic_full))
  cat(sprintf("  ΔAIC (m_dist − m_tier) = %.2f  [positive = tier better than distance alone]\n",
              aic_dist - aic_tier))
  cat(sprintf("  ΔAIC (m_tier − m_full) = %.2f  [positive = dist adds to tier; negative = tier absorbs dist]\n",
              aic_tier - aic_full))

  layer4_model_comparison <- data.frame(
    model     = c("m_tier: hp_abs_diff ~ sv_tier + (1|patient_id)",
                  "m_dist: hp_abs_diff ~ log10(dist+1) + (1|patient_id)",
                  "m_full: hp_abs_diff ~ sv_tier + log10(dist+1) + (1|patient_id)"),
    AIC           = round(c(aic_tier, aic_dist, aic_full), 2),
    delta_vs_tier = round(c(0, aic_dist - aic_tier, aic_full - aic_tier), 2),
    notes = c(
      "reference (tier only)",
      "negative delta = dist worse than tier",
      "negative delta = full model better than tier alone"
    )
  )

  layer4_gradient <- tier_dist_df %>%
    group_by(sv_tier) %>%
    summarise(
      n_admr           = n(),
      median_abs_delta = median(hp_abs_diff, na.rm = TRUE),
      q25              = quantile(hp_abs_diff, 0.25, na.rm = TRUE),
      q75              = quantile(hp_abs_diff, 0.75, na.rm = TRUE),
      .groups          = "drop"
    )

  cat("Layer 4 gradient (median_abs_delta by tier):\n")
  print(layer4_gradient)
}

} # end if has_tier Layer 4

# 5b. HBV-associated SV analysis (highest-priority tier; reported separately) ==
#     HP-specific Δβ for SVs with sv_tier == "HBV_associated"

cat("\n=== HBV-associated SVs: HP-specific Δβ ===\n")

hbv_hp_results <- data.frame()
hbv_hp_df      <- data.frame()

if (has_tier) {
  hbv_hp_results <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr) || length(sv_gr) == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    hbv_gr <- sv_gr[
      !is.na(mcols(sv_gr)$sv_tier) &
      as.character(mcols(sv_gr)$sv_tier) == "HBV_associated"
    ]
    if (length(hbv_gr) == 0) return(NULL)

    sv_hp_map <- infer_sv_hp_map(mcols(hbv_gr)$hp, mcols(hbv_gr)$block_id)
    hp_df <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
    if (is.null(hp_df) || nrow(hp_df) < 3L) {
      return(data.frame(
        patient_id    = pt,
        n_hbv_sv      = length(hbv_gr),
        n_admr        = 0L,
        wilcox_p      = NA_real_,
        rank_biserial = NA_real_,
        mean_sv_minus_wt = NA_real_
      ))
    }

    wt <- tryCatch(
      wilcox.test(hp_df$sv_hp_beta, hp_df$wt_hp_beta,
                  paired = TRUE, alternative = "two.sided", exact = FALSE),
      error = function(e) list(p.value = NA_real_, statistic = NA_real_)
    )
    n  <- nrow(hp_df)
    W  <- as.numeric(wt$statistic)
    rb <- if (!is.na(W)) (2 * W) / (n * (n + 1) / 2) - 1 else NA_real_

    data.frame(
      patient_id       = pt,
      n_hbv_sv         = length(hbv_gr),
      n_admr           = n,
      mean_sv_beta     = round(mean(hp_df$sv_hp_beta), 3),
      mean_wt_beta     = round(mean(hp_df$wt_hp_beta), 3),
      mean_sv_minus_wt = round(mean(hp_df$sv_minus_wt), 3),
      median_abs_delta = round(median(hp_df$hp_delta_abs), 3),
      wilcox_p         = round(wt$p.value, 5),
      rank_biserial    = round(rb, 3)
    )
  }))

  if (!is.null(hbv_hp_results) && nrow(hbv_hp_results) > 0) {
    hbv_hp_results$wilcox_fdr <- p.adjust(hbv_hp_results$wilcox_p, method = "BH")
    cat(sprintf("HBV-associated SVs found in %d / %d patients\n",
                sum(!is.na(hbv_hp_results$n_admr) & hbv_hp_results$n_admr > 0),
                length(PATIENT_IDS)))
    cat(sprintf("Total HBV SV count: %d\n", sum(hbv_hp_results$n_hbv_sv, na.rm = TRUE)))
    print(hbv_hp_results)

    if (sum(hbv_hp_results$n_hbv_sv, na.rm = TRUE) == 0) {
      message("WARNING: HBV SV count = 0 across all patients — check HBV contig naming in VCF.")
    }

    fwrite(hbv_hp_results, file.path(opt$outdir, "layer_hbv_hp_delta.csv"))
  }
}

# =============================================================================
# 5c. Methylation-competent SV classification
#     An SV is "methylation-competent" if ≥1 aDMR falls in the same phase block.
#     Logistic mixed model identifies which SV features predict competence.
# =============================================================================

cat("\n=== Section 5c: Methylation-competent SV classification ===\n")

sv_competence_df   <- data.frame()
sv_competence_summ <- data.frame()

if (has_tier) {
  sv_competence_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr)  || length(sv_gr)  == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    dmr_block_ids <- unique(na.omit(mcols(dmr_gr)$block_id))
    sv_blk        <- as.character(mcols(sv_gr)$block_id)

    data.frame(
      patient_id     = pt,
      bp_id          = names(sv_gr),
      sv_tier        = as.character(mcols(sv_gr)$sv_tier),
      sv_arch        = as.character(mcols(sv_gr)$sv_arch),
      sv_bio         = as.character(mcols(sv_gr)$sv_bio),
      arch_bio       = as.character(mcols(sv_gr)$arch_bio),
      sv_type        = as.character(mcols(sv_gr)$sv_type),
      svlen          = abs(as.numeric(mcols(sv_gr)$svlen)),
      meth_competent = (!is.na(sv_blk)) & (sv_blk %in% dmr_block_ids)
    )
  }))

  if (nrow(sv_competence_df) > 0) {
    sv_competence_df$sv_arch <- factor(sv_competence_df$sv_arch, levels = c("non_boundary","boundary"))
    sv_competence_df$sv_bio  <- factor(sv_competence_df$sv_bio,  levels = c("balanced","CN_changing"))
    sv_competence_df$sv_tier <- factor(sv_competence_df$sv_tier, levels = SV_TIER_LEVELS)
    sv_competence_df$log10_svlen <- log10(pmax(sv_competence_df$svlen, 1L))

    sv_competence_summ <- sv_competence_df %>%
      group_by(sv_tier, sv_arch, sv_bio) %>%
      summarise(
        n_sv          = n(),
        n_competent   = sum(meth_competent),
        pct_competent = round(mean(meth_competent) * 100, 1),
        .groups       = "drop"
      )
    cat("Methylation-competent SV fraction by tier:\n")
    print(sv_competence_summ)

    lr_competent <- tryCatch(
      lme4::glmer(meth_competent ~ sv_arch + sv_bio + log10_svlen + (1 | patient_id),
                  data   = sv_competence_df,
                  family = binomial,
                  control = glmerControl(optimizer = "bobyqa")),
      error = function(e) { message("glmer competence failed: ", e$message); NULL }
    )
    if (!is.null(lr_competent)) {
      cat("\nLogistic LME: meth_competent ~ sv_arch + sv_bio + log10_svlen + (1|patient_id)\n")
      coef_df <- as.data.frame(summary(lr_competent)$coefficients)
      coef_df$OR <- round(exp(coef_df$Estimate), 3)
      print(coef_df)
    }

    fwrite(sv_competence_df,  file.path(opt$outdir, "sv_meth_competence.csv.gz"))
    fwrite(sv_competence_summ, file.path(opt$outdir, "sv_meth_competence_summary.csv"))
  }
}


# =============================================================================
# 5d. Distance-bin × SV biology (CN_changing vs balanced) interaction
#     Bins: proximal (0–10 kb), near (10–50 kb), far (50–500 kb), distal (>500 kb)
#     Tests whether the distance-decay slope differs by SV biology class.
# =============================================================================

cat("\n=== Section 5d: Distance-bin × SV biology interaction ===\n")

dist_bio_df   <- data.frame()
dist_bin_summ <- data.frame()

if (has_tier) {
  dist_bio_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr)  || length(sv_gr)  == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    do.call(rbind, lapply(c("CN_changing", "balanced"), function(bio) {
      sv_bio_gr <- sv_gr[
        !is.na(mcols(sv_gr)$sv_bio) &
        as.character(mcols(sv_gr)$sv_bio) == bio
      ]
      if (length(sv_bio_gr) == 0) return(NULL)

      sv_hp_map <- infer_sv_hp_map(mcols(sv_bio_gr)$hp, mcols(sv_bio_gr)$block_id)
      hp_df     <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
      if (is.null(hp_df) || nrow(hp_df) == 0) return(NULL)

      dmr_blk   <- dmr_gr[mcols(dmr_gr)$block_id %in% names(sv_hp_map)]
      sv_in_blk <- sv_bio_gr[mcols(sv_bio_gr)$block_id %in% mcols(dmr_blk)$block_id]
      dist_vec  <- rep(NA_real_, length(dmr_blk))
      if (length(sv_in_blk) > 0) {
        dn <- distanceToNearest(dmr_blk, sv_in_blk)
        dist_vec[queryHits(dn)] <- mcols(dn)$distance
      }

      n <- nrow(hp_df)
      data.frame(
        patient_id      = pt,
        sv_bio          = bio,
        hp_abs_diff     = hp_df$hp_delta_abs,
        sv_minus_wt     = hp_df$sv_minus_wt,
        nearest_bp_dist = dist_vec[seq_len(n)]
      )
    }))
  }))

  if (nrow(dist_bio_df) > 0) {
    dist_bio_df <- dist_bio_df %>% filter(!is.na(nearest_bp_dist), !is.na(hp_abs_diff))
    dist_bio_df$log10_dist <- log10(pmax(dist_bio_df$nearest_bp_dist, 1L) + 1L)
    dist_bio_df$dist_bin   <- cut(
      dist_bio_df$nearest_bp_dist,
      breaks = c(0, 10e3, 50e3, 500e3, Inf),
      labels = c("0-10kb", "10-50kb", "50-500kb", ">500kb"),
      include.lowest = TRUE
    )
    dist_bio_df$sv_bio <- factor(dist_bio_df$sv_bio, levels = c("CN_changing", "balanced"))

    m_dist_bio_int  <- tryCatch(
      lme4::lmer(hp_abs_diff ~ log10_dist * sv_bio + (1 | patient_id),
                 data = dist_bio_df, REML = FALSE),
      error = function(e) { message("LME dist×bio interaction failed: ", e$message); NULL }
    )
    m_dist_bio_main <- tryCatch(
      lme4::lmer(hp_abs_diff ~ log10_dist + sv_bio + (1 | patient_id),
                 data = dist_bio_df, REML = FALSE),
      error = function(e) NULL
    )
    if (!is.null(m_dist_bio_int) && !is.null(m_dist_bio_main)) {
      cat(sprintf("LME dist×bio: AIC interaction=%.2f  additive=%.2f  ΔAIC=%.2f\n",
                  AIC(m_dist_bio_int), AIC(m_dist_bio_main),
                  AIC(m_dist_bio_main) - AIC(m_dist_bio_int)))
    }

    dist_bin_summ <- dist_bio_df %>%
      group_by(dist_bin, sv_bio) %>%
      summarise(
        n               = n(),
        median_abs_diff = median(hp_abs_diff, na.rm = TRUE),
        q25             = quantile(hp_abs_diff, 0.25, na.rm = TRUE),
        q75             = quantile(hp_abs_diff, 0.75, na.rm = TRUE),
        .groups         = "drop"
      )
    cat("\nDistance-bin × SV biology summary:\n")
    print(dist_bin_summ)

    fwrite(dist_bio_df  %>% dplyr::select(-log10_dist),
           file.path(opt$outdir, "dist_bin_svbio_hp_delta.csv.gz"))
    fwrite(dist_bin_summ, file.path(opt$outdir, "dist_bin_svbio_summary.csv"))
  }
}


# =============================================================================
# 5e. Signed Δβ direction: are SVs preferentially hyper- or hypo-methylating?
#     Stratified by SV type (DEL/DUP/INV/TRA/INS).
#     Expected (dosage model): DEL → SV-HP hypermethylation (>50%); DUP → variable.
# =============================================================================

cat("\n=== Section 5e: Signed Δβ direction by SV type ===\n")

dir_svtype_df   <- data.frame()
dir_svtype_summ <- data.frame()

if (has_tier) {
  dir_svtype_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr)  || length(sv_gr)  == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    do.call(rbind, lapply(c("DEL","DUP","INV","TRA","INS"), function(svt) {
      sv_t_gr <- sv_gr[!is.na(mcols(sv_gr)$sv_type) & as.character(mcols(sv_gr)$sv_type) == svt]
      if (length(sv_t_gr) == 0) return(NULL)

      sv_hp_map <- infer_sv_hp_map(mcols(sv_t_gr)$hp, mcols(sv_t_gr)$block_id)
      hp_df     <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
      if (is.null(hp_df) || nrow(hp_df) < 3L) return(NULL)

      n_valid <- sum(!is.na(hp_df$sv_minus_wt))
      n_hyper <- sum(hp_df$sv_minus_wt > 0, na.rm = TRUE)
      binom_p <- tryCatch(
        binom.test(n_hyper, n_valid, p = 0.5)$p.value,
        error = function(e) NA_real_
      )
      data.frame(
        patient_id       = pt,
        sv_type          = svt,
        n_admr           = n_valid,
        n_hyper          = n_hyper,
        pct_hyper        = round(n_hyper / n_valid * 100, 1),
        binom_p          = round(binom_p, 5),
        mean_sv_minus_wt = round(mean(hp_df$sv_minus_wt, na.rm = TRUE), 3)
      )
    }))
  }))

  if (nrow(dir_svtype_df) > 0) {
    dir_svtype_df$binom_fdr <- p.adjust(dir_svtype_df$binom_p, method = "BH")
    dir_svtype_df$sv_type   <- factor(dir_svtype_df$sv_type, levels = SV_TYPE_LEVELS)

    dir_svtype_summ <- dir_svtype_df %>%
      group_by(sv_type) %>%
      summarise(
        n_patients       = n(),
        total_admr       = sum(n_admr),
        pooled_n_hyper   = sum(n_hyper),
        pooled_pct_hyper = round(sum(n_hyper) / sum(n_admr) * 100, 1),
        # Meta-analysis of per-patient binomial p-values (Fisher's combined probability).
        # Pooling raw counts treats patient aDMRs as independent observations, inflating n.
        meta_binom_p = tryCatch({
          ps <- binom_p[!is.na(binom_p) & binom_p > 0 & binom_p <= 1]
          if (length(ps) < 2L) NA_real_
          else pchisq(-2 * sum(log(ps)), df = 2L * length(ps), lower.tail = FALSE)
        }, error = function(e) NA_real_),
        .groups = "drop"
      )
    dir_svtype_summ$meta_binom_fdr <- p.adjust(dir_svtype_summ$meta_binom_p, method = "BH")

    cat("Signed Δβ direction by SV type (>50% = SV-allele hypermethylation):\n")
    cat("Expected (dosage model): DEL → >50%; DUP → <50%\n")
    print(dir_svtype_summ)

    fwrite(dir_svtype_df,   file.path(opt$outdir, "direction_svtype_hp_delta.csv"))
    fwrite(dir_svtype_summ, file.path(opt$outdir, "direction_svtype_summary.csv"))
  }
}


# =============================================================================
# 5f. Revised Layer 2: 2D collapsed tier (sv_arch × sv_bio)
#     Tests whether boundary disruption adds explanatory power *within* each
#     SV biology class (CN_changing or balanced).
#     Key comparison: boundary_CN_changing vs non_boundary_CN_changing.
# =============================================================================

cat("\n=== Section 5f: Revised Layer 2 — 2D collapsed tier (arch × sv_bio) ===\n")

arch_bio_hp_df  <- data.frame()
arch_bio_result <- data.frame()

if (has_tier) {
  arch_bio_hp_df <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
    sv_gr  <- sv_phased[[pt]]
    dmr_gr <- admr_phased[[pt]]
    if (is.null(sv_gr)  || length(sv_gr)  == 0 ||
        is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

    do.call(rbind, lapply(ARCH_BIO_LEVELS, function(ab) {
      sv_ab_gr <- sv_gr[
        !is.na(mcols(sv_gr)$arch_bio) &
        as.character(mcols(sv_gr)$arch_bio) == ab
      ]
      if (length(sv_ab_gr) == 0) return(NULL)

      sv_hp_map <- infer_sv_hp_map(mcols(sv_ab_gr)$hp, mcols(sv_ab_gr)$block_id)
      hp_df     <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
      if (is.null(hp_df) || nrow(hp_df) == 0) return(NULL)

      hp_df$patient_id  <- pt
      hp_df$arch_bio    <- ab
      hp_df$hp_abs_diff <- hp_df$hp_delta_abs
      hp_df
    }))
  }))

  if (nrow(arch_bio_hp_df) > 0) {
    arch_bio_hp_df$arch_bio <- factor(arch_bio_hp_df$arch_bio, levels = ARCH_BIO_LEVELS)

    pt_ab_med <- arch_bio_hp_df %>%
      group_by(patient_id, arch_bio) %>%
      summarise(hp_abs_diff = median(hp_abs_diff, na.rm = TRUE), .groups = "drop")

    kw_ab <- tryCatch(
      kruskal.test(hp_abs_diff ~ arch_bio, data = pt_ab_med),
      error = function(e) list(p.value = NA_real_)
    )
    cat(sprintf("2D-tier KW (patient-median hp_abs_diff across arch×bio): p = %.4f\n", kw_ab$p.value))

    # Primary contrast: does boundary annotation add to CN-changing effect?
    bc <- pt_ab_med$hp_abs_diff[pt_ab_med$arch_bio == "boundary_CN_changing"]
    nc <- pt_ab_med$hp_abs_diff[pt_ab_med$arch_bio == "non_boundary_CN_changing"]
    if (length(bc) >= 3L && length(nc) >= 3L) {
      wt_cn_arch <- wilcox.test(bc, nc, alternative = "greater", exact = FALSE)
      cat(sprintf("Boundary vs non_boundary (CN_changing only): Wilcoxon p = %.4f\n",
                  wt_cn_arch$p.value))
    }

    # Secondary contrast: boundary balanced vs non_boundary balanced
    bb <- pt_ab_med$hp_abs_diff[pt_ab_med$arch_bio == "boundary_balanced"]
    nb <- pt_ab_med$hp_abs_diff[pt_ab_med$arch_bio == "non_boundary_balanced"]
    if (length(bb) >= 3L && length(nb) >= 3L) {
      wt_bal_arch <- wilcox.test(bb, nb, alternative = "greater", exact = FALSE)
      cat(sprintf("Boundary vs non_boundary (balanced only): Wilcoxon p = %.4f\n",
                  wt_bal_arch$p.value))
    }

    arch_bio_result <- arch_bio_hp_df %>%
      group_by(arch_bio) %>%
      summarise(
        n_admr           = n(),
        median_abs_delta = median(hp_abs_diff, na.rm = TRUE),
        q25              = quantile(hp_abs_diff, 0.25, na.rm = TRUE),
        q75              = quantile(hp_abs_diff, 0.75, na.rm = TRUE),
        .groups          = "drop"
      )
    cat("\n2D tier gradient (median |Δβ|):\n")
    print(arch_bio_result)

    fwrite(arch_bio_result, file.path(opt$outdir, "layer2_2dtier_gradient.csv"))
  }
}


# --- Visualization: Sections 5c–5f -------------------------------------------

## Fig 5c: % methylation-competent SVs by tier (bar chart)
if (nrow(sv_competence_summ) > 0) {
  p5c_df <- sv_competence_summ %>%
    mutate(sv_tier = factor(sv_tier, levels = SV_TIER_LEVELS))

  p5c <- ggplot(p5c_df, aes(x = sv_tier, y = pct_competent, fill = sv_tier)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_text(aes(label = sprintf("%d/%d", n_competent, n_sv),
                  y = pct_competent + 1.5), size = 3, color = "grey30") +
    scale_fill_manual(values = SV_TIER_COLORS, guide = "none") +
    scale_x_discrete(limits = SV_TIER_LEVELS) +
    labs(
      title    = "Section 5c. Methylation-competent SVs by tier",
      subtitle = "SVs with ≥1 aDMR in the same phase block",
      x        = "SV tier",
      y        = "% methylation-competent SVs",
      caption  = "Labels: n competent / n total SVs"
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(opt$outdir, "fig5c_meth_competent_by_tier.png"),
         p5c, width = 8, height = 5, device = png, dpi = 300)
  saveRDS(p5c, file.path(opt$outdir, "fig5c_meth_competent_by_tier.rds"))
}

## Fig 5d: Distance-decay × SV biology (line + ribbon)
if (nrow(dist_bin_summ) > 0) {
  p5d <- ggplot(dist_bin_summ,
                aes(x = dist_bin, y = median_abs_diff,
                    color = sv_bio, group = sv_bio)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = sv_bio), alpha = 0.15, color = NA) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(CN_changing = "#E24B4A", balanced = "#3B8BD4"),
                       name = "SV biology") +
    scale_fill_manual(values  = c(CN_changing = "#E24B4A", balanced = "#3B8BD4"), guide = "none") +
    labs(
      title    = "Section 5d. HP |Δβ| distance decay by SV biology class",
      subtitle = "CN_changing = DEL+DUP; balanced = INV+TRA+INS | ribbon = IQR",
      x        = "Distance to nearest SV breakpoint (bin)",
      y        = "Median |SV-HP β − WT-HP β|"
    ) +
    theme_hcc

  ggsave(file.path(opt$outdir, "fig5d_dist_decay_by_svbio.png"),
         p5d, width = 8, height = 5, device = png, dpi = 300)
  saveRDS(p5d, file.path(opt$outdir, "fig5d_dist_decay_by_svbio.rds"))
}

## Fig 5e: Signed Δβ direction by SV type (lollipop / diverging bar vs 50%)
if (nrow(dir_svtype_summ) > 0) {
  p5e_df <- dir_svtype_summ %>%
    filter(!is.na(pooled_pct_hyper)) %>%
    mutate(
      deviation = pooled_pct_hyper - 50,
      sig_label = ifelse(meta_binom_fdr < 0.05, "*", "")
    )

  p5e <- ggplot(p5e_df, aes(x = sv_type, y = deviation, fill = sv_type)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sprintf("%.0f%%%s", pooled_pct_hyper, sig_label),
                  y = deviation + sign(deviation) * 1.5),
              size = 3.2, color = "grey20") +
    scale_fill_manual(values = SV_TYPE_COLORS, guide = "none") +
    labs(
      title    = "Section 5e. SV-allele methylation directionality by SV type",
      subtitle = "Deviation from 50% (null = random); * = pooled binomial FDR < 0.05",
      x        = "SV type",
      y        = "% SV-HP hypermethylated − 50%",
      caption  = "Expected (dosage model): DEL → positive (hyper); DUP → negative (hypo)"
    ) +
    theme_hcc

  ggsave(file.path(opt$outdir, "fig5e_direction_by_svtype.png"),
         p5e, width = 7, height = 5, device = png, dpi = 300)
  saveRDS(p5e, file.path(opt$outdir, "fig5e_direction_by_svtype.rds"))
}

## Fig 5f: 2D tier violin (arch × sv_bio)
if (nrow(arch_bio_hp_df) > 0) {
  p5f <- ggplot(arch_bio_hp_df,
                aes(x = arch_bio, y = hp_abs_diff, fill = arch_bio)) +
    geom_violin(alpha = 0.7, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.4,
                 fill = "white", color = "grey30", alpha = 0.85) +
    scale_fill_manual(values = ARCH_BIO_COLORS, guide = "none") +
    scale_x_discrete(limits = ARCH_BIO_LEVELS) +
    labs(
      title    = "Section 5f. HP |Δβ| by 2D collapsed tier (architecture × SV biology)",
      subtitle = "Key contrast: boundary_CN_changing vs non_boundary_CN_changing",
      x        = "2D tier (sv_arch × sv_bio)",
      y        = "|SV-HP β − WT-HP β|",
      caption  = "boundary = TAD_CTCF|TAD_only|CTCF_only; CN_changing = DEL|DUP"
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  ggsave(file.path(opt$outdir, "fig5f_2dtier_hp_delta.png"),
         p5f, width = 9, height = 5.5, device = png, dpi = 300)
  saveRDS(p5f, file.path(opt$outdir, "fig5f_2dtier_hp_delta.rds"))
}


# Visualization ===============================================================
## 1. Figure 1 — Layer 1: enrichment_ratio by tier (bar + JT p-value) =========
if (has_tier && nrow(layer1_results) > 0) {
  tier_summ_l1 <- layer1_results %>%
    group_by(sv_tier) %>%
    summarise(
      med_ratio = median(enrichment_ratio, na.rm = TRUE),
      q25       = quantile(enrichment_ratio, 0.25, na.rm = TRUE),
      q75       = quantile(enrichment_ratio, 0.75, na.rm = TRUE),
      n_sig     = sum(p_perm_fdr < 0.05, na.rm = TRUE),
      n_total   = n(),
      .groups   = "drop"
    )

  p1 <- ggplot(tier_summ_l1, aes(x = sv_tier, y = med_ratio, fill = sv_tier)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, color = "grey40") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%d/%d", n_sig, n_total),
                  y     = pmax(q75, med_ratio) + 0.05),
              size = 3.2, color = "grey30") +
    scale_fill_manual(values = SV_TIER_COLORS, guide = "none") +
    scale_x_discrete(limits = SV_TIER_LEVELS) +
    labs(
      title    = sprintf("Layer 1. SV-block aDMR density enrichment by %s", STRAT_LABEL),
      subtitle = sprintf("Block-label permutation n=%d | IQR error bar | %s p = %.4f",
                         N_PERM_BLOCK,
                         if (opt$stratify_by == "tier") "JT" else "KW",
                         if (opt$stratify_by == "tier") jt_l1$p.value else kw_l1$p.value),
      x        = STRAT_LABEL,
      y        = "Enrichment ratio (SV-block / non-SV-block aDMR density)",
      caption  = "Width-corrected density (aDMRs / block Mb); dashed line = ratio 1; labels = n sig / n patients"
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  saveRDS(p1, file.path(opt$outdir, sprintf("hap_fig1_block_enrichment_by_%s.rds", opt$stratify_by)))
  ggsave(file.path(opt$outdir, sprintf("hap_fig1_block_enrichment_by_%s.png", opt$stratify_by)),
         p1, width = 9, height = 5.5, device = png, dpi = 300)
}

## 2. Figure 2 — Layer 2: hp_abs_diff by tier (violin + boxplot) ===============

if (has_tier && nrow(all_hp_df) > 0) {
  p2 <- ggplot(all_hp_df, aes(x = sv_tier, y = hp_abs_diff, fill = sv_tier)) +
    geom_violin(alpha = 0.7, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.4,
                 fill = "white", color = "grey30", alpha = 0.85) +
    scale_fill_manual(values = SV_TIER_COLORS, guide = "none") +
    scale_x_discrete(limits = SV_TIER_LEVELS) +
    labs(
      title    = sprintf("Layer 2. HP-specific |Δβ| by %s (all patients pooled)", STRAT_LABEL),
      subtitle = sprintf("KW p = %.4f%s",
                         kw_l2$p.value,
                         if (opt$stratify_by == "tier")
                           sprintf(" | JT p = %.4f (decreasing trend across tier hierarchy)", jt_l2$p.value)
                         else " | direction concordance test in console"),
      x        = STRAT_LABEL,
      y        = "|SV-HP β − WT-HP β|",
      caption  = sprintf("Per-patient paired Wilcoxon FDR in layer2_hp_by_%s.csv; Dunn post-hoc in console",
                         opt$stratify_by)
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  # Save figure
  saveRDS(p2, file.path(opt$outdir, sprintf("hap_fig2_hp_delta_by_%s.rds", opt$stratify_by)))
  ggsave(file.path(opt$outdir, sprintf("hap_fig2_hp_delta_by_%s.png", opt$stratify_by)),
         p2, width = 9, height = 5.5, device = png, dpi = 300)
}

## 3. Figure 3 — Layer 3: CRE overlap % by tier (grouped bar) ================
if (has_tier && nrow(layer3_results) > 0) {
  cre_long <- layer3_results %>%
    group_by(sv_tier) %>%
    summarise(
      Promoter = median(pct_promoter, na.rm = TRUE),
      CTCF     = median(pct_ctcf,     na.rm = TRUE),
      Enhancer = median(pct_enhancer, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    pivot_longer(c(Promoter, CTCF, Enhancer),
                 names_to = "cre_type", values_to = "pct")

  p3 <- ggplot(cre_long, aes(x = sv_tier, y = pct, fill = cre_type)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    scale_fill_manual(
      values = c(Promoter = "#E24B4A", CTCF = "#7F77DD", Enhancer = "#1D9E75"),
      name   = "CRE type"
    ) +
    scale_x_discrete(limits = SV_TIER_LEVELS) +
    labs(
      title    = sprintf("Layer 3. aDMR–CRE overlap fraction by %s", STRAT_LABEL),
      subtitle = if (opt$stratify_by == "tier")
        sprintf("TAD_CTCF vs non_boundary Wilcoxon p = %.4f | JT p = %.4f",
                wt_l3$p.value, jt_l3$p.value)
      else
        sprintf("KW p = %.4f | Dunn post-hoc in console", kw_l3$p.value),
      x        = STRAT_LABEL,
      y        = "Median % aDMRs overlapping CRE",
      caption  = "Enhancer: GeneHancer; CTCF: HepG2 ENCFF543WTP; Promoter: hg38_genes_promoters"
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  saveRDS(p3, file.path(opt$outdir, sprintf("hap_fig3_cre_overlap_by_%s.rds", opt$stratify_by)))
  ggsave(file.path(opt$outdir, sprintf("hap_fig3_cre_overlap_by_%s.png", opt$stratify_by)),
         p3, width = 9, height = 5.5, device = png, dpi = 300)
}

## 4. Figure 4 — Layer 4: tier gradient point + errorbar; AIC annotation ======
if (has_tier && nrow(layer4_gradient) > 0) {
  delta_aic <- if (!is.na(aic_tier) && !is.na(aic_dist)) round(aic_dist - aic_tier, 1) else NA_real_
  aic_subtitle <- if (!is.na(delta_aic)) {
    sprintf("Tier model ΔAIC = +%.1f vs distance model (positive = tier better fits hp_abs_diff)", delta_aic)
  } else {
    "LME comparison not available"
  }

  p4 <- ggplot(layer4_gradient, aes(x = sv_tier, y = median_abs_delta, color = sv_tier)) +
    geom_point(size = 4.5) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.25, linewidth = 1.1) +
    scale_color_manual(values = SV_TIER_COLORS, guide = "none") +
    scale_x_discrete(limits = SV_TIER_LEVELS) +
    labs(
      title    = "Layer 4. HP-specific |Δβ| gradient across SV tier hierarchy",
      subtitle = aic_subtitle,
      x        = "SV tier (high → low disruption)",
      y        = "Median |SV-HP β − WT-HP β| (IQR error bar)",
      caption  = "LME: hp_abs_diff ~ sv_tier + (1|patient_id) vs ~ log10(dist+1) + (1|patient_id); REML=FALSE"
    ) +
    theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  saveRDS(p4, file.path(opt$outdir, sprintf("hap_fig4_tier_gradient_%s.rds", opt$stratify_by)))
  ggsave(file.path(opt$outdir, sprintf("hap_fig4_tier_gradient_%s.png", opt$stratify_by)),
         p4, width = 8, height = 5, device = png, dpi = 300)
}

## 5. Sup1 — SV-admr distance distribution  =======================
# Distance distribution: aDMR midpoint → nearest SV breakpoint

### A. Per-patient nearest-SV-bp distance for every aDMR ===========
# For each aDMR, find: (a) nearest ANY SV bp, (b) nearest SV bp per tier
# (c) nearest SV bp per type.
# Compute (a) globally and record the tier/type of the nearest SV.
message("Computing aDMR → nearest-SV-bp distances …")

dist_df <- rbindlist(lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_phased[[pt]]
  dmr_gr <- admr_phased[[pt]]

  if (length(sv_gr) == 0 || length(dmr_gr) == 0) return(NULL)

  hits <- distanceToNearest(dmr_gr, sv_gr, ignore.strand = TRUE)

  if (length(hits) == 0) return(NULL)

  data.table(
    patient_id   = pt,
    dist_bp      = mcols(hits)$distance,
    nearest_tier = as.character(mcols(sv_gr[subjectHits(hits)])$sv_tier),
    nearest_type = as.character(mcols(sv_gr[subjectHits(hits)])$sv_type),
    hp_delta     = mcols(dmr_gr[queryHits(hits)])$hp_delta
  )
}))

message(sprintf("Total aDMR–SV pairs: %d", nrow(dist_df)))

cat("\n=== Global distance summary ===\n")
cat(sprintf("  N aDMRs with a nearest SV : %d\n", nrow(dist_df)))
cat(sprintf("  Median distance           : %s bp\n",
            format(median(dist_df$dist_bp), big.mark = ",")))
cat(sprintf("  Mean distance             : %s bp\n",
            format(round(mean(dist_df$dist_bp)), big.mark = ",")))
cat(sprintf("  %% within 50 kb           : %.1f%%\n",
            mean(dist_df$dist_bp <= 50000) * 100))
cat(sprintf("  %% within 100 kb          : %.1f%%\n",
            mean(dist_df$dist_bp <= 100000) * 100))
cat(sprintf("  %% within 500 kb          : %.1f%%\n",
            mean(dist_df$dist_bp <= 500000) * 100))
cat(sprintf("  %% > 1 Mb                 : %.1f%%\n",
            mean(dist_df$dist_bp > 1e6) * 100))

cat("\n=== Median distance by nearest SV tier ===\n")
print(dist_df %>%
  filter(!is.na(nearest_tier)) %>%
  group_by(nearest_tier) %>%
  summarise(
    n          = n(),
    median_bp  = median(dist_bp),
    pct_50kb   = round(mean(dist_bp <= 50000) * 100, 1),
    pct_100kb  = round(mean(dist_bp <= 100000) * 100, 1),
    .groups    = "drop"
  ) %>%
  arrange(median_bp))

cat("\n=== Median distance by nearest SV type ===\n")
print(dist_df %>%
  filter(!is.na(nearest_type)) %>%
  group_by(nearest_type) %>%
  summarise(
    n          = n(),
    median_bp  = median(dist_bp),
    pct_50kb   = round(mean(dist_bp <= 50000) * 100, 1),
    pct_100kb  = round(mean(dist_bp <= 100000) * 100, 1),
    .groups    = "drop"
  ) %>%
  arrange(median_bp))

### distance vs |HP Δβ| Spearman correlation ==================
rho_all <- cor(log10(dist_df$dist_bp + 1), dist_df$hp_delta,
               method = "spearman", use = "complete.obs")
cat(sprintf("\nSpearman ρ (log10 dist, |HP Δβ|) all aDMRs: %.4f\n", rho_all))

# Per-tier / per-type Spearman ρ 
cat("\n=== Spearman ρ (log10 dist vs |HP Δβ|) by nearest SV tier ===\n")
tier_rho <- dist_df %>%
  filter(!is.na(nearest_tier)) %>%
  group_by(nearest_tier) %>%
  summarise(
    n   = n(),
    rho = cor(log10(dist_bp + 1), hp_delta, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )
print(tier_rho)

cat("\n=== Spearman ρ (log10 dist vs |HP Δβ|) by nearest SV type ===\n")
type_rho <- dist_df %>%
  filter(!is.na(nearest_type)) %>%
  group_by(nearest_type) %>%
  summarise(
    n   = n(),
    rho = cor(log10(dist_bp + 1), hp_delta, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )
print(type_rho)

# Plots
dist_plot <- dist_df %>%
  mutate(log10_dist = log10(dist_bp + 1))

### 5a. Global density on log10 scale ============================
p_global <- ggplot(dist_plot, aes(x = log10_dist)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, fill = "#3B8BD4", alpha = 0.7, color = NA) +
  geom_density(color = "#BA7517", linewidth = 0.9) +
  geom_vline(xintercept = log10(c(50e3, 100e3, 500e3)),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = log10(50e3) + 0.05, y = Inf, label = "50 kb",
           angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
  annotate("text", x = log10(100e3) + 0.05, y = Inf, label = "100 kb",
           angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
  annotate("text", x = log10(500e3) + 0.05, y = Inf, label = "500 kb",
           angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
  scale_x_continuous(
    name   = "Distance to nearest SV breakpoint (log₁₀ bp)",
    breaks = 0:8,
    labels = c("0", "10", "100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  labs(
    title    = "aDMR → nearest SV breakpoint distance (all patients pooled)",
    subtitle = sprintf("N = %s aDMRs | median = %s bp | %%.within 50 kb = %.1f%%",
                       format(nrow(dist_plot), big.mark = ","),
                       format(median(dist_plot$dist_bp), big.mark = ","),
                       mean(dist_plot$dist_bp <= 50000) * 100),
    y        = "Density"
  ) +
  theme_hcc

### 5b. Density by SV tier (facet) =========================
tier_levels_present <- intersect(SV_TIER_LEVELS, unique(dist_plot$nearest_tier))
p_tier <- dist_plot %>%
  filter(nearest_tier %in% tier_levels_present) %>%
  mutate(nearest_tier = factor(nearest_tier, levels = tier_levels_present)) %>%
  ggplot(aes(x = log10_dist, fill = nearest_tier, color = nearest_tier)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = log10(c(50e3, 100e3)),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values  = SV_TIER_COLORS, name = "Nearest SV tier") +
  scale_color_manual(values = SV_TIER_COLORS, name = "Nearest SV tier") +
  scale_x_continuous(
    name   = "Distance to nearest SV breakpoint (log₁₀ bp)",
    breaks = 2:8,
    labels = c("100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  labs(
    title    = "aDMR distance distribution by nearest SV tier",
    subtitle = sprintf("Spearman ρ (all tiers pooled, log10 dist vs |HP Δβ|) = %.4f", rho_all), #nolint
    y        = "Density"
  ) +
  theme_hcc

### 5c. Density by SV type (facet) ========================
type_levels_present <- intersect(SV_TYPE_LEVELS, unique(dist_plot$nearest_type))
p_type <- dist_plot %>%
  filter(nearest_type %in% type_levels_present) %>%
  mutate(nearest_type = factor(nearest_type, levels = type_levels_present)) %>%
  ggplot(aes(x = log10_dist, fill = nearest_type, color = nearest_type)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = log10(c(50e3, 100e3)),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values  = SV_TYPE_COLORS, name = "Nearest SV type") +
  scale_color_manual(values = SV_TYPE_COLORS, name = "Nearest SV type") +
  scale_x_continuous(
    name   = "Distance to nearest SV breakpoint (log₁₀ bp)",
    breaks = 2:8,
    labels = c("100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  labs(
    title    = "aDMR distance distribution by nearest SV type",
    subtitle = "Each curve: aDMRs whose nearest SV belongs to that type",
    y        = "Density"
  ) +
  theme_hcc

### 5d. Scatter: log10 dist vs |HP Δβ|  =================
# coloured by SV tier, hexbin + loess
p_scatter_tier <- dist_plot %>%
  filter(!is.na(nearest_tier), nearest_tier %in% tier_levels_present) %>%
  mutate(nearest_tier = factor(nearest_tier, levels = tier_levels_present)) %>%
  ggplot(aes(x = log10_dist, y = hp_delta, color = nearest_tier)) +
  geom_smooth(method = "gam", se = TRUE, span = 0.4, linewidth = 0.8, alpha = 0.2) +
  scale_color_manual(values = SV_TIER_COLORS, name = "Nearest SV tier") +
  scale_x_continuous(
    name   = "Distance to nearest SV bp (log₁₀)",
    breaks = 2:8,
    labels = c("100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  labs(
    title    = "Distance-decay: |HP Δβ| vs distance to nearest SV bp (by tier)",
    subtitle = "LOESS smooths; flat curves = no distance-decay (locus-specific mechanism)",
    y        = "|HP1.Methy − HP2.Methy|"
  ) +
  theme_hcc

### 5e. Scatter by SV type =============================
p_scatter_type <- dist_plot %>%
  filter(!is.na(nearest_type), nearest_type %in% type_levels_present) %>%
  mutate(nearest_type = factor(nearest_type, levels = type_levels_present)) %>%
  ggplot(aes(x = log10_dist, y = hp_delta, color = nearest_type)) +
  geom_smooth(method = "gam", se = TRUE, span = 0.4, linewidth = 0.8, alpha = 0.2) +
  scale_color_manual(values = SV_TYPE_COLORS, name = "Nearest SV type") +
  scale_x_continuous(
    name   = "Distance to nearest SV bp (log₁₀)",
    breaks = 2:8,
    labels = c("100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  labs(
    title    = "Distance-decay: |HP Δβ| vs distance to nearest SV bp (by SV type)",
    subtitle = "LOESS smooths per SV type",
    y        = "|HP1.Methy − HP2.Methy|"
  ) +
  theme_hcc

### 5f. ECDF plot =======================================
# cumulative % of aDMRs within each distance threshold
ecdf_tiers <- dist_plot %>%
  filter(!is.na(nearest_tier), nearest_tier %in% tier_levels_present) %>%
  mutate(nearest_tier = factor(nearest_tier, levels = tier_levels_present))

p_ecdf <- ggplot(ecdf_tiers, aes(x = log10_dist, color = nearest_tier)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = log10(c(50e3, 100e3, 500e3)),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = SV_TIER_COLORS, name = "Nearest SV tier") +
  scale_x_continuous(
    name   = "Distance to nearest SV bp (log₁₀)",
    breaks = 2:8,
    labels = c("100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
  ) +
  scale_y_continuous(name = "Cumulative fraction of aDMRs", labels = scales::percent) +
  labs(
    title    = "ECDF: cumulative aDMR fraction by distance threshold",
    subtitle = "Curves closely overlapping = distance distribution independent of SV tier"
  ) +
  theme_hcc

# Save combined plots
combined <- (p_global / p_tier / p_type) +
  plot_annotation(title = "aDMR–SV breakpoint distance distributions",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(opt$outdir, "sv_admr_distance_distribution.png"), combined,
       width = 11, height = 18, device = "png", dpi = 300)
saveRDS(list(global = p_global, by_tier = p_tier, by_type = p_type), file.path(opt$outdir, "sv_admr_distance_plots.rds"))

decay_combined <- (p_scatter_tier / p_scatter_type / p_ecdf) +
  plot_annotation(title = "Distance-decay and ECDF analyses",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave(file.path(opt$outdir, "sv_admr_distance_decay.png"), decay_combined,
       width = 11, height = 18, device = "png", dpi = 300)
saveRDS(list(scatter_tier = p_scatter_tier, scatter_type = p_scatter_type, ecdf = p_ecdf), file.path(opt$outdir, "sv_admr_distance_decay_plots.rds"))

fwrite(dist_df, file.path(opt$outdir, "sv_admr_distance_per_admr.csv.gz"))

rm(dist_df, dist_plot, p_global, p_tier, p_type, p_scatter_tier, p_scatter_type, p_ecdf, decay_combined, combined)
gc()

## 6. Sup2 — SV-admr distance distribution  =======================
# Distance distribution within phase block

message("Computing within-block aDMR → nearest SV bp distances …")

block_dist_df <- data.table::rbindlist(lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_phased[[pt]]
  dmr_gr <- admr_phased[[pt]]
  if (length(sv_gr) == 0 || length(dmr_gr) == 0) return(NULL)

  sv_gr  <- sv_gr[!is.na(mcols(sv_gr)$block_id)]
  dmr_gr <- dmr_gr[!is.na(mcols(dmr_gr)$block_id)]
  if (length(sv_gr) == 0 || length(dmr_gr) == 0) return(NULL)

  blocks <- intersect(unique(mcols(dmr_gr)$block_id), unique(mcols(sv_gr)$block_id))
  if (length(blocks) == 0) return(NULL)

  data.table::rbindlist(lapply(blocks, function(b) {
    dmr_blk <- dmr_gr[mcols(dmr_gr)$block_id == b]
    sv_blk  <- sv_gr[mcols(sv_gr)$block_id == b]
    if (length(dmr_blk) == 0 || length(sv_blk) == 0) return(NULL)

    hits <- distanceToNearest(dmr_blk, sv_blk, ignore.strand = TRUE)
    if (length(hits) == 0) return(NULL)

    data.table(
      patient_id   = pt,
      block_id     = as.character(b),
      dist_bp      = mcols(hits)$distance,
      nearest_tier = as.character(mcols(sv_blk[subjectHits(hits)])$sv_tier),
      nearest_type = as.character(mcols(sv_blk[subjectHits(hits)])$sv_type)
    )
  }), fill = TRUE)
}), fill = TRUE)

if (!is.null(block_dist_df) && nrow(block_dist_df) > 0L) {
  block_dist_df <- block_dist_df[!is.na(dist_bp)]
  block_dist_plot <- block_dist_df %>%
    mutate(log10_dist = log10(dist_bp + 1))

  cat("\n=== Within-block distance summary ===\n")
  cat(sprintf("  N aDMRs with SVs in same block : %d\n", nrow(block_dist_df)))
  cat(sprintf("  Median distance               : %s bp\n",
              format(median(block_dist_df$dist_bp), big.mark = ",")))
  cat(sprintf("  Mean distance                 : %s bp\n",
              format(round(mean(block_dist_df$dist_bp)), big.mark = ",")))

  p_block <- ggplot(block_dist_plot, aes(x = log10_dist)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, fill = "#3B8BD4", alpha = 0.7, color = NA) +
    geom_density(color = "#BA7517", linewidth = 0.9) +
    geom_vline(xintercept = log10(c(50e3, 100e3, 500e3)),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    annotate("text", x = log10(50e3) + 0.05, y = Inf, label = "50 kb",
             angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
    annotate("text", x = log10(100e3) + 0.05, y = Inf, label = "100 kb",
             angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
    annotate("text", x = log10(500e3) + 0.05, y = Inf, label = "500 kb",
             angle = 90, vjust = 1.3, hjust = 1.1, size = 3, color = "grey40") +
    scale_x_continuous(
      name   = "Distance to nearest SV breakpoint (within block, log₁₀ bp)",
      breaks = 0:8,
      labels = c("0", "10", "100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
    ) +
    labs(
      title    = "aDMR → nearest SV breakpoint distance (within same phase block)",
      subtitle = sprintf("N = %s aDMRs | median = %s bp",
                         format(nrow(block_dist_df), big.mark = ","),
                         format(median(block_dist_df$dist_bp), big.mark = ",")),
      y        = "Density"
    ) +
    theme_hcc

  p_block_ecdf <- ggplot(block_dist_plot, aes(x = dist_bp)) +
    stat_ecdf(geom = "step", linewidth = 0.9, color = "#3B8BD4") +
    scale_x_continuous(
      trans  = "log10",
      breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
      labels = c("1", "10", "100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb")
    ) +
    labs(
      title    = "ECDF: within-block distance to nearest SV breakpoint",
      x        = "Distance (bp, log₁₀ scale)",
      y        = "Cumulative fraction of aDMRs"
    ) +
    theme_hcc

  block_combined <- (p_block / p_block_ecdf) +
    plot_annotation(title = "Within-phase-block aDMR–SV breakpoint distance distribution",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))

  ggsave(file.path(opt$outdir, "sv_admr_distance_within_block.png"), block_combined,
         width = 10, height = 8, device = "png", dpi = 300)
  saveRDS(list(global = p_block, ecdf = p_block_ecdf),
          file.path(opt$outdir, "sv_admr_distance_within_block_plots.rds"))
  fwrite(block_dist_df, file.path(opt$outdir, "sv_admr_distance_within_block.csv.gz"))

  rm(block_dist_df, block_dist_plot, p_block, p_block_ecdf, block_combined)
  gc()
} else {
  message("Sup2: no within-block SV–aDMR pairs available; skipping.")
}


# 10. Evidence summary + output files ========================================

cat("\n=== Evidence summary ===\n")

evidence_summary <- data.frame(
  layer   = paste0("Layer ", 1:4),
  test    = c(
    sprintf("Block-label permutation (n=%d) + Wilcoxon; Jonckheere-Terpstra across tiers", N_PERM_BLOCK),
    "Paired Wilcoxon per patient×tier (BH FDR); pooled KW + JT + Dunn post-hoc (Bonferroni)",
    "aDMR–CRE overlap % by tier; Wilcoxon TAD_CTCF vs non_boundary + JT",
    "LME AIC comparison: sv_tier vs log10(dist+1) as predictor of hp_abs_diff"
  ),
  null_hypo = c(
    "SV-block and non-SV-block aDMR density are equal; no monotonic decrease across tier hierarchy",
    "SV-HP β = WT-HP β; no tier-ordered magnitude gradient",
    "CRE overlap fraction does not differ by SV tier",
    "Breakpoint distance explains hp_abs_diff as well as functional tier annotation"
  ),
  what_it_shows = c(
    "Tier hierarchy predicts block-level SV–DMR co-localization enrichment",
    "TAD/CTCF-disrupting SVs elevate haplotype-specific methylation change",
    "CRE co-localization supports locus-specific regulatory disruption mechanism",
    "Functional annotation (tier) outperforms distance as cis-effect predictor"
  ),
  jt_p = round(c(
    jt_l1$p.value,
    jt_l2$p.value,
    jt_l3$p.value,
    NA_real_
  ), 4),
  aic_tier_better = c(NA, NA, NA,
                      if (!is.na(aic_tier) && !is.na(aic_dist)) aic_tier < aic_dist else NA)
)
print(evidence_summary, row.names = FALSE)

# Save outputs
S <- opt$stratify_by
if (has_tier) {
  if (nrow(layer1_results) > 0)
    fwrite(layer1_results,          file.path(opt$outdir, sprintf("layer1_coloc_by_%s.csv",    S))) #nolint
  if (nrow(layer2_results) > 0)
    fwrite(layer2_results,          file.path(opt$outdir, sprintf("layer2_hp_by_%s.csv",       S))) #nolint
  if (nrow(layer3_results) > 0)
    fwrite(layer3_results,          file.path(opt$outdir, sprintf("layer3_cre_overlap_by_%s.csv", S))) #nolint
  if (nrow(layer4_gradient) > 0)
    fwrite(layer4_gradient,         file.path(opt$outdir, sprintf("layer4_gradient_by_%s.csv", S))) #nolint
  if (nrow(layer4_model_comparison) > 0)
    fwrite(layer4_model_comparison, file.path(opt$outdir, sprintf("layer4_model_comparison_%s.csv", S))) #nolint
}
fwrite(evidence_summary, file.path(opt$outdir, sprintf("evidence_summary_%s.csv", S))) #nolint

# Sections 5c–5f outputs (written inside each section; log here for reference)
cat("Additional output files (sections 5c-5f):\n")
cat("  sv_meth_competence.csv.gz          : Per-SV methylation-competence flag\n")
cat("  sv_meth_competence_summary.csv     : % competent by tier × sv_bio\n")
cat("  dist_bin_svbio_hp_delta.csv.gz     : Per-aDMR hp_abs_diff with distance bin + sv_bio\n")
cat("  dist_bin_svbio_summary.csv         : Median |Δβ| by dist_bin × sv_bio\n")
cat("  direction_svtype_hp_delta.csv      : Per-patient signed Δβ direction by SV type\n")
cat("  direction_svtype_summary.csv       : Pooled directionality summary by SV type\n")
cat("  layer2_2dtier_gradient.csv         : Median |Δβ| by 2D collapsed tier\n")
cat("  fig5c_meth_competent_by_tier.png\n")
cat("  fig5d_dist_decay_by_svbio.png\n")
cat("  fig5e_direction_by_svtype.png\n")
cat("  fig5f_2dtier_hp_delta.png\n")

cat("\n=== Analysis complete ===\n")
cat("Output files:\n")
if (has_tier) {
  cat("  layer1_coloc_by_tier.csv          : Block co-localization enrichment by SV tier\n")
  cat("  layer2_hp_by_tier.csv             : HP-specific Δβ by SV tier\n")
  cat("  layer3_cre_overlap_by_tier.csv    : CRE overlap fraction by SV tier\n")
  cat("  layer4_tier_gradient.csv          : Tier gradient (median ± IQR)\n")
  cat("  layer4_model_comparison.csv       : LME AIC: tier vs dist vs full (tier+dist)\n")
  cat("  layer_hbv_hp_delta.csv            : HBV-associated SV HP-specific Δβ (per patient)\n")
  cat("  hap_fig1_block_enrichment_by_tier.pdf\n")
  cat("  hap_fig2_hp_delta_by_tier.pdf\n")
  cat("  hap_fig3_cre_overlap_by_tier.pdf\n")
  cat("  hap_fig4_tier_gradient.pdf\n")
  cat("  sv_admr_distance_distribution.png\n")
  cat("  sv_admr_distance_decay.png\n")
  cat("  sv_admr_distance_per_admr.csv.gz\n")
  cat("  sv_admr_distance_within_block.png\n")
  cat("  sv_admr_distance_within_block.csv.gz\n")
}
cat("  evidence_summary.csv\n")
rm()
gc()
# =============================================================================
# Methods language
# =============================================================================
# SV functional annotation (sv_tier) was used as the primary stratification axis.
# Each SV was assigned to one of five tiers based on overlap with HepG2 Micro-C
# TAD boundaries and ENCODE CTCF peaks (ENCFF543WTP), ordered by expected
# regulatory impact: TAD_CTCF > TAD_only > CTCF_only > copy_neutral > non_boundary.
# A monotonic decrease across this hierarchy was tested with the Jonckheere-Terpstra
# trend test (nperm=2000; clinfun). Layer 1 quantified width-corrected aDMR density
# (aDMRs/Mb) in SV-containing vs. SV-free phase blocks per tier using block-label
# permutation (n=500) and Wilcoxon rank-sum test, with BH FDR across patient-tier
# combinations. Layer 2 compared SV-bearing vs. WT haplotype methylation per aDMR
# using paired Wilcoxon (per patient per tier, BH FDR), pooled Kruskal-Wallis, and
# Dunn post-hoc (Bonferroni). Layer 3 measured the fraction of aDMRs in SV-containing
# blocks overlapping annotated CREs (GeneHancer enhancers, HepG2 CTCF peaks, gene
# promoters). Layer 4 compared two linear mixed-effects models (lme4, REML=FALSE)
# by AIC: hp_abs_diff ~ sv_tier + (1|patient_id) versus ~ log10(dist+1) + (1|patient_id),
# testing whether functional tier annotation explains haplotype-specific methylation
# disruption better than proximity alone. Absence of distance-decay is interpreted as
# evidence for locus-specific, not globally diffuse, cis-regulatory disruption.
