# conda env: renv
# =============================================================================
# Analysis D: Insulation score validation
# =============================================================================
# Maps each SV breakpoint to its reference (HepG2 Micro-C) insulation score and:
#   (i)  tests whether the breakpoint's reference insulation score differs by
#        SV tier (TAD+CTCF disrupting expected to be at low-insulation, i.e.
#        strong boundary, positions);
#   (ii) tests whether nearby boundary strength correlates with DMR enrichment
#        in the same patient (Spearman ρ).
#
# Note: With reference-only insulation, this characterises *positional context*
#   of SV breakpoints, not patient-specific 3D reorganisation. Patient-specific
#   Δinsulation requires patient Micro-C/Hi-C which is not available.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(stringr)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--sv_strat_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz"),
  make_option("--dmr_file", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/consensus_dmrs_per_patient.csv.gz"),
  make_option("--insulation_tsv", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/tad_ctcf_validation/insulation_8kb.tsv.gz"),
  make_option("--window_bp", type = "integer", default = 240000L,
              help = "Which insulation window col to use (must match a window_bp in --windows passed to 06b) [default: %default]"),
  make_option("--outdir", type = "character",
              default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/tad_ctcf_validation"),
  make_option("--run_id", type = "character", default = "tier_v2"),
  make_option("--dmr_window", type = "integer", default = 50000L,
              help = "±window for DMR count near bp [bp]"),
  make_option("--no_plot", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# ── Load data ───────────────────────────────────────────────────────────────
message("Reading SV: ", opt$sv_strat_file)
sv_df <- fread(opt$sv_strat_file)
sv_gr <- makeGRangesFromDataFrame(sv_df, keep.extra.columns = TRUE)

message("Reading insulation: ", opt$insulation_tsv)
ins_df <- fread(opt$insulation_tsv)
ins_cols <- colnames(ins_df)
log2_col <- grep(paste0("^log2_insulation_score_", opt$window_bp, "$"), ins_cols, value = TRUE)
bnd_col  <- grep(paste0("^boundary_strength_",     opt$window_bp, "$"), ins_cols, value = TRUE)
is_bnd_col <- grep(paste0("^is_boundary_",         opt$window_bp, "$"), ins_cols, value = TRUE)
if (length(log2_col) == 0) {
  cat("Available insulation columns:\n")
  print(ins_cols)
  stop(sprintf("No log2_insulation_score_%d column found", opt$window_bp))
}
message("Using: ", log2_col)

ins_df <- ins_df %>%
  dplyr::select(chrom, start, end, is_bad_bin,
                log2_ins = !!log2_col,
                boundary_strength = any_of(bnd_col),
                is_boundary       = any_of(is_bnd_col))

ins_gr <- makeGRangesFromDataFrame(ins_df, keep.extra.columns = TRUE,
                                    seqnames.field = "chrom")

# ── Annotate SV bp with insulation ──────────────────────────────────────────
ov <- findOverlaps(sv_gr, ins_gr, select = "first")
sv_df$log2_ins          <- ifelse(is.na(ov), NA_real_, mcols(ins_gr)$log2_ins[ov])
sv_df$boundary_strength <- ifelse(is.na(ov), NA_real_, mcols(ins_gr)$boundary_strength[ov])
sv_df$is_boundary       <- ifelse(is.na(ov), NA, mcols(ins_gr)$is_boundary[ov])

cat(sprintf("Annotated %d / %d SV bp with insulation\n",
            sum(!is.na(sv_df$log2_ins)), nrow(sv_df)))

# ── D1. Tier-wise insulation score distribution ─────────────────────────────
sv_df$stratification <- factor(sv_df$stratification, levels = STRAT_LEVELS)
D1_df <- sv_df %>% dplyr::filter(!is.na(log2_ins))

kw_test <- kruskal.test(log2_ins ~ stratification, data = D1_df)
cat(sprintf("[D1] Kruskal-Wallis log2_ins ~ tier: H = %.3f  df = %d  p = %.3g\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

# Pairwise Wilcoxon: TAD+CTCF disrupting vs Non-boundary
D1_wilcox <- tryCatch(
  wilcox.test(D1_df$log2_ins[D1_df$stratification == "TAD+CTCF disrupting"],
              D1_df$log2_ins[D1_df$stratification == "Non-boundary"],
              alternative = "less", exact = FALSE),
  error = function(e) list(p.value = NA_real_))
cat(sprintf("[D1] Wilcoxon (TAD+CTCF disrupting < Non-boundary) p = %.3g\n",
            D1_wilcox$p.value))

D1_summary <- D1_df %>%
  dplyr::group_by(stratification) %>%
  dplyr::summarise(
    n           = dplyr::n(),
    median_ins  = round(median(log2_ins, na.rm = TRUE), 3),
    mean_ins    = round(mean(log2_ins, na.rm = TRUE), 3),
    pct_boundary = round(mean(is_boundary, na.rm = TRUE) * 100, 1),
    median_bnd_strength = round(median(boundary_strength, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("[D1] Tier-wise insulation summary:\n"); print(D1_summary)
fwrite(D1_summary, file.path(opt$outdir, sprintf("%s_D1_tier_insulation_summary.csv", opt$run_id)))

if (!opt$no_plot) {
  pD1 <- ggplot(D1_df, aes(x = stratification, y = log2_ins, fill = stratification)) +
    geom_violin(alpha = 0.7, color = "grey40", linewidth = 0.4,
                trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.6,
                 fill = "white", color = "grey30", alpha = 0.8) +
    scale_fill_manual(values = STRAT_COLORS, guide = "none") +
    labs(
      title    = "D1. Reference insulation score at SV breakpoints",
      subtitle = sprintf("Kruskal-Wallis p = %.3g | TAD+CTCF disrupting vs Non-boundary p = %.3g (one-sided)",
                          kw_test$p.value, D1_wilcox$p.value),
      x = NULL,
      y = sprintf("log2 insulation (window = %d kb)\nlower = stronger boundary",
                  opt$window_bp / 1000),
      caption = "HepG2 Micro-C 8kb bins"
    ) + theme_hcc +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(opt$outdir, sprintf("%s_D1_tier_insulation.png", opt$run_id)),
         pD1, width = 9, height = 5, device = "png", dpi = 300)
  saveRDS(pD1, file.path(opt$outdir, sprintf("%s_D1_tier_insulation.rds", opt$run_id)))
}

# ── D2. Correlate per-SV insulation with nearby DMR count ──────────────────
message("Loading DMRs for per-patient bp-DMR count")
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

# Count DMRs in ±dmr_window per bp, per patient
sv_list <- sv_df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$sample)

bp_dmr_counts <- lapply(names(sv_list), function(pt) {
  sv <- sv_list[[pt]]
  dmr <- consensus_dmr[[pt]]
  if (is.null(dmr) || length(dmr) == 0) return(NULL)
  wins <- suppressWarnings(GenomicRanges::resize(sv, width = opt$dmr_window * 2L,
                                                  fix = "center"))
  hits <- findOverlaps(wins, dmr, minoverlap = 1L)
  n_dmr <- tabulate(queryHits(hits), nbins = length(wins))
  data.frame(
    patient_id = pt,
    bp_id      = sv$bp_id,
    sourceId   = sv$sourceId,
    n_dmr      = n_dmr,
    log2_ins   = mcols(sv)$log2_ins,
    bnd_str    = mcols(sv)$boundary_strength,
    tier       = as.character(mcols(sv)$stratification)
  )
}) %>% bind_rows()

bp_dmr_counts <- bp_dmr_counts %>% dplyr::filter(!is.na(log2_ins))

D2_cor <- bp_dmr_counts %>%
  dplyr::group_by(tier) %>%
  dplyr::summarise(
    n_bp           = dplyr::n(),
    spearman_rho   = cor(log2_ins, n_dmr, method = "spearman", use = "complete.obs"),
    p_value        = tryCatch(
      cor.test(log2_ins, n_dmr, method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_),
    .groups        = "drop"
  ) %>% dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH"))
cat("[D2] Spearman ρ (log2_ins vs nearby DMR count) by tier:\n"); print(D2_cor)
fwrite(D2_cor, file.path(opt$outdir, sprintf("%s_D2_ins_dmr_correlation.csv", opt$run_id)))
fwrite(bp_dmr_counts, file.path(opt$outdir, sprintf("%s_D2_bp_dmr_insulation.csv.gz", opt$run_id)))

if (!opt$no_plot) {
  bp_plot <- bp_dmr_counts %>%
    dplyr::filter(tier %in% STRAT_LEVELS) %>%
    dplyr::mutate(tier = factor(tier, levels = STRAT_LEVELS))

  pD2 <- ggplot(bp_plot, aes(x = log2_ins, y = n_dmr, color = tier)) +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8,
                aes(group = tier), span = 0.8) +
    scale_color_manual(values = STRAT_COLORS, name = NULL) +
    facet_wrap(~tier, nrow = 1, scales = "free_x") +
    labs(
      title    = sprintf("D2. Reference insulation vs nearby DMR count (±%dkb)",
                          opt$dmr_window / 1000),
      subtitle = "Negative log2_ins = stronger boundary",
      x        = sprintf("log2 insulation (%d kb window)", opt$window_bp / 1000),
      y        = "DMR # near bp"
    ) + theme_hcc + theme(legend.position = "none")

  ggsave(file.path(opt$outdir, sprintf("%s_D2_ins_vs_dmr.png", opt$run_id)),
         pD2, width = 14, height = 4, device = "png", dpi = 300)
  saveRDS(pD2, file.path(opt$outdir, sprintf("%s_D2_ins_vs_dmr.rds", opt$run_id)))
}

# ── D3. Within "TAD+CTCF disrupting" tier — does boundary_strength predict DMRs? ──
# Stratify by boundary_strength quartiles within that tier
D3_within <- bp_dmr_counts %>%
  dplyr::filter(tier == "TAD+CTCF disrupting", !is.na(bnd_str)) %>%
  dplyr::mutate(
    bnd_q = cut(bnd_str,
                breaks = quantile(bnd_str, probs = c(0, .25, .5, .75, 1), na.rm = TRUE),
                labels = c("Q1 (weakest)", "Q2", "Q3", "Q4 (strongest)"),
                include.lowest = TRUE)
  )

D3_summary <- D3_within %>%
  dplyr::group_by(bnd_q) %>%
  dplyr::summarise(n_bp        = dplyr::n(),
                   median_dmr  = round(median(n_dmr, na.rm = TRUE), 3),
                   mean_dmr    = round(mean(n_dmr, na.rm = TRUE), 3),
                   .groups     = "drop")
cat("[D3] Within TAD+CTCF disrupting — DMR vs boundary strength quartile:\n")
print(D3_summary)
fwrite(D3_summary, file.path(opt$outdir, sprintf("%s_D3_bnd_strength_quartile.csv", opt$run_id)))

D3_kw <- tryCatch(
  kruskal.test(n_dmr ~ bnd_q, data = D3_within),
  error = function(e) list(p.value = NA_real_,
                           statistic = NA_real_,
                           parameter = NA_real_))
.h  <- if (is.null(D3_kw$statistic)) NA_real_ else as.numeric(D3_kw$statistic)
.df <- if (is.null(D3_kw$parameter)) NA       else as.character(D3_kw$parameter)
cat(sprintf("[D3] Kruskal-Wallis: H = %.3f  df = %s  p = %.3g\n",
            .h, .df, D3_kw$p.value))

if (!opt$no_plot && nrow(D3_within) > 0) {
  pD3 <- ggplot(D3_within, aes(x = bnd_q, y = n_dmr, fill = bnd_q)) +
    geom_boxplot(width = 0.55, alpha = 0.8, outlier.size = 0.6) +
    scale_fill_brewer(palette = "Reds", guide = "none") +
    labs(
      title    = "D3. Within TAD+CTCF disrupting — boundary strength vs DMR count",
      subtitle = sprintf("Kruskal-Wallis p = %.3g", D3_kw$p.value),
      x        = sprintf("Boundary strength quartile (window=%dkb)", opt$window_bp / 1000),
      y        = sprintf("DMR # in ±%dkb of breakpoint", opt$dmr_window / 1000)
    ) + theme_hcc

  ggsave(file.path(opt$outdir, sprintf("%s_D3_bnd_strength_quartile.png", opt$run_id)),
         pD3, width = 8, height = 5, device = "png", dpi = 300)
  saveRDS(pD3, file.path(opt$outdir, sprintf("%s_D3_bnd_strength_quartile.rds", opt$run_id)))
}

cat("\n=== Analysis D complete ===\n")
cat(sprintf("Outputs in: %s\n", opt$outdir))
cat("  - <run>_D1_tier_insulation.{csv,png,rds}\n")
cat("  - <run>_D2_ins_dmr_correlation.csv + D2_bp_dmr_insulation.csv.gz + D2 plot\n")
cat("  - <run>_D3_bnd_strength_quartile.{csv,png,rds}\n")
