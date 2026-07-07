#!/usr/bin/env Rscript
# trans_negative_control.R
# Two-part trans negative control analysis:
#
# Part 1 (PRIMARY): Window-scale enrichment decay
#   - Uses tier_50kb_v2_window_enrich_full.csv (10 kb → 1 Mb windows)
#   - Shows enrichment_ratio decays to ~1 (null) at 1 Mb → cis-specificity
#   - Tests: paired Wilcoxon across patients (50 kb vs 1 Mb enrichment)
#
# Part 2 (EXPLORATORY): HP |Δβ| by SV proximity zone in aDMR data
#   - Uses sv_admr_distance_per_admr.csv.gz
#   - NOTE: expects no gradient because ALL rows are aDMRs (pre-selected for
#     high HP divergence). Absence of gradient is interpretable: aDMR HP
#     divergence is constitutive, not exclusively SV-driven. Consistent with
#     the proximity enrichment model (C2) and "no tier magnitude gradient" (C3).
#
# Usage: mamba run -n renv Rscript post_processing/trans_negative_control.R \
#          2>&1 | tee logs/trans_neg_ctrl.log

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

OUTDIR    <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")
ENRICH_CSV <- file.path(OUTDIR, "02.sv_dmr_enrichment/tier_50kb_v2_window_enrich_full.csv")
DIST_CSV   <- file.path(OUTDIR, "03.haplotype_sv_admr_analysis/sv_admr_distance_per_admr.csv.gz")
ENRICH_OUT <- file.path(OUTDIR, "02.sv_dmr_enrichment")
FIG_PNG    <- file.path(OUTDIR, "figs/png/fig_trans_neg_ctrl.png")
FIG_PDF    <- file.path(OUTDIR, "figs/panels/fig_trans_neg_ctrl.pdf")

dir.create(file.path(OUTDIR, "figs/png"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "figs/panels"), showWarnings = FALSE, recursive = TRUE)

theme_hcc <- theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(colour = "grey50", size = 10),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom"
  )

# PART 1 — Window-scale enrichment decay (primary trans negative control) ======
message("=== Part 1: Window-scale enrichment decay ===")
we <- fread(ENRICH_CSV)
setnames(we, tolower(names(we)))
# Keep only cnv_class (non-tier) and all window sizes
message(sprintf("  Loaded: %d rows | windows: %s kb",
  nrow(we),
  paste(sort(unique(we$window_kb)), collapse = ", ")))

# Aggregate enrichment_ratio across cnv_class per patient per window
# (all-SV average for scale analysis; class-specific in Part 1b)
we_pt <- we[, .(
  mean_enrich = mean(enrichment_ratio, na.rm = TRUE),
  mean_pct    = mean(pct_sv_w_dmr,    na.rm = TRUE),
  n_classes   = .N
), by = .(patient_id, window_kb)][order(patient_id, window_kb)]

# Part 1a: Paired Wilcoxon: enrichment at 50 kb vs 1000 kb =====================
pts_both <- intersect(
  we_pt[window_kb == 50,   patient_id],
  we_pt[window_kb == 1000, patient_id]
)
cis_enrich   <- we_pt[window_kb == 50   & patient_id %in% pts_both,
                       mean_enrich[order(patient_id)]]
trans_enrich <- we_pt[window_kb == 1000 & patient_id %in% pts_both,
                       mean_enrich[order(patient_id)]]

pt1_wt <- suppressWarnings(
  wilcox.test(cis_enrich, trans_enrich, paired = TRUE,
              alternative = "greater", exact = FALSE)
)
message(sprintf(
  "  Paired Wilcoxon (cis 50kb > trans 1Mb): p = %.4f (n=%d patients)",
  pt1_wt$p.value, length(pts_both)
))
message(sprintf(
  "  Median enrichment: 50kb = %.3f | 1000kb = %.3f",
  median(cis_enrich), median(trans_enrich)
))

# Part 1b: Enrichment_ratio approaching null (1.0) at 1 Mb =====================
# Wilcoxon one-sample: enrichment_ratio at 1000 kb not different from 1 (null)
trans_all <- we[window_kb == 1000 & !is.na(enrichment_ratio), enrichment_ratio]
oneside_wt <- suppressWarnings(
  wilcox.test(trans_all, mu = 1.0, alternative = "greater", exact = FALSE)
)
message(sprintf(
  "  One-sample Wilcoxon at 1Mb (H0: ratio = 1): p = %.4f | median ratio = %.3f",
  oneside_wt$p.value, median(trans_all, na.rm = TRUE)
))

# Save Part 1 stats ============================================================
pt1_stats <- data.table(
  analysis     = c("paired_Wilcoxon_50kb_vs_1Mb", "onesample_1Mb_vs_null"),
  n            = c(length(pts_both), length(trans_all)),
  median_A     = c(median(cis_enrich), median(trans_all, na.rm = TRUE)),
  median_B     = c(median(trans_enrich), 1.0),
  wilcox_p     = c(pt1_wt$p.value, oneside_wt$p.value),
  alt          = c("50kb > 1Mb (paired)", "1Mb > 1.0 (one-sample)")
)
fwrite(pt1_stats, file.path(ENRICH_OUT, "trans_neg_ctrl_window_stats.csv"))

# Summary per window (for figure)
we_sum <- we[, .(
  median_enrich = median(enrichment_ratio, na.rm = TRUE),
  q25_enrich    = quantile(enrichment_ratio, 0.25, na.rm = TRUE),
  q75_enrich    = quantile(enrichment_ratio, 0.75, na.rm = TRUE),
  pct_sig       = mean(p_fdr < 0.05, na.rm = TRUE) * 100
), by = window_kb][order(window_kb)]
fwrite(we_sum, file.path(ENRICH_OUT, "trans_neg_ctrl_window_summary.csv"))

# PART 2 — HP |Δβ| by SV proximity zone (aDMR data; expected: NO gradient) =====
message("\n=== Part 2: HP |Δβ| by SV proximity zone (aDMR-level) ===")
dt <- fread(DIST_CSV)
setnames(dt, tolower(names(dt)))
dt[, abs_hp_delta := abs(hp_delta)]
dt[, dist_bin := fcase(
  dist_bp <= 50000L,                       "Cis (<=50kb)",
  dist_bp > 50000L & dist_bp <= 1000000L,  "Mid (50kb-1Mb)",
  dist_bp > 1000000L,                      "Trans (>1Mb)",
  default = NA_character_
)]
dt[, dist_bin := factor(dist_bin,
  levels = c("Cis (<=50kb)", "Mid (50kb-1Mb)", "Trans (>1Mb)"))]

bin_counts <- dt[!is.na(dist_bin), .N, by = dist_bin][order(dist_bin)]
bin_stats  <- dt[!is.na(dist_bin), .(
  n          = .N,
  median_abs = round(median(abs_hp_delta, na.rm = TRUE), 4),
  q25        = round(quantile(abs_hp_delta, 0.25, na.rm = TRUE), 4),
  q75        = round(quantile(abs_hp_delta, 0.75, na.rm = TRUE), 4)
), by = dist_bin][order(dist_bin)]

message("  aDMR bin counts + median |hp_delta|:")
print(bin_stats)

# KW test (expected: ns — supports constitutive aDMR model)
# ns alone is not evidence of absence; report eta-squared to bound effect size
kw <- kruskal.test(abs_hp_delta ~ dist_bin, data = dt[!is.na(dist_bin)])
kw_n <- nrow(dt[!is.na(dist_bin)])
kw_k <- nlevels(dt[!is.na(dist_bin)]$dist_bin)
kw_eta_sq <- (kw$statistic - kw_k + 1) / (kw_n - kw_k)  # epsilon-squared approx
message(sprintf(
  "  KW: H=%.3f, p=%.3e, eta2=%.4f  (expected ns = constitutive; eta2 bounds effect size)",
  kw$statistic, kw$p.value, max(kw_eta_sq, 0)))

# Pairwise (cis vs trans)
a_cis   <- dt[dist_bin == "Cis (<=50kb)",  abs_hp_delta]
a_trans <- dt[dist_bin == "Trans (>1Mb)",  abs_hp_delta]
cis_trans_wt <- suppressWarnings(
  wilcox.test(a_cis, a_trans, alternative = "greater", exact = FALSE))
# Rank-biserial correlation as effect size
rbc_cis_trans <- 1 - 2 * cis_trans_wt$statistic / (length(a_cis) * length(a_trans))
message(sprintf(
  "  Wilcoxon cis > trans: p = %.4f, r = %.3f (expected ns; |r| < 0.1 = negligible)",
  cis_trans_wt$p.value, rbc_cis_trans
))

# Per-patient (check consistency)
patient_bin <- dt[!is.na(dist_bin), .(
  median_abs = median(abs_hp_delta, na.rm = TRUE), n = .N
), by = .(patient_id, dist_bin)][order(patient_id, dist_bin)]

fwrite(bin_stats,   file.path(ENRICH_OUT, "trans_neg_ctrl_admr_bin_stats.csv"))
fwrite(patient_bin, file.path(ENRICH_OUT, "trans_neg_ctrl_admr_patient_bin.csv"))

# FIGURES ======================================================================
WIN_COLORS <- c(
  "10"   = "#D6604D",
  "50"   = "#E24B4A",
  "100"  = "#E67E22",
  "500"  = "#B0BEC5",
  "1000" = "#607D8B"
)
BIN_COLORS <- c(
  "Cis (<=50kb)"    = "#E24B4A",
  "Mid (50kb-1Mb)"  = "#E67E22",
  "Trans (>1Mb)"    = "#95A5A6"
)

# P1A: Enrichment ratio decay curve (per patient, all cnv_class combined) ======
p1a_df <- we_pt
p1a_df[, window_label := factor(window_kb, levels = sort(unique(window_kb)),
                                 labels = paste0("±", sort(unique(window_kb)), "kb"))]

p1a <- ggplot(p1a_df, aes(x = as.factor(window_kb), y = mean_enrich, group = patient_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.7) +
  geom_line(colour = "grey70", linewidth = 0.5, alpha = 0.8) +
  geom_point(aes(colour = as.factor(window_kb)), size = 2.5, alpha = 0.85) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               colour = "#333333", linewidth = 1.1) +
  stat_summary(aes(group = 1), fun = median, geom = "point",
               colour = "#333333", size = 3, shape = 18) +
  scale_x_discrete(labels = paste0("±", sort(unique(we$window_kb)), "kb")) +
  scale_colour_manual(values = WIN_COLORS, guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           size = 3.2, colour = "grey40",
           label = sprintf(
             "Paired Wilcoxon\n50kb > 1Mb: p=%.3f\n1Mb vs null(1.0): p=%.3f",
             pt1_wt$p.value, oneside_wt$p.value
           )) +
  labs(
    title    = "A. Enrichment ratio decays to null at 1 Mb (spatial cis-control)",
    subtitle = "Lines = patients | thick = cohort median | dashed = null (ratio=1)",
    x        = "SV breakpoint window radius",
    y        = "DMR enrichment ratio (obs / null mean)"
  ) +
  theme_hcc

# P1B: % enrichment significant (FDR<0.05) by window ===========================
p1b_df <- we[, .(pct_sig = mean(p_fdr < 0.05, na.rm = TRUE) * 100), by = window_kb]

p1b <- ggplot(p1b_df, aes(x = as.factor(window_kb), y = pct_sig,
                            fill = as.factor(window_kb))) +
  geom_col(alpha = 0.85, colour = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = WIN_COLORS, guide = "none") +
  scale_x_discrete(labels = paste0("±", sort(unique(we$window_kb)), "kb")) +
  scale_y_continuous(labels = function(x) paste0(round(x), "%"),
                     limits = c(0, 100)) +
  labs(
    title    = "B. Fraction significant (FDR<0.05) by window",
    subtitle = "Enrichment significance drops with distance",
    x        = "Window radius",
    y        = "% patient-class pairs with FDR < 0.05"
  ) +
  theme_hcc

# P2: aDMR HP |Δβ| — no gradient (expected) ====================================
cnt_str <- paste(
  sprintf("%s n=%s", bin_counts$dist_bin, format(bin_counts$N, big.mark = ",")),
  collapse = " | "
)
p2 <- ggplot(dt[!is.na(dist_bin)],
             aes(x = dist_bin, y = abs_hp_delta, fill = dist_bin)) +
  geom_violin(alpha = 0.65, colour = "grey40", linewidth = 0.4,
              trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.07, outlier.size = 0.2, fill = "white",
               colour = "grey30", alpha = 0.85) +
  scale_fill_manual(values = BIN_COLORS, guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           size = 3.0, colour = "grey40",
           label = sprintf(
             "KW p = %.2e (ns expected)\ncis vs trans p = %.3f",
             kw$p.value, cis_trans_wt$p.value
           )) +
  labs(
    title    = "C. aDMR HP |Δβ| by SV proximity (no gradient expected)",
    subtitle = paste0("ALL entries are aDMRs — constitutively high HP divergence regardless of SV\n", cnt_str),
    x = NULL, y = "|HP Δβ| (SV-HP − WT-HP)"
  ) +
  theme_hcc +
  theme(plot.subtitle = element_text(size = 8))

# P3: Per-patient median by zone (paired lines) ================================
p3 <- ggplot(patient_bin, aes(x = dist_bin, y = median_abs, colour = dist_bin)) +
  geom_line(aes(group = patient_id), colour = "grey65",
            linewidth = 0.5, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3,
               colour = "grey20", linewidth = 0.6) +
  scale_colour_manual(values = BIN_COLORS, guide = "none") +
  labs(
    title    = "D. Per-patient median |Δβ| by zone",
    subtitle = "No systematic rank order = aDMR HP divergence is patient-specific,\nnot proximity-driven (supports constitutive model)",
    x = NULL, y = "Median |HP Δβ|"
  ) +
  theme_hcc +
  theme(plot.subtitle = element_text(size = 8))

# Assembly =====================================================================
fig_main <- (p1a | p1b) / (p2 | p3) +
  plot_annotation(
    title   = "Trans Negative Control: SV-DMR enrichment is cis-specific (decays to null at 1 Mb)",
    caption = sprintf(
      "Part 1: window_enrich_full.csv (%d patient-class-window rows) | Part 2: sv_admr_distance_per_admr.csv.gz (%s aDMR pairs)",
      nrow(we), format(nrow(dt), big.mark = ",")
    ),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(FIG_PNG, fig_main, width = 14, height = 10, dpi = 200)
ggsave(FIG_PDF, fig_main, width = 14, height = 10)
message(sprintf("\nFigure saved: %s", FIG_PNG))

message("\nDone.")
