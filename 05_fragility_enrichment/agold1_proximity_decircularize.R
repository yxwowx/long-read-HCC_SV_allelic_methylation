#!/usr/bin/env Rscript
# A_GOLD1: De-circularization — does Gold OR=3.80 reflect proximity selection bias?
#
# Problem: Gold criteria = phased cis-concordance + n_patients>=3 + bp_dist<=50kb.
# bp_dist<=50kb is defined relative to the patient's own SV. SVs enrich at SegDup
# (OR=2.11), so Gold loci inherit SV SegDup co-localization by construction.
# ALL 99 Gold unique loci have bp_dist <= 50kb (100%).
#
# Design:
#   Gold*  = n_patients >= 3 ONLY (drop proximity; phase concordance requires
#            SV context so cannot be preserved independently).
#            Source: confident_dmr_per_patient.csv.gz (n=98,727 unique loci)
#   Gold*_strong = n_patients >= 3 AND median |HP delta| >= 0.20
#   Gold*_R5     = n_patients >= 5 (most stringent recurrence alone)
#
# Test: if Gold OR=3.80 is proximity-inflated, Gold* OR should collapse toward
#       all-aDMR OR (~2.25). If Gold* OR remains elevated, recurrence itself
#       (not proximity) drives the fragility enrichment.
#
# ⚠ NOTE (2026-07, verified against manuscript_v5.3.md): the "Gold OR=3.80"
# sensitivity question above is no longer in the current manuscript (grepped,
# no hits for "3.80"/"OR=3.8"). But this script's own byproduct output,
# `a_gold1_or_table.csv` (the Gold*/Gold*_R3/Gold*_R5 recurrence-only rows),
# is the CURRENT primary source of the manuscript's headline "Constitutional
# aDMR hotspots scale with recurrence" result (Fig 4B: All aDMR OR=1.32,
# Gold*_R3 OR=1.85, Gold*_R5 OR=4.22 [3.82-4.67], p<10^-174 — cited in
# Abstract/Results/Discussion/Conclusions), consumed by
# `viz/v4/fig2_segdup_coexistence.R` panel D. Read this script as "the
# recurrence-gradient script" first; the original proximity de-circularization
# framing above is secondary/historical context for why it was written this way.
#
# Run: mamba run -n renv Rscript post_processing/agold1_proximity_decircularize.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
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

set.seed(42)

# Paths ========================================================================
ADM_FILE  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz")
GOLD_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILV_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/silver_tier.csv")
SEGDUP    <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD       <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW    <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
FAI       <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_CTRL_MULT <- 10L
HP_STRONG   <- 0.20   # |HP delta| threshold for Gold*_strong

# 0. Load reference files ======================================================
message("Loading reference files …")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z")) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

segdup_gr <- import(SEGDUP, format = "BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr    <- import(LAD,    format = "BED"); seqlevelsStyle(lad_gr)    <- "UCSC"
bw        <- BigWigFile(PC1_BW)

# Helper: annotate GRanges with fragility features =============================
annotate_gr <- function(gr) {
  gr$segdup <- overlapsAny(gr, segdup_gr)
  gr$lad    <- overlapsAny(gr, lad_gr)
  pc1_vals  <- summary(bw, which = gr, type = "mean", defaultValue = NA_real_)
  gr$pc1    <- unlist(lapply(pc1_vals, function(x)
    if (length(x$score) == 0) NA_real_ else x$score[1]))
  gr$b_compartment <- !is.na(gr$pc1) & gr$pc1 < 0
  gr
}

# Helper: generate matched random controls =====================================
make_controls <- function(case_gr, n_mult = N_CTRL_MULT) {
  ctrl_list <- lapply(seq_len(length(case_gr)), function(i) {
    chr <- as.character(seqnames(case_gr)[i])
    w   <- width(case_gr)[i]
    len <- chrom_sizes$len[chrom_sizes$chr == chr]
    if (!length(len)) return(NULL)
    max_s <- len - w
    if (max_s < 1) return(NULL)
    ss <- sample.int(max_s, n_mult, replace = (n_mult > max_s))
    GRanges(chr, IRanges(ss, ss + w - 1L), is_case = 0L)
  })
  do.call(c, Filter(Negate(is.null), ctrl_list))
}

# Helper: logistic GLM and extract OR ==========================================
# Returns list(primary = <original segdup-focused row, unchanged for Wald tests
# and the existing forest plot>, multifeature = <long-format SegDup/LAD/B-comp
# table, adjusted + marginal, for the 3-population x fragility-feature figure>)
run_glm <- function(case_gr, label, n_ctrl_mult = N_CTRL_MULT) {
  ctrl_gr <- make_controls(case_gr, n_ctrl_mult)
  mcols(case_gr)$is_case <- 1L
  all_gr   <- annotate_gr(c(case_gr, ctrl_gr))

  df <- data.frame(
    is_case       = all_gr$is_case,
    segdup        = as.integer(all_gr$segdup),
    lad           = as.integer(all_gr$lad),
    b_compartment = as.integer(all_gr$b_compartment)
  )
  df$obs_weight <- ifelse(df$is_case == 1L, 1, 1 / n_ctrl_mult)

  m  <- glm(is_case ~ segdup + lad + b_compartment,
            data = df, weights = obs_weight, family = binomial())
  co <- summary(m)$coefficients
  ci <- confint.default(m)

  pct_segdup_case <- 100 * mean(df$segdup[df$is_case == 1], na.rm = TRUE)
  pct_segdup_ctrl <- 100 * mean(df$segdup[df$is_case == 0], na.rm = TRUE)

  primary <- data.frame(
    label         = label,
    n_case        = sum(df$is_case),
    n_ctrl        = sum(!df$is_case),
    pct_segdup_case = pct_segdup_case,
    pct_segdup_ctrl = pct_segdup_ctrl,
    OR            = exp(co["segdup", 1]),
    CI_lo         = exp(ci["segdup", 1]),
    CI_hi         = exp(ci["segdup", 2]),
    beta          = co["segdup", 1],
    se            = co["segdup", 2],
    p             = co["segdup", 4],
    sig           = as.character(cut(co["segdup",4],
                      c(-Inf,.001,.01,.05,Inf), labels=c("***","**","*","ns")))
  )

  # Multi-feature table: adjusted (from m above) + marginal (one feature at a
  # time, same case/control df) for SegDup, LAD, B-compartment
  FEATS <- c("segdup", "lad", "b_compartment")
  mf_rows <- do.call(rbind, lapply(FEATS, function(feat) {
    adj_OR <- exp(co[feat, 1]); adj_CI <- exp(ci[feat, ]); adj_p <- co[feat, 4]

    m_uni  <- glm(as.formula(sprintf("is_case ~ %s", feat)),
                  data = df, weights = obs_weight, family = binomial())
    co_uni <- summary(m_uni)$coefficients
    ci_uni <- confint.default(m_uni)

    rbind(
      data.frame(label = label, feature = feat, framework = "multivariate",
                 OR = adj_OR, CI_lo = adj_CI[1], CI_hi = adj_CI[2], p = adj_p,
                 n_case = sum(df$is_case), n_ctrl = sum(!df$is_case)),
      data.frame(label = label, feature = feat, framework = "marginal",
                 OR = exp(co_uni[feat, 1]), CI_lo = exp(ci_uni[feat, 1]),
                 CI_hi = exp(ci_uni[feat, 2]), p = co_uni[feat, 4],
                 n_case = sum(df$is_case), n_ctrl = sum(!df$is_case))
    )
  }))
  mf_rows$sig <- as.character(cut(mf_rows$p, c(-Inf,.001,.01,.05,Inf),
                                   labels = c("***","**","*","ns")))

  list(primary = primary, multifeature = mf_rows)
}

# 1. Build locus sets ==========================================================
message("Building locus sets …")

# A) Original Gold (with proximity)
gold_raw <- fread(GOLD_FILE)
gold_u   <- unique(gold_raw[, .(admr_chr, admr_start, admr_end, nCG)]) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr))
cat(sprintf("Gold (original, with proximity): %d unique loci\n", nrow(gold_u)))

# B) Original Silver
silv_raw <- fread(SILV_FILE)
silv_u   <- unique(silv_raw[, .(admr_chr, admr_start, admr_end)]) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr))
cat(sprintf("Silver (original): %d unique loci\n", nrow(silv_u)))

# C) Load all confident aDMRs for Gold* derivation
message("Loading all confident aDMRs …")
adm_all <- fread(ADM_FILE) |>
  filter(!is.na(admr_chr), grepl("^chr[0-9XY]+$", admr_chr)) |>
  mutate(hp_delta = abs(HP1.Methy - HP2.Methy))

# Collapse to unique loci: take max n_patients, median hp_delta
adm_loci <- adm_all |>
  group_by(admr_chr, admr_start, admr_end) |>
  summarise(
    n_patients_max  = max(n_patients),
    hp_delta_median = median(hp_delta, na.rm = TRUE),
    nCG             = dplyr::first(nCG),
    .groups = "drop"
  ) |>
  as.data.table()

cat(sprintf("All unique aDMR loci: %d\n", nrow(adm_loci)))

# D) Gold* sets (proximity-free)
goldstar_R3 <- adm_loci[n_patients_max >= 3]
cat(sprintf("Gold*_R3 (n_patients>=3, no proximity): %d unique loci\n", nrow(goldstar_R3)))

goldstar_R3_strong <- adm_loci[n_patients_max >= 3 & hp_delta_median >= HP_STRONG]
cat(sprintf("Gold*_R3_strong (n_patients>=3, |HP delta|>=%.2f): %d unique loci\n",
            HP_STRONG, nrow(goldstar_R3_strong)))

goldstar_R5 <- adm_loci[n_patients_max >= 5]
cat(sprintf("Gold*_R5 (n_patients>=5, no proximity): %d unique loci\n", nrow(goldstar_R5)))

# All aDMRs (unique loci with has_admr_support)
all_u <- adm_loci
cat(sprintf("All aDMR unique loci: %d\n", nrow(all_u)))

# 2. Run GLMs for each set =====================================================
message("\nRunning GLMs …")

make_gr <- function(dt) {
  GRanges(dt$admr_chr,
          IRanges(dt$admr_start, dt$admr_end),
          is_case = 1L)
}

message("  Gold (original) …")
glm_gold <- run_glm(make_gr(gold_u),           "Gold (original, bp<=50kb)")
res_gold <- glm_gold$primary

message("  All aDMR …")
glm_all  <- run_glm(make_gr(all_u),            "All aDMR (no filter)")
res_all  <- glm_all$primary

message("  Gold*_R3 (n>=3, no proximity) …")
glm_r3   <- run_glm(make_gr(goldstar_R3),      "Gold*_R3 (n>=3, no prox)")
res_r3   <- glm_r3$primary

message("  Gold*_R3_strong (n>=3, |HP|>=0.20) …")
glm_r3s  <- run_glm(make_gr(goldstar_R3_strong),"Gold*_R3+|HP|>=0.20")
res_r3s  <- glm_r3s$primary

message("  Gold*_R5 (n>=5, no proximity) …")
glm_r5   <- run_glm(make_gr(goldstar_R5),      "Gold*_R5 (n>=5, no prox)")
res_r5   <- glm_r5$primary

results <- rbind(res_gold, res_all, res_r3, res_r3s, res_r5)

# Multi-feature (SegDup/LAD/B-comp, adjusted + marginal) table — feeds the
# 3-population x fragility-feature figure (viz/v1/figS9). Only the recurrence-
# only sets (All / Gold*_R3 / Gold*_R5) are the ones actually cited in the
# manuscript's recurrence-gradient claim (Fig 4B); Gold(original)/Gold*_R3_strong
# are included too for completeness/context.
multifeature_results <- do.call(rbind, list(
  glm_gold$multifeature, glm_all$multifeature,
  glm_r3$multifeature,   glm_r3s$multifeature,  glm_r5$multifeature
))
cat("\n=== A_GOLD1: multi-feature OR table (adjusted + marginal) ===\n")
print(multifeature_results[multifeature_results$label %in%
        c("All aDMR (no filter)", "Gold*_R3 (n>=3, no prox)", "Gold*_R5 (n>=5, no prox)"), ] |>
  transform(OR = round(OR, 3), CI_lo = round(CI_lo, 3), CI_hi = round(CI_hi, 3)))
fwrite(multifeature_results, file.path(OUT_DIR, "a_gold1_multifeature_or_table.csv"))
message("Wrote: a_gold1_multifeature_or_table.csv")

cat("\n=== A_GOLD1: SegDup OR by locus set ===\n")
print(results |>
  mutate(across(c(OR, CI_lo, CI_hi), ~round(., 3)),
         p = signif(p, 3),
         pct_segdup_case = round(pct_segdup_case, 1)) |>
  select(label, n_case, pct_segdup_case, OR, CI_lo, CI_hi, p, sig))

# 3. Wald test: Gold vs Gold* (does OR differ?) ================================
cat("\n=== Wald tests: Gold vs Gold* ===\n")
wald_vs_gold <- function(res_star, label_star) {
  z  <- (res_gold$beta - res_star$beta) / sqrt(res_gold$se^2 + res_star$se^2)
  p  <- 2 * pnorm(-abs(z))
  cat(sprintf("Gold vs %-40s: Δβ=%.3f  z=%.3f  p=%.4f  %s\n",
              label_star, res_gold$beta - res_star$beta, z, p,
              ifelse(p < 0.05, "DIFFERENT", "equivalent")))
}
wald_vs_gold(res_all,  "All aDMR")
wald_vs_gold(res_r3,   "Gold*_R3")
wald_vs_gold(res_r3s,  "Gold*_R3+|HP|>=0.20")
wald_vs_gold(res_r5,   "Gold*_R5")

# Proximity-driven fraction: % of Gold OR that disappears in Gold*_R3
pct_drop <- 100 * (1 - (log(res_r3$OR) / log(res_gold$OR)))
cat(sprintf("\nβ attenuation in Gold*_R3 vs Gold: %.1f%%\n", pct_drop))
cat(sprintf("  Gold  β=%.3f  OR=%.2f\n", res_gold$beta, res_gold$OR))
cat(sprintf("  R3    β=%.3f  OR=%.2f\n", res_r3$beta,   res_r3$OR))
cat(sprintf("  All   β=%.3f  OR=%.2f\n", res_all$beta,  res_all$OR))

# 4. Save ======================================================================
fwrite(results, file.path(OUT_DIR, "a_gold1_or_table.csv"))
message("Wrote: a_gold1_or_table.csv")

# 5. Figure ====================================================================
label_order <- c(
  "Gold (original, bp<=50kb)",
  "Gold*_R5 (n>=5, no prox)",
  "Gold*_R3+|HP|>=0.20",
  "Gold*_R3 (n>=3, no prox)",
  "All aDMR (no filter)"
)

plot_df <- results |>
  mutate(
    label  = factor(label, levels = label_order),
    color  = case_when(
      grepl("original", label)   ~ "Gold (proximity-defined)",
      grepl("All aDMR",  label)  ~ "All aDMR baseline",
      TRUE                       ~ "Gold* (proximity-free)"
    )
  ) |>
  filter(!is.na(label))

proxy_OR_inflation <- sprintf(
  "Gold proximity-defined OR: %.2f\nGold*_R3 (recurrence only): %.2f\nβ attenuation: %.0f%%",
  res_gold$OR, res_r3$OR, pct_drop)

p_forest <- ggplot(plot_df,
                   aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                       y = label, colour = color)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), width = 0.3) +
  geom_point(aes(size = n_case)) +
  geom_text(aes(x = CI_hi * 1.08,
                label = sprintf("%.2f %s\n(n=%d)", OR, sig, n_case)),
            hjust = 0, size = 3) +
  scale_colour_manual(values = c(
    "Gold (proximity-defined)" = "#d4a017",
    "Gold* (proximity-free)"   = "#2171b5",
    "All aDMR baseline"        = "#666666"
  )) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  scale_x_log10() +
  expand_limits(x = max(plot_df$CI_hi, na.rm = TRUE) * 1.8) +
  labs(
    x      = "SegDup OR (log scale)",
    y      = NULL,
    colour = NULL,
    title  = "A_GOLD1: Does Gold OR=3.80 reflect proximity selection bias?",
    subtitle = proxy_OR_inflation,
    caption = paste(
      "Model: is_aDMR ~ segdup + lad + b_compartment (10x matched controls).",
      "Gold* removes bp_dist<=50kb proximity criterion.",
      "If Gold* OR ≈ All-aDMR OR, Gold enrichment was proximity-driven.",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, family = "mono"))

# Panel B: % SegDup by locus set (raw proportions)
p_bar <- ggplot(plot_df,
                aes(x = label, y = pct_segdup_case, fill = color)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = results$pct_segdup_ctrl[1],
             linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.1f%%", pct_segdup_case)),
            vjust = -0.4, size = 3.5) +
  annotate("text", x = 0.5, y = results$pct_segdup_ctrl[1] + 0.3,
           label = sprintf("Control: %.1f%%", results$pct_segdup_ctrl[1]),
           hjust = 0, colour = "grey40", size = 3) +
  scale_fill_manual(values = c(
    "Gold (proximity-defined)" = "#d4a017",
    "Gold* (proximity-free)"   = "#2171b5",
    "All aDMR baseline"        = "#666666"
  ), guide = "none") +
  scale_x_discrete(labels = function(x)
    gsub(" \\(", "\n(", x)) +
  labs(x = NULL, y = "% overlapping SegDup",
       title = "Raw SegDup overlap by locus set") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8))

combined <- p_forest / p_bar + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG_DIR, "fig_agold1_decircularize.png"),
       combined, width = 10, height = 9, dpi = 150)
message("Saved: fig_agold1_decircularize.png")

cat("\n=== A_GOLD1 DONE ===\n")
cat(sprintf("Gold (proximity-defined) OR = %.2f [%.2f-%.2f] %s\n",
            res_gold$OR, res_gold$CI_lo, res_gold$CI_hi, res_gold$sig))
cat(sprintf("Gold*_R3 (no proximity)  OR = %.2f [%.2f-%.2f] %s\n",
            res_r3$OR,   res_r3$CI_lo,   res_r3$CI_hi,   res_r3$sig))
cat(sprintf("All aDMR (baseline)      OR = %.2f [%.2f-%.2f] %s\n",
            res_all$OR,  res_all$CI_lo,  res_all$CI_hi,  res_all$sig))
cat(sprintf("Beta attenuation (Gold → Gold*_R3): %.1f%%\n", pct_drop))
