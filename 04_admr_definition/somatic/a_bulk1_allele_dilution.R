#!/usr/bin/env Rscript
# A_BULK1 — Allele dilution / long-read detection power
#
# Tests whether the SegDup co-enrichment signal is allele-specific (only
# detectable by long-read phased 5mC) or also visible in simulated bulk
# sequencing (HP1 and HP2 averaged). Directly resolves [Mo4] and supports
# the thesis title framing "Allele-Specific ... Using Long-Read Sequencing."
#
# Design: same phased aDMR locus universe, two instability definitions:
#   AS arm   : |HP1.Methy - HP2.Methy|_max >= tau  (long-read only)
#   Bulk arm : |(HP1+HP2)/2 - normal.Methy|_max >= tau  (bulk-simulatable)
# Threshold sweep tau in {0.10, 0.15, 0.20, 0.25, 0.30}; matched random controls.
#
# Usage: mamba run -n renv Rscript post_processing/a_bulk1_allele_dilution.R
#
# Outputs -> result/:
#   a_bulk1_or_comparison.csv      headline OR/CI/p at tau=0.15
#   a_bulk1_threshold_sweep.csv    OR_AS vs OR_bulk across tau sweep
#   a_bulk1_locus_overlap.csv      AS-only / both / bulk-only / neither (tau=0.15)
#   figures/fig_a_bulk1.png        barplot + ribbon sweep plot

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# ── Paths ─────────────────────────────────────────────────────────────────────
CONF_DMR <- "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz"
SEGDUP   <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
LAD      <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed"
PC1_BW   <- "/node200data/kachungk/reference/GRCh38/3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw"
FAI      <- "/node200data/kachungk/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
OUT_DIR  <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_CTRL_MULT  <- 10
THRESHOLDS   <- c(0.10, 0.15, 0.20, 0.25, 0.30)
TAU_HEADLINE <- 0.15
C17_OR       <- 1.05   # TCGA-LIHC 450K reference OR from C17

# ── 1. Load and aggregate phased aDMR loci ────────────────────────────────────
message("Loading confident_dmr_per_patient.csv.gz ...")
raw <- fread(CONF_DMR, data.table = FALSE)
cat(sprintf("Raw rows: %d  unique patients: %d\n",
            nrow(raw), length(unique(raw$patient_code))))

raw <- raw |>
  filter(
    !is.na(HP1.Methy), !is.na(HP2.Methy), !is.na(normal.Methy),
    grepl("^chr[0-9XY]+$", admr_chr),
    !is.na(admr_start), !is.na(admr_end),
    admr_end > admr_start
  )

# One row per unique phased aDMR locus; max instability across patients
loci <- raw |>
  mutate(
    as_metric   = abs(HP1.Methy - HP2.Methy),
    bulk_metric = abs((HP1.Methy + HP2.Methy) / 2 - normal.Methy)
  ) |>
  group_by(admr_chr, admr_start, admr_end) |>
  summarise(
    n_patients_obs    = max(n_patients, na.rm = TRUE),
    hp_abs_diff_max   = max(as_metric,   na.rm = TRUE),
    bulk_abs_diff_max = max(bulk_metric, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("Unique phased aDMR loci: %d\n", nrow(loci)))
cat(sprintf("  Median |HP1-HP2|_max:         %.3f\n", median(loci$hp_abs_diff_max)))
cat(sprintf("  Median |(sim_bulk-normal)|_max: %.3f\n", median(loci$bulk_abs_diff_max)))

# ── 2. Preload reference GRanges and chrom sizes ──────────────────────────────
message("Loading SegDup, LAD, and chrom sizes ...")
segdup_gr   <- import(SEGDUP, format = "BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr      <- import(LAD,    format = "BED"); seqlevelsStyle(lad_gr)    <- "UCSC"
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z"), data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

# ── 3. Helpers ────────────────────────────────────────────────────────────────

# Build GRanges for a locus subset + matched random controls, annotate SegDup+LAD
build_annotated <- function(loci_sub, seed_offset = 0) {
  set.seed(123 + seed_offset)
  case_gr <- GRanges(
    seqnames = loci_sub$admr_chr,
    ranges   = IRanges(loci_sub$admr_start, loci_sub$admr_end),
    is_case  = 1L
  )
  widths <- width(case_gr)
  ctrl_list <- lapply(seq_len(length(case_gr)), function(i) {
    chr <- as.character(seqnames(case_gr)[i])
    w   <- widths[i]
    len <- chrom_sizes$len[chrom_sizes$chr == chr]
    if (length(len) == 0) return(NULL)
    max_start <- len - w
    if (max_start < 1) return(NULL)
    starts <- sample.int(max_start, size = N_CTRL_MULT, replace = TRUE)
    GRanges(seqnames = chr, ranges = IRanges(starts, starts + w - 1L), is_case = 0L)
  })
  ctrl_gr <- do.call(c, Filter(Negate(is.null), ctrl_list))
  all_gr  <- c(case_gr, ctrl_gr)
  all_gr$segdup <- overlapsAny(all_gr, segdup_gr)
  all_gr$lad    <- overlapsAny(all_gr, lad_gr)
  all_gr
}

# Fisher's exact test for one binary feature
fisher_feature <- function(feat, label, df_in) {
  tab <- table(is_case = df_in$is_case, feat = df_in[[feat]])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  ft <- fisher.test(tab, simulate.p.value = TRUE, B = 5000)
  data.frame(
    feature  = label,
    n_case   = sum(df_in$is_case == 1),
    n_ctrl   = sum(df_in$is_case == 0),
    pct_case = 100 * mean(df_in[[feat]][df_in$is_case == 1], na.rm = TRUE),
    pct_ctrl = 100 * mean(df_in[[feat]][df_in$is_case == 0], na.rm = TRUE),
    OR       = ft$estimate,
    CI_lo    = ft$conf.int[1],
    CI_hi    = ft$conf.int[2],
    p        = ft$p.value
  )
}

run_enrichment <- function(all_gr, arm_label, tau_val) {
  df <- data.frame(
    is_case = all_gr$is_case,
    segdup  = as.integer(all_gr$segdup),
    lad     = as.integer(all_gr$lad)
  )
  df$obs_weight <- ifelse(df$is_case == 1L, 1, 1 / N_CTRL_MULT)
  res <- bind_rows(
    fisher_feature("segdup", "SegDup", df),
    fisher_feature("lad",    "LAD",    df)
  ) |> mutate(arm = arm_label, tau = tau_val,
              sig = cut(p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")))
  res
}

# ── 4. Threshold sweep ────────────────────────────────────────────────────────
message("Running threshold sweep (AS vs bulk) ...")
sweep_list <- list()

for (i_tau in seq_along(THRESHOLDS)) {
  tau <- THRESHOLDS[i_tau]
  message(sprintf("  tau=%.2f", tau))

  loci_as   <- loci |> filter(hp_abs_diff_max   >= tau)
  loci_bulk <- loci |> filter(bulk_abs_diff_max  >= tau)

  cat(sprintf("    AS loci = %d   Bulk loci = %d\n", nrow(loci_as), nrow(loci_bulk)))

  if (nrow(loci_as) >= 3) {
    gr_as <- build_annotated(loci_as, seed_offset = i_tau * 10)
    sweep_list[[paste0("as_",   tau)]] <- run_enrichment(gr_as, "AS (long-read)", tau)
  }
  if (nrow(loci_bulk) >= 3) {
    gr_bk <- build_annotated(loci_bulk, seed_offset = i_tau * 10 + 5)
    sweep_list[[paste0("bulk_", tau)]] <- run_enrichment(gr_bk, "Bulk-simulated", tau)
  }
}

sweep_df <- bind_rows(sweep_list) |>
  mutate(arm = factor(arm, levels = c("AS (long-read)", "Bulk-simulated")))

cat("\n=== SegDup Fisher: threshold sweep ===\n")
print(sweep_df |> filter(feature == "SegDup") |>
      select(arm, tau, n_case, pct_case, pct_ctrl, OR, CI_lo, CI_hi, p, sig))

# ── 5. Multivariate GLM at tau=0.15 (PC1 extracted once) ─────────────────────
message("Extracting PC1 for multivariate GLM at tau=0.15 ...")

run_glm <- function(loci_sub, arm_label, seed_offset) {
  set.seed(123 + seed_offset)
  all_gr <- build_annotated(loci_sub, seed_offset = seed_offset)
  bw       <- BigWigFile(PC1_BW)
  pc1_vals <- summary(bw, which = all_gr, type = "mean", defaultValue = NA_real_)
  all_gr$pc1 <- unlist(lapply(pc1_vals, function(x)
    if (length(x$score) == 0) NA_real_ else x$score[1]
  ))
  all_gr$b_compartment <- !is.na(all_gr$pc1) & all_gr$pc1 < 0

  df <- data.frame(
    is_case       = all_gr$is_case,
    segdup        = as.integer(all_gr$segdup),
    lad           = as.integer(all_gr$lad),
    b_compartment = as.integer(all_gr$b_compartment)
  )
  df$obs_weight <- ifelse(df$is_case == 1L, 1, 1 / N_CTRL_MULT)
  df_cc <- df |> filter(!is.na(b_compartment))

  m <- glm(is_case ~ segdup + lad + b_compartment,
            data = df_cc, weights = obs_weight, family = binomial())
  co <- summary(m)$coefficients
  ci <- confint.default(m)
  terms <- rownames(co)[-1]
  data.frame(
    arm     = arm_label, tau = TAU_HEADLINE,
    feature = terms,
    OR      = exp(co[terms, 1]),
    CI_lo   = exp(ci[terms, 1]),
    CI_hi   = exp(ci[terms, 2]),
    p       = co[terms, 4],
    model   = "Multivariate GLM",
    n_case  = sum(df_cc$is_case == 1),
    n_ctrl  = sum(df_cc$is_case == 0)
  ) |> mutate(sig = cut(p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")))
}

loci_as_hl   <- loci |> filter(hp_abs_diff_max   >= TAU_HEADLINE)
loci_bulk_hl <- loci |> filter(bulk_abs_diff_max  >= TAU_HEADLINE)

glm_as   <- run_glm(loci_as_hl,   "AS (long-read)", seed_offset = 99)
glm_bulk <- run_glm(loci_bulk_hl, "Bulk-simulated", seed_offset = 100)
glm_df   <- bind_rows(glm_as, glm_bulk)

cat("\n=== Multivariate GLM at tau=0.15 ===\n")
print(glm_df |> select(arm, feature, n_case, OR, CI_lo, CI_hi, p, sig))

# ── 6. Headline table ─────────────────────────────────────────────────────────
headline_df <- bind_rows(
  sweep_df |> filter(tau == TAU_HEADLINE) |>
    mutate(model = "Fisher", n_case = as.integer(n_case), n_ctrl = as.integer(n_ctrl)),
  glm_df |>
    mutate(pct_case = NA_real_, pct_ctrl = NA_real_,
           n_case = as.integer(n_case), n_ctrl = as.integer(n_ctrl))
)

cat("\n=== Headline (tau=0.15) ===\n")
print(headline_df |> select(arm, model, feature, n_case, OR, CI_lo, CI_hi, p, sig))

# ── 7. Locus overlap at tau=0.15 ─────────────────────────────────────────────
loci_hl <- loci |>
  mutate(
    is_as   = hp_abs_diff_max   >= TAU_HEADLINE,
    is_bulk = bulk_abs_diff_max >= TAU_HEADLINE
  )

overlap_tab <- data.frame(
  category = c("AS-only", "Both AS+Bulk", "Bulk-only", "Neither"),
  n        = c(
    sum( loci_hl$is_as & !loci_hl$is_bulk),
    sum( loci_hl$is_as &  loci_hl$is_bulk),
    sum(!loci_hl$is_as &  loci_hl$is_bulk),
    sum(!loci_hl$is_as & !loci_hl$is_bulk)
  )
) |>
  mutate(
    pct          = round(100 * n / sum(n), 1),
    tau          = TAU_HEADLINE,
    total_loci   = nrow(loci_hl)
  )

cat(sprintf("\n=== Locus overlap at tau=%.2f ===\n", TAU_HEADLINE))
print(overlap_tab)

pct_as_lost <- round(
  100 * overlap_tab$n[overlap_tab$category == "AS-only"] /
    (overlap_tab$n[overlap_tab$category == "AS-only"] +
     overlap_tab$n[overlap_tab$category == "Both AS+Bulk"]),
  1)
cat(sprintf("\n%% of AS-instability loci invisible to bulk: %.1f%%\n", pct_as_lost))

# ── 8. Save ───────────────────────────────────────────────────────────────────
fwrite(headline_df, file.path(OUT_DIR, "a_bulk1_or_comparison.csv"))
fwrite(sweep_df,    file.path(OUT_DIR, "a_bulk1_threshold_sweep.csv"))
fwrite(overlap_tab, file.path(OUT_DIR, "a_bulk1_locus_overlap.csv"))
fwrite(glm_df,      file.path(OUT_DIR, "a_bulk1_glm.csv"))
message("Wrote: a_bulk1_or_comparison.csv, a_bulk1_threshold_sweep.csv, a_bulk1_locus_overlap.csv, a_bulk1_glm.csv")

# ── 9. Figure ─────────────────────────────────────────────────────────────────
message("Generating figures ...")

seg_headline <- sweep_df |>
  filter(tau == TAU_HEADLINE, feature == "SegDup") |>
  mutate(sig_lab = sprintf("OR=%.2f\n[%.2f–%.2f]\n%s", OR, CI_lo, CI_hi, sig))

pA <- ggplot(seg_headline, aes(x = arm, y = OR, ymin = CI_lo, ymax = CI_hi, color = arm)) +
  geom_hline(yintercept = 1,       linetype = "dashed",  color = "grey50", linewidth = 0.8) +
  geom_hline(yintercept = C17_OR,  linetype = "dotted",  color = "#984ea3", linewidth = 0.9) +
  geom_errorbar(width = 0.20, linewidth = 1.2) +
  geom_point(size = 5) +
  geom_text(aes(label = sig_lab), vjust = -0.6, size = 3.3, color = "black") +
  annotate("text", x = 2.48, y = C17_OR + 0.06,
           label = sprintf("TCGA 450K (C17)\nOR=%.2f", C17_OR),
           hjust = 1, size = 2.9, color = "#984ea3") +
  scale_color_manual(values = c("AS (long-read)" = "#d73027",
                                 "Bulk-simulated"  = "#4292c6")) +
  scale_y_continuous(limits = c(0.8, NA), expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    title    = "SegDup enrichment: allele-specific vs bulk-simulated",
    subtitle = sprintf("Phased aDMR loci | tau=%.2f | %d× matched random controls", TAU_HEADLINE, N_CTRL_MULT),
    x = NULL, y = "Odds Ratio (SegDup overlap)",
    caption  = "Fisher's exact test (B=5000); dotted = TCGA-LIHC 450K ref (C17)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

sweep_seg <- sweep_df |>
  filter(feature == "SegDup") |>
  mutate(arm = factor(arm, levels = c("AS (long-read)", "Bulk-simulated")))

pB <- ggplot(sweep_seg, aes(x = tau, y = OR, ymin = CI_lo, ymax = CI_hi,
                              color = arm, group = arm)) +
  geom_hline(yintercept = 1,      linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = C17_OR, linetype = "dotted", color = "#984ea3") +
  geom_ribbon(aes(fill = arm), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("AS (long-read)" = "#d73027", "Bulk-simulated" = "#4292c6")) +
  scale_fill_manual( values = c("AS (long-read)" = "#d73027", "Bulk-simulated" = "#4292c6")) +
  scale_x_continuous(breaks = THRESHOLDS, labels = sprintf("%.2f", THRESHOLDS)) +
  labs(
    title   = "Threshold sweep: AS vs bulk SegDup OR",
    x = "Threshold tau", y = "OR (SegDup)", color = NULL, fill = NULL,
    caption = "Ribbon = 95% CI; dotted = TCGA 450K (C17)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

combined <- pA | pB
ggsave(file.path(FIG_DIR, "fig_a_bulk1.png"), combined,
       width = 12, height = 5, dpi = 150)
message("Saved: figures/fig_a_bulk1.png")

cat("\nDone.\n")
