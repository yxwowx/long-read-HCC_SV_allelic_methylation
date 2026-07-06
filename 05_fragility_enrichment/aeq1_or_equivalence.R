#!/usr/bin/env Rscript
# A_EQ1: Formal OR-equivalence test for SegDup enrichment: SV vs. aDMR
#
# Core claim: OR(SV→SegDup) ≈ OR(aDMR→SegDup) supports shared fragility.
# Existing C13/C14 use non-parallel covariate sets. This script builds
# fully parallel models at three levels of adjustment, then runs a
# Wald equivalence test (H0: β_SV = β_aDMR) at each level.
#
# Levels:
#   L1 (base):      segdup + lad + b_compartment           [C13 binary = C14 base]
#   L2 (+ repeat):  segdup + lad + b_compartment + repeat_density   [C13 full, new aDMR]
#   L3 (+ CpG):     segdup + lad + b_compartment + repeat_density + log10_cpg [aDMR only]
#
# Wald test: z = (β_SV - β_aDMR) / sqrt(SE_SV² + SE_aDMR²); p = 2*pnorm(-|z|)
# H0: β_SV = β_aDMR (equivalence); H1: β_SV ≠ β_aDMR
# If p > 0.05: ORs are statistically equivalent → shared fragility supported.
#
# Run: mamba run -n renv Rscript post_processing/aeq1_or_equivalence.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ── Paths ──────────────────────────────────────────────────────────────────────
SV_FRAG   <- "/node200data/kachungk/hcc_data/DMR_SVs/result/sv_fragility_annotation.csv"
GOLD_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/gold_tier_final.csv"
SILV_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/silver_tier.csv"
SEGDUP    <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
LAD       <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed"
RMSK      <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/rmsk.bed"
PC1_BW    <- "/node200data/kachungk/reference/GRCh38/3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw"
FAI       <- "/node200data/kachungk/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
OUT_DIR   <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

SV_CTRL_MULT   <- 5L   # random regions per SV breakpoint  (matches C13)
ADM_CTRL_MULT  <- 10L  # matched random regions per aDMR    (matches C14/C22)

WINDOW_BP <- 5000L     # repeat-density window radius (matches fragility_multivariate.R)

# TOST equivalence margin on log-OR scale: ±log(1.5) treats any OR ratio within
# [1/1.5, 1.5] as practically equivalent (a defensible "small effect" bound).
DELTA <- log(1.5)  # ≈ 0.405

# ── Helper: extract OR + SE from GLM ──────────────────────────────────────────
extract_or_se <- function(model, term = "segdup") {
  co <- summary(model)$coefficients
  ci <- confint.default(model)
  list(
    beta  = co[term, 1],
    se    = co[term, 2],
    OR    = exp(co[term, 1]),
    CI_lo = exp(ci[term, 1]),
    CI_hi = exp(ci[term, 2]),
    p     = co[term, 4]
  )
}

# ── 0. Load shared reference files ────────────────────────────────────────────
message("Loading reference annotation files …")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z")) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

segdup_gr <- import(SEGDUP, format = "BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr    <- import(LAD,    format = "BED"); seqlevelsStyle(lad_gr)    <- "UCSC"
rmsk_gr   <- import(RMSK,   format = "BED"); seqlevelsStyle(rmsk_gr)   <- "UCSC"
bw        <- BigWigFile(PC1_BW)
bsgenome  <- BSgenome.Hsapiens.UCSC.hg38

annotate_gr <- function(gr) {
  gr$segdup        <- overlapsAny(gr, segdup_gr)
  gr$lad           <- overlapsAny(gr, lad_gr)
  windows          <- GRanges(seqnames(gr),
                              IRanges(pmax(1L, start(gr) - WINDOW_BP),
                                      end(gr) + WINDOW_BP))
  gr$repeat_density <- countOverlaps(windows, rmsk_gr)
  pc1_vals <- summary(bw, which = gr, type = "mean", defaultValue = NA_real_)
  gr$pc1   <- unlist(lapply(pc1_vals, function(x)
    if (length(x$score) == 0) NA_real_ else x$score[1]))
  gr$b_compartment <- !is.na(gr$pc1) & gr$pc1 < 0
  gr
}

# ── 1. SV dataset ──────────────────────────────────────────────────────────────
message("Building SV case–control dataset …")
sv_raw <- fread(SV_FRAG) |>
  filter(!is.na(seqnames), grepl("^chr[0-9XY]+$", seqnames)) |>
  distinct(bp_id, .keep_all = TRUE)
message(sprintf("  SV breakpoints: %d", nrow(sv_raw)))

sv_case_gr <- GRanges(sv_raw$seqnames, IRanges(sv_raw$start, sv_raw$start),
                      is_sv = 1L)

sv_ctrl_list <- lapply(names(table(seqnames(sv_case_gr))), function(chr) {
  n   <- table(seqnames(sv_case_gr))[[chr]] * SV_CTRL_MULT
  len <- chrom_sizes$len[chrom_sizes$chr == chr]
  if (!length(len) || !n) return(NULL)
  pos <- sample.int(len - 1L, n, replace = TRUE)
  GRanges(chr, IRanges(pos, pos), is_sv = 0L)
})
sv_ctrl_gr <- do.call(c, Filter(Negate(is.null), sv_ctrl_list))

sv_all_gr <- annotate_gr(c(sv_case_gr, sv_ctrl_gr))

# CpG density for SV breakpoints (±500 bp window)
message("  Counting CpGs around SV breakpoints (BSgenome) …")
sv_windows    <- GRanges(seqnames(sv_all_gr),
                          IRanges(pmax(1L, start(sv_all_gr) - 500L),
                                  end(sv_all_gr) + 500L))
sv_seqs       <- getSeq(bsgenome, sv_windows)
sv_cpg_cnt    <- vcountPattern("CG", sv_seqs)
sv_all_gr$nCG_dens  <- sv_cpg_cnt / 1.0  # per 1 kb window → per 1000bp
sv_all_gr$log10_cpg <- log10(sv_all_gr$nCG_dens + 1)

sv_df <- data.frame(
  event         = "SV",
  is_event      = sv_all_gr$is_sv,
  segdup        = as.integer(sv_all_gr$segdup),
  lad           = as.integer(sv_all_gr$lad),
  b_compartment = as.integer(sv_all_gr$b_compartment),
  repeat_density = sv_all_gr$repeat_density,
  log10_cpg     = sv_all_gr$log10_cpg
)
sv_df$obs_weight <- ifelse(sv_df$is_event == 1L, 1, 1 / SV_CTRL_MULT)
message(sprintf("  SV df: %d rows (case=%d, ctrl=%d)",
                nrow(sv_df), sum(sv_df$is_event), sum(!sv_df$is_event)))

# ── 2. aDMR dataset ────────────────────────────────────────────────────────────
message("Building aDMR case–control dataset …")
admr_cols <- c("tier_class","admr_chr","admr_start","admr_end","nCG")
admr_raw  <- bind_rows(
  fread(GOLD_FILE) |> select(all_of(admr_cols)),
  fread(SILV_FILE) |> select(all_of(admr_cols))
) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr)) |>
  distinct(admr_chr, admr_start, admr_end, .keep_all = TRUE)

admr_raw$width_kb  <- (admr_raw$admr_end - admr_raw$admr_start + 1) / 1000
admr_raw$nCG_dens  <- admr_raw$nCG / pmax(admr_raw$width_kb, 0.1)
admr_raw$log10_cpg <- log10(admr_raw$nCG_dens + 1)

message(sprintf("  aDMR loci: %d (Gold=%d, Silver=%d)",
                nrow(admr_raw),
                sum(admr_raw$tier_class == "Gold"),
                sum(admr_raw$tier_class == "Silver")))

admr_case_gr <- GRanges(
  admr_raw$admr_chr,
  IRanges(admr_raw$admr_start, admr_raw$admr_end),
  is_admr   = 1L,
  nCG_dens  = admr_raw$nCG_dens,
  log10_cpg = admr_raw$log10_cpg
)

admr_ctrl_list <- lapply(seq_len(nrow(admr_raw)), function(i) {
  chr <- admr_raw$admr_chr[i]
  w   <- admr_raw$admr_end[i] - admr_raw$admr_start[i] + 1L
  len <- chrom_sizes$len[chrom_sizes$chr == chr]
  if (!length(len)) return(NULL)
  max_s <- len - w
  if (max_s < 1) return(NULL)
  ss <- sample.int(max_s, ADM_CTRL_MULT, replace = FALSE)
  GRanges(chr, IRanges(ss, ss + w - 1L), is_admr = 0L,
          nCG_dens = NA_real_, log10_cpg = NA_real_)
})
admr_ctrl_gr <- do.call(c, Filter(Negate(is.null), admr_ctrl_list))

# CpG density for controls via BSgenome
message("  Counting CpGs in control regions (BSgenome) …")
ctrl_seqs            <- getSeq(bsgenome, admr_ctrl_gr)
ctrl_cpg_cnt         <- vcountPattern("CG", ctrl_seqs)
admr_ctrl_gr$nCG_dens  <- ctrl_cpg_cnt / (width(admr_ctrl_gr) / 1000)
admr_ctrl_gr$log10_cpg <- log10(admr_ctrl_gr$nCG_dens + 1)

admr_all_gr <- annotate_gr(c(admr_case_gr, admr_ctrl_gr))

admr_df <- data.frame(
  event         = "aDMR",
  is_event      = admr_all_gr$is_admr,
  segdup        = as.integer(admr_all_gr$segdup),
  lad           = as.integer(admr_all_gr$lad),
  b_compartment = as.integer(admr_all_gr$b_compartment),
  repeat_density = admr_all_gr$repeat_density,
  log10_cpg     = admr_all_gr$log10_cpg
)
admr_df$obs_weight <- ifelse(admr_df$is_event == 1L, 1, 1 / ADM_CTRL_MULT)
message(sprintf("  aDMR df: %d rows (case=%d, ctrl=%d)",
                nrow(admr_df), sum(admr_df$is_event), sum(!admr_df$is_event)))

# ── 3. Fit models at each level ────────────────────────────────────────────────
message("Fitting parallel GLMs …")

fit_model <- function(df, formula_str, label) {
  m <- glm(as.formula(formula_str), data = df,
           weights = obs_weight, family = binomial())
  r <- extract_or_se(m, "segdup")
  r$label   <- label
  r$n_case  <- sum(df$is_event)
  r$n_ctrl  <- sum(!df$is_event)
  r$formula <- formula_str
  r
}

# Level 1: base (same as C13 binary / C14)
sv_L1   <- fit_model(sv_df,   "is_event ~ segdup + lad + b_compartment", "SV — L1 base")
adm_L1  <- fit_model(admr_df, "is_event ~ segdup + lad + b_compartment", "aDMR — L1 base")

# Level 2: + repeat_density
sv_L2   <- fit_model(sv_df,   "is_event ~ segdup + lad + b_compartment + repeat_density",
                     "SV — L2 +repeat")
adm_L2  <- fit_model(admr_df, "is_event ~ segdup + lad + b_compartment + repeat_density",
                     "aDMR — L2 +repeat")

# Level 3: + log10_cpg (aDMR only; SV includes it for symmetry)
sv_L3   <- fit_model(sv_df   |> filter(!is.na(log10_cpg)),
                     "is_event ~ segdup + lad + b_compartment + repeat_density + log10_cpg",
                     "SV — L3 +repeat+CpG")
adm_L3  <- fit_model(admr_df |> filter(!is.na(log10_cpg)),
                     "is_event ~ segdup + lad + b_compartment + repeat_density + log10_cpg",
                     "aDMR — L3 +repeat+CpG")

# ── 4. Wald equivalence tests ─────────────────────────────────────────────────
message("\nRunning Wald equivalence tests …")

wald_test <- function(sv_res, adm_res, level_label) {
  beta_diff <- sv_res$beta - adm_res$beta
  SE_diff   <- sqrt(sv_res$se^2 + adm_res$se^2)

  # Wald difference test (H0: β_SV = β_aDMR)
  z <- beta_diff / SE_diff
  p <- 2 * pnorm(-abs(z))

  # TOST equivalence test (H0: |β_SV − β_aDMR| ≥ DELTA; H1: difference < DELTA)
  # Equivalence established (reject H0) if both one-sided p-values < 0.05.
  tost_p_lower <- pnorm((beta_diff + DELTA) / SE_diff, lower.tail = FALSE)
  tost_p_upper <- pnorm((beta_diff - DELTA) / SE_diff, lower.tail = TRUE)
  tost_p       <- max(tost_p_lower, tost_p_upper)

  data.frame(
    level          = level_label,
    sv_OR          = sv_res$OR,
    sv_CI_lo       = sv_res$CI_lo,
    sv_CI_hi       = sv_res$CI_hi,
    sv_p           = sv_res$p,
    admr_OR        = adm_res$OR,
    admr_CI_lo     = adm_res$CI_lo,
    admr_CI_hi     = adm_res$CI_hi,
    admr_p         = adm_res$p,
    beta_diff      = beta_diff,
    z_equiv        = z,
    p_equiv        = p,
    verdict        = ifelse(p > 0.05, "no significant difference (p>0.05)",
                            "significant difference (p≤0.05)"),
    tost_margin    = DELTA,
    tost_p_lower   = tost_p_lower,
    tost_p_upper   = tost_p_upper,
    tost_p         = tost_p,
    tost_verdict   = ifelse(tost_p < 0.05,
                            sprintf("EQUIVALENT within ±OR%.1f margin (TOST p=%.3f)",
                                    exp(DELTA), tost_p),
                            sprintf("not established within ±OR%.1f margin (TOST p=%.3f)",
                                    exp(DELTA), tost_p))
  )
}

eq_tests <- bind_rows(
  wald_test(sv_L1, adm_L1, "L1: +lad+b_comp"),
  wald_test(sv_L2, adm_L2, "L2: +lad+b_comp+repeat"),
  wald_test(sv_L3, adm_L3, "L3: +lad+b_comp+repeat+CpG")
)

cat("\n=== A_EQ1: Wald difference + TOST equivalence tests for SegDup OR ===\n")
print(eq_tests |>
  mutate(across(c(sv_OR, admr_OR), ~round(., 3)),
         across(c(beta_diff, z_equiv), ~round(., 4)),
         p_equiv = signif(p_equiv, 3),
         tost_p  = signif(tost_p,  3)) |>
  select(level, sv_OR, admr_OR, beta_diff, z_equiv, p_equiv, verdict,
         tost_p, tost_verdict))

# ── 5. Summary OR table ───────────────────────────────────────────────────────
cat("\n=== Full OR table ===\n")
all_results <- bind_rows(
  lapply(list(sv_L1, adm_L1, sv_L2, adm_L2, sv_L3, adm_L3), function(r) {
    data.frame(
      label    = r$label,
      OR       = r$OR,
      CI_lo    = r$CI_lo,
      CI_hi    = r$CI_hi,
      p        = r$p,
      sig      = cut(r$p, c(-Inf,.001,.01,.05,Inf), labels=c("***","**","*","ns")),
      n_case   = r$n_case,
      formula  = r$formula
    )
  })
)
print(all_results |> mutate(across(c(OR,CI_lo,CI_hi), ~round(.,3)),
                             p = signif(p,3)) |>
  select(label, OR, CI_lo, CI_hi, p, sig))

# ── 6. Save ───────────────────────────────────────────────────────────────────
fwrite(eq_tests,    file.path(OUT_DIR, "a_eq1_wald_tests.csv"))
fwrite(all_results, file.path(OUT_DIR, "a_eq1_or_table.csv"))
message("Wrote: a_eq1_wald_tests.csv, a_eq1_or_table.csv")

# ── 7. Figure: forest + equivalence summary ───────────────────────────────────
level_order <- c("L1: +lad+b_comp", "L2: +lad+b_comp+repeat", "L3: +lad+b_comp+repeat+CpG")

plot_df2 <- rbind(
  data.frame(level     = factor(eq_tests$level, levels = level_order),
             event_lab = "SV breakpoints",
             OR        = eq_tests$sv_OR,
             CI_lo     = eq_tests$sv_CI_lo,
             CI_hi     = eq_tests$sv_CI_hi,
             sig       = cut(eq_tests$sv_p, c(-Inf,.001,.01,.05,Inf),
                             labels = c("***","**","*","ns"))),
  data.frame(level     = factor(eq_tests$level, levels = level_order),
             event_lab = "aDMR loci",
             OR        = eq_tests$admr_OR,
             CI_lo     = eq_tests$admr_CI_lo,
             CI_hi     = eq_tests$admr_CI_hi,
             sig       = cut(eq_tests$admr_p, c(-Inf,.001,.01,.05,Inf),
                             labels = c("***","**","*","ns")))
)

p_forest <- ggplot(plot_df2,
                   aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                       y = level, colour = event_lab, shape = event_lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi),
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3.5, position = position_dodge(0.5)) +
  geom_text(aes(x = CI_hi * 1.08, label = paste0(round(OR,2), " ", sig)),
            position = position_dodge(0.5), hjust = 0, size = 3) +
  scale_colour_manual(values = c("SV breakpoints" = "#2171b5",
                                  "aDMR loci"       = "#cb181d")) +
  scale_x_log10() +
  expand_limits(x = max(plot_df2$CI_hi, na.rm=TRUE) * 1.6) +
  labs(
    x = "SegDup OR (log scale)", y = "Model (covariates added)",
    colour = NULL, shape = NULL,
    title = "A_EQ1: SegDup OR — SV vs. aDMR (parallel models)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Equivalence panel: Wald difference p + TOST p side-by-side
eq_long <- rbind(
  data.frame(level = factor(eq_tests$level, levels = level_order),
             test  = "Wald difference\n(H0: β_SV = β_aDMR)",
             p     = eq_tests$p_equiv,
             sig   = eq_tests$p_equiv > 0.05),   # TRUE = green (no diff)
  data.frame(level = factor(eq_tests$level, levels = level_order),
             test  = sprintf("TOST equivalence\n(margin ±OR%.1f)", exp(DELTA)),
             p     = eq_tests$tost_p,
             sig   = eq_tests$tost_p  < 0.05)    # TRUE = green (equivalent)
)
p_equiv <- ggplot(eq_long,
                  aes(x = factor(level, levels = level_order), y = p,
                      fill = sig)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
  geom_text(aes(label = sprintf("p=%.3f", p)),
            position = position_dodge(0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("FALSE" = "#fb6a4a", "TRUE" = "#74c476"),
                    guide = "none") +
  scale_x_discrete(labels = c("L1", "L2", "L3")) +
  scale_y_continuous(limits = c(0, max(eq_long$p, 0.06) * 1.5)) +
  facet_wrap(~test, nrow = 1) +
  labs(
    x = "Model level", y = "p-value",
    title = "Difference test (left) and TOST equivalence test (right)",
    caption = paste0("Wald: green = no significant difference (p>0.05).\n",
                     "TOST: green = equivalence established (p<0.05) ",
                     sprintf("within ±OR%.1f (Δβ=%.3f) margin.", exp(DELTA), DELTA))
  ) +
  theme_bw(base_size = 12)

combined_fig <- p_forest / p_equiv + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG_DIR, "fig_aeq1_equivalence.png"),
       combined_fig, width = 9, height = 8, dpi = 150)
message("Saved: fig_aeq1_equivalence.png")

cat("\n=== A_EQ1 DONE ===\n")
for (i in seq_len(nrow(eq_tests))) {
  cat(sprintf("  %s: SV OR=%.2f  aDMR OR=%.2f  Δβ=%.4f  z=%.3f  p=%.4f  → %s\n",
              eq_tests$level[i], eq_tests$sv_OR[i], eq_tests$admr_OR[i],
              eq_tests$beta_diff[i], eq_tests$z_equiv[i],
              eq_tests$p_equiv[i], eq_tests$verdict[i]))
}
