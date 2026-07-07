#!/usr/bin/env Rscript
# P0-D: aDMR-own fragility enrichment (SV-independent)
#
# Design:
#   Test whether Gold+Silver aDMR loci are THEMSELVES enriched at fragility-prone
#   regions (SegDup/LAD/B-comp), independent of SV proximity.
#   This completes the "shared fragility" two-axis story:
#     Axis 1 (C7, C11): SV breakpoints are enriched at fragility loci
#     Axis 2 (this): aDMR loci are enriched at the SAME fragility features
#
#   Control: matched random regions (same chr, same width distribution × 1000 permutations)
#   Test: Fisher's exact test per feature + multivariate logistic regression
#
# Output:
#   result/admr_fragility_enrichment.csv    per-feature OR, CI, p
#   result/figures/fig_admr_fragility.png   grouped bar + forest plot

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

set.seed(123)

GOLD_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILV_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/silver_tier.csv")
SEGDUP    <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD       <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW    <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
FAI       <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_CTRL_MULT <- 10   # matched random regions per aDMR

# 1. Load unique Gold+Silver aDMR coordinates ==================================
message("Loading Gold + Silver aDMR coordinates...")
admr_cols <- c("tier_class", "admr_chr", "admr_start", "admr_end")
gold <- fread(GOLD_FILE, data.table = FALSE) |> select(all_of(admr_cols))
silv <- fread(SILV_FILE, data.table = FALSE) |> select(all_of(admr_cols))

admr <- bind_rows(gold, silv) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr)) |>
  distinct(admr_chr, admr_start, admr_end, .keep_all = TRUE)

cat(sprintf("Unique aDMR loci: %d (Gold=%d, Silver=%d)\n",
            nrow(admr),
            sum(admr$tier_class == "Gold"),
            sum(admr$tier_class == "Silver")))

admr_gr <- GRanges(
  seqnames = admr$admr_chr,
  ranges   = IRanges(admr$admr_start, admr$admr_end),
  tier     = admr$tier_class,
  is_admr  = 1L
)
admr_widths <- width(admr_gr)

# 2. Generate matched random regions (same chr + size distribution) ============
message("Generating matched random control regions...")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z"), data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

ctrl_list <- lapply(seq_len(nrow(admr)), function(i) {
  chr <- as.character(seqnames(admr_gr)[i])
  w   <- admr_widths[i]
  len <- chrom_sizes$len[chrom_sizes$chr == chr]
  if (length(len) == 0) return(NULL)
  n_ctrl <- N_CTRL_MULT
  max_start <- len - w
  if (max_start < 1) return(NULL)
  starts <- sample.int(max_start, size = n_ctrl, replace = TRUE)
  GRanges(seqnames = chr,
          ranges   = IRanges(starts, starts + w - 1L),
          tier     = "Control",
          is_admr  = 0L)
})
ctrl_gr <- do.call(c, Filter(Negate(is.null), ctrl_list))
cat(sprintf("Controls: %d (%.0f× aDMRs)\n", length(ctrl_gr), length(ctrl_gr)/length(admr_gr)))

all_gr <- c(admr_gr, ctrl_gr)

# 3. Annotate with fragility features ==========================================
message("Loading and annotating fragility features...")

segdup_gr <- import(SEGDUP, format = "BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr    <- import(LAD,    format = "BED"); seqlevelsStyle(lad_gr)    <- "UCSC"

all_gr$segdup <- overlapsAny(all_gr, segdup_gr)
all_gr$lad    <- overlapsAny(all_gr, lad_gr)

message("Extracting PC1 from BigWig...")
bw       <- BigWigFile(PC1_BW)
pc1_vals <- summary(bw, which = all_gr, type = "mean", defaultValue = NA_real_)
all_gr$pc1 <- unlist(lapply(pc1_vals, function(x) {
  if (length(x$score) == 0) NA_real_ else x$score[1]
}))
all_gr$b_compartment <- !is.na(all_gr$pc1) & all_gr$pc1 < 0

cat(sprintf("Annotated %d regions (%d aDMR + %d ctrl)\n",
            length(all_gr), sum(all_gr$is_admr), sum(!all_gr$is_admr)))

# 4. Per-feature Fisher test ===================================================
df <- data.frame(
  is_admr       = all_gr$is_admr,
  tier          = all_gr$tier,
  segdup        = as.integer(all_gr$segdup),
  lad           = as.integer(all_gr$lad),
  b_compartment = as.integer(all_gr$b_compartment),
  pc1           = all_gr$pc1
)
df$obs_weight <- ifelse(df$is_admr == 1L, 1, 1 / N_CTRL_MULT)

fisher_feature <- function(feat, label, df_in = df) {
  tab <- table(is_admr = df_in$is_admr, feat = df_in[[feat]])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  ft  <- fisher.test(tab, simulate.p.value = TRUE, B = 5000)
  p   <- ft$p.value
  data.frame(
    feature  = label,
    n_admr   = sum(df_in$is_admr == 1),
    n_ctrl   = sum(df_in$is_admr == 0),
    pct_admr = 100 * mean(df_in[[feat]][df_in$is_admr == 1], na.rm = TRUE),
    pct_ctrl = 100 * mean(df_in[[feat]][df_in$is_admr == 0], na.rm = TRUE),
    OR       = ft$estimate,
    CI_lo    = ft$conf.int[1],
    CI_hi    = ft$conf.int[2],
    p        = p,
    sig      = cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","ns"))
  )
}

fisher_res <- bind_rows(
  fisher_feature("segdup",        "SegDup overlap"),
  fisher_feature("lad",           "LAD overlap"),
  fisher_feature("b_compartment", "B-compartment")
) |> mutate(
  p_fdr  = p.adjust(p, method = "BH"),
  sig    = cut(p_fdr, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","ns")),
  OR_lab = sprintf("OR=%.2f [%.2f–%.2f] %s", OR, CI_lo, CI_hi, sig)
)

# Per-tier breakdown (Gold vs Silver)
cat(sprintf("df tier distribution: %s\n",
            paste(capture.output(table(df$tier, df$is_admr)), collapse=" | ")))

tier_res <- bind_rows(
  lapply(c("Gold", "Silver"), function(t) {
    df_admr_t <- df[df$is_admr == 1L & !is.na(df$tier) & df$tier == t, ]
    df_ctrl   <- df[df$is_admr == 0L, ]
    message(sprintf("Tier %s: %d aDMR, %d ctrl", t, nrow(df_admr_t), nrow(df_ctrl)))
    df_t2 <- rbind(df_admr_t, df_ctrl)
    res <- bind_rows(
      fisher_feature("segdup",        "SegDup", df_t2),
      fisher_feature("lad",           "LAD",    df_t2),
      fisher_feature("b_compartment", "B-comp", df_t2)
    )
    if (nrow(res) == 0) return(NULL)
    res |> mutate(tier = t)
  })
)
if (nrow(tier_res) > 0) {
  tier_res <- tier_res |>
    mutate(
      p_fdr = p.adjust(p, method = "BH"),
      sig   = cut(p_fdr, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","ns"))
    )
}

cat("\n=== Per-feature Fisher (All aDMR vs random) ===\n")
print(fisher_res |> select(feature, pct_admr, pct_ctrl, OR, CI_lo, CI_hi, p, p_fdr, sig))

cat("\n=== Per-tier breakdown ===\n")
if (nrow(tier_res) > 0) {
  print(tier_res |> select(tier, feature, pct_admr, pct_ctrl, OR, p, p_fdr, sig))
} else {
  cat("(no tier results — check df$tier values above)\n")
}

# 5. Multivariate logistic regression on aDMRs =================================
message("Fitting logistic regression...")
df_cc <- df |> filter(!is.na(pc1))
m_admr <- glm(is_admr ~ segdup + lad + b_compartment,
              data = df, weights = obs_weight, family = binomial())
m_full <- glm(is_admr ~ segdup + lad + b_compartment,
              data = df_cc, weights = obs_weight, family = binomial())

multi_res <- lapply(list(list(m_admr, df, "Multivariate"),
                          list(m_full, df_cc, "Multivariate (cc)")), function(x) {
  m <- x[[1]]; lab <- x[[3]]
  co <- summary(m)$coefficients
  ci <- confint.default(m)
  terms <- rownames(co)[-1]
  data.frame(feature = terms, model = lab,
             OR = exp(co[terms,1]), CI_lo = exp(ci[terms,1]), CI_hi = exp(ci[terms,2]),
             p  = co[terms,4])
}) |> bind_rows() |>
  mutate(sig = cut(p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")))

cat("\n=== Multivariate logistic (aDMR vs random) ===\n")
print(multi_res |> select(feature, model, OR, CI_lo, CI_hi, p, sig))

# 6. Save results ==============================================================
out <- bind_rows(
  fisher_res |> mutate(model = "Univariate Fisher", tier = "All"),
  tier_res   |> mutate(model = "Univariate Fisher"),
  multi_res  |>
    mutate(feature = recode(feature, segdup="SegDup overlap",
                            lad="LAD overlap", b_compartment="B-compartment"),
           tier = "All", n_admr = NA_real_, n_ctrl = NA_real_,
           pct_admr = NA_real_, pct_ctrl = NA_real_)
)
fwrite(out, file.path(OUT_DIR, "admr_fragility_enrichment.csv"))
message("Wrote: admr_fragility_enrichment.csv")

# 7. Figures ===================================================================
feat_levels <- c("SegDup overlap", "LAD overlap", "B-compartment")

# Panel A: grouped bar — % overlap for aDMR vs control
bar_df <- bind_rows(
  fisher_res |> select(feature, pct = pct_admr, sig) |> mutate(group = "Gold+Silver aDMR"),
  fisher_res |> select(feature, pct = pct_ctrl, sig) |> mutate(group = "Random control")
) |> mutate(feature = factor(feature, levels = feat_levels))

pA <- ggplot(bar_df, aes(x = feature, y = pct, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.65, alpha = 0.85) +
  geom_text(
    data = fisher_res |> mutate(feature = factor(feature, levels = feat_levels),
                                 y_pos = pmax(pct_admr, pct_ctrl) + 2),
    aes(x = feature, y = y_pos,
        label = paste0("OR=", round(OR,2), "\n", sig),
        fill = NULL),
    size = 3, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Gold+Silver aDMR" = "#d73027",
                                "Random control"   = "#6baed6")) +
  labs(
    title    = "aDMR Fragility Enrichment vs Random Genome",
    subtitle = sprintf("Gold+Silver aDMR (n=%d) vs matched random regions (n=%d)",
                       sum(df$is_admr), sum(!df$is_admr)),
    x = NULL, y = "% overlapping feature", fill = NULL,
    caption = "Fisher's exact test (simulated p, B=5000)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Panel B: forest (multivariate aDMR)
forest_df <- multi_res |>
  filter(model == "Multivariate") |>
  mutate(feature_lab = recode(feature, segdup="SegDup overlap",
                               lad="LAD overlap", b_compartment="B-compartment"),
         feature_lab = factor(feature_lab, levels = rev(feat_levels)))

pB <- ggplot(forest_df, aes(x = OR, y = feature_lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, color = "#d73027") +
  geom_point(size = 4, color = "#d73027") +
  geom_text(aes(x = CI_hi * 1.08,
                label = sprintf("OR=%.2f\n%s", OR, sig)),
            hjust = 0, size = 3.2) +
  scale_x_log10() +
  expand_limits(x = max(forest_df$CI_hi) * 1.5) +
  labs(
    title    = "Multivariate: Independent Fragility Effects\non aDMR loci",
    x        = "Odds Ratio (log scale)",
    y        = NULL,
    caption  = "Logistic regression: is_aDMR ~ segdup + lad + b_compartment"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel C: Gold vs Silver tier comparison (SegDup OR)
tier_plot <- tier_res |>
  mutate(feature = factor(feature, levels = c("SegDup","LAD","B-comp")))

pC <- ggplot(tier_plot, aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                              y = feature, color = tier, shape = tier)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3.5, position = position_dodge(0.5)) +
  geom_text(aes(x = CI_hi * 1.1, label = sig),
            position = position_dodge(0.5), hjust = 0, size = 3.2) +
  scale_color_manual(values = c("Gold" = "#f4a460", "Silver" = "#708090")) +
  scale_x_log10() +
  labs(
    title   = "Fragility Enrichment by Tier",
    x       = "Odds Ratio (log scale)",
    y       = NULL, color = "Tier", shape = "Tier"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

combined <- (pA | pB) / pC + plot_layout(heights = c(1.4, 1))
ggsave(file.path(FIG_DIR, "fig_admr_fragility.png"), combined,
       width = 12, height = 9, dpi = 150)
message("Saved: fig_admr_fragility.png")

cat("\nDone.\n")
