#!/usr/bin/env Rscript
# P1-D: Normal liver inter-sample β-variance at Gold/Silver aDMR loci
#
# Purpose:
#   Test if Gold/Silver aDMR loci show elevated CpG methylation variability
#   in NORMAL liver tissue (TCGA-LIHC paired normals, n=50).
#   If TRUE → constitutive methylation instability (C8) is inherent to these
#   loci, not cancer-induced.
#
# Data: TCGA-LIHC 450K normal samples (n=50) cached from V1 analysis.
# Design:
#   For each aDMR locus: find overlapping HM450 probes (±2kb window)
#   Compute inter-sample variance across 50 normals
#   Compare: Gold/Silver aDMR probes vs matched-background probes (Wilcoxon)
#
# Run: mamba run -n renv Rscript post_processing/admr_normal_tissue_variance.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(SummarizedExperiment)
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
METH_CACHE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/external_validation_cache/LIHC_meth450.rds")
GOLD_CSV   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILVER_CSV <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/silver_tier.csv")
OUT_DIR    <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR    <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

PROBE_WINDOW_BP <- 2000L   # ±2kb around aDMR center to match probes

# 1. Load cached TCGA normal methylation =======================================
message("Loading TCGA-LIHC 450K (cached)...")
meth_se  <- readRDS(METH_CACHE)
meth_mat <- assay(meth_se)          # probes × samples
col_data <- colData(meth_se)

# Extract normal samples using sample_type column (populated by TCGAbiolinks)
if ("sample_type" %in% names(col_data)) {
  normal_idx <- which(grepl("Normal|normal", col_data$sample_type))
} else if ("shortLetterCode" %in% names(col_data)) {
  normal_idx <- which(col_data$shortLetterCode == "NT")
} else {
  # Fall back: barcodes with type code 10/11
  normal_idx <- which(substr(rownames(col_data), 14, 15) %in% c("10", "11"))
}
cat(sprintf("Normal samples: %d / %d total\n", length(normal_idx), ncol(meth_mat)))

if (length(normal_idx) < 5) stop("Too few normal samples found")

norm_mat <- meth_mat[, normal_idx, drop = FALSE]
cat(sprintf("Normal matrix: %d probes × %d samples\n", nrow(norm_mat), ncol(norm_mat)))

# Filter probes with >20% missingness across normals
miss_frac <- rowMeans(is.na(norm_mat))
keep      <- miss_frac < 0.2
norm_mat  <- norm_mat[keep, ]
cat(sprintf("After missingness filter: %d probes\n", nrow(norm_mat)))

# Per-probe variance across normals
probe_var <- matrixStats::rowVars(norm_mat, na.rm = TRUE)
probe_sd  <- sqrt(probe_var)

# Probe genomic coordinates
probe_gr_full <- rowRanges(meth_se)
probe_gr      <- probe_gr_full[keep]
seqlevelsStyle(probe_gr) <- "UCSC"

cat(sprintf("Probe variance computed: median=%.4f, IQR=[%.4f, %.4f]\n",
            median(probe_sd, na.rm = TRUE),
            quantile(probe_sd, 0.25, na.rm = TRUE),
            quantile(probe_sd, 0.75, na.rm = TRUE)))

# 2. Load aDMR loci ============================================================
message("\nLoading aDMR loci...")
load_admr <- function(csv, tier_name) {
  dt  <- fread(csv, data.table = FALSE)
  chr <- if ("admr_chr" %in% names(dt)) dt$admr_chr else dt$chr
  s   <- if ("admr_start" %in% names(dt)) dt$admr_start else dt$start
  e   <- if ("admr_end" %in% names(dt)) dt$admr_end else dt$end
  gr  <- GRanges(chr, IRanges(s, e), tier = tier_name)
  gr[grepl("^chr[0-9XY]+$", as.character(seqnames(gr)))]
}

gold_gr   <- load_admr(GOLD_CSV,   "Gold")
silver_gr <- load_admr(SILVER_CSV, "Silver")

# Deduplicate by locus (same admr_chr:admr_start-admr_end across patients)
gold_gr   <- unique(gold_gr)
silver_gr <- unique(silver_gr)
cat(sprintf("Unique Gold aDMR: %d | Silver: %d\n",
            length(gold_gr), length(silver_gr)))

# 3. Map probes to aDMR loci (±2kb window) =====================================
message("Matching probes to aDMR loci...")

# Expand aDMR loci by ±PROBE_WINDOW_BP
expand_gr <- function(gr, bp) {
  resize(gr, width = width(gr) + 2 * bp, fix = "center") |> trim()
}

gold_win   <- expand_gr(gold_gr,   PROBE_WINDOW_BP)
silver_win <- expand_gr(silver_gr, PROBE_WINDOW_BP)

# Probe annotation: is_gold, is_silver, or background
probe_is_gold   <- overlapsAny(probe_gr, gold_win)
probe_is_silver <- overlapsAny(probe_gr, silver_win) & !probe_is_gold

cat(sprintf("Probes overlapping Gold aDMR (±%dkb): %d\n",
            PROBE_WINDOW_BP / 1000, sum(probe_is_gold)))
cat(sprintf("Probes overlapping Silver aDMR (±%dkb): %d\n",
            PROBE_WINDOW_BP / 1000, sum(probe_is_silver)))

# 4. Per-aDMR locus: mean variance of matched probes ===========================
# For each Gold/Silver locus: collect probe SDs → use median as locus-level estimate
locus_variance <- function(admr_win_gr, label) {
  lapply(seq_len(length(admr_win_gr)), function(i) {
    hits  <- subjectHits(findOverlaps(admr_win_gr[i], probe_gr))
    if (length(hits) == 0) return(NULL)
    data.frame(
      locus = i,
      tier  = label,
      n_probes   = length(hits),
      median_sd  = median(probe_sd[hits], na.rm = TRUE),
      mean_sd    = mean(probe_sd[hits], na.rm = TRUE),
      max_sd     = max(probe_sd[hits], na.rm = TRUE)
    )
  }) |> bind_rows()
}

message("Computing per-locus variance for Gold aDMR...")
gold_var   <- locus_variance(gold_win,   "Gold")
message("Computing per-locus variance for Silver aDMR...")
silver_var <- locus_variance(silver_win, "Silver")

cat(sprintf("Gold loci with probes: %d / %d\n",
            nrow(gold_var), length(gold_gr)))
cat(sprintf("Silver loci with probes: %d / %d\n",
            nrow(silver_var), length(silver_gr)))

# 5. Background: matched random loci ===========================================
message("Generating matched background loci...")
# Sample random probes (same number as Gold, with replacement for 1000 permutations)
bg_probes_idx <- which(!probe_is_gold & !probe_is_silver)

# Null distribution: draw n_gold random probes, compute their SD
n_gold   <- nrow(gold_var)
n_perm   <- 1000L
null_sds <- replicate(n_perm, {
  idx <- sample(bg_probes_idx, n_gold, replace = TRUE)
  median(probe_sd[idx], na.rm = TRUE)
})

null_df <- data.frame(tier = "Background (permuted)", median_sd = null_sds)

# 6. Statistical tests =========================================================
message("Running statistical tests...")

all_probe_df <- data.frame(
  probe_sd = probe_sd,
  group    = case_when(
    probe_is_gold   ~ "Gold aDMR",
    probe_is_silver ~ "Silver aDMR",
    TRUE            ~ "Background"
  )
)

wt_gold   <- wilcox.test(
  probe_sd[probe_is_gold],
  probe_sd[bg_probes_idx],
  alternative = "greater"
)
wt_silver <- wilcox.test(
  probe_sd[probe_is_silver],
  probe_sd[bg_probes_idx],
  alternative = "greater"
)

# Effect size: median fold-change
med_gold   <- median(probe_sd[probe_is_gold],   na.rm = TRUE)
med_silver <- median(probe_sd[probe_is_silver], na.rm = TRUE)
med_bg     <- median(probe_sd[bg_probes_idx],   na.rm = TRUE)

cat(sprintf("\n=== P1-D Results ===\n"))
cat(sprintf("Median SD — Gold: %.4f | Silver: %.4f | Background: %.4f\n",
            med_gold, med_silver, med_bg))
cat(sprintf("Gold vs Background: FC=%.2fx, Wilcoxon p=%.3g\n",
            med_gold / med_bg, wt_gold$p.value))
cat(sprintf("Silver vs Background: FC=%.2fx, Wilcoxon p=%.3g\n",
            med_silver / med_bg, wt_silver$p.value))

# Empirical p from permutation
perm_p_gold <- mean(null_sds >= med_gold)
cat(sprintf("Gold permutation p (1000 draws): %.4f\n", perm_p_gold))

# Save results
result_df <- data.frame(
  group         = c("Gold aDMR", "Silver aDMR", "Background"),
  n_probes      = c(sum(probe_is_gold), sum(probe_is_silver), length(bg_probes_idx)),
  median_sd     = c(med_gold, med_silver, med_bg),
  fc_vs_bg      = c(med_gold / med_bg, med_silver / med_bg, 1.0),
  wilcox_p      = c(wt_gold$p.value, wt_silver$p.value, NA_real_),
  perm_p        = c(perm_p_gold, NA_real_, NA_real_)
)
fwrite(result_df, file.path(OUT_DIR, "admr_normal_variance.csv"))
message("Saved: admr_normal_variance.csv")

# 7. Figures ===================================================================
message("\nGenerating figures...")
theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

sig_label <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# Panel A: probe-level SD density by group (subsample background for speed)
bg_sample_idx <- sample(bg_probes_idx, min(50000L, length(bg_probes_idx)))
dens_df <- bind_rows(
  data.frame(sd = probe_sd[probe_is_gold],   group = "Gold aDMR"),
  data.frame(sd = probe_sd[probe_is_silver], group = "Silver aDMR"),
  data.frame(sd = probe_sd[bg_sample_idx],   group = "Background")
) |>
  filter(!is.na(sd)) |>
  mutate(group = factor(group, levels = c("Gold aDMR", "Silver aDMR", "Background")))

p_density <- ggplot(dens_df, aes(x = sd, color = group, fill = group)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(
    data = result_df |> filter(group != "Background"),
    aes(xintercept = median_sd, color = group),
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_color_manual(values = c("Gold aDMR" = "#E24B4A",
                                "Silver aDMR" = "#E07B39",
                                "Background"  = "#888780")) +
  scale_fill_manual(values  = c("Gold aDMR" = "#E24B4A",
                                "Silver aDMR" = "#E07B39",
                                "Background"  = "#888780")) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("Gold vs BG: FC=%.2fx, p=%s\nSilver vs BG: FC=%.2fx, p=%s",
                           med_gold / med_bg,   sig_label(wt_gold$p.value),
                           med_silver / med_bg, sig_label(wt_silver$p.value))) +
  labs(
    title    = "P1-D: Normal Liver CpG Methylation Variability",
    subtitle = sprintf("TCGA-LIHC normal liver (n=%d); probes ±%dkb of aDMR loci",
                       length(normal_idx), PROBE_WINDOW_BP / 1000),
    x        = "Inter-sample SD (β-value, 50 normal livers)",
    y        = "Density",
    color    = NULL, fill = NULL
  ) +
  theme_hcc

# Panel B: per-locus median SD boxplot
locus_df <- bind_rows(gold_var, silver_var) |>
  filter(!is.na(median_sd)) |>
  mutate(tier = factor(tier, levels = c("Gold", "Silver")))

bg_locus_df <- data.frame(
  locus = seq_len(n_perm),
  tier  = "Background",
  median_sd = null_sds
)

all_locus_df <- bind_rows(
  locus_df |> select(tier, median_sd),
  bg_locus_df |> select(tier, median_sd)
) |> mutate(tier = factor(tier, levels = c("Gold", "Silver", "Background")))

p_box <- ggplot(all_locus_df, aes(x = tier, y = median_sd, fill = tier)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1) +
  geom_jitter(data = locus_df, width = 0.15, alpha = 0.4, size = 1.5,
              aes(color = tier)) +
  scale_fill_manual(values  = c("Gold" = "#E24B4A", "Silver" = "#E07B39",
                                "Background" = "#888780")) +
  scale_color_manual(values = c("Gold" = "#C0392B", "Silver" = "#C0612E")) +
  labs(
    title = "Per-locus Median SD across 50 Normal Livers",
    x     = NULL,
    y     = "Median SD (β-value)",
    fill  = NULL
  ) +
  theme_hcc + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_s9a_admr_normal_variance_density.png"),
       p_density, width = 8, height = 5, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_s9b_admr_normal_variance_locus.png"),
       p_box,     width = 6, height = 5, dpi = 150)
message("Saved: fig_s9a + fig_s9b")

message("\nDone.")
