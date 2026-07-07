#!/usr/bin/env Rscript
# Replication timing analysis: PC1 (MicroC) vs ENCODE Repli-seq, aDMR/SV loci
#
# Sections:
#   1. Load query loci (Gold/Silver aDMR, SV breakpoints, random background)
#   2. Extract PC1 scores at all loci (hg38 bigWig, already in sv_fragility_annotation)
#   3. Load & liftOver ENCODE Repli-seq S1-S4 (hg19→hg38), compute RT score
#   4. Extract Repli-seq RT scores at query loci
#   5. Group-level comparisons (violin) — PC1 and RT
#   6. PC1 vs Repli-seq concordance (scatter + correlation)
#   7. B2 probe-level interaction: beta_SD ~ PC1/RT × is_admr
#   8. Output tables and figures
#
# Run: mamba run -n renv Rscript post_processing/replication_timing_admr.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(rtracklayer)
  library(matrixStats)
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
BASE       <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")
REPLISEQ_DIR <- file.path(BASE, "external_validation_cache/repliseq")
PC1_BW     <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
CHAIN_FILE <- Sys.getenv("LIFTOVER_CHAIN")
GOLD_CSV   <- file.path(BASE, "04.final_candidate/gold_tier_final.csv")
SILVER_CSV <- file.path(BASE, "04.final_candidate/silver_tier.csv")
SV_ANN     <- file.path(BASE, "result/sv_fragility_annotation.csv")
METH_CACHE <- file.path(BASE, "external_validation_cache/LIHC_meth450.rds")
OUT_DIR    <- file.path(BASE, "result")
FIG_DIR    <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

PHASE_FILES <- list(
  S1  = file.path(REPLISEQ_DIR, "ENCFF001GPK_S1_hg19.bigWig"),
  S2  = file.path(REPLISEQ_DIR, "ENCFF001GPO_S2_hg19.bigWig"),
  S3  = file.path(REPLISEQ_DIR, "ENCFF001GPU_S3_hg19.bigWig"),
  S4  = file.path(REPLISEQ_DIR, "ENCFF001GPX_S4_hg19.bigWig")
)

missing_bw <- Filter(function(f) !file.exists(f), unlist(PHASE_FILES))
if (length(missing_bw) > 0) {
  stop("Missing Repli-seq bigWig(s):\n",
       paste(missing_bw, collapse = "\n"),
       "\nRun: bash post_processing/download_encode_repliseq.sh first")
}

theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))
sig_label <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# 1. Load query loci ===========================================================
message("=== Section 1: Loading query loci ===")

load_admr_gr <- function(csv, tier_name) {
  dt  <- fread(csv, data.table = FALSE)
  chr <- if ("admr_chr" %in% names(dt)) dt$admr_chr else dt$chr
  s   <- if ("admr_start" %in% names(dt)) dt$admr_start else dt$start
  e   <- if ("admr_end" %in% names(dt)) dt$admr_end else dt$end
  gr  <- GRanges(chr, IRanges(s, e), tier = tier_name)
  gr  <- gr[grepl("^chr[0-9XY]+$", as.character(seqnames(gr)))]
  unique(gr)
}

gold_gr   <- load_admr_gr(GOLD_CSV,   "Gold aDMR")
silver_gr <- load_admr_gr(SILVER_CSV, "Silver aDMR")
cat(sprintf("Gold aDMR loci: %d | Silver: %d\n", length(gold_gr), length(silver_gr)))

sv_dt <- fread(SV_ANN)
sv_nb <- sv_dt[sv_dt$tier_group == "Non-boundary", ]
sv_nb <- sv_nb[!is.na(sv_nb$seqnames) & !is.na(sv_nb$start), ]
sv_gr <- GRanges(sv_nb$seqnames, IRanges(sv_nb$start, sv_nb$start + 1L),
                 tier = "SV breakpoint", pc1_score = sv_nb$pc1_score)
cat(sprintf("Non-boundary SV breakpoints: %d\n", length(sv_gr)))

# Random background: 500 windows sampled from UCSC hg38 simple repeats as genome proxy
# Use Silver locus sizes as window sizes for matched sampling
set.seed(42)
hg38_autosomes <- paste0("chr", c(1:22))
window_size    <- median(width(silver_gr))
rand_gr_list   <- lapply(hg38_autosomes, function(chr) {
  n <- round(500 / 22)
  chr_len <- seqlengths(silver_gr)[chr]
  if (is.na(chr_len)) chr_len <- 200000000L
  starts <- sample(1:(chr_len - window_size - 1L), n, replace = FALSE)
  GRanges(chr, IRanges(starts, starts + window_size - 1L), tier = "Random")
})
rand_gr <- do.call(c, Filter(Negate(is.null), rand_gr_list))
cat(sprintf("Random background windows: %d\n", length(rand_gr)))

all_query <- c(gold_gr, silver_gr, sv_gr[, "tier"], rand_gr)
seqlevelsStyle(all_query) <- "UCSC"

# 2. Extract PC1 scores ========================================================
message("\n=== Section 2: Extracting PC1 scores ===")

extract_bw_mean <- function(query_gr, bw_path) {
  scores <- numeric(length(query_gr))
  for (i in seq_along(query_gr)) {
    vals <- tryCatch(
      import.bw(bw_path, which = query_gr[i])$score,
      error = function(e) numeric(0)
    )
    scores[i] <- if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA_real_
  }
  scores
}

message("Extracting PC1 at all query loci (may take a few minutes)...")
all_query$pc1 <- extract_bw_mean(all_query, PC1_BW)
cat(sprintf("PC1 extraction done. NAs: %d / %d\n",
            sum(is.na(all_query$pc1)), length(all_query)))

# 3. Load & liftOver ENCODE Repli-seq (hg19 -> hg38) ===========================
message("\n=== Section 3: Loading Repli-seq bigWigs and lifting to hg38 ===")

chain <- import.chain(CHAIN_FILE)

load_phase_hg38 <- function(bw_path, phase_name) {
  message("  Loading phase: ", phase_name)
  gr <- import.bw(bw_path)
  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- gr[as.character(seqnames(gr)) %in% paste0("chr", c(1:22, "X", "Y"))]
  lifted_list <- liftOver(gr, chain)
  lifted <- unlist(lifted_list)
  lifted <- lifted[!is.na(lifted$score)]
  message(sprintf("    hg19: %d bins → hg38: %d bins (%.1f%% mapped)",
                  length(gr), length(lifted), 100 * length(lifted) / length(gr)))
  lifted
}

phase_gr <- lapply(names(PHASE_FILES), function(p) load_phase_hg38(PHASE_FILES[[p]], p))
names(phase_gr) <- names(PHASE_FILES)

# 4. Extract Repli-seq RT score at query loci ==================================
message("\n=== Section 4: Extracting Repli-seq scores at query loci ===")

score_at_loci <- function(query_gr, phase_gr_list) {
  # Returns matrix: n_loci × n_phases
  mat <- matrix(NA_real_, nrow = length(query_gr), ncol = length(phase_gr_list),
                dimnames = list(NULL, names(phase_gr_list)))
  for (ph in names(phase_gr_list)) {
    hits <- findOverlaps(query_gr, phase_gr_list[[ph]], ignore.strand = TRUE)
    if (length(hits) == 0) next
    phase_scores <- tapply(
      phase_gr_list[[ph]]$score[subjectHits(hits)],
      queryHits(hits),
      mean,
      na.rm = TRUE
    )
    mat[as.integer(names(phase_scores)), ph] <- phase_scores
  }
  mat
}

message("Scoring query loci against 4 Repli-seq phases...")
score_mat <- score_at_loci(all_query, phase_gr)

# RT score: log2((S1+S2+0.01)/(S3+S4+0.01)); positive = early, negative = late
all_query$rt_early <- rowMeans(score_mat[, c("S1","S2"), drop=FALSE], na.rm=TRUE)
all_query$rt_late  <- rowMeans(score_mat[, c("S3","S4"), drop=FALSE], na.rm=TRUE)
all_query$rt_score <- log2((all_query$rt_early + 0.01) / (all_query$rt_late + 0.01))

cat(sprintf("RT scores computed. NAs: %d / %d\n",
            sum(is.na(all_query$rt_score)), length(all_query)))

# Build summary data.frame
loci_df <- data.frame(
  tier    = all_query$tier,
  pc1     = all_query$pc1,
  rt_score = all_query$rt_score,
  stringsAsFactors = FALSE
) |>
  mutate(
    tier = factor(tier, levels = c("Gold aDMR", "Silver aDMR", "SV breakpoint", "Random"))
  ) |>
  filter(!is.na(pc1) | !is.na(rt_score))

cat(sprintf("Summary table: %d loci\n", nrow(loci_df)))

# 5. Group-level comparisons ===================================================
message("\n=== Section 5: Group-level comparisons ===")

rand_pc1 <- loci_df$pc1[loci_df$tier == "Random"]
rand_rt  <- loci_df$rt_score[loci_df$tier == "Random"]

stat_summary <- loci_df |>
  group_by(tier) |>
  summarise(
    n          = n(),
    pc1_median = median(pc1, na.rm = TRUE),
    pc1_q25    = quantile(pc1, 0.25, na.rm = TRUE),
    pc1_q75    = quantile(pc1, 0.75, na.rm = TRUE),
    rt_median  = median(rt_score, na.rm = TRUE),
    rt_q25     = quantile(rt_score, 0.25, na.rm = TRUE),
    rt_q75     = quantile(rt_score, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

wilcox_results <- lapply(c("Gold aDMR", "Silver aDMR", "SV breakpoint"), function(g) {
  grp_pc1 <- loci_df$pc1[loci_df$tier == g]
  grp_rt  <- loci_df$rt_score[loci_df$tier == g]
  wt_pc1  <- suppressWarnings(wilcox.test(grp_pc1, rand_pc1, alternative = "less"))
  wt_rt   <- suppressWarnings(wilcox.test(grp_rt,  rand_rt,  alternative = "less"))
  data.frame(
    tier       = g,
    pc1_vs_rand_p  = wt_pc1$p.value,
    rt_vs_rand_p   = wt_rt$p.value
  )
}) |> bind_rows()

stat_table <- left_join(stat_summary, wilcox_results, by = "tier")
fwrite(stat_table, file.path(OUT_DIR, "replication_timing_group.csv"))
message("Saved: replication_timing_group.csv")

cat("\n=== Group-Level Results ===\n")
print(stat_table)

# Violin: PC1 by group =========================================================
tier_colors <- c("Gold aDMR" = "#E24B4A", "Silver aDMR" = "#E07B39",
                 "SV breakpoint" = "#4A90D9", "Random" = "#888780")

p_pc1_violin <- ggplot(
    loci_df |> filter(!is.na(pc1)),
    aes(x = tier, y = pc1, fill = tier)
  ) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = tier_colors) +
  labs(
    title    = "MicroC PC1 at aDMR / SV Loci vs Random",
    subtitle = "Negative PC1 = B-compartment / late-replicating",
    x        = NULL, y = "HepG2 PC1 score",
    fill     = NULL
  ) +
  theme_hcc + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_b2a_pc1_violin.png"),
       p_pc1_violin, width = 7, height = 5, dpi = 150)
message("Saved: fig_b2a_pc1_violin.png")

# Violin: Repli-seq RT score by group ==========================================
p_rt_violin <- ggplot(
    loci_df |> filter(!is.na(rt_score)),
    aes(x = tier, y = rt_score, fill = tier)
  ) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = tier_colors) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
    label = paste(
      apply(wilcox_results, 1, function(r)
        sprintf("%s vs Rand: p=%s (RT)", r[["tier"]], sig_label(as.numeric(r[["rt_vs_rand_p"]])))),
      collapse = "\n"
    )
  ) +
  labs(
    title    = "ENCODE Repli-seq RT Score at aDMR / SV Loci vs Random",
    subtitle = "HepG2 log2(S1+S2 / S3+S4); positive = earlier replication",
    x        = NULL, y = "RT score [log2 early/late]",
    fill     = NULL
  ) +
  theme_hcc + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_b2b_rt_repliseq_violin.png"),
       p_rt_violin, width = 7, height = 5, dpi = 150)
message("Saved: fig_b2b_rt_repliseq_violin.png")

# 6. PC1 vs Repli-seq concordance ==============================================
message("\n=== Section 6: PC1 vs Repli-seq concordance ===")

concordance_df <- loci_df |> filter(!is.na(pc1) & !is.na(rt_score))
cor_all <- cor.test(concordance_df$pc1, concordance_df$rt_score, method = "spearman")
cat(sprintf("PC1 vs RT score: Spearman rho=%.3f, p=%.3g (n=%d)\n",
            cor_all$estimate, cor_all$p.value, nrow(concordance_df)))

# Save concordance table
fwrite(concordance_df, file.path(OUT_DIR, "replication_timing_loci.csv"))
message("Saved: replication_timing_loci.csv")

p_concordance <- ggplot(
    concordance_df |> sample_n(min(5000, nrow(concordance_df))),
    aes(x = pc1, y = rt_score, color = tier)
  ) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = tier_colors) +
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
    label = sprintf("Spearman ρ=%.3f, p=%.2g", cor_all$estimate, cor_all$p.value)
  ) +
  labs(
    title    = "PC1 vs Repli-seq RT Score Concordance",
    subtitle = "HepG2 MicroC PC1 vs ENCODE2 Repli-seq at aDMR/SV/Random loci",
    x        = "MicroC PC1 (positive = A-comp/early)",
    y        = "Repli-seq RT score [log2 early/late]",
    color    = NULL
  ) +
  theme_hcc

ggsave(file.path(FIG_DIR, "fig_b2c_pc1_repliseq_concordance.png"),
       p_concordance, width = 7, height = 5, dpi = 150)
message("Saved: fig_b2c_pc1_repliseq_concordance.png")

# 7. B2 probe-level interaction model ==========================================
message("\n=== Section 7: B2 probe-level interaction (beta_SD ~ RT × is_admr) ===")

message("Loading TCGA 450K methylation SE...")
meth_se  <- readRDS(METH_CACHE)
meth_mat <- assay(meth_se)
col_data <- colData(meth_se)

if ("sample_type" %in% names(col_data)) {
  normal_idx <- which(grepl("Normal|normal", col_data$sample_type))
} else if ("shortLetterCode" %in% names(col_data)) {
  normal_idx <- which(col_data$shortLetterCode == "NT")
} else {
  normal_idx <- which(substr(rownames(col_data), 14, 15) %in% c("10", "11"))
}
cat(sprintf("Normal samples: %d\n", length(normal_idx)))

norm_mat  <- meth_mat[, normal_idx, drop = FALSE]
miss_frac <- rowMeans(is.na(norm_mat))
keep      <- miss_frac < 0.2
norm_mat  <- norm_mat[keep, ]
probe_var <- rowVars(norm_mat, na.rm = TRUE)
probe_sd  <- sqrt(probe_var)

probe_gr_full <- rowRanges(meth_se)
probe_gr      <- probe_gr_full[keep]
seqlevelsStyle(probe_gr) <- "UCSC"
cat(sprintf("Probes after missingness filter: %d\n", length(probe_gr)))

# Label probes by aDMR group
probe_is_gold   <- overlapsAny(probe_gr, gold_gr + 2000L)
probe_is_silver <- overlapsAny(probe_gr, silver_gr + 2000L) & !probe_is_gold
bg_idx          <- which(!probe_is_gold & !probe_is_silver)

# Sample: Gold (all), Silver (all), Background (10k)
bg_sample <- sample(bg_idx, min(10000L, length(bg_idx)))
use_idx   <- c(which(probe_is_gold), which(probe_is_silver), bg_sample)
cat(sprintf("Probe sample: Gold=%d, Silver=%d, Background=%d (total=%d)\n",
            sum(probe_is_gold), sum(probe_is_silver),
            length(bg_sample), length(use_idx)))

probe_df <- data.frame(
  beta_sd    = probe_sd[use_idx],
  is_admr    = as.integer(probe_is_gold[use_idx] | probe_is_silver[use_idx]),
  group      = dplyr::case_when(
    probe_is_gold[use_idx]   ~ "Gold aDMR",
    probe_is_silver[use_idx] ~ "Silver aDMR",
    TRUE                     ~ "Background"
  )
)

# Extract PC1 and RT at probe locations
query_probe_gr <- probe_gr[use_idx]

message("Extracting PC1 at probe loci...")
probe_df$pc1 <- extract_bw_mean(query_probe_gr, PC1_BW)

message("Extracting Repli-seq scores at probe loci...")
probe_scores <- score_at_loci(query_probe_gr, phase_gr)
probe_df$rt_early <- rowMeans(probe_scores[, c("S1","S2"), drop=FALSE], na.rm=TRUE)
probe_df$rt_late  <- rowMeans(probe_scores[, c("S3","S4"), drop=FALSE], na.rm=TRUE)
probe_df$rt_score <- log2((probe_df$rt_early + 0.01) / (probe_df$rt_late + 0.01))

# Probe metadata: nCG and GC content from rowData if available
rd <- as.data.frame(rowData(meth_se))
if ("nCpG" %in% names(rd)) {
  probe_df$log10_cpg <- log10(rd[keep, ][use_idx, "nCpG"] + 1)
} else {
  probe_df$log10_cpg <- NA_real_
}

fwrite(probe_df, file.path(OUT_DIR, "replication_timing_probe.csv"))
message("Saved: replication_timing_probe.csv")

# Fit interaction models
probe_clean <- probe_df |>
  filter(!is.na(beta_sd) & !is.na(pc1) & !is.na(rt_score) & beta_sd > 0)

# Model A: PC1
if (!all(is.na(probe_clean$log10_cpg))) {
  lm_pc1 <- lm(log10(beta_sd + 0.001) ~ pc1 * is_admr + log10_cpg, data = probe_clean)
  lm_rt  <- lm(log10(beta_sd + 0.001) ~ rt_score * is_admr + log10_cpg, data = probe_clean)
} else {
  lm_pc1 <- lm(log10(beta_sd + 0.001) ~ pc1 * is_admr, data = probe_clean)
  lm_rt  <- lm(log10(beta_sd + 0.001) ~ rt_score * is_admr, data = probe_clean)
}

cat("\n--- Model A: beta_SD ~ PC1 × is_admr ---\n")
print(summary(lm_pc1)$coefficients)
cat(sprintf("Adjusted R²: %.4f\n", summary(lm_pc1)$adj.r.squared))

cat("\n--- Model B: beta_SD ~ RT_score × is_admr ---\n")
print(summary(lm_rt)$coefficients)
cat(sprintf("Adjusted R²: %.4f\n", summary(lm_rt)$adj.r.squared))

# Interaction coefficients
int_pc1 <- coef(summary(lm_pc1))["pc1:is_admr", ]
int_rt  <- coef(summary(lm_rt))["rt_score:is_admr", ]

cat(sprintf("\nInteraction pc1:is_admr — β=%.4f, p=%.3g\n",
            int_pc1["Estimate"], int_pc1["Pr(>|t|)"]))
cat(sprintf("Interaction rt_score:is_admr — β=%.4f, p=%.3g\n",
            int_rt["Estimate"], int_rt["Pr(>|t|)"]))

# Interaction figure ===========================================================
# Binned plot showing the interaction: RT bin × is_admr vs mean log10(beta_SD)
probe_clean$rt_bin  <- cut(probe_clean$rt_score, breaks = 5, labels = FALSE)
probe_clean$pc1_bin <- cut(probe_clean$pc1,      breaks = 5, labels = FALSE)

bin_summary <- probe_clean |>
  mutate(admr_grp = ifelse(is_admr == 1, "aDMR", "Background")) |>
  group_by(rt_bin, admr_grp) |>
  summarise(
    mean_log_sd = mean(log10(beta_sd + 0.001), na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  filter(!is.na(rt_bin))

p_interaction_rt <- ggplot(bin_summary, aes(x = rt_bin, y = mean_log_sd,
                                             color = admr_grp, group = admr_grp)) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(size = n), alpha = 0.8) +
  scale_color_manual(values = c("aDMR" = "#E24B4A", "Background" = "#888780")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
           label = sprintf("Interaction β=%.4f, p=%s",
                           int_rt["Estimate"], sig_label(int_rt["Pr(>|t|)"]))) +
  labs(
    title    = "B2: Replication Timing × aDMR Status Interaction",
    subtitle = "Repli-seq RT score (bin 1=late, 5=early) vs inter-sample β-SD",
    x        = "Replication timing quintile (1=late → 5=early)",
    y        = "Mean log10(β-SD + 0.001)",
    color    = NULL, size    = "n probes"
  ) +
  theme_hcc

bin_summary_pc1 <- probe_clean |>
  mutate(admr_grp = ifelse(is_admr == 1, "aDMR", "Background")) |>
  group_by(pc1_bin, admr_grp) |>
  summarise(
    mean_log_sd = mean(log10(beta_sd + 0.001), na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  filter(!is.na(pc1_bin))

int_pc1_coef <- coef(summary(lm_pc1))["pc1:is_admr", ]
p_interaction_pc1 <- ggplot(bin_summary_pc1, aes(x = pc1_bin, y = mean_log_sd,
                                                   color = admr_grp, group = admr_grp)) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(size = n), alpha = 0.8) +
  scale_color_manual(values = c("aDMR" = "#E24B4A", "Background" = "#888780")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
           label = sprintf("Interaction β=%.4f, p=%s",
                           int_pc1_coef["Estimate"], sig_label(int_pc1_coef["Pr(>|t|)"]))) +
  labs(
    title    = "B2: MicroC PC1 × aDMR Status Interaction",
    subtitle = "PC1 quintile (1=B-comp/late → 5=A-comp/early) vs β-SD",
    x        = "PC1 quintile (1=late/B → 5=early/A)",
    y        = "Mean log10(β-SD + 0.001)",
    color    = NULL, size    = "n probes"
  ) +
  theme_hcc

ggsave(file.path(FIG_DIR, "fig_b2d_interaction_pc1.png"),
       p_interaction_pc1, width = 7, height = 5, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_b2e_interaction_rt.png"),
       p_interaction_rt,  width = 7, height = 5, dpi = 150)
message("Saved: fig_b2d_interaction_pc1.png + fig_b2e_interaction_rt.png")

# 8. Final comparison summary ==================================================
message("\n=== Section 8: Final comparison summary ===")

comparison_df <- data.frame(
  proxy      = c("MicroC PC1", "ENCODE Repli-seq"),
  source     = c("HepG2 MicroC GSE278978", "ENCODE2 HepG2 S1-S4"),
  assembly   = c("hg38 (native)", "hg19 → liftOver hg38"),
  cor_rho    = c(1.0, cor_all$estimate),
  int_beta   = c(int_pc1_coef["Estimate"], int_rt["Estimate"]),
  int_p      = c(int_pc1_coef["Pr(>|t|)"], int_rt["Pr(>|t|)"]),
  gold_vs_rand_p  = c(
    wilcox_results$pc1_vs_rand_p[wilcox_results$tier == "Gold aDMR"],
    wilcox_results$rt_vs_rand_p[wilcox_results$tier == "Gold aDMR"]
  )
)
fwrite(comparison_df, file.path(OUT_DIR, "replication_timing_comparison.csv"))
message("Saved: replication_timing_comparison.csv")

cat("\n=== COMPARISON: PC1 vs Repli-seq ===\n")
print(comparison_df)

message("\nDone. Output files:")
message("  result/replication_timing_group.csv")
message("  result/replication_timing_loci.csv")
message("  result/replication_timing_probe.csv")
message("  result/replication_timing_comparison.csv")
message("  figures/fig_b2a_pc1_violin.png")
message("  figures/fig_b2b_rt_repliseq_violin.png")
message("  figures/fig_b2c_pc1_repliseq_concordance.png")
message("  figures/fig_b2d_interaction_pc1.png")
message("  figures/fig_b2e_interaction_rt.png")
