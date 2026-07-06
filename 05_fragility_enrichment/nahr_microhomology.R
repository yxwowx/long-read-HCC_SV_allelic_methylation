#!/usr/bin/env Rscript
# P0-B: NAHR breakpoint signature at SegDup-overlapping non-boundary SVs
#
# NOTE: Severus (long-read) reports HOMLEN=0 for all SVs in this cohort.
# Alternative NAHR proxies used:
#   (1) SV size distribution for DEL/DUP (NAHR events cluster at inter-repeat distances)
#   (2) INSIDE_VNTR flag (tandem repeat context)
#   (3) Breakpoint PC1 / LAD enrichment across SegDup vs non-SegDup
#
# Output: result/nahr_microhomology.csv + result/figures/fig_nahr_svsize.png

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

SV_ANN   <- "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz"
FRAG_ANN <- "/node200data/kachungk/hcc_data/DMR_SVs/result/sv_fragility_annotation.csv"
OUT_DIR  <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load SV annotation ---
message("Reading SV annotation...")
sv_ann <- fread(cmd = paste("zcat", SV_ANN), data.table = FALSE) |>
  select(bp_id, HOMLEN, svtype, svLen, INSIDE_VNTR, sample, sv_tier, is_hbv)

cat(sprintf("HOMLEN range: %d – %d (all-zero: %s)\n",
            min(sv_ann$HOMLEN), max(sv_ann$HOMLEN),
            all(sv_ann$HOMLEN == 0)))

# --- 2. Load fragility annotation; deduplicate bp_id ---
frag <- fread(FRAG_ANN, data.table = FALSE) |>
  select(bp_id, segdup_overlap, tier_group, lad_overlap, b_compartment, pc1_score) |>
  distinct(bp_id, .keep_all = TRUE)

# --- 3. Merge; restrict to non-boundary SVs ---
df <- left_join(sv_ann, frag, by = "bp_id") |>
  filter(tier_group == "Non-boundary")

cat(sprintf("\nNon-boundary SVs: %d total\n", nrow(df)))
cat(sprintf("SegDup-overlapping: %d | Non-SegDup: %d\n",
            sum(df$segdup_overlap, na.rm = TRUE),
            sum(!df$segdup_overlap, na.rm = TRUE)))
cat(sprintf("SV type counts:\n"))
print(table(df$svtype))

# --- 4. HOMLEN check (document limitation) ---
cat(sprintf("\nHOMLEN uniformly 0: %s — Severus long-read limitation\n",
            all(df$HOMLEN == 0)))
cat("Switching to SV size (svLen) and VNTR context as NAHR proxies.\n\n")

# --- 5a. DEL/DUP only: SV size by SegDup overlap ---
del_dup <- df |>
  filter(svtype %in% c("DEL", "DUP"), !is.na(svLen), svLen > 0) |>
  mutate(
    group     = ifelse(segdup_overlap, "SegDup-overlapping", "Non-SegDup"),
    group     = factor(group, levels = c("Non-SegDup", "SegDup-overlapping")),
    log10_len = log10(abs(svLen))
  )

cat(sprintf("DEL/DUP with svLen > 0: %d (SegDup=%d, non-SegDup=%d)\n",
            nrow(del_dup),
            sum(del_dup$segdup_overlap),
            sum(!del_dup$segdup_overlap)))

wt_size <- wilcox.test(log10_len ~ group, data = del_dup)
size_stats <- del_dup |>
  group_by(group) |>
  summarise(n = n(), median_bp = median(abs(svLen)),
            mean_bp = mean(abs(svLen)), .groups = "drop")
cat("\nSV size stats (DEL+DUP):\n")
print(size_stats)
cat(sprintf("Wilcoxon log10(svLen) SegDup vs non-SegDup: W=%.0f, p=%.4f\n",
            wt_size$statistic, wt_size$p.value))

# --- 5b. VNTR context: INSIDE_VNTR by SegDup overlap ---
vntr_tab <- table(segdup = df$segdup_overlap, inside_vntr = df$INSIDE_VNTR)
cat("\nContingency table (SegDup × INSIDE_VNTR):\n")
print(vntr_tab)
# Only run Fisher test if the table has variation
if (nrow(vntr_tab) >= 2 && ncol(vntr_tab) >= 2 && all(vntr_tab > 0)) {
  ft_vntr <- fisher.test(vntr_tab)
  cat(sprintf("Fisher VNTR: OR=%.3f (95%% CI %.3f–%.3f), p=%.4f\n",
              ft_vntr$estimate, ft_vntr$conf.int[1], ft_vntr$conf.int[2], ft_vntr$p.value))
} else {
  ft_vntr <- list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  cat("Fisher test skipped (insufficient VNTR variation).\n")
}

# --- 5c. PC1 (compartment score) by SegDup overlap ---
wt_pc1 <- wilcox.test(pc1_score ~ segdup_overlap, data = df[!is.na(df$pc1_score), ])
cat(sprintf("Wilcoxon PC1 (SegDup vs non-SegDup): W=%.0f, p=%.4f\n",
            wt_pc1$statistic, wt_pc1$p.value))

pc1_stats <- df |>
  filter(!is.na(pc1_score)) |>
  group_by(segdup_overlap) |>
  summarise(n = n(), median_pc1 = median(pc1_score), mean_pc1 = mean(pc1_score),
            .groups = "drop")
cat("\nPC1 stats by SegDup:\n")
print(pc1_stats)

# --- 6. Save result CSV ---
result_df <- data.frame(
  metric              = c("svLen_DEL_DUP", "INSIDE_VNTR_Fisher", "PC1_Wilcoxon"),
  n_segdup            = c(sum(del_dup$segdup_overlap),
                          sum(df$segdup_overlap, na.rm = TRUE),
                          sum(df$segdup_overlap & !is.na(df$pc1_score), na.rm = TRUE)),
  n_non_segdup        = c(sum(!del_dup$segdup_overlap),
                          sum(!df$segdup_overlap, na.rm = TRUE),
                          sum(!df$segdup_overlap & !is.na(df$pc1_score), na.rm = TRUE)),
  stat                = c(wt_size$statistic, ft_vntr$estimate, wt_pc1$statistic),
  p_value             = c(wt_size$p.value, ft_vntr$p.value, wt_pc1$p.value),
  interpretation      = c(
    ifelse(wt_size$p.value < 0.05, "SegDup SVs differ in size (NAHR size signature)", "No size difference"),
    ifelse(!is.na(ft_vntr$p.value) && ft_vntr$p.value < 0.05, "SegDup SVs enriched at VNTR loci", "No VNTR enrichment"),
    ifelse(wt_pc1$p.value < 0.05, "SegDup SVs in distinct chromatin compartment", "No compartment difference")
  ),
  homlen_note         = "HOMLEN=0 for all SVs; Severus long-read limitation — NAHR direct signature not testable"
)
fwrite(result_df, file.path(OUT_DIR, "nahr_microhomology.csv"))
message("Wrote: ", file.path(OUT_DIR, "nahr_microhomology.csv"))

# --- 7. Figures ---
p_label_size <- ifelse(wt_size$p.value < 0.001, "p<0.001",
                sprintf("p=%.3f", wt_size$p.value))

# Panel A: SV size violin
pA <- ggplot(del_dup, aes(x = group, y = log10_len, fill = group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Non-SegDup" = "#6baed6", "SegDup-overlapping" = "#d73027")) +
  annotate("text", x = 1.5, y = max(del_dup$log10_len, na.rm = TRUE),
           label = p_label_size, size = 3.5, hjust = 0.5) +
  labs(
    title    = "SV Size at Non-Boundary Breakpoints",
    subtitle = "DEL + DUP only",
    x        = NULL,
    y        = "log10(SV size, bp)",
    caption  = "Wilcoxon rank-sum; NAHR proxy (HOMLEN=0 in long-read data)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# Panel B: PC1 violin (compartment context)
p_label_pc1 <- ifelse(wt_pc1$p.value < 0.001, "p<0.001",
               sprintf("p=%.3f", wt_pc1$p.value))

pB <- df |>
  filter(!is.na(pc1_score)) |>
  mutate(group = ifelse(segdup_overlap, "SegDup-overlapping", "Non-SegDup"),
         group = factor(group, levels = c("Non-SegDup", "SegDup-overlapping"))) |>
  ggplot(aes(x = group, y = pc1_score, fill = group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Non-SegDup" = "#6baed6", "SegDup-overlapping" = "#d73027")) +
  annotate("text", x = 1.5, y = max(df$pc1_score, na.rm = TRUE),
           label = p_label_pc1, size = 3.5, hjust = 0.5) +
  labs(
    title    = "HepG2 Compartment Score at Breakpoints",
    subtitle = "Non-boundary SVs",
    x        = NULL,
    y        = "PC1 score (negative = B-compartment)",
    caption  = "Wilcoxon rank-sum; HepG2 Micro-C PC1"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

combined <- pA | pB
ggsave(file.path(FIG_DIR, "fig_nahr_svsize.png"), combined,
       width = 10, height = 5, dpi = 150)
message("Saved: ", file.path(FIG_DIR, "fig_nahr_svsize.png"))

cat("\nDone.\n")
