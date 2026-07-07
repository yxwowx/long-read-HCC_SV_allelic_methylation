#!/usr/bin/env Rscript
# Compare primary_50kb vs primary_100kb window enrichment results
# Usage: Rscript compare_window_runs.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

OUTDIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/window_enrichment/nperm_1000")

theme_hcc <- theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey50", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom"
  )

RUN_COLORS <- c("primary_50kb" = "#E24B4A", "primary_100kb" = "#3B8BD4")

# 1. Load and combine files ====================================================
f50  <- file.path(OUTDIR, "primary_50kb_window_enrich_full.csv")
f100 <- file.path(OUTDIR, "primary_100kb_window_enrich_full.csv")

if (!file.exists(f50))  stop("File not found: ", f50)
if (!file.exists(f100)) stop("File not found: ", f100)

d50  <- fread(f50)  |> mutate(run_id = "primary_50kb",  primary_window_kb = 50)
d100 <- fread(f100) |> mutate(run_id = "primary_100kb", primary_window_kb = 100)

combined <- bind_rows(d50, d100) |>
  mutate(run_id = factor(run_id, levels = c("primary_50kb", "primary_100kb")))

out_csv <- file.path(OUTDIR, "combined_window_enrich.csv")
fwrite(combined, out_csv, row.names = FALSE, quote = FALSE)
cat("Combined CSV saved:", out_csv, "\n")

# 2. Graph A: enrichment_ratio curve ===========================================

CNV_ORDER <- c("copy_neutral", "copy_gaining", "copy_losing", "insertion", "COM")
CNV_LABEL <- c(
  copy_neutral = "Copy-neutral\n(INV/TRA)",
  copy_gaining = "Copy-gaining\n(DUP)",
  copy_losing  = "Copy-losing\n(DEL)",
  insertion    = "Insertion\n(INS)",
  COM          = "Complex\n(COM)"
)

plot_A_df <- combined |>
  filter(cnv_class %in% CNV_ORDER) |>
  mutate(cnv_class = factor(cnv_class, levels = CNV_ORDER,
                             labels = CNV_LABEL[CNV_ORDER])) |>
  group_by(run_id, cnv_class, window_kb) |>
  summarise(
    mean_ratio = mean(enrichment_ratio, na.rm = TRUE),
    se_ratio   = sd(enrichment_ratio, na.rm = TRUE) /
                   sqrt(sum(!is.na(enrichment_ratio))),
    .groups    = "drop"
  )

p_A <- ggplot(plot_A_df,
              aes(x = window_kb, y = mean_ratio,
                  color = run_id, fill = run_id, group = run_id)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_ribbon(aes(ymin = mean_ratio - se_ratio,
                  ymax = mean_ratio + se_ratio),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_x_log10(
    limits = c(9, 1100),
    breaks = c(10, 50, 100, 500, 1000),
    labels = paste0(c(10, 50, 100, 500, 1000), "kb")
  ) +
  scale_color_manual(values = RUN_COLORS, name = "Run") +
  scale_fill_manual(values  = RUN_COLORS, name = "Run") +
  facet_wrap(~cnv_class, nrow = 1) +
  labs(
    title    = "Enrichment ratio curve: primary_50kb vs primary_100kb",
    subtitle = "Mean ± SE across patients | dashed line: ratio = 1",
    x        = "Window size (log scale)", y = "Enrichment ratio (obs / null)",
    caption  = "FDR: BH | n_perm = 1000"
  ) +
  theme_hcc +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

out_A <- file.path(OUTDIR, "compare_50kb_vs_100kb_ratio.pdf")
ggsave(out_A, p_A, width = 16, height = 5, device = cairo_pdf)
cat("Graph A saved:", out_A, "\n")

# 3. Graph B: p_fdr distribution (primary_window rows only) ====================

plot_B_df <- combined |>
  filter(
    (run_id == "primary_50kb"  & window_kb == 50) |
    (run_id == "primary_100kb" & window_kb == 100)
  ) |>
  filter(cnv_class %in% CNV_ORDER) |>
  mutate(
    cnv_class   = factor(cnv_class, levels = CNV_ORDER,
                          labels = CNV_LABEL[CNV_ORDER]),
    neg_log_fdr = -log10(pmax(p_fdr, 1e-10))
  )

sig_line <- -log10(0.05)

p_B <- ggplot(plot_B_df,
              aes(x = cnv_class, y = neg_log_fdr, color = run_id)) +
  geom_hline(yintercept = sig_line, linetype = "dashed",
             color = "grey40", linewidth = 0.7) +
  annotate("text", x = "Copy-neutral\n(INV/TRA)", y = sig_line + 0.08,
           label = "-log10(0.05)", size = 3, color = "grey40", hjust = 0.5) +
  geom_jitter(size = 2.2, alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width  = 0.5)) +
  stat_summary(fun = median, geom = "crossbar",
               aes(group = run_id),
               width = 0.35, linewidth = 0.5,
               position = position_dodge(width = 0.5)) +
  scale_color_manual(values = RUN_COLORS, name = "Run") +
  labs(
    title    = "p_fdr distribution, at the primary window",
    subtitle = "50kb run uses window_kb=50 rows, 100kb run uses window_kb=100 rows | horizontal line: median",
    x        = NULL, y = "-log10(FDR-adjusted p)",
    caption  = "Dashed line: -log10(0.05) | each point = one patient"
  ) +
  theme_hcc +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

out_B <- file.path(OUTDIR, "compare_50kb_vs_100kb_pvalue.pdf")
ggsave(out_B, p_B, width = 10, height = 5, device = cairo_pdf)
cat("Graph B saved:", out_B, "\n")

# 4. Compare significant cases (p_fdr < 0.05) ==================================
cat("\n=== Comparing significant cases (p_fdr < 0.05), at the primary window ===\n\n")

sig_cases <- combined |>
  filter(
    (run_id == "primary_50kb"  & window_kb == 50) |
    (run_id == "primary_100kb" & window_kb == 100)
  ) |>
  filter(p_fdr < 0.05)

sig_summary <- sig_cases |>
  group_by(run_id) |>
  summarise(
    n_sig_total = n(),
    cnv_classes = paste(sort(unique(cnv_class)), collapse = ", "),
    .groups = "drop"
  )

cat("Significant case count by run:\n")
print(as.data.frame(sig_summary))

cat("\nDetail by cnv_class x patient:\n")
sig_detail <- sig_cases |>
  select(run_id, cnv_class, patient_id, window_kb, enrichment_ratio, p_fdr) |>
  arrange(run_id, cnv_class, patient_id)
print(as.data.frame(sig_detail))

# Cases shared between both runs vs unique to one run
key50  <- sig_cases |> filter(run_id == "primary_50kb") |>
            mutate(key = paste(cnv_class, patient_id, sep = ":")) |>
            pull(key)
key100 <- sig_cases |> filter(run_id == "primary_100kb") |>
            mutate(key = paste(cnv_class, patient_id, sep = ":")) |>
            pull(key)

cat("\nComparison summary:\n")
cat(sprintf("primary_50kb  significant: %d\n", length(key50)))
cat(sprintf("primary_100kb significant: %d\n", length(key100)))
cat(sprintf("Shared (both runs): %d  — %s\n",
            length(intersect(key50, key100)),
            if (length(intersect(key50, key100)) > 0)
              paste(intersect(key50, key100), collapse = ", ")
            else "none"))
cat(sprintf("primary_50kb only:  %d  — %s\n",
            length(setdiff(key50, key100)),
            if (length(setdiff(key50, key100)) > 0)
              paste(setdiff(key50, key100), collapse = ", ")
            else "none"))
cat(sprintf("primary_100kb only: %d  — %s\n",
            length(setdiff(key100, key50)),
            if (length(setdiff(key100, key50)) > 0)
              paste(setdiff(key100, key50), collapse = ", ")
            else "none"))

winner <- if (length(key50) > length(key100)) "primary_50kb" else
          if (length(key100) > length(key50)) "primary_100kb" else "tied"
cat(sprintf("\nRun with more significant cases: %s\n", winner))

cat("\n=== Done ===\n")
