#!/usr/bin/env -S conda run -n renv Rscript --vanilla
# Plot enrichment_ratio ~ window_kb for all PRIMARY_WINDOW runs
# Usage: Rscript plot_window_enrichment.R --outdir DIR [--nperm 100] [--output FILE]

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option("--outdir", type = "character", default = ".",
              help = "Directory containing *_window_enrich_full.csv files"),
  make_option("--nperm", type = "integer", default = 100L,
              help = "n_perm tag used in run_id (for title annotation) [default: %default]"),
  make_option("--output", type = "character", default = NULL,
              help = "Output PDF path [default: <outdir>/combined_enrichment_curves.pdf]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$output)) {
  opt$output <- file.path(opt$outdir, sprintf("combined_enrichment_curves_nperm%d.pdf", opt$nperm))
}

# Load all *_window_enrich_full.csv files ======================================
csv_files <- list.files(opt$outdir, pattern = "_window_enrich_full\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No *_window_enrich_full.csv files found in: ", opt$outdir)

message("Found ", length(csv_files), " CSV file(s):")
message(paste(" -", csv_files, collapse = "\n"))

df_list <- lapply(csv_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  # Extract primary_window label from filename (e.g. primary_25kb_window_enrich_full.csv)
  run_id <- sub("_window_enrich_full\\.csv$", "", basename(f))
  d$run_id <- run_id
  # Extract kb value for ordering
  kb_val <- as.integer(sub(".*_(\\d+)kb$", "\\1", run_id))
  d$primary_window_kb <- ifelse(is.na(kb_val), NA_integer_, kb_val)
  d
})

df <- bind_rows(df_list)

# Summarise: mean enrichment_ratio across patients per (run_id, cnv_class, window_kb)
df_sum <- df %>%
  group_by(run_id, primary_window_kb, cnv_class, window_kb) %>%
  summarise(
    mean_ratio = mean(enrichment_ratio, na.rm = TRUE),
    se_ratio   = sd(enrichment_ratio, na.rm = TRUE) / sqrt(sum(!is.na(enrichment_ratio))),
    .groups = "drop"
  ) %>%
  mutate(
    run_label = sprintf("Primary %dkb", primary_window_kb),
    run_label = factor(run_label, levels = sprintf("Primary %dkb", sort(unique(primary_window_kb))))
  )

CNV_COLORS <- c(
  "copy_neutral" = "#3B8BD4",
  "copy_gaining" = "#E24B4A",
  "copy_losing"  = "#BA7517",
  "insertion"    = "#7F77DD"
)

p <- ggplot(df_sum, aes(x = window_kb, y = mean_ratio,
                         color = run_label, group = run_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_ribbon(aes(ymin = mean_ratio - se_ratio, ymax = mean_ratio + se_ratio,
                  fill = run_label), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_log10(
    breaks = c(10, 50, 100, 500, 1000),
    labels = function(x) paste0(x, "kb")
  ) +
  scale_color_brewer(palette = "Set1", name = "Primary window") +
  scale_fill_brewer(palette = "Set1", name = "Primary window") +
  facet_wrap(~cnv_class, nrow = 2, scales = "free_y") +
  labs(
    title    = "DMR enrichment near SV breakpoints — window sensitivity",
    subtitle = sprintf("n_perm = %d | mean ± SE across patients", opt$nperm),
    x        = "Window size (log scale)",
    y        = "Enrichment ratio (obs / null)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey50", size = 10),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

ggsave(opt$output, p, width = 10, height = 7)
message("Plot saved: ", opt$output)

# Also save the summarised data
out_csv <- sub("\\.pdf$", ".csv", opt$output)
fwrite(df_sum, out_csv, row.names = FALSE, quote = FALSE)
message("Summary CSV saved: ", out_csv)
