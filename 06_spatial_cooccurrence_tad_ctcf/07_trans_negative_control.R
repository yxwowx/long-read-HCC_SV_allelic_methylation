#!/usr/bin/env Rscript
# P0-1: Trans negative control
# Tests whether HP|Δβ| decays with distance from SV (cis ≤50kb vs mid 50kb–1Mb vs trans >1Mb)
# proving that the distance-decay signal is spatial, not just window-based.
#
# Run: mamba run -n renv Rscript pipeline/07_trans_negative_control.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(clinfun)
  library(optparse)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--phased_ov", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/03.haplotype_sv_admr_analysis/all_hp_admr_tier.csv.gz",
    help = "Phased SV-aDMR overlap file (all_hp_admr_tier.csv.gz)"),
  make_option("--outdir", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/03.haplotype_sv_admr_analysis",
    help = "Output directory"),
  make_option("--run_id", type = "character", default = "tier_v2",
    help = "Run identifier prefix for output files"),
  make_option("--min_abs_db", type = "double", default = 0.1,
    help = "Minimum |Δβ| to include aDMR in analysis (default 0.1)")
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR   <- opt$outdir
RUN_ID   <- opt$run_id
LOG_FILE <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

DIST_BREAKS <- c(-Inf, 50e3, 1e6, Inf)
DIST_LABELS <- c("cis (≤50kb)", "mid (50kb–1Mb)", "trans (>1Mb)")
DIST_COLORS <- c("cis (≤50kb)" = "#E24B4A", "mid (50kb–1Mb)" = "#BA7517", "trans (>1Mb)" = "#888780")

# ── 1. Load ──────────────────────────────────────────────────────────────────
message("Reading phased overlap: ", opt$phased_ov)
dat <- fread(opt$phased_ov)

# Flexible column resolution
find_col <- function(df, candidates) {
  m <- intersect(candidates, names(df))
  if (length(m) == 0) stop("None of [", paste(candidates, collapse=", "), "] found in file")
  m[[1]]
}

col_dist   <- find_col(dat, c("dist", "distance", "abs_dist", "dist_to_sv"))
col_abs_db <- find_col(dat, c("abs_hp_delta", "abs_delta_beta", "abs_db", "hp_abs_diff"))
col_tier   <- find_col(dat, c("sv_tier", "stratification", "tier"))
col_pt     <- find_col(dat, c("sample", "patient_code", "patient_id"))
col_type   <- find_col(dat, c("sv_type", "geom_type", "cnv_class"))

message(sprintf("Using cols: dist=%s  abs_db=%s  tier=%s  patient=%s",
                col_dist, col_abs_db, col_tier, col_pt))

dat <- dat %>%
  dplyr::rename(
    raw_dist = !!col_dist,
    abs_db   = !!col_abs_db,
    sv_tier  = !!col_tier,
    sample   = !!col_pt,
    sv_type  = !!col_type
  ) %>%
  dplyr::mutate(
    abs_dist = abs(raw_dist),
    dist_class = cut(abs_dist,
      breaks = DIST_BREAKS,
      labels = DIST_LABELS,
      include.lowest = TRUE, right = TRUE)
  ) %>%
  dplyr::filter(!is.na(abs_db), abs_db >= opt$min_abs_db, !is.na(dist_class))

message(sprintf("Rows after filter (abs_db >= %.2f): %d", opt$min_abs_db, nrow(dat)))
print(table(dat$dist_class))

# ── 2. KS tests: cis vs mid, cis vs trans, mid vs trans ──────────────────────
ks_pairs <- list(
  c("cis (≤50kb)", "mid (50kb–1Mb)"),
  c("cis (≤50kb)", "trans (>1Mb)"),
  c("mid (50kb–1Mb)", "trans (>1Mb)")
)

ks_res <- lapply(ks_pairs, function(pair) {
  x <- dat$abs_db[dat$dist_class == pair[1]]
  y <- dat$abs_db[dat$dist_class == pair[2]]
  if (length(x) < 5 || length(y) < 5) return(NULL)
  t <- ks.test(x, y, alternative = "greater")
  data.frame(group1 = pair[1], group2 = pair[2],
             n1 = length(x), n2 = length(y),
             KS_D = round(t$statistic, 4),
             p_value = round(t$p.value, 6),
             stringsAsFactors = FALSE)
}) %>% bind_rows()
cat("\n=== KS test (group1 > group2 in |Δβ|) ===\n")
print(ks_res)

# ── 3. JT trend test (cis > mid > trans) ─────────────────────────────────────
jt_df <- dat %>% filter(!is.na(dist_class))
jt_order <- as.integer(factor(jt_df$dist_class,
                                levels = rev(DIST_LABELS)))  # cis=3, mid=2, trans=1 → decreasing
jt_res <- tryCatch(
  clinfun::jonckheere.test(jt_df$abs_db, jt_order,
                           alternative = "decreasing", nperm = 2000L),
  error = function(e) { message("JT failed: ", e$message); list(p.value = NA) }
)
cat(sprintf("\nJT trend (cis>mid>trans in |Δβ|): p = %.4f\n", jt_res$p.value))

# ── 4. Per-patient median summary ────────────────────────────────────────────
pt_summary <- dat %>%
  group_by(sample, dist_class) %>%
  summarise(
    n_dmr       = n(),
    median_abs_db = median(abs_db, na.rm = TRUE),
    mean_abs_db   = mean(abs_db, na.rm = TRUE),
    .groups = "drop"
  )

# ── 5. Plots ──────────────────────────────────────────────────────────────────
p_violin <- ggplot(dat, aes(x = dist_class, y = abs_db, fill = dist_class)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = DIST_COLORS) +
  scale_y_continuous(limits = c(0, quantile(dat$abs_db, 0.99))) +
  labs(title = "P0-1: Trans negative control",
       subtitle = sprintf("JT p = %.4f  |  n = %d DMR pairs", jt_res$p.value, nrow(dat)),
       x = "Distance class", y = "HP |Δβ|") +
  theme_hcc

p_pt <- ggplot(pt_summary, aes(x = dist_class, y = median_abs_db,
                                group = sample, color = dist_class)) +
  geom_line(color = "grey70", linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = DIST_COLORS) +
  labs(title = "Per-patient median |Δβ|",
       x = "Distance class", y = "Median HP |Δβ|") +
  theme_hcc

p_tier <- ggplot(dat, aes(x = dist_class, y = abs_db, fill = dist_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~sv_tier) +
  scale_fill_manual(values = DIST_COLORS) +
  scale_y_continuous(limits = c(0, quantile(dat$abs_db, 0.99))) +
  labs(title = "By SV tier", x = NULL, y = "HP |Δβ|") +
  theme_hcc

p_combined <- p_violin | p_pt | p_tier
ggsave(file.path(OUTDIR, paste0(RUN_ID, "_P01_trans_neg_control.png")),
       p_combined, width = 16, height = 5, dpi = 150)

# ── 6. Save CSVs ──────────────────────────────────────────────────────────────
fwrite(ks_res,     file.path(OUTDIR, paste0(RUN_ID, "_P01_ks_test.csv")))
fwrite(pt_summary, file.path(OUTDIR, paste0(RUN_ID, "_P01_pt_summary.csv")))
fwrite(data.frame(JT_p = jt_res$p.value, n_obs = nrow(jt_df)),
       file.path(OUTDIR, paste0(RUN_ID, "_P01_jt_result.csv")))

cat(append = TRUE,
    text   = sprintf("[%s] P0-1 trans_negative_control: JT p=%.4f; KS cis>trans D=%.3f p=%.4f; n=%d\n",
                     Sys.Date(), jt_res$p.value,
                     ks_res$KS_D[ks_res$group2 == "trans (>1Mb)"],
                     ks_res$p_value[ks_res$group2 == "trans (>1Mb)"],
                     nrow(dat)),
    file   = LOG_FILE)

message("Done: P0-1 outputs in ", OUTDIR)
