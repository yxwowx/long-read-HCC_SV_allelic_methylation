#!/usr/bin/env Rscript
# C17 probe-density audit: HM450 probe depletion in SegDup regions
#
# Quantifies how under-represented HM450 probes are in SegDup regions,
# explaining the attenuated C17 OR=1.05 vs in-house aDMR OR=2.22.
#
# Strategy:
#   1. Extract all probe coordinates from cached SE (all probes, not just DMRs)
#   2. Calculate probes/kb in SegDup vs non-SegDup genome
#   3. Compute depletion ratio and expected OR attenuation
#   4. Secondary: among SegDup-overlapping probes, what % come from DMR set?
#
# Output:
#   result/c17_hm450_segdup_probe_density.csv   — summary table
#   result/figures/fig_c17_probe_density.png    — depletion bar + scatter

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
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

METH_RDS <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/tcga_cache/tcga_lihc_meth450k_se.rds")
SEGDUP   <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
FAI      <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

# 1. All HM450 probe coordinates ===============================================
message("Loading HM450 probes from cached SE...")
meth_se   <- readRDS(METH_RDS)
probe_gr  <- rowRanges(meth_se)
seqlevelsStyle(probe_gr) <- "UCSC"
probe_gr  <- keepStandardChromosomes(probe_gr, pruning.mode = "coarse")
probe_gr  <- probe_gr[grepl("^chr[0-9XY]+$", as.character(seqnames(probe_gr)))]
cat(sprintf("Total HM450 probes (autosomes + chrX/Y): %d\n", length(probe_gr)))

# DMR probe flag (probes used in C17 analysis: FDR<0.05, |Δβ|≥0.1)
# Re-derive from the proximity result: n_bg_probes = 415372, n_dmr_probes = 124153
# We tag by overlap with DMR β result if available; otherwise use all probes as background
dmr_csv <- file.path(OUT_DIR, "tcga_lihc_dmr_segdup_proximity.csv")
has_dmr_list <- file.exists(dmr_csv)

# 2. SegDup regions ============================================================
message("Loading SegDup...")
segdup_raw <- import(SEGDUP, format = "BED")
seqlevelsStyle(segdup_raw) <- "UCSC"
segdup_gr  <- keepStandardChromosomes(segdup_raw, pruning.mode = "coarse")
segdup_gr  <- segdup_gr[grepl("^chr[0-9XY]+$", as.character(seqnames(segdup_gr)))]
segdup_gr  <- reduce(segdup_gr)
cat(sprintf("SegDup regions (merged): %d, total bp: %s\n",
            length(segdup_gr),
            format(sum(as.numeric(width(segdup_gr))), big.mark = ",")))

# 3. Chrom sizes ===============================================================
chrom_df <- fread(FAI,
                  col.names = c("chr", "len", "x", "y", "z"),
                  data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr))
genome_bp    <- sum(as.numeric(chrom_df$len))
segdup_bp    <- sum(as.numeric(width(segdup_gr)))
nonsegdup_bp <- genome_bp - segdup_bp
cat(sprintf("Genome: %.2f Gb | SegDup: %.0f Mb (%.1f%%) | non-SegDup: %.2f Gb\n",
            genome_bp / 1e9, segdup_bp / 1e6,
            100 * segdup_bp / genome_bp, nonsegdup_bp / 1e9))

# 4. Probe overlap with SegDup =================================================
message("Overlapping probes with SegDup...")
probe_in_segdup  <- overlapsAny(probe_gr, segdup_gr)
n_probe_segdup   <- sum(probe_in_segdup)
n_probe_nonseg   <- sum(!probe_in_segdup)
n_probe_total    <- length(probe_gr)

density_segdup   <- n_probe_segdup   / (segdup_bp    / 1e3)   # probes per kb
density_nonseg   <- n_probe_nonseg   / (nonsegdup_bp / 1e3)
density_genome   <- n_probe_total    / (genome_bp    / 1e3)

depletion_ratio  <- density_segdup / density_nonseg  # < 1 = depleted

cat(sprintf("\nProbe counts:\n"))
cat(sprintf("  SegDup    : %d probes / %.0f Mb = %.4f probes/kb\n",
            n_probe_segdup, segdup_bp / 1e6, density_segdup))
cat(sprintf("  Non-SegDup: %d probes / %.2f Gb = %.4f probes/kb\n",
            n_probe_nonseg, nonsegdup_bp / 1e9, density_nonseg))
cat(sprintf("  Genome    : %d probes / %.2f Gb = %.4f probes/kb\n",
            n_probe_total, genome_bp / 1e9, density_genome))
cat(sprintf("  Depletion ratio (SegDup/non-SegDup): %.3f (%.1f%% of expected)\n",
            depletion_ratio, 100 * depletion_ratio))

# 5. Expected OR attenuation ===================================================
# If true SegDup enrichment OR = 2.22 (aDMR), and HM450 probes are depleted
# in SegDup by factor D, the observable OR is attenuated to:
#   OR_obs ≈ OR_true * D   (first-order approximation for small enrichment)
# More precisely: use the probe-density-corrected expected OR calculation.
true_or_admr <- 2.22
true_or_sv   <- 2.11
obs_or_c17   <- 1.05

# Probe-count based correction: fraction of SegDup "accessible" by HM450
pct_segdup_accessible <- n_probe_segdup / (segdup_bp / 1e3) /
                         (n_probe_total / (genome_bp / 1e3))

expected_or_under_depletion_admr <- 1 + (true_or_admr - 1) * depletion_ratio
expected_or_under_depletion_sv   <- 1 + (true_or_sv   - 1) * depletion_ratio

cat(sprintf("\n--- Expected OR attenuation under probe depletion ---\n"))
cat(sprintf("True aDMR OR (in-house)      : %.2f\n", true_or_admr))
cat(sprintf("Probe depletion factor       : %.3f\n", depletion_ratio))
cat(sprintf("Expected attenuated OR (aDMR): %.3f\n", expected_or_under_depletion_admr))
cat(sprintf("Observed C17 OR              : %.3f\n", obs_or_c17))
cat(sprintf("Residual unexplained attenuation: %.3f\n",
            obs_or_c17 - expected_or_under_depletion_admr))

# 6. Per-chromosome breakdown ==================================================
chr_stats <- lapply(chrom_df$chr, function(ch) {
  pr_chr <- probe_gr[seqnames(probe_gr) == ch]
  sd_chr <- segdup_gr[seqnames(segdup_gr) == ch]
  clen   <- chrom_df$len[chrom_df$chr == ch]
  sd_bp  <- if (length(sd_chr) > 0) sum(as.numeric(width(sd_chr))) else 0
  n_sd   <- if (length(sd_chr) > 0) sum(overlapsAny(pr_chr, sd_chr)) else 0
  n_ns   <- length(pr_chr) - n_sd
  data.frame(
    chr           = ch,
    chrom_len     = clen,
    segdup_bp     = sd_bp,
    n_probes      = length(pr_chr),
    n_probe_segdup = n_sd,
    n_probe_nonseg = n_ns,
    pct_probe_in_segdup = ifelse(length(pr_chr) > 0, 100 * n_sd / length(pr_chr), NA),
    pct_bp_segdup = 100 * sd_bp / clen,
    probe_density_segdup = ifelse(sd_bp > 0, n_sd / (sd_bp / 1e3), NA),
    probe_density_nonseg = ifelse(clen - sd_bp > 0, n_ns / ((clen - sd_bp) / 1e3), NA)
  )
}) |> bind_rows()

chr_stats <- chr_stats |>
  mutate(depletion = probe_density_segdup / probe_density_nonseg)

# 7. Save summary ==============================================================
summary_df <- data.frame(
  metric = c(
    "n_probes_total", "n_probes_in_segdup", "n_probes_not_segdup",
    "segdup_bp", "nonsegdup_bp", "genome_bp",
    "probe_density_segdup_per_kb", "probe_density_nonsegdup_per_kb",
    "probe_density_genome_per_kb", "depletion_ratio",
    "pct_probes_in_segdup", "pct_genome_segdup",
    "true_OR_admr_inhouse", "obs_OR_c17",
    "expected_OR_under_depletion", "residual_attenuation"
  ),
  value = c(
    n_probe_total, n_probe_segdup, n_probe_nonseg,
    segdup_bp, nonsegdup_bp, genome_bp,
    density_segdup, density_nonseg,
    density_genome, depletion_ratio,
    100 * n_probe_segdup / n_probe_total,
    100 * segdup_bp / genome_bp,
    true_or_admr, obs_or_c17,
    expected_or_under_depletion_admr,
    obs_or_c17 - expected_or_under_depletion_admr
  )
)

fwrite(summary_df, file.path(OUT_DIR, "c17_hm450_segdup_probe_density.csv"))
fwrite(chr_stats,  file.path(OUT_DIR, "c17_hm450_segdup_probe_density_per_chr.csv"))
message("Wrote: c17_hm450_segdup_probe_density.csv")

# 8. Figure ====================================================================
# Panel A: probe density comparison bar
bar_df <- data.frame(
  region  = c("SegDup", "Non-SegDup", "Genome"),
  density = c(density_segdup, density_nonseg, density_genome)
)
bar_df$region <- factor(bar_df$region, levels = c("Non-SegDup", "Genome", "SegDup"))

pA <- ggplot(bar_df, aes(x = region, y = density, fill = region)) +
  geom_col(width = 0.6, color = "white") +
  scale_fill_manual(values = c("SegDup" = "#d73027",
                               "Non-SegDup" = "#4575b4",
                               "Genome" = "#74add1")) +
  labs(title = "HM450 probe density by genomic region",
       x = NULL, y = "Probes per kb",
       caption = sprintf("Depletion ratio (SegDup/Non-SegDup) = %.3f", depletion_ratio)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# Panel B: per-chromosome depletion scatter
chr_plot <- chr_stats |>
  filter(!is.na(depletion), pct_bp_segdup > 0.1) |>
  mutate(label = ifelse(depletion < 0.3 | depletion > 1.2, chr, ""))

pB <- ggplot(chr_plot, aes(x = pct_bp_segdup, y = depletion, label = label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = depletion_ratio, linetype = "dotted",
             color = "#d73027", linewidth = 0.8) +
  geom_point(aes(size = n_probes), color = "#4575b4", alpha = 0.7) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 10) +
  scale_size_continuous(range = c(2, 8), name = "N probes") +
  annotate("text", x = max(chr_plot$pct_bp_segdup) * 0.7,
           y = depletion_ratio + 0.04,
           label = sprintf("Overall = %.3f", depletion_ratio),
           color = "#d73027", size = 3.5) +
  labs(title = "Per-chromosome probe depletion in SegDup",
       x = "% chromosome bp in SegDup",
       y = "Probe density ratio (SegDup / Non-SegDup)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Panel C: expected vs observed OR
or_df <- data.frame(
  label = c("True OR\n(in-house aDMR)", "Expected OR\n(under depletion)", "Observed OR\n(TCGA-LIHC C17)"),
  OR    = c(true_or_admr, expected_or_under_depletion_admr, obs_or_c17),
  type  = c("true", "expected", "observed")
)
or_df$label <- factor(or_df$label, levels = or_df$label)

pC <- ggplot(or_df, aes(x = label, y = OR, fill = type)) +
  geom_col(width = 0.5, color = "white") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("true" = "#1a9641",
                               "expected" = "#fdae61",
                               "observed" = "#d73027")) +
  labs(title = "OR attenuation due to probe depletion",
       x = NULL, y = "Odds Ratio (SegDup enrichment)",
       caption = sprintf("Depletion factor = %.3f; expected OR = 1 + (OR_true-1)*D = %.3f",
                         depletion_ratio, expected_or_under_depletion_admr)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

combined <- (pA | pC) / pB + plot_layout(heights = c(1, 1.5))
ggsave(file.path(FIG_DIR, "fig_c17_probe_density.png"), combined,
       width = 11, height = 9, dpi = 150)
message("Saved: fig_c17_probe_density.png")

cat("\n=== C17 PROBE DENSITY SUMMARY ===\n")
print(summary_df[summary_df$metric %in% c(
  "n_probes_total", "n_probes_in_segdup",
  "probe_density_segdup_per_kb", "probe_density_nonsegdup_per_kb",
  "depletion_ratio", "pct_probes_in_segdup", "pct_genome_segdup",
  "true_OR_admr_inhouse", "obs_OR_c17",
  "expected_OR_under_depletion", "residual_attenuation"
), ], row.names = FALSE)
cat("Done.\n")
