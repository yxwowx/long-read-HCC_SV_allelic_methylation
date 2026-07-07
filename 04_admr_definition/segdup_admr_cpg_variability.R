#!/usr/bin/env Rscript
# P1-B: SegDup-overlapping vs non-SegDup Gold aDMR
#   - CpG density (nCG per locus bp)
#   - Allelic methylation variability (hp_abs_diff = |HP Δβ|)
# Wilcoxon test + boxplot (per unique locus)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(gridExtra)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

# Paths ========================================================================
GOLD_CSV <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SEGDUP   <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
OUT_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result/figures")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Load data ====================================================================
cat("Loading Gold aDMR loci...\n")
gold <- fread(GOLD_CSV)

# Deduplicate to unique loci; for hp_abs_diff take median across patients
loci <- gold[, .(
  nCG         = nCG[1],
  hp_abs_diff = median(hp_abs_diff, na.rm = TRUE),
  n_patients  = .N
), by = .(admr_chr, admr_start, admr_end)]

loci[, locus_width := admr_end - admr_start]
loci[, cpg_density := nCG / locus_width * 1000]  # CpGs per kb

cat(sprintf("Unique Gold aDMR loci: %d\n", nrow(loci)))

# SegDup annotation ============================================================
cat("Annotating SegDup overlap...\n")
sd_gr <- import(SEGDUP, format = "BED")
seqlevelsStyle(sd_gr) <- "UCSC"

loci_gr <- GRanges(
  seqnames = loci$admr_chr,
  ranges   = IRanges(loci$admr_start + 1L, loci$admr_end)
)
loci[, segdup := factor(countOverlaps(loci_gr, sd_gr) > 0,
                        levels = c(FALSE, TRUE),
                        labels = c("non-SegDup", "SegDup"))]

n_sd  <- sum(loci$segdup == "SegDup")
n_nsd <- sum(loci$segdup == "non-SegDup")
cat(sprintf("SegDup-overlapping: %d | non-SegDup: %d\n", n_sd, n_nsd))

# Wilcoxon tests ===============================================================
run_wilcox <- function(var, label) {
  x <- loci[segdup == "SegDup",     get(var)]
  y <- loci[segdup == "non-SegDup", get(var)]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  w <- wilcox.test(x, y, exact = FALSE)
  cat(sprintf("\n%s:\n  SegDup median=%.4f (n=%d) | non-SegDup median=%.4f (n=%d)\n  Wilcoxon W=%.0f, p=%.4g\n",
              label, median(x), length(x), median(y), length(y), w$statistic, w$p.value))
  list(p = w$p.value, W = w$statistic,
       med_sd = median(x), med_nsd = median(y),
       n_sd = length(x), n_nsd = length(y))
}

cat("\n=== P1-B: SegDup vs non-SegDup Gold aDMR ===")
res_ncg  <- run_wilcox("cpg_density",  "CpG density (CpGs/kb)")
res_meth <- run_wilcox("hp_abs_diff",  "|HP Δβ| (allelic methylation variability)")

# Theme ========================================================================
COL_SD  <- "#C0392B"
COL_NSD <- "#7F8C8D"
theme_hcc <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        plot.title  = element_text(face = "bold", size = 11),
        axis.text   = element_text(size = 10))

fmt_p <- function(p) {
  if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
}

# Panel A: CpG density =========================================================
pA <- ggplot(loci, aes(x = segdup, y = cpg_density, fill = segdup)) +
  geom_boxplot(width = 0.5, outlier.size = 1.2, outlier.alpha = 0.5) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.4) +
  scale_x_discrete(labels = c("non-SegDup" = sprintf("non-SegDup\n(n=%d)", n_nsd),
                               "SegDup"     = sprintf("SegDup\n(n=%d)",    n_sd))) +
  scale_fill_manual(values = c("non-SegDup" = COL_NSD, "SegDup" = COL_SD)) +
  annotate("text", x = 1.5, y = max(loci$cpg_density, na.rm = TRUE) * 1.05,
           label = fmt_p(res_ncg$p), size = 3.5) +
  labs(title = "A", subtitle = "CpG density at Gold aDMR loci",
       x = NULL, y = "CpGs per kb") +
  theme_hcc + theme(legend.position = "none")

# Panel B: |HP delta-beta| =====================================================
pB <- ggplot(loci, aes(x = segdup, y = hp_abs_diff, fill = segdup)) +
  geom_boxplot(width = 0.5, outlier.size = 1.2, outlier.alpha = 0.5) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.4) +
  scale_x_discrete(labels = c("non-SegDup" = sprintf("non-SegDup\n(n=%d)", n_nsd),
                               "SegDup"     = sprintf("SegDup\n(n=%d)",    n_sd))) +
  scale_fill_manual(values = c("non-SegDup" = COL_NSD, "SegDup" = COL_SD)) +
  annotate("text", x = 1.5, y = max(loci$hp_abs_diff, na.rm = TRUE) * 1.05,
           label = fmt_p(res_meth$p), size = 3.5) +
  labs(title = "B", subtitle = "Allelic methylation variability at Gold aDMR loci",
       x = NULL, y = "|HP1 - HP2| (median across patients)") +
  theme_hcc + theme(legend.position = "none")

# Save =========================================================================
fig_path <- file.path(FIG_DIR, "fig_p1b_segdup_admr_cpg_variability.png")
fig <- arrangeGrob(pA, pB, ncol = 2)
ggsave(fig_path, fig, width = 9, height = 5, dpi = 300)
cat(sprintf("\nSaved: %s\n", fig_path))

# CSV summary ==================================================================
out <- data.table(
  metric    = c("CpG density (CpGs/kb)", "|HP abs_diff|"),
  n_segdup  = c(res_ncg$n_sd,   res_meth$n_sd),
  n_nsd     = c(res_ncg$n_nsd,  res_meth$n_nsd),
  med_sd    = c(res_ncg$med_sd, res_meth$med_sd),
  med_nsd   = c(res_ncg$med_nsd,res_meth$med_nsd),
  fold_sd_nsd = c(res_ncg$med_sd / res_ncg$med_nsd,
                  res_meth$med_sd / res_meth$med_nsd),
  wilcox_p  = c(res_ncg$p, res_meth$p)
)
fwrite(out, file.path(OUT_DIR, "p1b_segdup_admr_cpg_variability.csv"))
cat("Summary:\n"); print(out)
