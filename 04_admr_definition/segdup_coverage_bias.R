#!/usr/bin/env Rscript
# segdup_coverage_bias.R
# Quantify mapping bias at SegDup regions using mosdepth 500bp-bin coverage.
# For each tumor (and matched normal) sample: compare median coverage at
# SegDup bins vs. autosomal non-SegDup bins, and report fold-ratio.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

# Paths ========================================================================
COV_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "minimap2.out_hg38/coverage")
SEGDUP   <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
OUT_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
OUT_CSV  <- file.path(OUT_DIR, "segdup_coverage_bias.csv")

# Load SegDup regions ==========================================================
cat("Loading SegDup BED...\n")
sd_gr <- import(SEGDUP, format = "BED")
seqlevelsStyle(sd_gr) <- "UCSC"
# keep autosomes only
sd_gr <- sd_gr[seqnames(sd_gr) %in% paste0("chr", 1:22)]
cat(sprintf("  SegDup: %d intervals, %.1f Mb (autosomes)\n",
            length(sd_gr), sum(width(sd_gr)) / 1e6))

# Helper: compute per-sample stats =============================================
process_sample <- function(bed_gz, label) {
  cat(sprintf("\nProcessing %s ...\n", label))
  dt <- fread(cmd = paste("zcat", bed_gz),
              col.names = c("chr", "start", "end", "cov"),
              showProgress = FALSE)

  # autosomes only
  dt <- dt[chr %in% paste0("chr", 1:22)]

  # make GRanges for bins
  bins_gr <- GRanges(seqnames = dt$chr,
                     ranges   = IRanges(dt$start + 1L, dt$end))

  # overlap with SegDup
  hits     <- findOverlaps(bins_gr, sd_gr)
  sd_idx   <- unique(queryHits(hits))
  nsd_idx  <- setdiff(seq_len(nrow(dt)), sd_idx)

  sd_cov   <- dt$cov[sd_idx]
  nsd_cov  <- dt$cov[nsd_idx]
  gw_cov   <- dt$cov

  # coverage ratio: normalise SegDup median by genome-wide median
  gw_med   <- median(gw_cov)
  sd_med   <- median(sd_cov)
  nsd_med  <- median(nsd_cov)
  ratio    <- ifelse(gw_med > 0, sd_med / gw_med, NA_real_)

  # proportion of zero-coverage bins
  zero_sd  <- mean(sd_cov  == 0)
  zero_nsd <- mean(nsd_cov == 0)

  cat(sprintf("  GW median: %.1f  |  SegDup median: %.1f  |  non-SegDup: %.1f\n",
              gw_med, sd_med, nsd_med))
  cat(sprintf("  SegDup/GW ratio: %.3f  |  zero-cov SegDup: %.1f%%  non-SegDup: %.1f%%\n",
              ratio, zero_sd * 100, zero_nsd * 100))

  data.table(
    sample      = label,
    n_bins_gw   = nrow(dt),
    n_bins_sd   = length(sd_idx),
    n_bins_nsd  = length(nsd_idx),
    gw_median   = gw_med,
    sd_median   = sd_med,
    nsd_median  = nsd_med,
    sd_gw_ratio = ratio,
    zero_pct_sd = round(zero_sd  * 100, 2),
    zero_pct_nsd= round(zero_nsd * 100, 2)
  )
}

# Enumerate tumor + normal samples =============================================
bed_files <- list.files(COV_DIR, pattern = "\\.regions\\.bed\\.gz$", full.names = TRUE)
bed_files <- bed_files[!grepl("\\.csi$", bed_files)]

results <- lapply(bed_files, function(f) {
  label <- sub("\\.regions\\.bed\\.gz$", "", basename(f))
  process_sample(f, label)
})
res_dt <- rbindlist(results)

# Summary across samples =======================================================
cat("\n\n=== SegDup Coverage Bias Summary ===\n")
cat(sprintf("Samples analysed: %d\n", nrow(res_dt)))
cat(sprintf("Median SegDup/GW coverage ratio: %.3f [%.3f, %.3f]\n",
            median(res_dt$sd_gw_ratio),
            quantile(res_dt$sd_gw_ratio, 0.25),
            quantile(res_dt$sd_gw_ratio, 0.75)))
cat(sprintf("Median zero-cov rate  SegDup: %.1f%%  non-SegDup: %.1f%%\n",
            median(res_dt$zero_pct_sd), median(res_dt$zero_pct_nsd)))

# Tumor vs Normal breakdown
res_dt[, type := ifelse(grepl("^T_", sample), "Tumor", "Normal")]
cat("\nBy sample type:\n")
print(res_dt[, .(median_ratio = round(median(sd_gw_ratio), 3),
                  IQR_ratio    = paste0("[", round(quantile(sd_gw_ratio, 0.25), 3),
                                        ", ", round(quantile(sd_gw_ratio, 0.75), 3), "]")),
             by = type])

# Save =========================================================================
fwrite(res_dt, OUT_CSV)
cat(sprintf("\nSaved: %s\n", OUT_CSV))
