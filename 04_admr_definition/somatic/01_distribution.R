#!/usr/bin/env Rscript
# 01_distribution.R
# Compute 1Mb bin density for SV breakpoints and somatic aDMR,
# SV type x CNV class cross-table, somatic aDMR width/nCG summary.
# Output: SV_aDMR/distribution_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
})

UTILS_DIR <- "/home/kachungk/script/SV-DMR/shared_file/pipeline"
source(file.path(UTILS_DIR, "shared_utils.R"))

DATA_ROOT    <- "/node200data/kachungk/hcc_data"
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")

BIN_SIZE <- 1e6L

MAIN_CHROMS <- paste0("chr", c(1:22, "X"))

# ── Load data ─────────────────────────────────────────────────────────────────
cat("[1] Loading annotated somatic aDMR...\n")
admr_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
cat(sprintf("  %d rows\n", nrow(admr_dt)))

cat("[2] Loading SV annotation...\n")
sv_dt <- fread(file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz"))
sv_dt <- sv_dt[seqnames %in% MAIN_CHROMS]
cat(sprintf("  %d SV breakpoints, %d patients\n", nrow(sv_dt), uniqueN(sv_dt$sample)))

# ── 1Mb bin density ───────────────────────────────────────────────────────────
cat("[3] Computing 1Mb bin density...\n")

make_bins <- function(chrom_lens, bin_size) {
  lapply(names(chrom_lens), function(chr) {
    starts <- seq(1L, chrom_lens[chr], by = bin_size)
    data.table(seqnames = chr, start = starts,
               end = pmin(starts + bin_size - 1L, chrom_lens[chr]),
               bin_id = paste0(chr, ":", starts))
  }) |> rbindlist()
}

bins_dt <- make_bins(CHROM_LENS[MAIN_CHROMS], BIN_SIZE)
bins_gr <- makeGRangesFromDataFrame(bins_dt, keep.extra.columns = TRUE)

# aDMR density
admr_gr <- makeGRangesFromDataFrame(admr_dt[seqnames %in% MAIN_CHROMS],
                                     seqnames.field = "seqnames")
bins_dt[, n_admr := countOverlaps(bins_gr, admr_gr)]

# SV breakpoint density (each row = one breakpoint)
sv_gr <- makeGRangesFromDataFrame(sv_dt[, .(seqnames, start, end = start + 1L)])
bins_dt[, n_sv_bp := countOverlaps(bins_gr, sv_gr)]

# SegDup overlap per bin (for track comparison)
segdup_gr <- rtracklayer::import.bed(
  "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
bins_dt[, pct_segdup := countOverlaps(bins_gr, segdup_gr) / width(bins_gr) * BIN_SIZE]

# ── SV type × CNV class cross-table ──────────────────────────────────────────
cat("[4] SV type x CNV class cross-table...\n")
sv_cross <- sv_dt[, .N, by = .(geom_type, cnv_class)]
setorder(sv_cross, geom_type, cnv_class)
cat("\nSV type x CNV class:\n")
print(dcast(sv_cross, geom_type ~ cnv_class, value.var = "N", fill = 0L))

# ── Somatic aDMR width and nCG summary ────────────────────────────────────────
cat("[5] Somatic aDMR width/nCG summary...\n")
admr_dt[, width := end - start + 1L]

width_summ <- data.frame(
  metric = c("median_width_bp","mean_width_bp","Q1_width","Q3_width",
             "p90_width","median_nCG","mean_nCG","pct_ov_bulk",
             "pct_segdup","pct_lad","pct_recur_3plus"),
  value  = c(
    median(admr_dt$width), mean(admr_dt$width),
    quantile(admr_dt$width, 0.25), quantile(admr_dt$width, 0.75),
    quantile(admr_dt$width, 0.90),
    median(admr_dt$nCG), mean(admr_dt$nCG),
    mean(admr_dt$ov_bulk, na.rm = TRUE) * 100,
    mean(admr_dt$segdup) * 100,
    mean(admr_dt$lad) * 100,
    mean(admr_dt$n_patients >= 3) * 100
  )
)
cat("\nSomatic aDMR width/nCG:\n")
print(width_summ)

# ── Save ──────────────────────────────────────────────────────────────────────
out_list <- list(
  bin_density  = bins_dt,
  sv_cross     = sv_cross,
  admr_summary = as.data.table(width_summ)
)

fwrite(bins_dt,     file.path(OUT_ROOT, "distribution_bin_density.csv.gz"))
fwrite(sv_cross,    file.path(OUT_ROOT, "distribution_sv_cross.csv"))
fwrite(as.data.table(width_summ), file.path(OUT_ROOT, "distribution_summary.csv"))

cat("\nSaved: distribution_bin_density.csv.gz, distribution_sv_cross.csv, distribution_summary.csv\n")

log_decision("01_distribution: 1Mb bin density (SV/aDMR/segdup), SV type x CNV cross-table, aDMR width/nCG summary computed")

cat("=== Done: 01_distribution.R ===\n")
