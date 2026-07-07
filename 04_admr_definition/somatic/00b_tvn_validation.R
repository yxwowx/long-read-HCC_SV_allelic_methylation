#!/usr/bin/env Rscript
# 00b_tvn_validation.R
# Generate Tumor-vs-Normal |HP Δβ| table for Fig 2A
#
# For each somatic aDMR locus (tumor-acquired, constitutional AMR removed):
#   Tumor:  |diff.Methy| = |HP1.Methy − HP2.Methy| from somatic_admr_annotated.csv.gz
#   Normal: |diff.Methy| at the same locus from the per-patient DSS normal aDMR files
#           (*.normal_aDMR.sorted.txt from DMR_minimap2.out_hg38/DSS/)
#           Loci with no normal coverage get the genome-wide normal background β-SD (~0.05).
#
# Output: SV_aDMR/tvn_hp_delta.csv
#   Columns: locus_id, patient_code, tissue_type (Tumor/Normal), hp_abs_diff
#
# Run:
#   mamba run -n renv Rscript 00b_tvn_validation.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(stringr)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

DATA_ROOT    <- Sys.getenv("HCC_DATA_DIR")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
NORMAL_ADMr_DIR <- file.path(DATA_ROOT, "DMR_minimap2.out_hg38/DSS")

MAIN_CHROMS   <- paste0("chr", c(1:22, "X"))
# Background |HP Δβ| assigned to somatic aDMR loci with no matching normal aDMR
# (reflects normal tissue inter-HP variability ≈ 0.05 from TCGA-LIHC data C21)
NORMAL_BG_DELTA <- 0.05

cat("=== 00b_tvn_validation.R ===\n")

# Load patient code mapping ====================================================
pmap <- fread(PATIENT_MAP_PATH)
# columns: Samples_ID (e.g. JJT_HCC), patient_code (e.g. P1)
cat(sprintf("[1] Patient map: %d entries\n", nrow(pmap)))

# Load somatic aDMR (tumor allelic methylation) ================================
cat("[2] Loading somatic aDMR (tumor)...\n")
somatic_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
somatic_dt <- somatic_dt[seqnames %in% MAIN_CHROMS]
somatic_dt[, locus_id := paste0(seqnames, ":", start, "-", end)]
somatic_dt[, hp_abs_diff := abs(diff.Methy)]
cat(sprintf("  %d somatic aDMR loci across %d patients\n",
            nrow(somatic_dt), uniqueN(somatic_dt$patient_code)))

tumor_dt <- somatic_dt[, .(locus_id, seqnames, start, end, patient_code,
                             tissue_type = "Tumor",
                             hp_abs_diff)]

# Load normal aDMR files; look up |HP Δβ| at somatic loci ======================
cat("[3] Loading normal aDMR files and matching to somatic loci...\n")
normal_files <- list.files(NORMAL_ADMr_DIR, pattern = "\\.normal_aDMR\\.sorted\\.txt$",
                            full.names = TRUE)
cat(sprintf("  Found %d normal aDMR files\n", length(normal_files)))

# Map sample name → patient_code via pmap
sample_to_pt <- setNames(pmap$patient_code, pmap$Samples_ID)

# Build a combined normal aDMR GRanges (all patients)
normal_list <- lapply(normal_files, function(f) {
  sample_name <- str_remove(basename(f), "\\.normal_aDMR\\.sorted\\.txt$")
  pt <- sample_to_pt[sample_name]
  if (is.na(pt)) {
    cat(sprintf("  Skipping %s (no patient_code match)\n", sample_name))
    return(NULL)
  }
  dt <- fread(f, header = TRUE)
  setnames(dt, c("meanMethy1","meanMethy2"), c("HP1.Methy","HP2.Methy"), skip_absent = TRUE)
  setnames(dt, c("chr","seqnames")[c("chr","seqnames") %in% names(dt)][1], "seqnames", skip_absent = TRUE)
  dt <- dt[seqnames %in% MAIN_CHROMS]
  dt[, `:=`(patient_code  = pt,
             hp_abs_diff  = abs(diff.Methy),
             tissue_type  = "Normal")]
  dt[, .(seqnames, start, end, patient_code, hp_abs_diff, tissue_type)]
})
normal_dt_all <- rbindlist(Filter(Negate(is.null), normal_list))
cat(sprintf("  %d normal aDMR loci loaded across %d patients\n",
            nrow(normal_dt_all), uniqueN(normal_dt_all$patient_code)))

# For each somatic aDMR locus, find normal |HP Δβ| at the same coordinates =====
# Strategy: per-patient nearest-overlap. If the somatic aDMR locus overlaps a
# normal aDMR, take its |HP Δβ|. Otherwise assign the background (NORMAL_BG_DELTA).
cat("[4] Matching normal |HP Δβ| at somatic loci...\n")
normal_rows <- list()
patients    <- unique(somatic_dt$patient_code)

for (pt in patients) {
  som_pt  <- somatic_dt[patient_code == pt]
  norm_pt <- normal_dt_all[patient_code == pt]

  if (nrow(som_pt) == 0) next

  # Build somatic GRanges for this patient
  som_gr <- makeGRangesFromDataFrame(
    som_pt[, .(seqnames, start, end)],
    seqnames.field = "seqnames"
  )

  if (nrow(norm_pt) == 0) {
    # No normal aDMR data for this patient → assign background
    cat(sprintf("  %s: no normal aDMR — assigning background %.2f to %d loci\n",
                pt, NORMAL_BG_DELTA, nrow(som_pt)))
    normal_rows[[pt]] <- data.table(
      locus_id     = som_pt$locus_id,
      patient_code = pt,
      tissue_type  = "Normal",
      hp_abs_diff  = NORMAL_BG_DELTA
    )
    next
  }

  norm_gr <- makeGRangesFromDataFrame(
    norm_pt[, .(seqnames, start, end, hp_abs_diff)],
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames"
  )

  # Find overlaps: somatic locus → normal aDMR (any overlap)
  hits   <- findOverlaps(som_gr, norm_gr, select = "first")
  n_matched <- sum(!is.na(hits))
  n_total   <- length(hits)

  hp_normal <- rep(NORMAL_BG_DELTA, n_total)   # default = background
  hp_normal[!is.na(hits)] <- norm_gr$hp_abs_diff[hits[!is.na(hits)]]

  cat(sprintf("  %s: %d / %d somatic loci matched to normal aDMR (%.1f%% assigned normal value)\n",
              pt, n_matched, n_total, 100*n_matched/n_total))

  normal_rows[[pt]] <- data.table(
    locus_id     = som_pt$locus_id,
    patient_code = pt,
    tissue_type  = "Normal",
    hp_abs_diff  = hp_normal
  )
}

normal_final_dt <- rbindlist(normal_rows)

# Combine tumor and normal =====================================================
cat("[5] Combining tumor and normal rows...\n")
combined_dt <- rbind(
  tumor_dt[, .(locus_id, patient_code, tissue_type, hp_abs_diff)],
  normal_final_dt
)
combined_dt[, tissue_type := factor(tissue_type, levels = c("Normal","Tumor"))]

cat(sprintf("  Combined: %d rows (%d Tumor, %d Normal)\n",
            nrow(combined_dt),
            sum(combined_dt$tissue_type == "Tumor"),
            sum(combined_dt$tissue_type == "Normal")))

# Quick sanity: median per tissue
summ <- combined_dt[, .(median_abs = median(hp_abs_diff, na.rm=TRUE),
                         mean_abs   = mean(hp_abs_diff, na.rm=TRUE),
                         n          = .N), by = tissue_type]
cat("\n=== T-vs-N summary ===\n")
print(summ)

# Save =========================================================================
out_path <- file.path(OUT_ROOT, "tvn_hp_delta.csv")
fwrite(combined_dt, out_path)
cat(sprintf("\nSaved: tvn_hp_delta.csv (%d rows)\n", nrow(combined_dt)))

cat("=== Done: 00b_tvn_validation.R ===\n")
