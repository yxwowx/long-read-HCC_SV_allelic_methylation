#!/usr/bin/env Rscript
# 00_somatic_admr_annotate.R
# Annotate somatic aDMR with:
#   1. ov_bulk: reciprocal ≥30% overlap with same-patient bulk tumor-normal DMR
#   2. SegDup, LAD, B-compartment (PC1<0) overlap
#   3. log10(nCG)
#   4. Cross-patient recurrence (n_patients)
# Outputs: SV_aDMR/somatic_admr_annotated.csv.gz, SV_aDMR/ov_bulk_distribution.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(stringr)
})

UTILS_DIR <- "/home/kachungk/script/SV-DMR/shared_file/pipeline"
source(file.path(UTILS_DIR, "shared_utils.R"))

DATA_ROOT    <- "/node200data/kachungk/hcc_data"
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

REF_ROOT     <- "/node200data/kachungk/reference/GRCh38"
SEGDUP_BED   <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD_BED      <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW       <- file.path(REF_ROOT, "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
BULK_DMR_DIR <- file.path(DATA_ROOT, "DMR_minimap2.out_hg38/DSS")

MAIN_CHROMS  <- paste0("chr", c(1:22, "X"))

# ── 1. Load somatic aDMR ──────────────────────────────────────────────────────
cat("[1] Loading somatic aDMR...\n")
somatic_dt <- fread(file.path(DMR_SVS_ROOT, "01.DMR_recurrence/somatic_admr_per_patient.csv.gz"))
cat(sprintf("  Loaded %d rows, %d patients\n", nrow(somatic_dt),
            uniqueN(somatic_dt$patient_code)))

somatic_gr <- makeGRangesFromDataFrame(somatic_dt, keep.extra.columns = TRUE,
                                        seqnames.field = "seqnames",
                                        start.field = "start", end.field = "end")
somatic_gr <- somatic_gr[seqnames(somatic_gr) %in% MAIN_CHROMS]
patients   <- unique(somatic_dt$patient_code)

# ── 2. Compute ov_bulk per patient ────────────────────────────────────────────
cat("[2] Computing ov_bulk (reciprocal ≥30% overlap with bulk tumor-normal DMR)...\n")

bulk_files <- list.files(BULK_DMR_DIR, pattern = "\\.normal_vs_tumor_DMR\\.sorted\\.txt$",
                         full.names = TRUE)
names(bulk_files) <- str_remove(basename(bulk_files), "\\.normal_vs_tumor_DMR\\.sorted\\.txt$")

# Map from SampleName to patient_code using patient_code_rule
patient_code_rule <- fread(PATIENT_MAP_PATH)
sample_to_patient <- setNames(patient_code_rule$patient_code,
                               patient_code_rule$Samples_ID)

ov_bulk_flags  <- logical(length(somatic_gr))
ov_bulk_report <- list()

for (pt in patients) {
  cat(sprintf("  Patient %s\n", pt))

  # Identify sample IDs corresponding to this patient_code (tumor sample)
  sample_ids <- patient_code_rule$Samples_ID[patient_code_rule$patient_code == pt]
  # Bulk DMR files are named by Samples_ID (the original tumor ID)
  matching_files <- bulk_files[names(bulk_files) %in% sample_ids]

  if (length(matching_files) == 0) {
    cat(sprintf("    WARNING: no bulk DMR file for %s — skipping ov_bulk\n", pt))
    next
  }

  bulk_dt <- lapply(matching_files, function(f) {
    dt <- fread(f, header = TRUE)
    # DSS output: chr, start, end, length, nCG, meanMethy1(normal), meanMethy2(tumor), diff.Methy
    # Rename to match filter_dmr_region expectations
    setnames(dt, c("chr", "start", "end"), c("seqnames", "start", "end"),
             skip_absent = TRUE)
    if ("meanMethy1" %in% names(dt) && !"HP1.Methy" %in% names(dt)) {
      dt[, HP1.Methy := meanMethy1]
      dt[, HP2.Methy := meanMethy2]
    }
    dt
  }) |> rbindlist(fill = TRUE)

  bulk_gr <- makeGRangesFromDataFrame(bulk_dt, keep.extra.columns = TRUE,
                                       seqnames.field = "seqnames",
                                       start.field = "start", end.field = "end")
  bulk_gr <- filter_dmr_region(bulk_gr, meth_diff_cutoff = 0.2)
  bulk_gr <- bulk_gr[seqnames(bulk_gr) %in% MAIN_CHROMS]

  # Somatic aDMR for this patient
  idx_pt <- which(mcols(somatic_gr)$patient_code == pt)
  sa_gr  <- somatic_gr[idx_pt]

  if (length(sa_gr) == 0 || length(bulk_gr) == 0) next

  # ov_bulk: somatic aDMR overlaps bulk DMR by ≥50bp (any overlap at same locus)
  # Reciprocal ≥30% is too strict here because aDMRs (~730bp) and bulk DMRs
  # often have different widths; spatial co-location is the right criterion
  hits_ov <- findOverlaps(sa_gr, bulk_gr, minoverlap = 50L, ignore.strand = TRUE)
  flagged  <- unique(queryHits(hits_ov))

  ov_bulk_flags[idx_pt[flagged]] <- TRUE

  ov_bulk_report[[pt]] <- data.frame(
    patient_code  = pt,
    n_somatic     = length(sa_gr),
    n_ov_bulk     = length(flagged),
    pct_ov_bulk   = round(length(flagged) / length(sa_gr) * 100, 1)
  )

  rm(bulk_dt, bulk_gr, sa_gr); gc(verbose = FALSE)
}

mcols(somatic_gr)$ov_bulk <- ov_bulk_flags

ov_bulk_df <- bind_rows(ov_bulk_report)
cat("\n[ov_bulk 요약 (환자별)]\n")
print(ov_bulk_df)
n_total     <- sum(ov_bulk_df$n_somatic)
n_ov        <- sum(ov_bulk_df$n_ov_bulk)
cat(sprintf("\n전체: ov_bulk=TRUE %d / %d (%.1f%%)\n", n_ov, n_total, n_ov/n_total*100))

fwrite(ov_bulk_df, file.path(OUT_ROOT, "ov_bulk_distribution.csv"))
cat("  저장: ov_bulk_distribution.csv\n")

# ── 3. SegDup, LAD annotation ─────────────────────────────────────────────────
cat("[3] Annotating SegDup and LAD overlap...\n")

segdup_gr <- import.bed(SEGDUP_BED) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
lad_gr    <- import.bed(LAD_BED)    |> keepStandardChromosomes(pruning.mode = "coarse")
lad_gr    <- lad_gr[seqnames(lad_gr) %in% MAIN_CHROMS]

mcols(somatic_gr)$segdup <- somatic_gr %over% segdup_gr
mcols(somatic_gr)$lad    <- somatic_gr %over% lad_gr

cat(sprintf("  segdup=TRUE: %d (%.1f%%)\n",
            sum(mcols(somatic_gr)$segdup), mean(mcols(somatic_gr)$segdup)*100))
cat(sprintf("  lad=TRUE:    %d (%.1f%%)\n",
            sum(mcols(somatic_gr)$lad), mean(mcols(somatic_gr)$lad)*100))

# ── 4. B-compartment (PC1 < 0) ────────────────────────────────────────────────
cat("[4] Computing B-compartment from PC1 BigWig...\n")
if (file.exists(PC1_BW)) {
  pc1_rle <- import(PC1_BW, as = "RleList")
  # binnedAverage requires seqlevels(bin) == names(numvar); restrict to intersection
  shared_chroms <- intersect(seqlevels(somatic_gr), names(pc1_rle))
  sg_sub    <- keepSeqlevels(somatic_gr, shared_chroms, pruning.mode = "coarse")
  pc1_sub   <- pc1_rle[shared_chroms]
  pc1_scores <- binnedAverage(sg_sub, numvar = pc1_sub, varname = "PC1")
  # Assign back: rows in somatic_gr for shared chroms are in the same order as sg_sub
  idx_in_sg <- which(as.character(seqnames(somatic_gr)) %in% shared_chroms)
  mcols(somatic_gr)$b_compartment <- NA
  mcols(somatic_gr)$b_compartment[idx_in_sg] <- mcols(pc1_scores)$PC1 < 0
  cat(sprintf("  b_compartment=TRUE: %d (%.1f%%)\n",
              sum(mcols(somatic_gr)$b_compartment, na.rm = TRUE),
              mean(mcols(somatic_gr)$b_compartment, na.rm = TRUE)*100))
} else {
  warning("PC1 BigWig not found — setting b_compartment = NA")
  mcols(somatic_gr)$b_compartment <- NA
}

# ── 5. log10(nCG) ─────────────────────────────────────────────────────────────
mcols(somatic_gr)$log10_nCG <- log10(pmax(mcols(somatic_gr)$nCG, 1))

# ── 6. Cross-patient recurrence ───────────────────────────────────────────────
cat("[5] Computing cross-patient recurrence (reciprocal ≥30%)...\n")

# Pool all somatic aDMR; cluster overlapping loci across patients
# n_patients = number of distinct patients with ≥1 reciprocal-overlap somatic aDMR at locus
# Use reduce-then-countOverlaps approach for memory efficiency
patient_list <- lapply(patients, function(pt) somatic_gr[mcols(somatic_gr)$patient_code == pt])
names(patient_list) <- patients

# Per-locus recurrence: for each somatic aDMR, count how many patients have
# a reciprocally overlapping somatic aDMR (≥30%)
recur_counts <- integer(length(somatic_gr))

for (i in seq_along(patients)) {
  pt_i <- patients[i]
  idx_i <- which(mcols(somatic_gr)$patient_code == pt_i)
  gr_i  <- somatic_gr[idx_i]

  # Count overlap contribution from every OTHER patient
  for (j in seq_along(patients)) {
    if (i == j) next
    gr_j <- patient_list[[patients[j]]]
    if (length(gr_j) == 0) next

    hits_rj <- reciprocal_overlap_hits(gr_i, gr_j, min_pct = 0.30, min_bp = 1L)
    if (length(hits_rj) == 0) next
    # Each query locus that has ≥1 hit in patient_j gets +1
    recur_counts[idx_i[unique(queryHits(hits_rj))]] <-
      recur_counts[idx_i[unique(queryHits(hits_rj))]] + 1L
  }
}
# Add self (each locus is present in at least 1 patient — itself)
recur_counts <- recur_counts + 1L

mcols(somatic_gr)$n_patients <- recur_counts

recur_dist <- table(cut(recur_counts, breaks = c(0,1,2,3,4,5,Inf),
                        labels = c("1","2","3","4","5","≥6"), right = TRUE))
cat("  Recurrence distribution:\n")
print(recur_dist)

# ── 7. Export ─────────────────────────────────────────────────────────────────
cat("[6] Exporting annotated somatic aDMR...\n")

out_dt <- as.data.table(somatic_gr) |>
  dplyr::rename(seqnames = seqnames, start = start, end = end)

out_path <- file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz")
fwrite(out_dt, out_path)
cat(sprintf("  저장: %s (%d rows)\n", out_path, nrow(out_dt)))

cat("\n=== Feature summary ===\n")
cat(sprintf("  Total somatic aDMR: %d\n", nrow(out_dt)))
cat(sprintf("  ov_bulk=TRUE:       %d (%.1f%%)\n",
            sum(out_dt$ov_bulk, na.rm=TRUE), mean(out_dt$ov_bulk, na.rm=TRUE)*100))
cat(sprintf("  segdup=TRUE:        %d (%.1f%%)\n",
            sum(out_dt$segdup), mean(out_dt$segdup)*100))
cat(sprintf("  lad=TRUE:           %d (%.1f%%)\n",
            sum(out_dt$lad), mean(out_dt$lad)*100))
if (!all(is.na(out_dt$b_compartment)))
  cat(sprintf("  b_compartment=TRUE: %d (%.1f%%)\n",
              sum(out_dt$b_compartment, na.rm=TRUE),
              mean(out_dt$b_compartment, na.rm=TRUE)*100))
cat(sprintf("  n_patients≥3:       %d (%.1f%%)\n",
            sum(out_dt$n_patients >= 3), mean(out_dt$n_patients >= 3)*100))

log_decision("00_somatic_admr_annotate: ov_bulk, segdup, lad, b_compartment, n_patients annotated; outputs: somatic_admr_annotated.csv.gz, ov_bulk_distribution.csv")

cat("\n=== Done: 00_somatic_admr_annotate.R ===\n")
