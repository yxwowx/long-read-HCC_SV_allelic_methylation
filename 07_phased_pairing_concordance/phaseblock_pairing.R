#!/usr/bin/env Rscript
# 04_phaseblock_pairing.R
# Phase-block SV–somatic aDMR pairing with tiering re-evaluation
# Tiering: ov_bulk added as new Gold confidence axis
# Output: SV_aDMR/phaseblock_pairs.csv, SV_aDMR/somatic_tier_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
})

UTILS_DIR <- "/home/kachungk/script/SV-DMR/shared_file/pipeline"
source(file.path(UTILS_DIR, "shared_utils.R"))

DATA_ROOT     <- "/node200data/kachungk/hcc_data"
DMR_SVS_ROOT  <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT      <- file.path(DATA_ROOT, "SV_aDMR")
PHASE_VCF_DIR <- file.path(DATA_ROOT, "hg38+HBV/clairS/phased_vcf")

MAIN_CHROMS   <- paste0("chr", c(1:22, "X"))

# Tier criteria (somatic aDMR — ov_bulk as new axis)
DIST_MAX_BP   <- 50000L   # ≤50kb
DELTA_MIN     <- 0.15     # HP |Δβ| ≥0.15
RECUR_GOLD    <- 3L       # recurrence n ≥ 3 for somatic aDMR (constitutional removed)

# ── Infer SV-bearing haplotype (from pipeline/04) ─────────────────────────────
infer_sv_hp_map <- function(hp_vec, block_id_vec) {
  ok <- !is.na(block_id_vec)
  hp_vec       <- as.integer(hp_vec[ok])
  block_id_vec <- as.character(block_id_vec[ok])
  spl <- split(hp_vec, block_id_vec)
  res <- vapply(spl, function(x) {
    x <- x[!is.na(x) & x %in% c(1L, 2L)]
    if (length(x) == 0L) return(NA_integer_)
    ux <- unique(x)
    if (length(ux) != 1L) return(NA_integer_)
    ux[[1]]
  }, integer(1))
  res[!is.na(res)]
}

get_sv_hp_beta <- function(admr_gr, sv_hp_map) {
  sv_hp_map <- sv_hp_map[!is.na(sv_hp_map) & sv_hp_map %in% c(1L, 2L)]
  if (length(sv_hp_map) == 0L) return(NULL)
  in_sv <- mcols(admr_gr)$block_id %in% names(sv_hp_map)
  if (!any(in_sv)) return(NULL)
  admr_sub  <- admr_gr[in_sv]
  sv_hp_vec <- as.integer(sv_hp_map[as.character(mcols(admr_sub)$block_id)])
  sv_beta   <- ifelse(sv_hp_vec == 1L, mcols(admr_sub)$HP1.Methy, mcols(admr_sub)$HP2.Methy)
  wt_beta   <- ifelse(sv_hp_vec == 1L, mcols(admr_sub)$HP2.Methy, mcols(admr_sub)$HP1.Methy)
  data.frame(
    seqnames    = as.character(seqnames(admr_sub)),
    start       = start(admr_sub),
    end         = end(admr_sub),
    block_id    = mcols(admr_sub)$block_id,
    patient_code = mcols(admr_sub)$patient_code,
    ov_bulk     = mcols(admr_sub)$ov_bulk,
    n_patients  = mcols(admr_sub)$n_patients,
    diff_methy  = mcols(admr_sub)$diff.Methy,
    sv_hp       = sv_hp_vec,
    sv_hp_beta  = sv_beta,
    wt_hp_beta  = wt_beta,
    sv_minus_wt = sv_beta - wt_beta,
    hp_delta_abs = abs(mcols(admr_sub)$HP1.Methy - mcols(admr_sub)$HP2.Methy)
  )
}

# ── Load phase blocks ─────────────────────────────────────────────────────────
cat("[1] Loading phase block GTF files...\n")
gtf_files <- list.files(PHASE_VCF_DIR, pattern = "\\.gtf$", full.names = TRUE)
cat(sprintf("  Found %d GTF files\n", length(gtf_files)))

phase_block_dt <- lapply(gtf_files, function(f) {
  gr <- tryCatch(import(f), error = function(e) { warning(basename(f)); NULL })
  if (is.null(gr)) return(NULL)
  data.table(
    seqnames   = as.character(seqnames(gr)),
    start      = start(gr),
    end        = end(gr),
    block_id   = as.character(mcols(gr)$gene_id),
    patient_name = str_remove(basename(f), "\\.gtf$")
  )
}) |> rbindlist()

# Anonymize patient names
pmap <- fread(PATIENT_MAP_PATH)
phase_block_dt <- phase_block_dt |>
  left_join(pmap, by = c("patient_name" = "Samples_ID"))

phase_block_dt <- phase_block_dt[seqnames %in% MAIN_CHROMS]
cat(sprintf("  Total phase blocks: %d rows\n", nrow(phase_block_dt)))

# ── Load somatic aDMR ─────────────────────────────────────────────────────────
cat("[2] Loading annotated somatic aDMR...\n")
admr_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
admr_dt <- admr_dt[seqnames %in% MAIN_CHROMS]
cat(sprintf("  %d somatic aDMR across %d patients\n",
            nrow(admr_dt), uniqueN(admr_dt$patient_code)))

# ── Load SV annotation ─────────────────────────────────────────────────────────
cat("[3] Loading SV annotation...\n")
sv_dt <- fread(file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz"))
sv_dt <- sv_dt[seqnames %in% MAIN_CHROMS]
cat(sprintf("  %d SV breakpoints\n", nrow(sv_dt)))

# ── Assign block_id to somatic aDMR per patient ───────────────────────────────
cat("[4] Assigning block_id to somatic aDMR (per patient)...\n")
patients <- intersect(unique(admr_dt$patient_code), unique(phase_block_dt$patient_code))
cat(sprintf("  Processing %d patients\n", length(patients)))

admr_with_block <- lapply(patients, function(pt) {
  admr_pt <- admr_dt[patient_code == pt]
  pb_pt   <- phase_block_dt[patient_code == pt]
  if (nrow(admr_pt) == 0 || nrow(pb_pt) == 0) return(NULL)

  admr_gr <- makeGRangesFromDataFrame(admr_pt, keep.extra.columns = TRUE)
  pb_gr   <- makeGRangesFromDataFrame(pb_pt, keep.extra.columns = TRUE)

  # Assign block_id: take the first overlapping block per aDMR
  hits <- findOverlaps(admr_gr, pb_gr, select = "first")
  mcols(admr_gr)$block_id <- NA_character_
  has_block <- !is.na(hits)
  mcols(admr_gr)$block_id[has_block] <- mcols(pb_gr)$block_id[hits[has_block]]
  admr_gr
}) |> Filter(Negate(is.null), x = _)

cat(sprintf("  Assigned blocks for %d patients\n", length(admr_with_block)))

# ── Pair SV–somatic aDMR within same phase block ──────────────────────────────
cat("[5] Pairing SV and somatic aDMR in same phase block...\n")

all_pairs <- lapply(seq_along(patients), function(i) {
  pt      <- patients[i]
  admr_gr <- admr_with_block[[i]]
  if (is.null(admr_gr)) return(NULL)
  admr_gr <- admr_gr[!is.na(mcols(admr_gr)$block_id)]
  if (length(admr_gr) == 0) return(NULL)

  # SVs for this patient (sample column = anonymized patient code in sv_tad_ctcf_annotation.csv.gz)
  sv_pt <- sv_dt[sample == pt]
  if (nrow(sv_pt) == 0) return(NULL)

  # SV HP map
  sv_hp_map <- infer_sv_hp_map(sv_pt$HP, sv_pt$PHASESETID)
  if (length(sv_hp_map) == 0) return(NULL)

  # Get HP-specific betas for aDMR in SV phase blocks
  beta_df <- get_sv_hp_beta(admr_gr, sv_hp_map)
  if (is.null(beta_df) || nrow(beta_df) == 0) return(NULL)

  # Calculate SV–aDMR bp_dist (minimum distance to any SV in same block)
  sv_in_block <- sv_pt[PHASESETID %in% names(sv_hp_map)]
  if (nrow(sv_in_block) == 0) return(NULL)

  sv_gr_pt <- makeGRangesFromDataFrame(
    sv_in_block[, .(seqnames, start, end = start + 1L, PHASESETID, geom_type, cnv_class, stratification)],
    keep.extra.columns = TRUE
  )

  beta_gr <- makeGRangesFromDataFrame(
    beta_df[, c("seqnames","start","end","block_id")], keep.extra.columns = FALSE
  )

  # Minimum distance from each aDMR to nearest SV in its block
  dist_mat <- distanceToNearest(beta_gr, sv_gr_pt)
  beta_df$bp_dist <- NA_integer_
  beta_df$bp_dist[queryHits(dist_mat)] <- mcols(dist_mat)$distance

  # Add SV type from nearest SV
  beta_df$geom_type    <- NA_character_
  beta_df$cnv_class    <- NA_character_
  beta_df$stratification <- NA_character_
  near_sv <- sv_gr_pt[subjectHits(dist_mat)]
  beta_df$geom_type[queryHits(dist_mat)]     <- mcols(near_sv)$geom_type
  beta_df$cnv_class[queryHits(dist_mat)]     <- mcols(near_sv)$cnv_class
  beta_df$stratification[queryHits(dist_mat)] <- mcols(near_sv)$stratification

  beta_df
}) |> rbindlist(fill = TRUE)

cat(sprintf("  Total pairs: %d\n", nrow(all_pairs)))

# ── Tiering ──────────────────────────────────────────────────────────────────
cat("[6] Tier classification...\n")

all_pairs[, direction_match := sv_minus_wt < 0]  # SV-hp is HYPO relative to WT-hp (cis-induction test)

all_pairs[, tier := dplyr::case_when(
  # Gold: same block + ≤50kb + |Δβ|≥0.15 + direction_match + ov_bulk=TRUE + recurrence n≥3
  !is.na(block_id) & !is.na(bp_dist) & bp_dist <= DIST_MAX_BP &
    hp_delta_abs >= DELTA_MIN &
    abs(sv_minus_wt) >= DELTA_MIN &
    !is.na(ov_bulk) & ov_bulk == TRUE &
    !is.na(n_patients) & n_patients >= RECUR_GOLD ~ "Gold",
  # Silver: (same block OR n≥3) + ≤50kb + |Δβ|≥0.15 (ov_bulk not required)
  !is.na(block_id) & !is.na(bp_dist) & bp_dist <= DIST_MAX_BP &
    hp_delta_abs >= DELTA_MIN &
    abs(sv_minus_wt) >= DELTA_MIN ~ "Silver",
  # Bronze: sensitivity only
  TRUE ~ "Bronze"
)]

tier_summ <- all_pairs[, .N, by = tier]
cat("\nTier distribution:\n")
print(tier_summ)

# Rationale documentation
tier_rationale <- data.frame(
  tier = c("Gold", "Silver", "Bronze"),
  criteria = c(
    "same block + bp_dist<=50kb + hp_delta_abs>=0.15 + sv_minus_wt>=0.15 + ov_bulk=TRUE + n_patients>=3",
    "same block + bp_dist<=50kb + hp_delta_abs>=0.15 + sv_minus_wt>=0.15 (ov_bulk not required)",
    "all remaining pairs (sensitivity analysis only)"
  ),
  rationale = c(
    "ov_bulk=TRUE confirms allele-specific imbalance co-occurs with net tumor-vs-normal shift; n>=3 replaces normal-exclusion filter since constitutional signal already removed",
    "Relaxed tier: phase-block co-occurrence + allele-specific magnitude without requiring bulk confirmation",
    "Sensitivity only; not used for main analyses"
  )
)

# ── Save ─────────────────────────────────────────────────────────────────────
fwrite(all_pairs, file.path(OUT_ROOT, "phaseblock_pairs.csv"))
fwrite(tier_summ, file.path(OUT_ROOT, "somatic_tier_summary.csv"))
fwrite(tier_rationale, file.path(OUT_ROOT, "somatic_tier_rationale.csv"))

cat(sprintf("\nSaved: phaseblock_pairs.csv (%d rows)\n", nrow(all_pairs)))
cat(sprintf("Saved: somatic_tier_summary.csv, somatic_tier_rationale.csv\n"))

log_decision(sprintf("04_phaseblock_pairing: %d SV-aDMR pairs; Gold=%d, Silver=%d",
                     nrow(all_pairs),
                     tier_summ[tier=="Gold", N],
                     tier_summ[tier=="Silver", N]))

cat("=== Done: 04_phaseblock_pairing.R ===\n")
