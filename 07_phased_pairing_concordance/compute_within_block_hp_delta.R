#!/usr/bin/env Rscript
# Compute within-phase-block SV-aDMR distances WITH hp_delta
# Supplements sv_admr_distance_within_block.csv.gz (which lacks hp_delta)
#
# Output: sv_admr_within_block_hp_delta.csv.gz
#   cols: patient_id, block_id, dist_bp, nearest_tier, nearest_type, hp_delta, hp_abs_diff
#
# Run: mamba run -n renv Rscript post_processing/compute_within_block_hp_delta.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
})

SV_FILE   <- "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/sv_tad_ctcf_annotation.v2.csv.gz"
ADM_FILE  <- "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz"
GTF_DIR   <- "/node200data/kachungk/hcc_data/hg38+HBV/clairS/phased_vcf"
PMAP_FILE <- "/home/kachungk/patient_code_mapping.csv"
OUTDIR    <- "/node200data/kachungk/hcc_data/DMR_SVs/03.haplotype_sv_admr_analysis"
OUT_FILE  <- file.path(OUTDIR, "sv_admr_within_block_hp_delta.csv.gz")

TIER_RECODE <- c(
  "TAD+CTCF disrupting" = "TAD_CTCF",
  "TAD-only"            = "TAD_only",
  "CTCF-only"           = "CTCF_only",
  "Copy-neutral"        = "copy_neutral",
  "Non-boundary"        = "non_boundary",
  "HBV-associated"      = "HBV_associated"
)

# ── Patient code mapping ──────────────────────────────────────────────────────
pmap <- fread(PMAP_FILE)  # cols: Samples_ID, patient_code
name2code <- setNames(pmap$patient_code, pmap$Samples_ID)

# ── Load phase blocks (one GRangesList per patient_code) ─────────────────────
message("Loading phase blocks …")
gtf_files <- list.files(GTF_DIR, pattern = "\\.gtf$", full.names = TRUE)
phase_blocks <- lapply(gtf_files, function(f) {
  gr <- import(f)
  pt_name <- str_remove(basename(f), "\\.gtf$")
  pt_code <- name2code[pt_name]
  if (is.na(pt_code)) return(NULL)
  mcols(gr)$patient_id <- pt_code
  mcols(gr)$block_id   <- mcols(gr)$gene_id
  gr[, c("patient_id", "block_id")]
}) %>%
  Filter(Negate(is.null), .) %>%
  do.call(c, .) %>%
  split(mcols(.)$patient_id)

PATIENT_IDS <- names(phase_blocks)
message(sprintf("Patients with phase blocks: %d — %s",
                length(PATIENT_IDS), paste(PATIENT_IDS, collapse = ", ")))

# ── Load aDMRs ────────────────────────────────────────────────────────────────
message("Loading aDMR data …")
admr_raw <- fread(ADM_FILE) %>%
  filter(!is.na(admr_chr), !is.na(admr_start), !is.na(admr_end)) %>%
  mutate(
    hp_delta   = HP1.Methy - HP2.Methy,
    patient_id = patient_code
  ) %>%
  filter(patient_id %in% PATIENT_IDS)

admr_gr_list <- admr_raw %>%
  dplyr::select(admr_chr, admr_start, admr_end, patient_id, hp_delta) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "admr_chr", start.field = "admr_start", end.field = "admr_end",
    keep.extra.columns = TRUE
  ) %>%
  split(mcols(.)$patient_id)

# ── Assign block_id to aDMRs via phase block overlap ─────────────────────────
message("Assigning block_id to aDMRs …")
admr_gr_list <- lapply(PATIENT_IDS, function(pt) {
  dmr_gr <- admr_gr_list[[pt]]
  blk_gr <- phase_blocks[[pt]]
  if (is.null(dmr_gr) || is.null(blk_gr) || length(dmr_gr) == 0 || length(blk_gr) == 0) return(NULL)
  hits <- findOverlaps(dmr_gr, blk_gr, select = "first")
  mcols(dmr_gr)$block_id <- blk_gr$block_id[hits]
  dmr_gr[!is.na(mcols(dmr_gr)$block_id)]
}) %>% setNames(PATIENT_IDS)

# ── Load SVs (already have PHASESETID = block_id) ────────────────────────────
message("Loading SV data …")
sv_dt <- fread(SV_FILE, select = c("seqnames","start","end","PHASESETID","sample","geom_type","stratification"))
setnames(sv_dt, c("PHASESETID","sample","geom_type"), c("block_id","patient_id","sv_type"))
sv_dt[, sv_tier := TIER_RECODE[stratification]]
sv_dt <- sv_dt[!is.na(block_id) & patient_id %in% PATIENT_IDS]
sv_raw <- sv_dt %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$patient_id)

# ── Compute within-block distances + hp_delta ─────────────────────────────────
message("Computing within-block distances …")
result <- rbindlist(lapply(PATIENT_IDS, function(pt) {
  sv_gr  <- sv_raw[[pt]]
  dmr_gr <- admr_gr_list[[pt]]
  if (is.null(sv_gr) || is.null(dmr_gr) || length(sv_gr) == 0 || length(dmr_gr) == 0) return(NULL)

  sv_gr  <- sv_gr[!is.na(mcols(sv_gr)$block_id)]
  dmr_gr <- dmr_gr[!is.na(mcols(dmr_gr)$block_id)]

  shared_blocks <- intersect(unique(mcols(dmr_gr)$block_id), unique(mcols(sv_gr)$block_id))
  if (length(shared_blocks) == 0) return(NULL)

  rbindlist(lapply(shared_blocks, function(b) {
    dmr_b <- dmr_gr[mcols(dmr_gr)$block_id == b]
    sv_b  <- sv_gr[mcols(sv_gr)$block_id == b]
    if (length(dmr_b) == 0 || length(sv_b) == 0) return(NULL)

    hits <- distanceToNearest(dmr_b, sv_b, ignore.strand = TRUE)
    if (length(hits) == 0) return(NULL)

    sv_hit <- sv_b[subjectHits(hits)]
    dm_hit <- dmr_b[queryHits(hits)]

    data.table(
      patient_id   = pt,
      block_id     = as.character(b),
      dist_bp      = mcols(hits)$distance,
      nearest_tier = as.character(mcols(sv_hit)$sv_tier),
      nearest_type = as.character(mcols(sv_hit)$sv_type),
      hp_delta     = mcols(dm_hit)$hp_delta,
      hp_abs_diff  = abs(mcols(dm_hit)$hp_delta)
    )
  }), fill = TRUE)
}), fill = TRUE)

result <- result[!is.na(dist_bp) & !is.na(hp_delta)]
message(sprintf("Done. N = %d within-block pairs across %d patients",
                nrow(result), length(unique(result$patient_id))))

cat("\n=== Summary ===\n")
cat(sprintf("  Median dist_bp    : %s bp\n", format(median(result$dist_bp), big.mark = ",")))
cat(sprintf("  Median hp_abs_diff: %.4f\n",  median(result$hp_abs_diff, na.rm = TRUE)))
cat(sprintf("  N by SV type:\n"))
print(result[, .N, by = nearest_type][order(-N)])
cat(sprintf("  N by SV tier:\n"))
print(result[, .N, by = nearest_tier][order(-N)])

# Spearman rho per type (pooled)
cat("\n=== Pooled Spearman rho (dist vs |HP delta|) ===\n")
for (t in c("DEL","DUP","INS","INV","TRA")) {
  sub <- result[nearest_type == t & !is.na(hp_abs_diff)]
  if (nrow(sub) < 10) next
  rho <- cor(sub$dist_bp, sub$hp_abs_diff, method = "spearman", use = "complete.obs")
  cat(sprintf("  %s: n=%d, rho=%.4f\n", t, nrow(sub), rho))
}

fwrite(result, OUT_FILE)
message("Saved: ", OUT_FILE)
