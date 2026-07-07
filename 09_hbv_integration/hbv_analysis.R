#!/usr/bin/env Rscript
# 12_hbv_analysis.R
# HBV integration characterization + allele-specific cis-methylation disruption
#
# Inputs:
#   --hbv_chimeric_dir : HBV/breakpoint/  (*_chimeric_breakpoint_hm_hbv.bed, T_ + N_)
#   --hbv_ins_dir      : HBV/ins/breakpoint/ (*_INS_*  + *_clipped_*, T_ + N_)
#   --sv_file          : sv_tad_ctcf_annotation.csv.gz
#   --admr_file        : 01.DMR_recurrence/confident_dmr_per_patient.csv.gz
#   --phase_vcf_dir    : hg38+HBV/clairS/phased_vcf/  (*.gtf phase blocks)
#   --patient_map      : patient_code_mapping.csv
#   --bg_hp_file       : (optional) 03.haplotype_sv_admr_analysis/all_hp_admr_tier.csv.gz
#
# Analyses:
#   A. HBV BED clustering — T_ and N_ separately; somatic filter (T_ minus N_)
#   B. SV–HBV positional match (±match_bp) with correct P-code mapping
#   C. Phase-block HP-specific |Δβ| for HBV-proximal SVs
#   D. HBV genome integration map (hbv_loc position distribution)
#   E. SV tier enrichment near HBV BK (Fisher's exact test)
#   F. dist_to_HBV_BK vs HP|Δβ| Spearman per patient
#
# Output : /node200data/kachungk/hcc_data/DMR_SVs/12.HBV_analysis/
# Run    : mamba run -n renv Rscript pipeline/12_hbv_analysis.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(optparse)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(stringr)
  library(scales)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--hbv_chimeric_dir", type = "character",
    default = "/node200data/kachungk/hcc_data/hg38+HBV/HBV/breakpoint/",
    help    = "Directory with *_chimeric_breakpoint_hm_hbv.bed files (T_ and N_)"),
  make_option("--hbv_ins_dir", type = "character",
    default = "/node200data/kachungk/hcc_data/hg38+HBV/HBV/ins/breakpoint/",
    help    = "Directory with *_INS_breakpoint_hm_hbv.bed and *_clipped_breakpoint_hm_hbv.bed"),
  make_option("--sv_file", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.csv.gz"),
  make_option("--admr_file", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz"),
  make_option("--phase_vcf_dir", type = "character",
    default = "/node200data/kachungk/hcc_data/hg38+HBV/clairS/phased_vcf",
    help    = "Directory containing per-patient *.gtf phase-block files"),
  make_option("--patient_map", type = "character",
    default = "/home/kachungk/patient_code_mapping.csv",
    help    = "CSV: Samples_ID (e.g. JJT_HCC) → patient_code (e.g. P1)"),
  make_option("--bg_hp_file", type = "character", default = NULL,
    help    = "Optional all_hp_admr_tier.csv.gz from pipeline 04 (non-HBV background)"),
  make_option("--cluster_bp", type = "integer", default = 100L,
    help    = "Gap window (bp) for greedy read clustering into a single integration locus"),
  make_option("--match_bp", type = "integer", default = 500L,
    help    = "Positional window (bp) for matching somatic HBV BK to SV records"),
  make_option("--somatic_bp", type = "integer", default = 1000L,
    help    = "T_ BK within somatic_bp of any N_ BK in same patient → germline, excluded"),
  make_option("--outdir", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/12.HBV_analysis"),
  make_option("--run_id", type = "character", default = "hbv_v1")
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR   <- opt$outdir
RUN_ID   <- opt$run_id
LOG_FILE <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

COLORS_HBV <- c(HBV_BND = "#E24B4A", non_HBV = "#3B8BD4")

# ── Helper functions ───────────────────────────────────────────────────────────

patient_map_dt <- fread(opt$patient_map)  # Samples_ID, patient_code

normalize_sampleid <- function(sid) sub("^[TN]_", "", sid)

sid_to_pcode <- function(sid_norm) {
  patient_map_dt$patient_code[match(sid_norm, patient_map_dt$Samples_ID)]
}

BED_COLS <- c("chrom", "start", "end", "hbv_loc", "sampleid", "readname")

load_bed_dir <- function(dir_path, patterns) {
  if (is.null(dir_path) || !dir.exists(dir_path)) return(NULL)
  beds <- unlist(lapply(patterns, function(p) Sys.glob(file.path(dir_path, p))))
  if (length(beds) == 0L) return(NULL)
  rbindlist(lapply(beds, function(f) {
    tryCatch(fread(f, header = FALSE, col.names = BED_COLS, sep = "\t"),
             error = function(e) { warning("Cannot read ", basename(f), ": ", e$message); NULL })
  }), fill = TRUE)
}

cluster_bed <- function(raw, window) {
  if (is.null(raw) || nrow(raw) == 0L) return(NULL)
  raw <- copy(raw)
  raw[, pos      := as.integer((start + end) %/% 2L)]
  raw[, norm_id  := normalize_sampleid(sampleid)]
  raw[, pcode    := sid_to_pcode(norm_id)]
  raw[, is_tumor := startsWith(sampleid, "T_")]
  raw <- raw[!is.na(pcode)]  # drop rows without mapping
  if (nrow(raw) == 0L) return(NULL)
  setorder(raw, pcode, chrom, pos)
  raw[, grp := {
    is_new <- c(TRUE,
                chrom[-1]    != chrom[-.N] |
                pcode[-1]    != pcode[-.N] |
                (pos[-1] - pos[-.N]) > window)
    cumsum(is_new)
  }]
  raw[, .(chrom    = chrom[1],
          pcode    = pcode[1],
          is_tumor = is_tumor[1],
          hbv_loc  = hbv_loc[1],
          pos      = as.integer(median(pos)),
          n_reads  = .N),
      by = grp][, grp := NULL][]
}

somatic_filter <- function(all_bk, somatic_bp) {
  tumor  <- all_bk[is_tumor == TRUE]
  normal <- all_bk[is_tumor == FALSE]
  if (nrow(normal) == 0L) {
    message("No N_ BED files found — all T_ loci treated as somatic.")
    return(tumor)
  }
  tumor[, is_somatic := mapply(function(pt, chr, p) {
    sub <- normal[pcode == pt & chrom == chr, pos]
    length(sub) == 0L || all(abs(p - sub) > somatic_bp)
  }, pcode, chrom, pos)]
  tumor[is_somatic == TRUE][, is_somatic := NULL][]
}

# Copied from pipeline 04 (infer SV haplotype per phase block)
infer_sv_hp_map <- function(hp_vec, block_id_vec) {
  ok           <- !is.na(block_id_vec)
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

# Copied from pipeline 04 (extract HP beta values for aDMRs in SV blocks)
get_sv_hp_beta <- function(admr_gr, sv_block_ids, sv_hp_map) {
  sv_hp_map <- sv_hp_map[!is.na(sv_hp_map) & sv_hp_map %in% c(1L, 2L)]
  if (length(sv_hp_map) == 0L) return(NULL)
  in_sv <- admr_gr$block_id %in% names(sv_hp_map)
  if (!any(in_sv)) return(NULL)
  admr_sub  <- admr_gr[in_sv]
  sv_hp_vec <- as.integer(sv_hp_map[as.character(admr_sub$block_id)])
  sv_beta   <- ifelse(sv_hp_vec == 1L, admr_sub$hp1_beta, admr_sub$hp2_beta)
  wt_beta   <- ifelse(sv_hp_vec == 1L, admr_sub$hp2_beta, admr_sub$hp1_beta)
  data.frame(
    block_id     = admr_sub$block_id,
    sv_hp        = sv_hp_vec,
    sv_hp_beta   = sv_beta,
    wt_hp_beta   = wt_beta,
    sv_minus_wt  = sv_beta - wt_beta,
    hp_delta_abs = abs(admr_sub$hp1_beta - admr_sub$hp2_beta)
  )
}


# ── A. Load + cluster HBV BED files (T_ and N_ from both directories) ─────────
message("Loading HBV BED files...")

raw_chimeric <- load_bed_dir(opt$hbv_chimeric_dir, c("*_chimeric_breakpoint_hm_hbv.bed"))
raw_ins      <- load_bed_dir(opt$hbv_ins_dir,
                             c("*_INS_breakpoint_hm_hbv.bed",
                               "*_clipped_breakpoint_hm_hbv.bed"))

raw_all <- rbindlist(list(raw_chimeric, raw_ins), fill = TRUE)
if (is.null(raw_all) || nrow(raw_all) == 0L)
  stop("No HBV BED files found. Check --hbv_chimeric_dir and --hbv_ins_dir.")

message(sprintf("Raw reads: %d chimeric, %d INS/clipped",
                if (is.null(raw_chimeric)) 0L else nrow(raw_chimeric),
                if (is.null(raw_ins))      0L else nrow(raw_ins)))

all_bk <- cluster_bed(raw_all, opt$cluster_bp)

cat("\n=== HBV breakpoints (all, pre-somatic-filter) ===\n")
print(all_bk[, .(n_loci = .N, n_reads = sum(n_reads)), by = .(pcode, is_tumor)][order(pcode)])


# ── B. Somatic filter (T_ clusters not within somatic_bp of any N_ cluster) ───
somatic_bk <- somatic_filter(all_bk, opt$somatic_bp)

cat(sprintf("\n=== Somatic HBV breakpoints (T_ minus N_ within %d bp) ===\n", opt$somatic_bp))
print(somatic_bk[, .(n_somatic_loci = .N, n_reads = sum(n_reads)), by = pcode][order(pcode)])

total_somatic <- nrow(somatic_bk)
if (total_somatic == 0L) {
  cat(append = TRUE, file = LOG_FILE,
      text = sprintf("[%s] 12 hbv_analysis: STOPPED — somatic HBV loci = 0\n", Sys.Date()))
  stop("No somatic HBV integration loci found after subtracting germline (N_) sites.")
}

fwrite(somatic_bk, file.path(OUTDIR, paste0(RUN_ID, "_somatic_hbv_loci.csv")))


# ── C. Match somatic HBV BK to SV records ─────────────────────────────────────
message("Reading SV annotation: ", opt$sv_file)
sv <- fread(opt$sv_file)

# sv$sample is already P-code (P1–P12); somatic_bk$pcode is also P-code
sv[, is_hbv_bnd := mapply(function(pt, chr, pos) {
  sub <- somatic_bk[pcode == pt & chrom == chr, pos]
  length(sub) > 0L && any(abs(pos - sub) <= opt$match_bp, na.rm = TRUE)
}, sample, seqnames, start)]

n_hbv_bnd <- sum(sv$is_hbv_bnd, na.rm = TRUE)
message(sprintf("HBV-proximal SVs: %d (±%d bp positional match)", n_hbv_bnd, opt$match_bp))

sv[, sv_tier_clean := TIER_RECODE[stratification]]

if (n_hbv_bnd == 0L) {
  cat(append = TRUE, file = LOG_FILE,
      text = sprintf("[%s] 12 hbv_analysis: STOPPED — 0 SVs matched HBV BK\n", Sys.Date()))
  stop("No SVs matched to somatic HBV integration sites. Check --match_bp or SV coordinate columns.")
}

cat("\n=== HBV-proximal SVs per patient and tier ===\n")
print(sv[is_hbv_bnd == TRUE, .N, by = .(sample, sv_tier_clean)])

hbv_bnd_ids <- sv[is_hbv_bnd == TRUE, bp_id]


# ── D. Load phase blocks ───────────────────────────────────────────────────────
message("Loading phase blocks from: ", opt$phase_vcf_dir)
phase_gtfs <- list.files(opt$phase_vcf_dir, pattern = "\\.gtf$", full.names = TRUE)
if (length(phase_gtfs) == 0L)
  stop("No GTF phase-block files found in: ", opt$phase_vcf_dir)

phase_blocks <- lapply(phase_gtfs, function(x) {
  gr <- tryCatch(import(x), error = function(e) { warning(basename(x), ": ", e$message); NULL })
  if (is.null(gr)) return(NULL)
  mcols(gr)$patient_name <- str_remove(basename(x), "\\.gtf$")
  as.data.frame(gr)
}) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows() %>%
  dplyr::select(seqnames, start, end, gene_id, patient_name) %>%
  dplyr::left_join(patient_map_dt, by = c("patient_name" = "Samples_ID")) %>%
  dplyr::rename(sample = patient_code, block_id = gene_id) %>%
  dplyr::filter(!is.na(sample)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$sample)

PATIENT_IDS <- sort(unique(sv$sample))


# ── E. Load aDMRs + assign block_id via phase block overlap ───────────────────
message("Reading aDMR file: ", opt$admr_file)
admr_raw <- fread(opt$admr_file)

# Detect coordinate and HP column names (handle raw vs renamed variants)
seq_col   <- intersect(c("admr_chr", "seqnames"),  names(admr_raw))[1]
start_col <- intersect(c("admr_start", "start"),   names(admr_raw))[1]
end_col   <- intersect(c("admr_end",   "end"),      names(admr_raw))[1]
hp1_col   <- intersect(c("hp1_beta", "HP1.Methy"), names(admr_raw))[1]
hp2_col   <- intersect(c("hp2_beta", "HP2.Methy"), names(admr_raw))[1]
pt_col    <- intersect(c("patient_code", "sample"), names(admr_raw))[1]

admr_dt <- as.data.table(admr_raw)

# Rename coordinate + HP cols to working names
setnames(admr_dt, seq_col,   "seqnames_gr")
setnames(admr_dt, start_col, "start_gr")
setnames(admr_dt, end_col,   "end_gr")
setnames(admr_dt, hp1_col,   "hp1_beta")
setnames(admr_dt, hp2_col,   "hp2_beta")
setnames(admr_dt, pt_col,    "pcode")

# Rename any remaining GRanges reserved names that are now metadata columns
# (e.g. "seqnames", "start", "end" left over when admr_chr/admr_start/admr_end were used above)
for (.rc in c("seqnames", "start", "end", "width", "strand", "element")) {
  if (.rc %in% names(admr_dt)) setnames(admr_dt, .rc, paste0("meta_", .rc))
}

admr_phased <- admr_dt %>%
  makeGRangesFromDataFrame(
    seqnames.field     = "seqnames_gr",
    start.field        = "start_gr",
    end.field          = "end_gr",
    keep.extra.columns = TRUE
  ) %>%
  split(mcols(.)$pcode)

admr_phased <- lapply(PATIENT_IDS, function(pt) {
  dmr_gr <- admr_phased[[pt]]
  blk_gr <- phase_blocks[[pt]]
  if (is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)
  if (is.null(blk_gr) || length(blk_gr) == 0) {
    mcols(dmr_gr)$block_id <- NA_character_
    return(dmr_gr)
  }
  hits <- findOverlaps(dmr_gr, blk_gr, select = "first")
  mcols(dmr_gr)$block_id <- blk_gr$block_id[hits]
  dmr_gr[!is.na(mcols(dmr_gr)$block_id)]
}) %>% setNames(PATIENT_IDS)


# ── F. HP-specific |Δβ| for HBV-proximal SVs ──────────────────────────────────
message("Computing HP-specific |Δβ| for HBV-proximal SVs...")

sv_hbv_gr_list <- sv[is_hbv_bnd == TRUE] %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames", start.field = "start", end.field = "end",
    keep.extra.columns = TRUE
  ) %>%
  split(mcols(.)$sample)

hbv_hp_results <- do.call(rbind, lapply(PATIENT_IDS, function(pt) {
  sv_hbv <- sv_hbv_gr_list[[pt]]
  dmr_gr  <- admr_phased[[pt]]
  n_sv    <- if (is.null(sv_hbv)) 0L else length(sv_hbv)

  if (n_sv == 0L || is.null(dmr_gr) || length(dmr_gr) == 0) {
    return(data.frame(patient_id = pt, n_hbv_sv = n_sv, n_admr = 0L,
                      mean_sv_beta = NA_real_, mean_wt_beta = NA_real_,
                      mean_sv_minus_wt = NA_real_, median_hp_abs = NA_real_,
                      wilcox_p = NA_real_, rank_biserial = NA_real_))
  }

  sv_hp_map <- infer_sv_hp_map(mcols(sv_hbv)$HP, mcols(sv_hbv)$PHASESETID)
  hp_df     <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)

  if (is.null(hp_df) || nrow(hp_df) < 3L) {
    return(data.frame(patient_id = pt, n_hbv_sv = n_sv, n_admr = 0L,
                      mean_sv_beta = NA_real_, mean_wt_beta = NA_real_,
                      mean_sv_minus_wt = NA_real_, median_hp_abs = NA_real_,
                      wilcox_p = NA_real_, rank_biserial = NA_real_))
  }

  wt <- tryCatch(
    wilcox.test(hp_df$sv_hp_beta, hp_df$wt_hp_beta,
                paired = TRUE, alternative = "two.sided", exact = FALSE),
    error = function(e) list(p.value = NA_real_, statistic = NA_real_)
  )
  n  <- nrow(hp_df)
  W  <- as.numeric(wt$statistic)
  rb <- if (!is.na(W)) (2 * W) / (n * (n + 1) / 2) - 1 else NA_real_

  data.frame(
    patient_id       = pt,
    n_hbv_sv         = n_sv,
    n_admr           = n,
    mean_sv_beta     = round(mean(hp_df$sv_hp_beta),      3),
    mean_wt_beta     = round(mean(hp_df$wt_hp_beta),      3),
    mean_sv_minus_wt = round(mean(hp_df$sv_minus_wt),     3),
    median_hp_abs    = round(median(hp_df$hp_delta_abs),  3),
    wilcox_p         = round(wt$p.value, 5),
    rank_biserial    = round(rb, 3)
  )
}))

if (!is.null(hbv_hp_results) && nrow(hbv_hp_results) > 0) {
  hbv_hp_results$wilcox_fdr <- p.adjust(hbv_hp_results$wilcox_p, method = "BH")
  cat("\n=== HBV-proximal SV HP-specific |Δβ| per patient ===\n")
  print(hbv_hp_results)
  fwrite(hbv_hp_results, file.path(OUTDIR, paste0(RUN_ID, "_hbv_hp_delta.csv")))
}

# Overall Wilcoxon (HBV_BND vs non-HBV background)
# Use pipeline 04 all_hp_admr_tier.csv.gz if available; otherwise skip comparison
bg_hp <- NULL
bg_hp_file <- if (!is.null(opt$bg_hp_file)) opt$bg_hp_file else
  file.path(dirname(opt$outdir), "03.haplotype_sv_admr_analysis", "all_hp_admr_tier.csv.gz")

if (file.exists(bg_hp_file)) {
  bg_hp <- fread(bg_hp_file)
  bg_hp <- bg_hp[sv_tier != "HBV_associated"]
  message(sprintf("Background HP loaded: %d non-HBV aDMR entries", nrow(bg_hp)))
}

# HBV vs non-HBV overall Wilcoxon (pooled aDMR level)
hbv_abs_db <- unlist(lapply(PATIENT_IDS, function(pt) {
  sv_hbv <- sv_hbv_gr_list[[pt]]
  dmr_gr  <- admr_phased[[pt]]
  if (is.null(sv_hbv) || length(sv_hbv) == 0 || is.null(dmr_gr)) return(NULL)
  sv_hp_map <- infer_sv_hp_map(mcols(sv_hbv)$HP, mcols(sv_hbv)$PHASESETID)
  hp_df <- get_sv_hp_beta(dmr_gr, names(sv_hp_map), sv_hp_map)
  if (!is.null(hp_df)) hp_df$hp_delta_abs else NULL
}))

bg_abs_db <- if (!is.null(bg_hp)) bg_hp$hp_abs_diff else NULL

wt_overall <- if (length(hbv_abs_db) >= 3L && length(bg_abs_db) >= 5L) {
  wilcox.test(hbv_abs_db, bg_abs_db, alternative = "two.sided", exact = FALSE)
} else NULL

if (!is.null(wt_overall))
  cat(sprintf("\nWilcoxon HBV-proximal vs non-HBV |Δβ|: p = %.4f  (n_hbv=%d, n_bg=%d)\n",
              wt_overall$p.value, length(hbv_abs_db), length(bg_abs_db)))

# Export per-aDMR HP |Δβ| for figS10 Panel C violin
if (length(hbv_abs_db) > 0L) {
  admr_export <- data.frame(hp_abs = hbv_abs_db, group = "HBV_proximal")
  if (!is.null(bg_abs_db) && length(bg_abs_db) > 0L) {
    set.seed(42L)
    bg_sample <- sample(bg_abs_db, min(5000L, length(bg_abs_db)))
    admr_export <- rbind(admr_export,
                         data.frame(hp_abs = bg_sample, group = "background"))
  }
  fwrite(admr_export, file.path(OUTDIR, paste0(RUN_ID, "_hp_delta_per_admr.csv")))
  message(sprintf("Per-aDMR HP |Δβ| export: n_hbv=%d, n_bg_sample=%d",
                  length(hbv_abs_db), sum(admr_export$group == "background")))
}

# ── G. HBV genome integration map ─────────────────────────────────────────────
# Parse hbv_loc: "chrHBV.C2:5036-5037" → genotype=chrHBV.C2, pos=5036
message("Analyzing HBV genome integration positions...")

somatic_bk[, hbv_genotype := sub(":.*", "", hbv_loc)]
somatic_bk[, hbv_pos      := suppressWarnings(
  as.integer(sub(".*:(\\d+)-.*", "\\1", hbv_loc)))]

cat("\n=== HBV integration loci by genotype ===\n")
print(somatic_bk[!is.na(hbv_pos), .(n_loci = .N, n_reads = sum(n_reads)),
                 by = .(hbv_genotype, pcode)][order(hbv_genotype, pcode)])

fwrite(somatic_bk[!is.na(hbv_pos),
                  .(hbv_genotype, hbv_pos, pcode, chrom, pos, n_reads)],
       file.path(OUTDIR, paste0(RUN_ID, "_hbv_genome_positions.csv")))


# ── H. SV tier enrichment near HBV BK (Fisher's exact test) ───────────────────
message("Testing SV tier enrichment near HBV breakpoints...")

tier_counts <- sv[!is.na(sv_tier_clean) & sv_tier_clean != "HBV_associated",
  .(n_total = .N, n_hbv = sum(is_hbv_bnd, na.rm = TRUE)), by = sv_tier_clean]
tier_counts[, n_non_hbv := n_total - n_hbv]

total_hbv_sv  <- sum(tier_counts$n_hbv)
total_nhbv_sv <- sum(tier_counts$n_non_hbv)

tier_counts[, `:=`(
  pct_hbv  = round(n_hbv / n_total * 100, 1),
  fisher_p = mapply(function(nh, nn) {
    tryCatch(
      fisher.test(matrix(c(nh, nn,
                            total_hbv_sv  - nh,
                            total_nhbv_sv - nn), 2L, 2L))$p.value,
      error = function(e) NA_real_
    )
  }, n_hbv, n_non_hbv)
)]
tier_counts[, fisher_fdr := p.adjust(fisher_p, method = "BH")]

cat("\n=== SV tier enrichment near HBV BK ===\n")
print(tier_counts[order(-pct_hbv)])
fwrite(tier_counts, file.path(OUTDIR, paste0(RUN_ID, "_tier_enrichment.csv")))


# ── I. Distance-decay: dist_to_nearest_HBV_BK vs HP|Δβ| ──────────────────────
message("Computing distance to nearest HBV BK per aDMR...")

hbv_patients <- intersect(unique(somatic_bk$pcode), PATIENT_IDS)

hbv_dist_df <- rbindlist(lapply(hbv_patients, function(pt) {
  dmr_gr <- admr_phased[[pt]]
  if (is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  bk_sub <- somatic_bk[pcode == pt & !is.na(pos)]
  if (nrow(bk_sub) == 0L) return(NULL)

  hbv_gr <- GRanges(seqnames = bk_sub$chrom,
                    ranges   = IRanges(bk_sub$pos, bk_sub$pos))

  hits <- tryCatch(
    distanceToNearest(dmr_gr, hbv_gr, ignore.strand = TRUE),
    error = function(e) NULL
  )
  if (is.null(hits) || length(hits) == 0L) return(NULL)

  dmr_sub <- dmr_gr[queryHits(hits)]
  data.table(
    pcode    = pt,
    dist_hbv = as.integer(mcols(hits)$distance),
    hp_abs   = abs(mcols(dmr_sub)$hp1_beta - mcols(dmr_sub)$hp2_beta)
  )
}), fill = TRUE)

hbv_dist_cor <- hbv_dist_df[!is.na(dist_hbv) & !is.na(hp_abs),
  .(n          = .N,
    spearman_r = if (.N >= 10L) cor(dist_hbv, hp_abs, method = "spearman") else NA_real_,
    spearman_p = if (.N >= 10L) {
      tryCatch(cor.test(dist_hbv, hp_abs, method = "spearman",
                        exact = FALSE)$p.value, error = function(e) NA_real_)
    } else NA_real_
  ), by = pcode]

hbv_dist_cor[, spearman_fdr := p.adjust(spearman_p, method = "BH")]

cat("\n=== Spearman(dist_to_HBV_BK, HP|Δβ|) per patient ===\n")
print(hbv_dist_cor)
fwrite(hbv_dist_cor, file.path(OUTDIR, paste0(RUN_ID, "_dist_spearman.csv")))


# ── Visualizations ─────────────────────────────────────────────────────────────

## Panel A — HBV genome integration map
genome_dat <- somatic_bk[!is.na(hbv_pos)]
if (nrow(genome_dat) > 0L) {
  p_genome <- ggplot(genome_dat, aes(x = hbv_pos, fill = hbv_genotype)) +
    geom_histogram(bins = 64, alpha = 0.75, position = "stack", colour = NA) +
    scale_fill_manual(
      values = c("chrHBV.C2" = "#E24B4A", "chrHBV.D" = "#3B8BD4"),
      name = "HBV genotype"
    ) +
    labs(title = "A. HBV integration site distribution",
         subtitle = sprintf("N = %d somatic loci across %d patients",
                            total_somatic, length(hbv_patients)),
         x = "HBV genome position (bp)", y = "n somatic loci") +
    theme_hcc
} else {
  p_genome <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No parsed HBV positions") + theme_void()
}

## Panel B — SV tier enrichment
p_tier <- ggplot(tier_counts, aes(x = reorder(sv_tier_clean, pct_hbv), y = pct_hbv,
                                   fill = sv_tier_clean)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%d\nFDR=%.2f", n_hbv, fisher_fdr),
                y = pct_hbv + 0.2), size = 2.8, colour = "grey30") +
  scale_fill_manual(values = SV_TIER_COLORS, guide = "none") +
  coord_flip() +
  labs(title = "B. % SVs proximal to HBV BK by tier",
       subtitle = sprintf("Match window ±%d bp; Fisher's FDR", opt$match_bp),
       x = NULL, y = "% SVs within match window of HBV BK") +
  theme_hcc

## Panel C — HP |Δβ| per patient (HBV-proximal SVs)
hp_plot_df <- hbv_hp_results[!is.na(hbv_hp_results$median_hp_abs) &
                              hbv_hp_results$n_admr >= 3L, ]
if (nrow(hp_plot_df) > 0L) {
  p_hp <- ggplot(hp_plot_df,
                 aes(x = reorder(patient_id, median_hp_abs), y = median_hp_abs,
                     fill = wilcox_fdr < 0.05)) +
    geom_col(alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#E24B4A", "FALSE" = "#888780"),
                      name = "Wilcoxon FDR < 0.05") +
    labs(title = "C. Median HP|Δβ| for HBV-proximal SV blocks",
         subtitle = "Paired Wilcoxon (SV-HP vs WT-HP); red = FDR < 0.05",
         x = "Patient", y = "Median |SV-HP β − WT-HP β|") +
    coord_flip() + theme_hcc +
    theme(legend.position = "right")
} else {
  p_hp <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Insufficient phased aDMRs in HBV SV blocks\n(n_admr < 3)") +
    theme_void()
}

## Panel D — Distance-decay scatter
if (!is.null(hbv_dist_df) && nrow(hbv_dist_df) > 50L) {
  rho_overall <- cor(hbv_dist_df$dist_hbv, hbv_dist_df$hp_abs,
                     method = "spearman", use = "complete.obs")
  sdat <- if (nrow(hbv_dist_df) > 5000L)
    hbv_dist_df[sample(.N, 5000L)] else hbv_dist_df

  p_decay <- ggplot(sdat, aes(x = dist_hbv / 1e3, y = hp_abs, colour = pcode)) +
    geom_point(alpha = 0.2, size = 0.6, shape = 16) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                colour = "black", linewidth = 0.9) +
    scale_x_log10(labels = comma_format()) +
    scale_colour_discrete(guide = "none") +
    labs(title = "D. HP|Δβ| vs distance to nearest HBV BK",
         subtitle = sprintf("Overall Spearman ρ = %.3f (all HBV+ patients pooled)",
                            rho_overall),
         x = "Distance to nearest HBV BK (kb, log₁₀)",
         y = "HP |Δβ|") +
    theme_hcc
} else {
  p_decay <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Insufficient aDMR–HBV-BK pairs for distance-decay") +
    theme_void()
}

p_combined <- (p_genome | p_tier) / (p_hp | p_decay)
ggsave(file.path(OUTDIR, paste0(RUN_ID, "_hbv_analysis.png")),
       p_combined + plot_annotation(
         title    = "HBV integration: somatic loci, SV co-localisation, and cis-methylation",
         theme    = theme(plot.title = element_text(face = "bold", size = 14))
       ),
       width = 14, height = 10, dpi = 150)


# ── Log ────────────────────────────────────────────────────────────────────────
wt_p_log <- if (!is.null(wt_overall)) round(wt_overall$p.value, 4) else NA_real_
cat(append = TRUE, file = LOG_FILE,
    text = sprintf(
      "[%s] 12 hbv_analysis: somatic_loci=%d; hbv_svs=%d; patients_hbv=%d; wilcox_p=%s; spearman_n=%d\n",
      Sys.Date(), total_somatic, n_hbv_bnd, length(hbv_patients),
      ifelse(is.na(wt_p_log), "NA", as.character(wt_p_log)),
      nrow(hbv_dist_cor[!is.na(spearman_r)])))

message("Done — outputs in: ", OUTDIR)
