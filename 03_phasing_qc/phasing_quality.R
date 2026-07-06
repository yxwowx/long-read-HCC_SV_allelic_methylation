#!/usr/bin/env Rscript
# phasing_quality.R — A-I: No-trio phasing QC
#
# (a) whatshap stats  — block contiguity (N50, phased fraction) on the plain
#     hg38 phasing used for all aDMR HP betas (primary phasing).
# (b) whatshap compare — pairwise reference-frame stability between the plain
#     hg38 phasing and the hg38+HBV phasing (lower bound on discordance; both
#     share the clairS caller — see Methods caveat).
#
# Without a trio, (a)+(b) provide contiguity + reference-robustness; A-II
# imprinting validation provides the orthogonal biological fidelity check.
#
# Usage: mamba run -n renv Rscript post_processing/phasing_quality.R
# (whatshap is called via full path from the hifiasm env)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
VCF_HG38_DIR <- "/node200data/kachungk/hcc_data/clairS_minimap2.out_hg38/phased_vcf"
VCF_HBV_DIR  <- "/node200data/kachungk/hcc_data/hg38+HBV/clairS/phased_vcf"
MAPPING_CSV  <- path.expand("~/patient_code_mapping.csv")
OUT_DIR      <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
LOG_FILE     <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
TMP_DIR      <- tempdir()
WHATSHAP     <- "/home/kachungk/miniforge3/envs/hifiasm/bin/whatshap"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Patient mapping ────────────────────────────────────────────────────────────
pm <- fread(MAPPING_CSV)
# Samples_ID = "JJT_HCC"; VCF naming: hg38 = {Samples_ID}_normal.phased.vcf.gz
#                                       HBV  = {name}_HCC.phased.vcf.gz
pm[, name := sub("_HCC$", "", Samples_ID)]
pm[, vcf_hg38 := file.path(VCF_HG38_DIR, paste0(Samples_ID, "_normal.phased.vcf.gz"))]
pm[, vcf_hbv  := file.path(VCF_HBV_DIR,  paste0(name, "_HCC.phased.vcf.gz"))]

pm_ok <- pm[file.exists(vcf_hg38) & file.exists(vcf_hbv)]
cat(sprintf("Patients with both VCFs present: %d / %d\n", nrow(pm_ok), nrow(pm)))
if (nrow(pm_ok) == 0) stop("No matching VCF pairs found.")

# ── A. whatshap stats (hg38 primary phasing) ─────────────────────────────────
message("\nRunning whatshap stats...")
stats_list <- lapply(seq_len(nrow(pm_ok)), function(i) {
  pt  <- pm_ok$patient_code[i]
  vcf <- pm_ok$vcf_hg38[i]
  tsv <- file.path(TMP_DIR, paste0(pt, "_stats.tsv"))

  ret    <- system2(WHATSHAP, args = c("stats", "--tsv", tsv, vcf),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(ret, "status") %||% 0L

  if (status != 0L || !file.exists(tsv)) {
    warning(pt, ": whatshap stats failed (exit=", status, ")")
    return(NULL)
  }
  dt <- tryCatch(fread(tsv, skip = "#"), error = function(e) {
    # fallback: try without skip
    fread(tsv, header = TRUE)
  })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  # Rename leading '#' from header if present
  if (names(dt)[1] == "#sample") setnames(dt, "#sample", "sample_id")
  dt[, patient_code := pt]
  dt
})
stats_all <- rbindlist(Filter(Negate(is.null), stats_list), fill = TRUE)

if (nrow(stats_all) == 0) stop("whatshap stats produced no output.")

# Extract genome-wide row (chromosome == "ALL" or "all")
chr_col <- if ("chromosome" %in% names(stats_all)) "chromosome" else
           if ("#chromosome" %in% names(stats_all)) "#chromosome" else NA_character_
if (!is.na(chr_col)) {
  stats_genome <- stats_all[get(chr_col) %in% c("ALL", "all")]
  if (nrow(stats_genome) == 0)
    stats_genome <- stats_all[, .SD[.N], by = patient_code]
} else {
  stats_genome <- stats_all
}

fwrite(stats_genome, file.path(OUT_DIR, "phasing_quality_stats.csv"))
message("Wrote: phasing_quality_stats.csv (", nrow(stats_genome), " patients)")

key_cols <- intersect(c("patient_code", chr_col, "heterozygous",
                         "phased", "phased_fraction", "blocks", "block_n50"),
                      names(stats_genome))
print(stats_genome[, ..key_cols])

# ── B. whatshap compare (hg38 vs hg38+HBV) ───────────────────────────────────
# The two phasings share the same clairS caller + sample but differ in reference
# (HBV contig is an extra decoy). Autosomal het-SNP sets overlap; HBV contig
# variants simply drop out. --ignore-sample-name handles header mismatch.
message("\nRunning whatshap compare (hg38 vs hg38+HBV)...")
cmp_list <- lapply(seq_len(nrow(pm_ok)), function(i) {
  pt   <- pm_ok$patient_code[i]
  vcf1 <- pm_ok$vcf_hg38[i]
  vcf2 <- pm_ok$vcf_hbv[i]
  tsv  <- file.path(TMP_DIR, paste0(pt, "_compare.tsv"))

  ret    <- system2(WHATSHAP,
                    args = c("compare",
                             "--names", "hg38,hg38HBV",
                             "--ignore-sample-name",
                             "--tsv-pairwise", tsv,
                             vcf1, vcf2),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(ret, "status") %||% 0L

  if (status != 0L || !file.exists(tsv)) {
    warning(pt, ": whatshap compare failed (exit=", status, ")\n",
            paste(head(ret, 5), collapse = "\n"))
    return(NULL)
  }
  dt <- tryCatch(fread(tsv), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (names(dt)[1] == "#dataset_name0") setnames(dt, "#dataset_name0", "dataset_name0")
  dt[, patient_code := pt]
  dt
})
cmp_all <- rbindlist(Filter(Negate(is.null), cmp_list), fill = TRUE)

if (nrow(cmp_all) > 0) {
  fwrite(cmp_all, file.path(OUT_DIR, "phasing_quality_compare.csv"))
  message("Wrote: phasing_quality_compare.csv (", nrow(cmp_all), " patients)")
  key_cmp <- intersect(c("patient_code", "all_switch_rate", "all_switchflip_rate",
                           "blockwise_hamming_rate", "all_assessed_pairs"),
                        names(cmp_all))
  print(cmp_all[, ..key_cmp])
} else {
  message("WARNING: whatshap compare produced no output — check sample-name compatibility.")
}

# ── Cohort summaries for log ──────────────────────────────────────────────────
med_n50 <- if ("block_n50" %in% names(stats_genome))
  median(as.numeric(stats_genome$block_n50), na.rm = TRUE) else NA_real_
med_sw  <- if (nrow(cmp_all) > 0 && "all_switch_rate" %in% names(cmp_all))
  median(as.numeric(cmp_all$all_switch_rate), na.rm = TRUE) else NA_real_
med_pf  <- if ("phased_fraction" %in% names(stats_genome))
  median(as.numeric(stats_genome$phased_fraction), na.rm = TRUE) else NA_real_

cat(sprintf("\nCohort: median block_N50=%.0f bp, phased_fraction=%.3f, compare switch_rate=%.4f\n",
            med_n50, med_pf, med_sw))

cat(append = TRUE, file = LOG_FILE,
    sprintf("[%s] phasing_quality.R (A-I): %d patients; median block_N50=%.0f bp; ",
            Sys.Date(), nrow(stats_genome), med_n50))
cat(append = TRUE, file = LOG_FILE,
    sprintf("phased_frac=%.3f; compare switch_rate=%.4f\n", med_pf, med_sw))

cat("Done.\n")
