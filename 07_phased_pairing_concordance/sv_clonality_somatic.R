#!/usr/bin/env Rscript
# sv_clonality_somatic.R — Somatic-context CCF analysis using sv_admr_analysis data
#
# Re-computes:
#   A. SegDup SV CCF comparison using SVclone CCF (somatic context)
#   B. Clonal SV (CCF >= 0.8) HYPO-concordance from phaseblock_pairs.csv
#
# Note: direction_match column in phaseblock_pairs.csv is buggy (all TRUE).
#       HYPO-concordance is derived directly from sv_minus_wt < 0.
#
# Outputs:
#   SV_aDMR/result/sv_clonality_somatic_ccf.csv
#   SV_aDMR/result/sv_clonality_somatic_stats.txt

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

set.seed(42)

BASE_DMR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")
BASE_ADM  <- file.path(Sys.getenv("HCC_DATA_DIR"), "SV_aDMR")
OUT_DIR   <- file.path(BASE_ADM, "result")
dir.create(OUT_DIR, showWarnings = FALSE)

# 1. Load phaseblock_pairs =====================================================
message("Loading phaseblock_pairs ...")
pairs <- fread(file.path(BASE_ADM, "phaseblock_pairs.csv"))
# Fix: derive HYPO-concordance from sv_minus_wt (direction_match column is buggy)
pairs[, hypo_concordant := sv_minus_wt < 0]
cat(sprintf("pairs: n=%d, HYPO=%.1f%%\n", nrow(pairs), 100*mean(pairs$hypo_concordant)))

# 2. Load SVclone CCF (sv_clonality_fragility.csv) =============================
message("Loading SVclone CCF ...")
ccf_dt <- fread(file.path(BASE_DMR, "result/sv_clonality_fragility.csv"))
ccf_dt <- ccf_dt[!is.na(CCF)]
cat(sprintf("SVclone CCF: %d SVs with valid CCF\n", nrow(ccf_dt)))

# 3. Load SV annotation to get PHASESETID -> bp_id linkage =====================
message("Loading SV annotation for PHASESETID linkage ...")
sv_annot <- fread(
  file.path(BASE_DMR, "sv_tad_ctcf_annotation.csv.gz"),
  select = c("bp_id", "PHASESETID", "sample", "seqnames", "start")
)
# Keep only SVs that have CCF in sv_clonality_fragility
sv_annot_ccf <- merge(
  sv_annot[!is.na(PHASESETID)],
  ccf_dt[, .(bp_id, CCF, timing, segdup_overlap)],
  by = "bp_id", all = FALSE
)
sv_annot_ccf[, block_id := as.character(PHASESETID)]
sv_annot_ccf[, patient_code := as.character(sample)]
cat(sprintf("SV annotation with CCF: %d rows\n", nrow(sv_annot_ccf)))

# 4. Join pairs to CCF via block_id + patient_code =============================
message("Joining pairs to SVclone CCF via phase block ...")
pairs[, block_id_chr := as.character(block_id)]

# For each pair, find matching SVs in same block × patient
# Then pick the SV whose position is closest to (admr_start ± bp_dist)
sv_lookup <- sv_annot_ccf[, .(
  block_id = block_id, patient_code = patient_code,
  sv_seqnames = seqnames, sv_start = as.integer(start),
  CCF = CCF, timing = timing, segdup_overlap = segdup_overlap
)]

pairs_joined <- sv_lookup[pairs, on = .(block_id = block_id_chr, patient_code),
                           allow.cartesian = TRUE]

# For same-chr SVs: pick SV nearest to admr_start at distance = bp_dist
pairs_joined[sv_seqnames == seqnames,
             dist_err := abs(abs(as.integer(start) - sv_start) - as.integer(bp_dist))]

# For each original pair (identified by row position), keep the best SV match
# Use seqnames + start + patient_code + block_id as pair key
pairs_joined[, pair_id := .GRP, by = .(patient_code, block_id, seqnames, start)]

best_match <- pairs_joined[
  !is.na(dist_err),
  .SD[which.min(dist_err)],
  by = pair_id
]

# For pairs with no same-chr SV match, use any SV in the block (TRA etc.)
unmatched_ids <- setdiff(pairs_joined[, unique(pair_id)],
                          best_match[, unique(pair_id)])
fallback <- pairs_joined[pair_id %in% unmatched_ids,
                          .SD[1], by = pair_id]

pairs_ccf <- rbind(best_match, fallback, fill = TRUE)
pairs_ccf <- pairs_ccf[order(pair_id)]

n_ccf <- sum(!is.na(pairs_ccf$CCF))
cat(sprintf("Pairs with SVclone CCF: %d / %d (%.1f%%)\n",
            n_ccf, nrow(pairs_ccf), 100 * n_ccf / nrow(pairs_ccf)))

# Sanity: should have same rows as input pairs
cat(sprintf("Output pairs: %d (input: %d)\n", nrow(pairs_ccf), nrow(pairs)))

# 5. Analysis A: SegDup SV CCF comparison ======================================
message("\n=== Analysis A: SegDup SV CCF comparison ===")
segdup_dt <- ccf_dt[, .(segdup_overlap, CCF)]
segdup_dt[, segdup_lab := ifelse(segdup_overlap, "SegDup SV", "Non-SegDup SV")]

wt_seg <- wilcox.test(CCF ~ segdup_overlap, data = segdup_dt, exact = FALSE)

sum_seg <- segdup_dt[!is.na(segdup_overlap), .(
  n          = .N,
  median_CCF = round(median(CCF, na.rm = TRUE), 3),
  pct_clonal = round(100 * mean(CCF >= 0.8, na.rm = TRUE), 1)
), by = segdup_lab]

cat("\nSegDup CCF summary (SVclone, somatic context):\n"); print(sum_seg)
med_seg    <- segdup_dt[segdup_lab == "SegDup SV",     median(CCF, na.rm=TRUE)]
med_nonseg <- segdup_dt[segdup_lab == "Non-SegDup SV", median(CCF, na.rm=TRUE)]
cat(sprintf("Median SegDup CCF:     %.3f\n", med_seg))
cat(sprintf("Median Non-SegDup CCF: %.3f\n", med_nonseg))
cat(sprintf("Wilcoxon p:            %.4g\n", wt_seg$p.value))

# 6. Analysis B: HYPO-concordance — all pairs ==================================
message("\n=== Analysis B: HYPO-concordance ===")
all_n    <- nrow(pairs)
all_hypo <- mean(pairs$hypo_concordant, na.rm = TRUE)
cat(sprintf("All pairs: n=%d, HYPO-concordance=%.1f%%\n", all_n, 100*all_hypo))

# Tier breakdown
tier_conc <- pairs[, .(n=.N, hypo=sum(hypo_concordant),
                        rate=round(100*mean(hypo_concordant),1)), by=tier]
cat("Tier concordance:\n"); print(tier_conc)

# Per distance bin (all pairs)
dist_bins   <- c(0, 50000, 200000, 500000, Inf)
bin_labels  <- c("<=50kb", "50-200kb", "200-500kb", ">500kb")
pairs[, dist_bin := cut(bp_dist, breaks=dist_bins, labels=bin_labels,
                         right=TRUE, include.lowest=TRUE)]
all_dist <- pairs[!is.na(dist_bin), .(
  n=.N, hypo=sum(hypo_concordant),
  rate=round(100*mean(hypo_concordant),1)
), by=dist_bin][order(dist_bin)]
cat("\nAll pairs, concordance by distance bin:\n"); print(all_dist)

# 7. Analysis B2: Clonal SV concordance (CCF >= 0.8) ===========================
message("\n=== Analysis B2: Clonal SV concordance (CCF >= 0.8) ===")
clonal_pairs <- pairs_ccf[!is.na(CCF) & CCF >= 0.8]
n_cl   <- nrow(clonal_pairs)
p_cl   <- mean(clonal_pairs$hypo_concordant, na.rm=TRUE)
cat(sprintf("Clonal pairs (CCF>=0.8): n=%d, HYPO-concordance=%.1f%%\n", n_cl, 100*p_cl))

# GLMER
or_cl <- NA; p_glmer <- NA
if (n_cl > 20 && clonal_pairs[, uniqueN(patient_code)] >= 3) {
  clonal_pairs[, hypo_int := as.integer(hypo_concordant)]
  m <- tryCatch(
    glmer(hypo_int ~ 1 + (1 | patient_code), data=clonal_pairs, family=binomial),
    error = function(e) { message("GLMER error: ", e$message); NULL }
  )
  if (!is.null(m)) {
    cf   <- fixef(m)["(Intercept)"]
    or_cl   <- exp(cf)
    su   <- summary(m)$coefficients
    z_val <- su["(Intercept)", "z value"]
    p_glmer <- 2 * pnorm(-abs(z_val))
    cat(sprintf("GLMER: OR=%.3f, p=%.4g\n", or_cl, p_glmer))
  }
} else {
  bt <- binom.test(sum(clonal_pairs$hypo_concordant, na.rm=TRUE), n_cl)
  p_glmer <- bt$p.value
  cat(sprintf("Binomial test: p=%.4g\n", p_glmer))
}

# Distance bin (clonal)
clonal_pairs[, dist_bin := cut(bp_dist, breaks=dist_bins, labels=bin_labels,
                                 right=TRUE, include.lowest=TRUE)]
clonal_dist <- clonal_pairs[!is.na(dist_bin), .(
  n=.N, hypo=sum(hypo_concordant),
  rate=round(100*mean(hypo_concordant),1)
), by=dist_bin][order(dist_bin)]
cat("\nClonal pairs by distance bin:\n"); print(clonal_dist)

p_50kb_cl <- clonal_pairs[bp_dist <= 50000,
                            round(100*mean(hypo_concordant, na.rm=TRUE), 1)]
n_50kb_cl <- clonal_pairs[bp_dist <= 50000, .N]
cat(sprintf("Clonal <=50kb: n=%d, concordance=%.1f%%\n", n_50kb_cl, p_50kb_cl))

# 8. Save outputs ==============================================================
fwrite(pairs_ccf, file.path(OUT_DIR, "sv_clonality_somatic_ccf.csv"))

sink(file.path(OUT_DIR, "sv_clonality_somatic_stats.txt"))
cat("=== sv_clonality_somatic.R — Statistics for Manuscript Update ===\n\n")
cat("--- A. SegDup SV CCF (SVclone, somatic context) ---\n")
print(sum_seg)
cat(sprintf("\nMedian SegDup CCF:     %.3f\n", med_seg))
cat(sprintf("Median Non-SegDup CCF: %.3f\n", med_nonseg))
cat(sprintf("Wilcoxon p:            %.4g\n", wt_seg$p.value))
cat("\n--- B. HYPO-concordance (all pairs) ---\n")
cat(sprintf("All pairs: n=%d, concordance=%.1f%%\n", all_n, 100*all_hypo))
print(all_dist)
cat("\n--- B2. Clonal SV concordance (CCF >= 0.8) ---\n")
cat(sprintf("Clonal (CCF>=0.8): n=%d, concordance=%.1f%%\n", n_cl, 100*p_cl))
if (!is.na(or_cl)) cat(sprintf("GLMER: OR=%.3f, p=%.4g\n", or_cl, p_glmer))
print(clonal_dist)
cat(sprintf("Clonal <=50kb: n=%d, concordance=%.1f%%\n", n_50kb_cl, p_50kb_cl))
sink()
message("Wrote: sv_clonality_somatic_stats.txt")

message("Done.")
