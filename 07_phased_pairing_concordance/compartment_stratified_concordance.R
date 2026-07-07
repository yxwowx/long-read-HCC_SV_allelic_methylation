#!/usr/bin/env Rscript
# compartment_stratified_concordance.R
#
# Re-derives the compartment-stratified haplotype concordance test against the
# CURRENT somatic aDMR phase-block pairing (phaseblock_pairs.csv, n=39,212;
# produced by 04_phaseblock_pairing.R), joining pairs to HepG2 Micro-C PC1
# compartment calls.
#
# This supersedes remodeled_constitutional_AMR/post_processing/
# concordance_distance_diagnostic.R (A_HP3), which performs the same join but
# against the older, explicitly-superseded n=8,616 pairing
# (sv_admr_hp_concordance_pairs.csv.gz) and therefore does not reproduce the
# manuscript's reported numbers. See docs/OPEN_ISSUES.md issue #1 for the
# verification that this script reproduces:
#   B-compartment: 53.9%, n=22,056, GLMER OR=1.172 [1.045-1.316], p=0.0068
#   A-compartment: 49.5%, n=14,781, GLMER OR=0.990 [0.917-1.068], p=0.79
#   (2,375 pairs unclassified / no PC1 bin overlap)
#
# Run: mamba run -n renv Rscript compartment_stratified_concordance.R

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(lme4)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

PAIRS_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "SV_aDMR/phaseblock_pairs.csv")
PC1_BW     <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
OUTFILE    <- file.path(Sys.getenv("HCC_DATA_DIR"), "SV_aDMR/result/compartment_stratified_concordance.csv")

pairs <- fread(PAIRS_FILE)
# direction_match in phaseblock_pairs.csv is buggy (always TRUE by construction
# of the tier-selection filter); recompute concordance directly from sv_minus_wt.
pairs[, concordant := as.integer(sv_minus_wt < 0)]
N <- nrow(pairs)
cat(sprintf("Loaded %d pairs across %d patients\n", N, uniqueN(pairs$patient_code)))

pc1_bw   <- import(PC1_BW, as = "GRanges")
pairs_gr <- makeGRangesFromDataFrame(pairs, keep.extra.columns = TRUE)
hits     <- findOverlaps(pairs_gr, pc1_bw, select = "first")
pairs[, pc1_score := pc1_bw$score[hits]]
pairs[, compartment := fifelse(is.na(pc1_score), NA_character_,
                                fifelse(pc1_score < 0, "B", "A"))]
cat(sprintf("Compartment-classified: %d / %d (unclassified: %d)\n",
            sum(!is.na(pairs$compartment)), N, sum(is.na(pairs$compartment))))

comp_dt <- pairs[!is.na(compartment), {
  bt <- binom.test(sum(concordant), .N, p = 0.5)
  list(n = .N, pct_concordant = round(100 * mean(concordant), 1), binom_p = bt$p.value)
}, by = compartment]
cat("\nNaive compartment concordance:\n"); print(comp_dt)

glmer_by_compartment <- rbindlist(lapply(c("A", "B"), function(comp) {
  sub <- pairs[compartment == comp]
  mm  <- glmer(concordant ~ 1 + (1 | patient_code), data = sub, family = binomial,
               control = glmerControl(optimizer = "bobyqa"))
  cf  <- fixef(mm); ci <- confint(mm, parm = "(Intercept)", method = "Wald")
  data.table(compartment = comp, n = nrow(sub),
             pct_concordant = round(100 * mean(sub$concordant), 1),
             OR = exp(cf["(Intercept)"]), OR_lo = exp(ci[1]), OR_hi = exp(ci[2]),
             p = summary(mm)$coefficients["(Intercept)", "Pr(>|z|)"])
}))
cat("\nPatient-blocked GLMER by compartment:\n"); print(glmer_by_compartment)

dir.create(dirname(OUTFILE), showWarnings = FALSE, recursive = TRUE)
fwrite(glmer_by_compartment, OUTFILE)
cat(sprintf("\nSaved: %s\n", OUTFILE))
