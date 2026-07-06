#!/usr/bin/env Rscript
# concordance_power_icc.R (A_HP2, re-pointed to the current pairing)
#
# Post-hoc power / ICC / design-effect analysis for the cis-concordance test,
# against the CURRENT phaseblock_pairs.csv (n=39,212). Supersedes
# remodeled_constitutional_AMR/post_processing/concordance_power_icc.R, whose
# PAIRS_FILE points at the older, superseded n=8,616 pairing
# (sv_admr_hp_concordance_pairs.csv.gz) and reports a stale OR=1.065/p=0.146
# in its own header comment.
#
# Verified to reproduce (docs/OPEN_ISSUES.md issue #1):
#   ICC (logit scale) = 0.0041, DEFF = 14.38, n_eff = 2727,
#   MDES @ 80% power = 52.7% concordance
#
# Run: mamba run -n renv Rscript concordance_power_icc.R

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(pwr)
})

PAIRS_FILE <- "/node200data/kachungk/hcc_data/SV_aDMR/phaseblock_pairs.csv"
OUTFILE    <- "/node200data/kachungk/hcc_data/SV_aDMR/result/concordance_power_icc.csv"

pairs <- fread(PAIRS_FILE)
pairs[, concordant := as.integer(sv_minus_wt < 0)]
setnames(pairs, "patient_code", "patient_id")

N     <- nrow(pairs)
n_cl  <- length(unique(pairs$patient_id))
m_bar <- N / n_cl
p_obs <- mean(pairs$concordant)
cat(sprintf("N=%d pairs, %d patients, m_bar=%.1f, p_obs=%.4f\n", N, n_cl, m_bar, p_obs))

mm   <- glmer(concordant ~ 1 + (1 | patient_id), data = pairs, family = binomial,
              control = glmerControl(optimizer = "bobyqa"))
tau2 <- VarCorr(mm)$patient_id[1, 1]
ICC  <- tau2 / (tau2 + pi^2 / 3)
DEFF <- 1 + (m_bar - 1) * ICC
n_eff <- N / DEFF
cat(sprintf("ICC=%.4f, DEFF=%.2f, n_eff=%.1f\n", ICC, DEFF, n_eff))

mdes   <- pwr.p.test(n = n_eff, sig.level = 0.05, power = 0.80, alternative = "two.sided")
p_mdes <- sin(mdes$h / 2 + asin(sqrt(0.5)))^2
cat(sprintf("MDES @ 80%% power: %.2f%% concordance (Cohen's h=%.4f)\n", 100 * p_mdes, mdes$h))

true_conc <- c(0.52, 0.525, 0.53, 0.55, 0.60)
cf_power  <- sapply(true_conc, function(p) {
  pwr.p.test(n = n_eff, h = abs(ES.h(p, 0.5)), sig.level = 0.05,
             alternative = "two.sided")$power
})
for (i in seq_along(true_conc))
  cat(sprintf("  power at true concordance=%.1f%%: %.1f%%\n", 100 * true_conc[i], 100 * cf_power[i]))

results <- data.table(
  metric = c("N_pairs", "n_patients", "m_bar", "p_observed", "ICC_logit_scale",
             "DEFF", "n_eff", "MDES_80pct_concordance",
             paste0("power_at_", round(100 * true_conc), "pct")),
  value  = c(N, n_cl, m_bar, p_obs, ICC, DEFF, n_eff, p_mdes, cf_power)
)
dir.create(dirname(OUTFILE), showWarnings = FALSE, recursive = TRUE)
fwrite(results, OUTFILE)
cat(sprintf("\nSaved: %s\n", OUTFILE))
