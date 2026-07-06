#!/usr/bin/env Rscript
# P0-A: Power calibration for non-significant patients (P3, P11, P12)
# Output: result/power_calibration_ns_patients.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

OUT_DIR  <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
PURITY_DIR <- "/node200data/kachungk/hcc_data/cnv_deepsomatic.out_hg38/purple"
SV_ANN   <- "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz"
COOCCUR  <- "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/tier_50kb_v2_cooccur_pct_results.csv"

patient_map <- data.frame(
  sample_id = c("JJT_HCC","KIS_HCC","KSJ_HCC","LHS_HCC","LSS_HCC",
                "PJS_HCC","PSY_HCC","WSY_HCC","YJS_HCC","NSH_HCC","MSB_HCC","YMS_HCC"),
  patient_id = paste0("P", c(1,2,3,4,5,6,7,8,9,10,11,12)),
  stringsAsFactors = FALSE
)

# --- 1. Read purity from PURPLE ---
purity_files <- list.files(PURITY_DIR, pattern = "\\.purple\\.purity\\.tsv$", full.names = TRUE)
purity_df <- lapply(purity_files, function(f) {
  sample_id <- sub("_tumor\\.purple\\.purity\\.tsv$", "", basename(f))
  d <- fread(f, nrows = 1)
  data.frame(sample_id = sample_id, tumor_purity = d$purity, ploidy = d$ploidy)
}) |> bind_rows() |>
  left_join(patient_map, by = "sample_id") |>
  select(patient_id, tumor_purity, ploidy)

# --- 2. Count HBV BND SVs per patient ---
message("Reading SV annotation for HBV BND counts...")
sv_ann <- fread(cmd = paste("zcat", SV_ANN), data.table = FALSE)
hbv_counts <- sv_ann |>
  filter(is_hbv == TRUE) |>
  count(sample, name = "n_hbv_bnd") |>
  rename(patient_id = sample)

# --- 3. Load enrichment results ---
cooccur <- fread(COOCCUR, data.table = FALSE)
cooccur <- cooccur |>
  mutate(
    fold_change = obs_pct / null_mean,
    sig_label   = factor(sig_label, levels = c("ns", "*", "**", "***"))
  )

# --- 4. Leave-one-out cohort enrichment fold-change ---
# For each patient: compare mean(obs_pct) excl. that patient vs full cohort mean
full_mean_fc <- mean(cooccur$fold_change)

loo_fc <- sapply(cooccur$patient_id, function(pid) {
  mean(cooccur$fold_change[cooccur$patient_id != pid])
})
loo_df <- data.frame(patient_id = cooccur$patient_id, loo_cohort_fc = loo_fc)

# --- 5. Assemble power calibration table ---
power_tbl <- cooccur |>
  select(patient_id, n_sv, n_dmr, obs_pct, null_mean, wilcox_p, wilcox_fdr, sig_label, fold_change) |>
  left_join(purity_df, by = "patient_id") |>
  left_join(hbv_counts, by = "patient_id") |>
  left_join(loo_df, by = "patient_id") |>
  mutate(
    n_hbv_bnd       = replace(n_hbv_bnd, is.na(n_hbv_bnd), 0),
    is_ns            = sig_label == "ns",
    power_note       = case_when(
      !is_ns                          ~ "Significant — no power concern",
      n_sv < 100                      ~ "Low SV count (<100) — likely insufficient power",
      tumor_purity < 0.4              ~ "Low purity (<0.4) — subclonal SVs may be missed",
      fold_change >= 0.95 & fold_change <= 1.05 ~ "True negative (fold-change ≈ 1.0)",
      TRUE                            ~ "Marginal — assess LOO impact"
    ),
    cohort_full_fc   = full_mean_fc
  ) |>
  arrange(desc(is_ns), patient_id)

fwrite(power_tbl, file.path(OUT_DIR, "power_calibration_ns_patients.csv"))
message("Wrote: ", file.path(OUT_DIR, "power_calibration_ns_patients.csv"))

# --- 6. Print summary ---
cat("\n=== Power Calibration Summary ===\n")
ns_pts <- power_tbl |> filter(is_ns)
cat(sprintf("Non-significant patients: %s\n\n", paste(ns_pts$patient_id, collapse = ", ")))
print(ns_pts |> select(patient_id, n_sv, n_dmr, tumor_purity, n_hbv_bnd,
                        obs_pct, null_mean, fold_change, loo_cohort_fc, power_note))

cat(sprintf("\nFull cohort mean fold-change (obs/null): %.3f\n", full_mean_fc))
cat("\nLOO cohort fold-change when each ns patient excluded:\n")
print(ns_pts |> select(patient_id, fold_change, loo_cohort_fc))
