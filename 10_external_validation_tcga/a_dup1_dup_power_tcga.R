#!/usr/bin/env Rscript
# A_DUP1: Formal power calculation + TCGA CN-gain breakpoint proximity vs. Δβ test
#
# Purpose:
#   1. Power: show our DUP (n=446 pairs, 10 pt) and INV (n=188, 8 pt) data are
#      underpowered to confirm ρ=−0.088/−0.112 at 80% power.
#   2. TCGA leverage: use 50 TCGA-LIHC paired tumor-normal 450K samples with
#      matched CNV segments (CN gain = Segment_Mean > 0.3, proxy for DUP).
#      Test: does distance to CN-gain breakpoint correlate with Δβ (tumor - normal)?
#      If null in n=50 TCGA pairs → the small ρ in our cohort is consistent with noise.
#
# Run: mamba run -n renv Rscript post_processing/a_dup1_dup_power_tcga.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(pwr)
  library(matrixStats)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

set.seed(42)

METH_RDS <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/tcga_cache/tcga_lihc_meth450k_se.rds")
CNV_RDS  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/tcga_cache/tcga_lihc_cnv_segment.rds")
OUT_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
dir.create(OUT_DIR, showWarnings = FALSE)

GAIN_THRESH  <- 0.3    # log2(CN/2) > 0.3 ≈ CN ≥ 2.46 (gain)
PROX_WIN     <- 5e5L   # 500 kb window around breakpoints
MISS_THRESH  <- 0.2    # max probe missingness

# Section 1: Power calculation =================================================
message("=== Section 1: Power calculation ===")

our_cohort <- data.frame(
  svtype   = c("DEL", "INS", "DUP", "INV"),
  n_pairs  = c(3120, 3120, 446, 188),
  n_pat    = c(12, 12, 10, 8),
  rho_obs  = c(-0.002, 0.003, -0.088, -0.112)
)

power_res <- lapply(seq_len(nrow(our_cohort)), function(i) {
  r <- abs(our_cohort$rho_obs[i])
  n <- our_cohort$n_pairs[i]
  pw <- if (r > 0) pwr.r.test(n = n, r = r, sig.level = 0.05)$power else NA_real_
  # n required for 80% power
  n80 <- if (r > 0) ceiling(pwr.r.test(r = r, sig.level = 0.05, power = 0.80)$n) else NA_integer_
  data.frame(
    svtype    = our_cohort$svtype[i],
    n_pairs   = n,
    n_patients = our_cohort$n_pat[i],
    rho_obs   = our_cohort$rho_obs[i],
    power_obs = round(pw, 3),
    n_for_80  = n80
  )
}) |> bind_rows()

cat("\n=== A_DUP1: Power summary ===\n")
print(power_res)
fwrite(power_res, file.path(OUT_DIR, "a_dup1_power_calc.csv"))
message("Saved: a_dup1_power_calc.csv")

# Section 2: TCGA paired CN-gain breakpoint proximity test =====================
message("\n=== Section 2: TCGA CN-gain breakpoint methylation test ===")

# 2a. Load methylation
message("Loading TCGA-LIHC 450K methylation...")
meth_se <- readRDS(METH_RDS)
cd <- as.data.frame(colData(meth_se))

# Identify paired patients (both tumor and normal)
cd$patient <- substr(cd$barcode, 1, 12)
cd$is_normal <- cd$sample_type == "Solid Tissue Normal"
cd$is_tumor  <- cd$sample_type == "Primary Tumor"

normal_patients <- unique(cd$patient[cd$is_normal])
tumor_patients  <- unique(cd$patient[cd$is_tumor])
paired_patients <- intersect(normal_patients, tumor_patients)
cat(sprintf("Paired patients (tumor + normal): %d\n", length(paired_patients)))

# For each patient: pick one tumor and one normal
pick_one <- function(pat, type_flag) {
  rows <- which(cd$patient == pat & type_flag)
  if (length(rows) == 0) return(NA_integer_)
  rows[1]
}

pair_idx <- lapply(paired_patients, function(p) {
  t_idx <- pick_one(p, cd$is_tumor)
  n_idx <- pick_one(p, cd$is_normal)
  if (anyNA(c(t_idx, n_idx))) return(NULL)
  list(patient = p, tumor_col = t_idx, normal_col = n_idx)
}) |> Filter(f = Negate(is.null))

cat(sprintf("Usable paired samples: %d\n", length(pair_idx)))

# 2b. Precompute Δβ matrix (tumor - normal) for paired patients
tumor_cols  <- sapply(pair_idx, `[[`, "tumor_col")
normal_cols <- sapply(pair_idx, `[[`, "normal_col")

meth_mat <- assay(meth_se)

# Probe missingness filter across all paired columns
used_cols   <- c(tumor_cols, normal_cols)
miss_frac   <- rowMeans(is.na(meth_mat[, used_cols, drop = FALSE]))
keep_probes <- which(miss_frac < MISS_THRESH)
cat(sprintf("Probes after missingness filter: %d / %d\n",
            length(keep_probes), nrow(meth_mat)))

# Δβ per patient
delta_mat <- meth_mat[keep_probes, tumor_cols, drop = FALSE] -
             meth_mat[keep_probes, normal_cols, drop = FALSE]
# Column names = patient IDs
colnames(delta_mat) <- sapply(pair_idx, `[[`, "patient")

# Probe genomic coordinates
probe_gr <- rowRanges(meth_se)[keep_probes]
seqlevelsStyle(probe_gr) <- "UCSC"

cat(sprintf("Δβ matrix: %d probes × %d patients\n",
            nrow(delta_mat), ncol(delta_mat)))

# 2c. Load CNV and define CN-gain breakpoints for paired patients
message("Loading CNV segments...")
cnv <- as.data.table(readRDS(CNV_RDS))
cnv[, patient := substr(Sample, 1, 12)]

# Filter to paired patients + gain segments
gain_cnv <- cnv[patient %in% paired_patients & Segment_Mean > GAIN_THRESH]
cat(sprintf("CN-gain segments in paired patients: %d\n", nrow(gain_cnv)))

# Each segment contributes two breakpoints: Start and End
bp_dt <- rbind(
  gain_cnv[, .(patient, chr = Chromosome, pos = Start)],
  gain_cnv[, .(patient, chr = Chromosome, pos = End)]
)
# Standardise chromosome names
bp_dt[, chr := paste0("chr", gsub("^chr", "", chr))]

# GRanges of breakpoints (1-bp intervals)
bp_gr <- GRanges(
  seqnames = bp_dt$chr,
  ranges   = IRanges(bp_dt$pos, bp_dt$pos),
  patient  = bp_dt$patient
)

cat(sprintf("Total CN-gain breakpoints: %d\n", length(bp_gr)))

# 2d. For each patient: test Spearman ρ(dist_to_nearest_breakpoint, Δβ)
message("Computing per-patient Spearman ρ(dist, Δβ)...")

compute_rho <- function(pat) {
  pat_bp <- bp_gr[bp_gr$patient == pat]
  if (length(pat_bp) == 0) return(NULL)

  # Find probes within PROX_WIN of any breakpoint
  # Extend breakpoints by PROX_WIN
  bp_win <- resize(pat_bp, width = 2 * PROX_WIN + 1, fix = "center") |> trim()
  hits   <- findOverlaps(probe_gr, bp_win)
  if (length(hits) == 0) return(NULL)

  probe_idx <- unique(queryHits(hits))

  # Distance to nearest breakpoint
  dist_to_nearest <- mcols(distanceToNearest(probe_gr[probe_idx], pat_bp))$distance

  # Δβ for this patient
  db <- delta_mat[probe_idx, pat]
  complete <- !is.na(db)
  if (sum(complete) < 30) return(NULL)

  rho_res <- cor.test(dist_to_nearest[complete], db[complete],
                      method = "spearman", exact = FALSE)
  data.frame(
    patient   = pat,
    n_probes  = sum(complete),
    n_breakpoints = length(pat_bp),
    rho       = rho_res$estimate,
    p_value   = rho_res$p.value
  )
}

rho_list <- lapply(paired_patients, compute_rho)
rho_dt   <- bind_rows(Filter(Negate(is.null), rho_list))
cat(sprintf("Patients with ≥30 probes near CN-gain breakpoints: %d / %d\n",
            nrow(rho_dt), length(paired_patients)))

cat("\n=== Per-patient Spearman ρ summary ===\n")
cat(sprintf("Median ρ(dist, Δβ): %.4f\n", median(rho_dt$rho, na.rm = TRUE)))
cat(sprintf("Mean ρ: %.4f ± %.4f (SD)\n",
            mean(rho_dt$rho, na.rm = TRUE),
            sd(rho_dt$rho, na.rm = TRUE)))

# One-sample t-test: is mean ρ different from 0?
t_res <- t.test(rho_dt$rho, mu = 0)
cat(sprintf("One-sample t-test (H0: mean ρ = 0): t=%.3f, df=%d, p=%.4f\n",
            t_res$statistic, t_res$parameter, t_res$p.value))

# Sign test: fraction with ρ < 0 (directionally consistent with cis-induction)
n_neg <- sum(rho_dt$rho < 0, na.rm = TRUE)
n_tot <- sum(!is.na(rho_dt$rho))
sign_binom <- binom.test(n_neg, n_tot, p = 0.5)
cat(sprintf("Sign test (# patients ρ < 0): %d/%d, p=%.4f\n",
            n_neg, n_tot, sign_binom$p.value))

fwrite(rho_dt, file.path(OUT_DIR, "a_dup1_tcga_rho_per_patient.csv"))
message("Saved: a_dup1_tcga_rho_per_patient.csv")

# 2e. Pooled Spearman across all patients (aggregate probe-patient records)
message("Computing pooled Spearman across all patients...")
pooled_list <- lapply(paired_patients, function(pat) {
  pat_bp <- bp_gr[bp_gr$patient == pat]
  if (length(pat_bp) == 0) return(NULL)
  bp_win   <- resize(pat_bp, width = 2 * PROX_WIN + 1, fix = "center") |> trim()
  hits     <- findOverlaps(probe_gr, bp_win)
  if (length(hits) == 0) return(NULL)
  probe_idx <- unique(queryHits(hits))
  dist_nn   <- mcols(distanceToNearest(probe_gr[probe_idx], pat_bp))$distance
  db        <- delta_mat[probe_idx, pat]
  complete  <- !is.na(db)
  if (sum(complete) < 10) return(NULL)
  data.frame(
    dist_bp = dist_nn[complete],
    delta_b = db[complete],
    patient = pat
  )
})
pooled_df <- bind_rows(Filter(Negate(is.null), pooled_list))
cat(sprintf("Pooled probe-patient records: %d\n", nrow(pooled_df)))

pooled_rho <- cor.test(pooled_df$dist_bp, pooled_df$delta_b,
                       method = "spearman", exact = FALSE)
cat(sprintf("\nPooled Spearman ρ(dist, Δβ) = %.4f, p = %.4g\n",
            pooled_rho$estimate, pooled_rho$p.value))

# Distance bin analysis
pooled_df$dist_bin <- cut(
  pooled_df$dist_bp,
  breaks = c(0, 50e3, 100e3, 200e3, 500e3),
  labels = c("0-50kb", "50-100kb", "100-200kb", "200-500kb"),
  include.lowest = TRUE
)

bin_summary <- pooled_df |>
  filter(!is.na(dist_bin)) |>
  group_by(dist_bin) |>
  summarise(
    n       = n(),
    mean_db = mean(delta_b, na.rm = TRUE),
    median_db = median(delta_b, na.rm = TRUE),
    se_db   = sd(delta_b, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

cat("\n=== Distance-bin Δβ summary ===\n")
print(bin_summary)
fwrite(bin_summary, file.path(OUT_DIR, "a_dup1_tcga_bin_summary.csv"))

# Section 3: Figures ===========================================================
message("\nGenerating figures...")
theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

sig_label <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# Panel A: Power curves for DUP and INV
n_seq <- seq(100, 2000, by = 50)
power_curves <- bind_rows(
  data.frame(n = n_seq, svtype = "DUP (ρ=0.088)",
             power = sapply(n_seq, function(n) pwr.r.test(n = n, r = 0.088, sig.level = 0.05)$power)),
  data.frame(n = n_seq, svtype = "INV (ρ=0.112)",
             power = sapply(n_seq, function(n) pwr.r.test(n = n, r = 0.112, sig.level = 0.05)$power))
)

pa <- ggplot(power_curves, aes(x = n, y = power, color = svtype)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey40") +
  geom_vline(
    data = data.frame(svtype = c("DUP (ρ=0.088)", "INV (ρ=0.112)"),
                      n_obs  = c(446L, 188L)),
    aes(xintercept = n_obs, color = svtype),
    linetype = "dotted", linewidth = 0.8
  ) +
  annotate("text", x = 446, y = 0.05, label = "DUP\nn=446", hjust = -0.1, size = 3) +
  annotate("text", x = 188, y = 0.05, label = "INV\nn=188", hjust = -0.1, size = 3) +
  scale_color_manual(values = c("DUP (ρ=0.088)" = "#E24B4A", "INV (ρ=0.112)" = "#5B9BD5")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title    = "A_DUP1 Panel A: Statistical Power for DUP/INV Distance-Decay",
    subtitle = sprintf("DUP power=%.1f%% (n=446); INV power=%.1f%% (n=188); 80%% threshold dashed",
                       power_res$power_obs[power_res$svtype == "DUP"] * 100,
                       power_res$power_obs[power_res$svtype == "INV"] * 100),
    x = "Number of within-block SV–aDMR pairs",
    y = "Power (α = 0.05, two-sided)",
    color = NULL
  ) +
  theme_hcc

# Panel B: Distance-bin Δβ in TCGA (CN gain breakpoints)
pb <- ggplot(bin_summary, aes(x = dist_bin, y = mean_db, fill = dist_bin)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_db - se_db, ymax = mean_db + se_db), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("Pooled ρ=%.4f\np=%s",
                           pooled_rho$estimate,
                           sig_label(pooled_rho$p.value))) +
  labs(
    title    = "A_DUP1 Panel B: TCGA-LIHC CN-Gain Breakpoint Proximity vs. Δβ",
    subtitle = sprintf("50 tumor-normal pairs; CN gain (log2R>0.3) breakpoints; %dk probe-patient records",
                       round(nrow(pooled_df) / 1000)),
    x = "Distance to nearest CN-gain breakpoint",
    y = "Mean Δβ (tumor − normal)",
    fill = NULL
  ) +
  theme_hcc + theme(legend.position = "none")

# Panel C: Per-patient ρ distribution
pc <- ggplot(rho_dt, aes(x = rho)) +
  geom_histogram(bins = 20, fill = "#5B9BD5", alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = mean(rho_dt$rho, na.rm = TRUE),
             color = "#E24B4A", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("Median ρ=%.4f\nt-test p=%s\nSign test p=%s",
                           median(rho_dt$rho, na.rm = TRUE),
                           sig_label(t_res$p.value),
                           sig_label(sign_binom$p.value))) +
  labs(
    title    = "A_DUP1 Panel C: Per-Patient Spearman ρ Distribution (TCGA)",
    subtitle = sprintf("%d patients with CN-gain breakpoints; red line = mean",
                       nrow(rho_dt)),
    x = "Spearman ρ (distance to CN-gain bp vs. Δβ)",
    y = "Count"
  ) +
  theme_hcc

# Save panels
FIG_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/figs/v2")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(FIG_DIR, "fig_a_dup1_power.png"),     pa, width = 8, height = 5, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_a_dup1_tcga_bin.png"),  pb, width = 7, height = 5, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_a_dup1_tcga_rho.png"),  pc, width = 7, height = 5, dpi = 150)
message("Saved: fig_a_dup1_*.png")

# Final summary ================================================================
cat("\n=== A_DUP1 FINAL SUMMARY ===\n")
dup_pow <- power_res[power_res$svtype == "DUP", ]
inv_pow <- power_res[power_res$svtype == "INV", ]
cat(sprintf("DUP: n=%d pairs, ρ=%.3f, power=%.1f%% (need n=%d for 80%%)\n",
            dup_pow$n_pairs, dup_pow$rho_obs, dup_pow$power_obs * 100, dup_pow$n_for_80))
cat(sprintf("INV: n=%d pairs, ρ=%.3f, power=%.1f%% (need n=%d for 80%%)\n",
            inv_pow$n_pairs, inv_pow$rho_obs, inv_pow$power_obs * 100, inv_pow$n_for_80))
cat(sprintf("TCGA pooled ρ=%.4f, p=%.4g (%s)\n",
            pooled_rho$estimate, pooled_rho$p.value,
            sig_label(pooled_rho$p.value)))
cat(sprintf("TCGA sign test: %d/%d ρ<0, p=%.4f (%s)\n",
            n_neg, n_tot, sign_binom$p.value,
            sig_label(sign_binom$p.value)))

message("\nA_DUP1 done.")
