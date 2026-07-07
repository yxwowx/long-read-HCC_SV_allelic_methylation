#!/usr/bin/env Rscript
# hbv_allele_anchored.R — Level B: allele-anchored normal deviation (C33)
#
# For aDMRs in phase blocks containing HBV-proximal SVs, computes:
#   dev_HBV = HBV-haplotype beta  − normal.Methy
#   dev_WT  = WT-haplotype beta   − normal.Methy
#
# Primary test: is |dev_HBV| > |dev_WT|? (HBV-hap deviates more from normal)
# Direction:    is dev_HBV < 0? (HBV-hap hypomethylated vs normal)
# Specificity:  same asymmetry at non-HBV SV blocks in same patients (internal ctrl)
#
# Reuses infer_sv_hp_map() and get_sv_hp_beta() copied from pipeline/04 and 12.
# Adds normal_beta from confident_dmr_per_patient.csv.gz to the HP-oriented frame.
#
# Usage: mamba run -n renv Rscript post_processing/hbv_allele_anchored.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

# Paths ========================================================================
HBV_LOCI_CSV  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/12.HBV_analysis/hbv_v1_somatic_hbv_loci.csv")
SV_FILE       <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/02.sv_dmr_enrichment/sv_tad_ctcf_annotation.v2.csv.gz")
CONF_DMR      <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz")
PHASE_VCF_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "hg38+HBV/clairS/phased_vcf")
MAPPING_CSV   <- Sys.getenv("PATIENT_CODE_MAP")
OUT_DIR       <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR       <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/figs/v2")

HBV_MATCH_BP  <- 500L   # positional tolerance for SV → HBV integration match

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Helpers (copied from pipeline/04 and 12) =====================================
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

# Extended get_sv_hp_beta: adds normal_beta column
get_sv_hp_beta_ext <- function(admr_gr, sv_hp_map) {
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
    normal_beta  = admr_sub$normal_beta,
    hp_delta_abs = abs(admr_sub$hp1_beta - admr_sub$hp2_beta),
    n_cpg        = admr_sub$n_cpg
  )
}

# 1. Load HBV somatic loci =====================================================
message("Loading somatic HBV loci...")
hbv_loci <- fread(HBV_LOCI_CSV)
hbv_loci <- hbv_loci[as.logical(is_tumor) == TRUE]
cat(sprintf("Somatic HBV loci: %d rows, %d patients\n",
            nrow(hbv_loci), length(unique(hbv_loci$pcode))))

# 2. Load SV table, flag HBV-proximal SVs ======================================
message("Loading SV annotation...")
sv <- fread(SV_FILE)

# Reconstruct is_hbv_bnd (replicates 12_hbv_analysis.R lines 221-224)
sv[, is_hbv_bnd := mapply(function(pt, chr, svpos) {
  sub <- hbv_loci[pcode == pt & chrom == chr, pos]
  length(sub) > 0L && any(abs(svpos - sub) <= HBV_MATCH_BP, na.rm = TRUE)
}, sample, seqnames, start)]

n_hbv_sv <- sum(sv$is_hbv_bnd, na.rm = TRUE)
cat(sprintf("HBV-proximal SVs: %d (±%d bp match)\n", n_hbv_sv, HBV_MATCH_BP))
if (n_hbv_sv == 0L) stop("No SVs matched HBV loci.")

hbv_patients <- unique(sv[is_hbv_bnd == TRUE, sample])
cat("HBV patients:", paste(sort(hbv_patients), collapse=", "), "\n")

# 3. Load patient mapping ======================================================
patient_map <- fread(MAPPING_CSV)
setnames(patient_map, "Samples_ID", "Samples_ID")  # keep as-is

# 4. Load phase blocks (GTF -> GRanges list, per patient) ======================
message("Loading phase block GTF files...")
gtf_files <- list.files(PHASE_VCF_DIR, pattern = "\\.gtf$", full.names = TRUE)
if (length(gtf_files) == 0L) stop("No GTF files found in: ", PHASE_VCF_DIR)

phase_blocks <- lapply(gtf_files, function(x) {
  gr <- tryCatch(import(x), error = function(e) {
    warning(basename(x), ": ", e$message); NULL
  })
  if (is.null(gr)) return(NULL)
  mcols(gr)$patient_name <- sub("\\.gtf$", "", basename(x))
  gr
})
phase_blocks <- Filter(Negate(is.null), phase_blocks)

# Flatten to data.frame, join patient code on Samples_ID (same convention as 12_hbv_analysis.R)
pb_df <- lapply(phase_blocks, function(gr) {
  df_gr <- as.data.frame(gr)
  data.frame(
    seqnames     = as.character(df_gr$seqnames),
    start        = df_gr$start,
    end          = df_gr$end,
    block_id     = df_gr$gene_id,
    patient_name = df_gr$patient_name,
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

pb_df <- left_join(
  pb_df,
  patient_map[, c("Samples_ID", "patient_code")],
  by = c("patient_name" = "Samples_ID")
) |> filter(!is.na(patient_code))

phase_gr_list <- pb_df |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
  (\(gr) split(gr, gr$patient_code))()

# 5. Load aDMRs, assign block_id, add normal_beta ==============================
message("Loading confident aDMR table...")
admr_raw <- fread(CONF_DMR)

# Detect column names (handle raw vs renamed variants)
seq_col   <- intersect(c("admr_chr",   "seqnames"), names(admr_raw))[1]
start_col <- intersect(c("admr_start", "start"),    names(admr_raw))[1]
end_col   <- intersect(c("admr_end",   "end"),       names(admr_raw))[1]
hp1_col   <- intersect(c("HP1.Methy",  "hp1_beta"), names(admr_raw))[1]
hp2_col   <- intersect(c("HP2.Methy",  "hp2_beta"), names(admr_raw))[1]
pt_col    <- intersect(c("patient_code","sample"),   names(admr_raw))[1]
norm_col  <- intersect(c("normal.Methy","normal_meth","normal_beta"), names(admr_raw))[1]
ncg_col   <- intersect(c("nCG", "n_cpg"),            names(admr_raw))[1]

admr_dt <- as.data.table(admr_raw)
setnames(admr_dt, seq_col,   "seqnames_gr")
setnames(admr_dt, start_col, "start_gr")
setnames(admr_dt, end_col,   "end_gr")
setnames(admr_dt, hp1_col,   "hp1_beta")
setnames(admr_dt, hp2_col,   "hp2_beta")
setnames(admr_dt, pt_col,    "pcode")
setnames(admr_dt, norm_col,  "normal_beta")
setnames(admr_dt, ncg_col,   "n_cpg")
admr_dt[, hp_delta := hp1_beta - hp2_beta]

# Protect reserved GRanges column names
for (.rc in c("seqnames","start","end","width","strand")) {
  if (.rc %in% names(admr_dt)) setnames(admr_dt, .rc, paste0("meta_", .rc))
}

admr_dt <- admr_dt[grepl("^chr[0-9XY]+$", seqnames_gr)]
admr_gr_all <- makeGRangesFromDataFrame(as.data.frame(admr_dt),
                                        seqnames.field = "seqnames_gr",
                                        start.field    = "start_gr",
                                        end.field      = "end_gr",
                                        keep.extra.columns = TRUE)
admr_gr_list <- split(admr_gr_all, admr_gr_all$pcode)

# Assign block_id via overlap with phase blocks
admr_phased_list <- lapply(hbv_patients, function(pt) {
  dmr_gr <- admr_gr_list[[pt]]
  blk_gr <- phase_gr_list[[pt]]
  if (is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)
  if (is.null(blk_gr) || length(blk_gr) == 0) {
    mcols(dmr_gr)$block_id <- NA_character_
    return(dmr_gr)
  }
  hits <- findOverlaps(dmr_gr, blk_gr, select = "first")
  mcols(dmr_gr)$block_id <- blk_gr$block_id[hits]
  dmr_gr[!is.na(mcols(dmr_gr)$block_id)]
}) |> setNames(hbv_patients)

# 6. HBV-proximal phase blocks: extract HP-oriented betas + normal =============
message("Extracting HP-oriented betas for HBV-proximal SV blocks...")

hbv_pairs_list <- lapply(hbv_patients, function(pt) {
  sv_hbv <- sv[sample == pt & is_hbv_bnd == TRUE]
  dmr_gr  <- admr_phased_list[[pt]]
  if (nrow(sv_hbv) == 0 || is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  sv_hp_map <- infer_sv_hp_map(sv_hbv$HP, sv_hbv$PHASESETID)
  if (length(sv_hp_map) == 0) return(NULL)

  hp_df <- get_sv_hp_beta_ext(dmr_gr, sv_hp_map)
  if (is.null(hp_df) || nrow(hp_df) == 0) return(NULL)

  hp_df$patient_id <- pt
  hp_df$group      <- "HBV-proximal"
  hp_df
}) |> Filter(Negate(is.null), x = _)

if (length(hbv_pairs_list) == 0) stop("No HBV-proximal aDMR pairs found.")
hbv_df <- bind_rows(hbv_pairs_list)

cat(sprintf("\nHBV-proximal aDMR pairs: %d (across %d patients)\n",
            nrow(hbv_df), length(unique(hbv_df$patient_id))))
print(hbv_df |> count(patient_id, name = "n_pairs"))

# 7. Specificity control: non-HBV SV blocks, same patients =====================
message("Extracting non-HBV SV blocks (specificity control)...")

ctrl_pairs_list <- lapply(hbv_patients, function(pt) {
  # Non-HBV phased SVs in same patient
  sv_ctrl <- sv[sample == pt & is_hbv_bnd == FALSE &
                !is.na(HP) & HP %in% c(1L, 2L) & !is.na(PHASESETID)]
  dmr_gr   <- admr_phased_list[[pt]]
  if (nrow(sv_ctrl) == 0 || is.null(dmr_gr) || length(dmr_gr) == 0) return(NULL)

  sv_hp_map <- infer_sv_hp_map(sv_ctrl$HP, sv_ctrl$PHASESETID)
  if (length(sv_hp_map) == 0) return(NULL)

  hp_df <- get_sv_hp_beta_ext(dmr_gr, sv_hp_map)
  if (is.null(hp_df) || nrow(hp_df) == 0) return(NULL)

  hp_df$patient_id <- pt
  hp_df$group      <- "Non-HBV SV"
  hp_df
}) |> Filter(Negate(is.null), x = _)

ctrl_df <- if (length(ctrl_pairs_list) > 0) bind_rows(ctrl_pairs_list) else data.frame()
cat(sprintf("Non-HBV control pairs: %d\n", nrow(ctrl_df)))

# 8. Compute deviation from normal =============================================
for (df in list(hbv_df, ctrl_df)) {
  if (nrow(df) == 0) next
  df$dev_hbv <- df$sv_hp_beta  - df$normal_beta
  df$dev_wt  <- df$wt_hp_beta  - df$normal_beta
  df$abs_dev_hbv <- abs(df$dev_hbv)
  df$abs_dev_wt  <- abs(df$dev_wt)
}
hbv_df$dev_hbv      <- hbv_df$sv_hp_beta  - hbv_df$normal_beta
hbv_df$dev_wt       <- hbv_df$wt_hp_beta  - hbv_df$normal_beta
hbv_df$abs_dev_hbv  <- abs(hbv_df$dev_hbv)
hbv_df$abs_dev_wt   <- abs(hbv_df$dev_wt)
if (nrow(ctrl_df) > 0) {
  ctrl_df$dev_hbv     <- ctrl_df$sv_hp_beta - ctrl_df$normal_beta
  ctrl_df$dev_wt      <- ctrl_df$wt_hp_beta - ctrl_df$normal_beta
  ctrl_df$abs_dev_hbv <- abs(ctrl_df$dev_hbv)
  ctrl_df$abs_dev_wt  <- abs(ctrl_df$dev_wt)
}

# 9. Statistical tests =========================================================
cat("\n=== Test 1: |dev_HBV| vs |dev_WT| — HBV-proximal pairs ===\n")
wt1 <- tryCatch(
  wilcox.test(hbv_df$abs_dev_hbv, hbv_df$abs_dev_wt,
              paired = TRUE, alternative = "greater", exact = FALSE),
  error = function(e) list(p.value = NA_real_, statistic = NA_real_)
)
cat(sprintf("Wilcoxon signed-rank (|dev_HBV| > |dev_WT|): W=%.1f, p=%.4f\n",
            as.numeric(wt1$statistic), wt1$p.value))
cat(sprintf("Median |dev_HBV|=%.3f  Median |dev_WT|=%.3f\n",
            median(hbv_df$abs_dev_hbv, na.rm=TRUE),
            median(hbv_df$abs_dev_wt,  na.rm=TRUE)))

cat("\n=== Test 2: Directionality — is HBV-hap hypomethylated vs normal? ===\n")
n_hypo  <- sum(hbv_df$dev_hbv < 0, na.rm = TRUE)
n_total <- sum(!is.na(hbv_df$dev_hbv))
binom_p <- binom.test(n_hypo, n_total, p = 0.5, alternative = "greater")$p.value
cat(sprintf("dev_HBV < 0: %d / %d (%.1f%%), binomial p=%.4f (H0=50%% hypo)\n",
            n_hypo, n_total, 100 * n_hypo / n_total, binom_p))

cat("\n=== Test 3: Per-patient sign test ===\n")
per_pat <- hbv_df |>
  group_by(patient_id) |>
  summarise(
    n_pairs      = n(),
    n_hypo       = sum(dev_hbv < 0, na.rm = TRUE),
    pct_hypo     = round(100 * mean(dev_hbv < 0, na.rm = TRUE), 1),
    mean_dev_hbv = round(mean(dev_hbv, na.rm = TRUE), 3),
    mean_dev_wt  = round(mean(dev_wt,  na.rm = TRUE), 3),
    .groups = "drop"
  )
print(per_pat)
n_pat_hypo <- sum(per_pat$pct_hypo > 50, na.rm = TRUE)
n_pat_total <- nrow(per_pat)
sign_p <- binom.test(n_pat_hypo, n_pat_total, p = 0.5)$p.value
cat(sprintf("Per-patient sign test: %d/%d patients majority-hypo, p=%.4f\n",
            n_pat_hypo, n_pat_total, sign_p))

# Specificity control: is asymmetry larger at HBV vs non-HBV blocks?
test_df <- list()
if (nrow(ctrl_df) > 0) {
  cat("\n=== Test 4: Specificity — |dev_HBV| − |dev_WT| at HBV vs non-HBV blocks ===\n")
  hbv_df$asym  <- hbv_df$abs_dev_hbv  - hbv_df$abs_dev_wt
  ctrl_df$asym <- ctrl_df$abs_dev_hbv - ctrl_df$abs_dev_wt
  wt4 <- tryCatch(
    wilcox.test(hbv_df$asym, ctrl_df$asym, alternative = "greater", exact = FALSE),
    error = function(e) list(p.value = NA_real_, statistic = NA_real_)
  )
  cat(sprintf("MWU (HBV asym > non-HBV asym): W=%.1f, p=%.4f\n",
              as.numeric(wt4$statistic), wt4$p.value))
  cat(sprintf("Median asym: HBV=%.3f  non-HBV=%.3f\n",
              median(hbv_df$asym, na.rm=TRUE),
              median(ctrl_df$asym, na.rm=TRUE)))
  test_df[["specificity"]] <- data.frame(
    test = "MWU asym HBV>nonHBV",
    stat = as.numeric(wt4$statistic), p = wt4$p.value,
    median_hbv  = median(hbv_df$asym,  na.rm=TRUE),
    median_ctrl = median(ctrl_df$asym, na.rm=TRUE)
  )
}

# 10. Compile + save ===========================================================
all_pairs <- bind_rows(hbv_df, ctrl_df)
fwrite(all_pairs,
       file.path(OUT_DIR, "hbv_allele_anchored_pairs.csv"))

tests_df <- data.frame(
  test            = c("Wilcoxon |dev_HBV|>|dev_WT| (paired)",
                      "Binomial dev_HBV<0 (hypo)"),
  stat            = c(as.numeric(wt1$statistic), NA_real_),
  p               = c(wt1$p.value, binom_p),
  n               = c(nrow(hbv_df), n_total),
  median_hbv_arm  = c(median(hbv_df$abs_dev_hbv, na.rm=TRUE),
                       n_hypo / n_total),
  median_wt_arm   = c(median(hbv_df$abs_dev_wt, na.rm=TRUE), 0.5)
)
if (length(test_df) > 0) tests_df <- bind_rows(tests_df, test_df[["specificity"]])
fwrite(tests_df, file.path(OUT_DIR, "hbv_allele_anchored_tests.csv"))
message("Wrote: hbv_allele_anchored_pairs.csv, hbv_allele_anchored_tests.csv")

# 11. Figure ===================================================================
message("Generating figure...")

plot_df <- all_pairs |>
  select(patient_id, group, dev_hbv, dev_wt, abs_dev_hbv, abs_dev_wt) |>
  tidyr::pivot_longer(
    cols = c(dev_hbv, dev_wt),
    names_to = "arm", values_to = "dev"
  ) |>
  mutate(arm = recode(arm, dev_hbv = "HBV-hap", dev_wt = "WT-hap"),
         arm = factor(arm, levels = c("HBV-hap", "WT-hap")))

pA <- ggplot(plot_df |> filter(group == "HBV-proximal"),
             aes(x = arm, y = dev, fill = arm)) +
  geom_violin(alpha = 0.65, trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("HBV-hap" = "#d73027", "WT-hap" = "#4292c6")) +
  facet_wrap(~patient_id, nrow = 1, scales = "free_y") +
  labs(title = "Deviation from normal methylation",
       subtitle = "HBV-proximal SV phase blocks",
       x = NULL, y = "β_hap − β_normal",
       caption = sprintf("Wilcoxon paired |dev_HBV|>|dev_WT|: p=%.4f  Binomial hypo: p=%.4f",
                         wt1$p.value, binom_p)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9),
        plot.title = element_text(face = "bold"))

pB <- ggplot(all_pairs |> filter(!is.na(abs_dev_hbv)) |>
               mutate(asym = abs_dev_hbv - abs_dev_wt),
             aes(x = group, y = asym, fill = group)) +
  geom_violin(alpha = 0.65, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("HBV-proximal" = "#d73027", "Non-HBV SV" = "#74c476")) +
  labs(title = "Specificity control",
       subtitle = "|dev_HBV| − |dev_WT| at HBV vs non-HBV blocks",
       x = NULL, y = "|dev_HBV| − |dev_WT|") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

combined <- pA / pB + plot_layout(heights = c(1.4, 1))
ggsave(file.path(FIG_DIR, "fig_hbv_allele_anchored.png"), combined,
       width = 11, height = 8, dpi = 150)
message("Saved: fig_hbv_allele_anchored.png")

cat("\nDone (C33 Level B).\n")

# A-III: Covariate-adjusted LMM (primary) + matched resample (sensitivity) =====
# Resolves C33 Test 4 marginal MWU p≈0.055 (unmatched controls).
cat("\n====== A-III: HBV specificity — covariate-adjusted LMM + matched resample ======\n")

suppressPackageStartupMessages({
  library(lme4)
})

# Ensure asym is present for all rows
if (!"asym" %in% names(all_pairs)) {
  all_pairs$asym <- all_pairs$abs_dev_hbv - all_pairs$abs_dev_wt
}
all_pairs$is_hbv <- all_pairs$group == "HBV-proximal"

# 1. Join SV covariates (cnv_class, svtype, VAF) via block_id == PHASESETID ====
# Use the sv data.table already loaded; take one representative SV per block
# (HBV blocks typically have a single proximal SV; non-HBV blocks take first).
sv_cov <- unique(sv[!is.na(PHASESETID) & !is.na(HP),
                    .(block_id   = as.character(PHASESETID),
                      patient_id = sample,
                      cnv_class,
                      svtype,
                      VAF)])
sv_cov_agg <- sv_cov[, .(cnv_class = cnv_class[1],
                          svtype    = svtype[1],
                          VAF       = mean(VAF, na.rm = TRUE)),
                     by = .(block_id, patient_id)]

lmm_input <- merge(
  as.data.table(all_pairs)[!is.na(asym),
                            .(block_id = as.character(block_id),
                              patient_id,
                              asym,
                              is_hbv)],
  sv_cov_agg,
  by  = c("block_id", "patient_id"),
  all.x = TRUE
)
lmm_clean <- lmm_input[!is.na(cnv_class) & !is.na(svtype) & !is.na(VAF)]

cat(sprintf("LMM data: %d pairs (%d HBV, %d non-HBV) after covariate join\n",
            nrow(lmm_clean),
            sum(lmm_clean$is_hbv, na.rm = TRUE),
            sum(!lmm_clean$is_hbv, na.rm = TRUE)))

lmm_result <- NULL
lmm_is_hbv_p <- NA_real_

if (nrow(lmm_clean) >= 10 && length(unique(lmm_clean$patient_id)) >= 3) {
  lmm_fit <- tryCatch(
    lme4::lmer(
      asym ~ is_hbv + cnv_class + svtype + VAF + (1 | patient_id),
      data = as.data.frame(lmm_clean),
      REML = TRUE
    ),
    error   = function(e) { message("LMM error: ", e$message); NULL },
    warning = function(w) {
      suppressWarnings(
        lme4::lmer(asym ~ is_hbv + cnv_class + svtype + VAF + (1 | patient_id),
                   data = as.data.frame(lmm_clean), REML = TRUE)
      )
    }
  )

  if (!is.null(lmm_fit)) {
    sc        <- as.data.frame(summary(lmm_fit)$coefficients)
    sc$term   <- rownames(sc)
    # lme4 (no lmerTest): columns are Estimate, Std. Error, t value
    # Use asymptotic z-approximation for p (valid for large n)
    names(sc)[names(sc) == "Estimate"]   <- "beta"
    names(sc)[names(sc) == "Std. Error"] <- "se"
    names(sc)[names(sc) == "t value"]    <- "t"
    sc$p <- 2 * pnorm(-abs(sc$t))  # asymptotic two-sided p
    lmm_result <- as.data.table(sc)

    hbv_row <- lmm_result[grepl("is_hbvTRUE", term)]
    if (nrow(hbv_row) > 0) {
      lmm_is_hbv_p <- hbv_row$p[1]
      cat(sprintf("LMM is_hbvTRUE: β=%.4f, SE=%.4f, t=%.3f, p=%.4f (asymptotic)\n",
                  hbv_row$beta[1], hbv_row$se[1], hbv_row$t[1], hbv_row$p[1]))
    }
    fwrite(lmm_result, file.path(OUT_DIR, "hbv_specificity_lmm.csv"))
    message("Wrote: hbv_specificity_lmm.csv")
  }
} else {
  cat(sprintf("Insufficient data for LMM (n=%d, n_pat=%d) — skipping.\n",
              nrow(lmm_clean), length(unique(lmm_clean$patient_id))))
}

# 2. Matched resample sensitivity ==============================================
# Restrict control pool to Non-boundary SVs (matching all HBV-proximal SVs).
# Match on cnv_class × svtype strata; draw 1000 times; compute empirical p.
set.seed(42L)
N_DRAWS <- 1000L

hbv_asym <- all_pairs[all_pairs$group == "HBV-proximal" & !is.na(all_pairs$asym),
                       "asym", drop = TRUE]

# Join stratification to all_pairs
ctrl_pool_sv_raw <- sv[is_hbv_bnd == FALSE & !is.na(HP) & HP %in% c(1L, 2L) &
                       !is.na(PHASESETID) & stratification == "Non-boundary",
                       .(block_id   = as.character(PHASESETID),
                         patient_id = sample,
                         cnv_class,
                         svtype)]
# One row per (block_id, patient_id) to avoid cartesian join
ctrl_pool_sv <- ctrl_pool_sv_raw[, .(cnv_class = cnv_class[1], svtype = svtype[1]),
                                  by = .(block_id, patient_id)]
ctrl_pool_pairs <- merge(
  as.data.table(all_pairs)[group == "Non-HBV SV" & !is.na(asym),
                            .(block_id = as.character(block_id), patient_id, asym)],
  ctrl_pool_sv,
  by = c("block_id", "patient_id"),
  all.x = FALSE
)
cat(sprintf("Non-boundary control pool for resample: %d pairs\n", nrow(ctrl_pool_pairs)))

# Strata from HBV SVs
hbv_sv_strata <- sv[is_hbv_bnd == TRUE, paste(cnv_class, svtype, sep = "|")]
strata_tbl    <- table(hbv_sv_strata)

# Observed W statistic (full non-HBV pool, no matching)
obs_W <- if (length(hbv_asym) >= 3 && nrow(ctrl_pool_pairs) >= 3) {
  wilcox.test(hbv_asym, ctrl_pool_pairs$asym,
              alternative = "greater", exact = FALSE)$statistic
} else {
  NA_real_
}

# Draw matched control sets
draw_matched_asym <- function(pool_dt, strata_tbl) {
  samp <- lapply(names(strata_tbl), function(s) {
    n   <- as.integer(strata_tbl[s])
    sub <- pool_dt[paste(cnv_class, svtype, sep = "|") == s, asym]
    if (length(sub) == 0L) return(rep(NA_real_, n))
    sample(sub, size = n, replace = TRUE)
  })
  unlist(samp)
}

null_W <- vapply(seq_len(N_DRAWS), function(i) {
  ctrl_asym <- draw_matched_asym(ctrl_pool_pairs, strata_tbl)
  ctrl_asym <- ctrl_asym[!is.na(ctrl_asym)]
  if (length(ctrl_asym) < 3 || length(hbv_asym) < 3) return(NA_real_)
  tryCatch(
    wilcox.test(hbv_asym, ctrl_asym, alternative = "greater", exact = FALSE)$statistic,
    error = function(e) NA_real_
  )
}, numeric(1))

emp_p <- mean(null_W >= obs_W, na.rm = TRUE)
cat(sprintf("Matched resample (%d draws): observed W=%.1f, empirical p=%.4f\n",
            sum(!is.na(null_W)), obs_W, emp_p))

resample_df <- data.frame(
  analysis         = "matched_nonboundary_resample",
  n_hbv_pairs      = length(hbv_asym),
  n_ctrl_pool      = nrow(ctrl_pool_pairs),
  n_draws          = sum(!is.na(null_W)),
  observed_W       = obs_W,
  empirical_p      = emp_p,
  null_W_median    = median(null_W, na.rm = TRUE),
  null_W_p95       = as.numeric(quantile(null_W, 0.95, na.rm = TRUE))
)
fwrite(resample_df, file.path(OUT_DIR, "hbv_specificity_matched_resample.csv"))
message("Wrote: hbv_specificity_matched_resample.csv")

cat("\nDone (A-III).\n")
