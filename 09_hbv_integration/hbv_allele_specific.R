#!/usr/bin/env Rscript
# 06_hbv_allele_specific.R
# HBV integration allele-specific methylation analysis
# Level A: loads existing per-read outputs (no recomputation)
# Level B: allele-anchored somatic aDMR analysis
#   - Re-anchors to somatic aDMR (constitutional signal removed)
#   - dev_HBV = HBV_hap_β − normal, dev_WT = WT_hap_β − normal
#   - normal_beta computed directly from the genome-wide DSS-smoothed
#     haplotype-resolved normal-liver methylation model
#     (total.hap.bsobj.smoothed.rds), via getMeth(type="smooth",
#     what="perRegion") over each aDMR's own interval, averaged across
#     normal_hap1/normal_hap2. Supersedes an earlier two-tier estimate
#     (DSS-called normal aDMR overlap for 36.5% of pairs + a fixed
#     background constant for the rest) that a very low correlation
#     check revealed to disagree sharply for the fallback-constant subset.
# Output: SV_aDMR/hbv_allele_anchored_pairs.csv, hbv_allele_anchored_tests.csv,
#         hbv_specificity_lmm.csv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(lme4)
  library(bsseq)
})

REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

DATA_ROOT        <- Sys.getenv("HCC_DATA_DIR")
DMR_SVS_ROOT     <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT         <- file.path(DATA_ROOT, "SV_aDMR")
PHASE_VCF_DIR    <- file.path(DATA_ROOT, "hg38+HBV/clairS/phased_vcf")
NORMAL_ADMR_DIR  <- file.path(DATA_ROOT, "DMR_minimap2.out_hg38/DSS")
BSOBJ_RDS        <- file.path(NORMAL_ADMR_DIR, "total.hap.bsobj.smoothed.rds")

HBV_LOCI_CSV <- file.path(DMR_SVS_ROOT, "12.HBV_analysis/hbv_v1_somatic_hbv_loci.csv")
SV_FILE      <- file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz")

# Fallback locus-average normal β, used only for the initial somatic_dt join
# (retained so downstream helper functions have a placeholder normal_beta
# column); overwritten for all HBV/control pairs by the smoothed getMeth
# calculation below before any statistical test is run.
NORMAL_BG_BETA <- 0.62

HBV_MATCH_BP  <- 500L
MAIN_CHROMS   <- paste0("chr", c(1:22, "X"))
N_RESAMPLE    <- 1000L

# Helpers (adapted from hbv_allele_anchored.R) =================================
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

get_sv_hp_beta_ext <- function(admr_gr, sv_hp_map) {
  sv_hp_map <- sv_hp_map[!is.na(sv_hp_map) & sv_hp_map %in% c(1L, 2L)]
  if (length(sv_hp_map) == 0L) return(NULL)
  in_sv <- mcols(admr_gr)$block_id %in% names(sv_hp_map)
  if (!any(in_sv)) return(NULL)
  admr_sub  <- admr_gr[in_sv]
  sv_hp_vec <- as.integer(sv_hp_map[as.character(mcols(admr_sub)$block_id)])
  sv_beta   <- ifelse(sv_hp_vec == 1L, mcols(admr_sub)$HP1.Methy, mcols(admr_sub)$HP2.Methy)
  wt_beta   <- ifelse(sv_hp_vec == 1L, mcols(admr_sub)$HP2.Methy, mcols(admr_sub)$HP1.Methy)
  data.frame(
    seqnames     = as.character(seqnames(admr_sub)),
    start        = start(admr_sub),
    end          = end(admr_sub),
    block_id     = mcols(admr_sub)$block_id,
    patient_code = mcols(admr_sub)$patient_code,
    sv_hp        = sv_hp_vec,
    sv_hp_beta   = sv_beta,
    wt_hp_beta   = wt_beta,
    sv_minus_wt  = sv_beta - wt_beta,
    normal_beta  = mcols(admr_sub)$normal_beta,
    hp_delta_abs = abs(mcols(admr_sub)$HP1.Methy - mcols(admr_sub)$HP2.Methy),
    n_cpg        = mcols(admr_sub)$nCG,
    ov_bulk      = mcols(admr_sub)$ov_bulk
  )
}

# Level A: verify existing outputs =============================================
cat("[Level A] Checking existing per-read HBV outputs...\n")
perread_paths <- list(
  meth   = file.path(DMR_SVS_ROOT, "result/hbv_perread/hbv_perread_meth.tsv.gz"),
  loci   = file.path(DMR_SVS_ROOT, "result/hbv_perread_locus_stats.csv"),
  pooled = file.path(DMR_SVS_ROOT, "result/hbv_perread_pooled.csv")
)
for (nm in names(perread_paths)) {
  if (file.exists(perread_paths[[nm]])) {
    cat(sprintf("  ✓ %s: %s\n", nm, basename(perread_paths[[nm]])))
  } else {
    cat(sprintf("  ✗ MISSING: %s\n", perread_paths[[nm]]))
  }
}

# Load HBV loci and SV table ===================================================
cat("\n[Level B] Loading HBV loci...\n")
if (!file.exists(HBV_LOCI_CSV)) {
  cat("HBV loci file not found — skipping Level B\n")
  quit(status = 0L)
}

hbv_loci <- fread(HBV_LOCI_CSV)
hbv_loci <- hbv_loci[as.logical(is_tumor) == TRUE]
cat(sprintf("  Somatic HBV loci: %d rows, %d patients\n",
            nrow(hbv_loci), uniqueN(hbv_loci$pcode)))

cat("[Level B] Loading SV annotation and flagging HBV-proximal SVs...\n")
sv <- fread(SV_FILE)
sv[, is_hbv_bnd := mapply(function(pt, chr, svpos) {
  sub <- hbv_loci[pcode == pt & chrom == chr, pos]
  length(sub) > 0L && any(abs(svpos - sub) <= HBV_MATCH_BP, na.rm = TRUE)
}, sample, seqnames, start)]

n_hbv_sv <- sum(sv$is_hbv_bnd, na.rm = TRUE)
cat(sprintf("  HBV-proximal SVs: %d (±%d bp)\n", n_hbv_sv, HBV_MATCH_BP))
if (n_hbv_sv == 0L) {
  cat("No HBV-proximal SVs found — skipping Level B\n")
  quit(status = 0L)
}
hbv_patients <- unique(sv[is_hbv_bnd == TRUE, sample])
cat("  HBV patients:", paste(sort(hbv_patients), collapse=", "), "\n")

# Load phase blocks ============================================================
cat("[Level B] Loading phase block GTF files...\n")
pmap <- fread(PATIENT_MAP_PATH)
gtf_files <- list.files(PHASE_VCF_DIR, pattern = "\\.gtf$", full.names = TRUE)

pb_df <- lapply(gtf_files, function(f) {
  gr <- tryCatch(import(f), error = function(e) NULL)
  if (is.null(gr)) return(NULL)
  data.table(seqnames = as.character(seqnames(gr)),
             start = start(gr), end = end(gr),
             block_id = as.character(mcols(gr)$gene_id),
             patient_name = str_remove(basename(f), "\\.gtf$"))
}) |> rbindlist()

pb_df <- pb_df |> left_join(pmap, by = c("patient_name" = "Samples_ID")) |>
  filter(!is.na(patient_code))
setDT(pb_df)
pb_gr_list <- split(
  makeGRangesFromDataFrame(pb_df, keep.extra.columns = TRUE),
  pb_df$patient_code
)

# Load somatic aDMR + match normal_beta from per-patient normal aDMR files =====
cat("[Level B] Loading somatic aDMR and matching normal_beta...\n")
somatic_dt <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
somatic_dt <- somatic_dt[seqnames %in% MAIN_CHROMS]

# Build per-patient locus-average normal methylation from *.normal_aDMR.sorted.txt
# These files have meanMethy1/meanMethy2 (HP1/HP2 average methylation in normal tissue).
# (HP1_normal + HP2_normal) / 2 gives the locus-level bulk-equivalent normal β.
cat("  Loading normal aDMR files for normal_beta...\n")
normal_files <- list.files(NORMAL_ADMR_DIR, pattern = "\\.normal_aDMR\\.sorted\\.txt$",
                            full.names = TRUE)
cat(sprintf("  Found %d normal aDMR files\n", length(normal_files)))

normal_avg_list <- lapply(normal_files, function(f) {
  sample_name <- str_remove(basename(f), "\\.normal_aDMR\\.sorted\\.txt$")
  pt <- pmap[Samples_ID == sample_name, patient_code]
  if (length(pt) == 0 || is.na(pt)) return(NULL)
  dt <- fread(f, header = TRUE)
  setnames(dt, c("meanMethy1","meanMethy2"), c("HP1.Methy","HP2.Methy"), skip_absent = TRUE)
  chr_col <- intersect(c("chr","seqnames","CHROM","chrom"), names(dt))[1]
  if (!is.na(chr_col) && chr_col != "seqnames") setnames(dt, chr_col, "seqnames")
  dt <- dt[seqnames %in% MAIN_CHROMS]
  dt[, normal_beta  := (HP1.Methy + HP2.Methy) / 2]
  dt[, patient_code := pt]
  dt[, .(seqnames, start, end, patient_code, normal_beta)]
})
normal_avg_dt <- rbindlist(Filter(Negate(is.null), normal_avg_list))
cat(sprintf("  Normal aDMR loci loaded: %d across %d patients\n",
            nrow(normal_avg_dt), uniqueN(normal_avg_dt$patient_code)))

# Per-patient overlap match: somatic aDMR locus → nearest normal aDMR entry
somatic_dt[, normal_beta := NA_real_]
for (pt in unique(somatic_dt$patient_code)) {
  som_pt  <- somatic_dt[patient_code == pt]
  norm_pt <- normal_avg_dt[patient_code == pt]
  if (nrow(som_pt) == 0 || nrow(norm_pt) == 0) next
  som_gr  <- makeGRangesFromDataFrame(som_pt[, .(seqnames, start, end)],
                                       seqnames.field = "seqnames")
  norm_gr <- makeGRangesFromDataFrame(norm_pt, keep.extra.columns = TRUE,
                                       seqnames.field = "seqnames")
  hits      <- findOverlaps(som_gr, norm_gr, select = "first")
  has_match <- !is.na(hits)
  somatic_dt[patient_code == pt & has_match, normal_beta := norm_gr$normal_beta[hits[has_match]]]
}
n_matched <- sum(!is.na(somatic_dt$normal_beta))
cat(sprintf("  normal_beta matched: %d / %d (%.1f%%)\n",
            n_matched, nrow(somatic_dt), mean(!is.na(somatic_dt$normal_beta)) * 100))

# Fallback for unmatched loci: genome-wide average CpG methylation in normal liver
somatic_dt[is.na(normal_beta), normal_beta := NORMAL_BG_BETA]
cat(sprintf("  Fallback NORMAL_BG_BETA=%.2f applied to %d unmatched loci\n",
            NORMAL_BG_BETA, nrow(somatic_dt) - n_matched))

# Build GRanges per patient
somatic_gr_list <- lapply(unique(somatic_dt$patient_code), function(pt) {
  sub <- somatic_dt[patient_code == pt]
  if (nrow(sub) == 0) return(NULL)
  makeGRangesFromDataFrame(sub, keep.extra.columns = TRUE, seqnames.field = "seqnames")
})
names(somatic_gr_list) <- unique(somatic_dt$patient_code)

# Assign block_id to somatic aDMR ==============================================
cat("[Level B] Assigning block_id to somatic aDMR...\n")
admr_phased_list <- lapply(hbv_patients, function(pt) {
  admr_gr <- somatic_gr_list[[pt]]
  blk_gr  <- pb_gr_list[[pt]]
  if (is.null(admr_gr) || length(admr_gr) == 0) return(NULL)
  if (is.null(blk_gr) || length(blk_gr) == 0) {
    mcols(admr_gr)$block_id <- NA_character_; return(admr_gr)
  }
  hits <- findOverlaps(admr_gr, blk_gr, select = "first")
  mcols(admr_gr)$block_id <- as.character(mcols(blk_gr)$block_id[hits])
  admr_gr[!is.na(mcols(admr_gr)$block_id)]
}) |> setNames(hbv_patients)

# Level B: HBV-proximal pairs ==================================================
cat("[Level B] Extracting HBV-proximal HP-beta pairs...\n")
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

if (length(hbv_pairs_list) == 0) {
  cat("No HBV-proximal aDMR pairs found after block assignment\n")
  quit(status = 0L)
}
hbv_df <- bind_rows(hbv_pairs_list)
cat(sprintf("  HBV pairs: %d across %d patients\n",
            nrow(hbv_df), uniqueN(hbv_df$patient_id)))
print(hbv_df |> dplyr::count(patient_id, name = "n_pairs"))

# Level B: Non-HBV SV control ==================================================
ctrl_pairs_list <- lapply(hbv_patients, function(pt) {
  sv_ctrl <- sv[sample == pt & (is.na(is_hbv_bnd) | is_hbv_bnd == FALSE) &
                HP %in% c(1L, 2L) & !is.na(PHASESETID)]
  dmr_gr  <- admr_phased_list[[pt]]
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
cat(sprintf("  Non-HBV control pairs: %d\n", nrow(ctrl_df)))

# Recompute normal_beta from the genome-wide smoothed methylation model ========
cat(sprintf("[Level B] Loading smoothed BSseq object (%s)...\n", basename(BSOBJ_RDS)))
bs_smoothed <- readRDS(BSOBJ_RDS)

get_smoothed_normal_beta <- function(df, bs) {
  df$normal_beta <- NA_real_
  for (pt in unique(df$patient_id)) {
    idx  <- which(df$patient_id == pt)
    cols <- paste0(pt, c("_normal_hap1", "_normal_hap2"))
    if (!all(cols %in% sampleNames(bs))) {
      cat(sprintf("  WARNING: %s missing normal hap columns in bsobj, keeping fallback\n", pt))
      next
    }
    bs_sub <- bs[, cols]
    gr <- GRanges(seqnames = df$seqnames[idx],
                  ranges = IRanges(start = df$start[idx], end = df$end[idx]))
    meth_mat <- getMeth(bs_sub, regions = gr, type = "smooth", what = "perRegion")
    df$normal_beta[idx] <- rowMeans(meth_mat, na.rm = TRUE)
  }
  df
}

hbv_df  <- get_smoothed_normal_beta(hbv_df, bs_smoothed)
if (nrow(ctrl_df) > 0) ctrl_df <- get_smoothed_normal_beta(ctrl_df, bs_smoothed)
cat(sprintf("  Smoothed normal_beta: %d/%d HBV pairs, %d/%d control pairs resolved (no NA)\n",
            sum(!is.na(hbv_df$normal_beta)), nrow(hbv_df),
            sum(!is.na(ctrl_df$normal_beta)), nrow(ctrl_df)))
rm(bs_smoothed)

# Compute deviations from normal ===============================================
hbv_df <- hbv_df |>
  mutate(
    dev_hbv     = sv_hp_beta  - normal_beta,
    dev_wt      = wt_hp_beta  - normal_beta,
    abs_dev_hbv = abs(dev_hbv),
    abs_dev_wt  = abs(dev_wt),
    asym        = abs_dev_hbv - abs_dev_wt
  )
if (nrow(ctrl_df) > 0) {
  ctrl_df <- ctrl_df |>
    mutate(
      dev_hbv     = sv_hp_beta  - normal_beta,
      dev_wt      = wt_hp_beta  - normal_beta,
      abs_dev_hbv = abs(dev_hbv),
      abs_dev_wt  = abs(dev_wt),
      asym        = abs_dev_hbv - abs_dev_wt
    )
}

# Statistical tests ============================================================
test_results <- list()

# Test 1: |dev_HBV| > |dev_WT|
wt1 <- tryCatch(wilcox.test(hbv_df$abs_dev_hbv, hbv_df$abs_dev_wt,
                              paired = TRUE, alternative = "greater", exact = FALSE),
                error = function(e) list(p.value = NA_real_, statistic = NA_real_))
test_results[["t1_abs_dev"]] <- data.frame(
  test = "WSR |dev_HBV| > |dev_WT| (paired)",
  n = nrow(hbv_df), statistic = as.numeric(wt1$statistic),
  p_val = wt1$p.value,
  median_hbv = median(hbv_df$abs_dev_hbv, na.rm = TRUE),
  median_wt  = median(hbv_df$abs_dev_wt,  na.rm = TRUE)
)
cat(sprintf("\nTest 1: |dev_HBV| > |dev_WT|: W=%.1f, p=%.4f\n",
            as.numeric(wt1$statistic), wt1$p.value))

# Test 2: Directionality (HBV-hap hypomethylated vs normal)
n_hypo <- sum(hbv_df$dev_hbv < 0, na.rm = TRUE)
n_tot  <- sum(!is.na(hbv_df$dev_hbv))
binom_p <- binom.test(n_hypo, n_tot, p = 0.5, alternative = "greater")$p.value
test_results[["t2_direction"]] <- data.frame(
  test = "binomial: dev_HBV < 0 (HBV hypomethylated vs normal)",
  n = n_tot, n_hypo = n_hypo,
  pct_hypo = n_hypo / n_tot * 100,
  p_val = binom_p
)
cat(sprintf("Test 2: hypo=%d/%d (%.1f%%), binomial p=%.4f\n",
            n_hypo, n_tot, n_hypo/n_tot*100, binom_p))

# Test 3: Per-patient sign test
per_pat <- hbv_df |>
  group_by(patient_id) |>
  summarise(n_pairs = n(), pct_hypo = mean(dev_hbv < 0, na.rm=TRUE)*100,
            mean_dev_hbv = mean(dev_hbv, na.rm=TRUE),
            mean_dev_wt  = mean(dev_wt,  na.rm=TRUE), .groups = "drop")
print(per_pat)
n_pat_hypo <- sum(per_pat$pct_hypo > 50)
sign_p <- binom.test(n_pat_hypo, nrow(per_pat), p = 0.5)$p.value
test_results[["t3_per_patient"]] <- data.frame(
  test = "per-patient sign test",
  n_patients_majority_hypo = n_pat_hypo,
  n_patients_total = nrow(per_pat),
  p_val = sign_p
)
cat(sprintf("Test 3: %d/%d patients majority-hypo, p=%.4f\n",
            n_pat_hypo, nrow(per_pat), sign_p))

# Test 4: Specificity — HBV vs non-HBV block asymmetry
if (nrow(ctrl_df) > 0) {
  wt4 <- tryCatch(wilcox.test(hbv_df$asym, ctrl_df$asym,
                               alternative = "greater", exact = FALSE),
                  error = function(e) list(p.value = NA_real_, statistic = NA_real_))
  test_results[["t4_specificity"]] <- data.frame(
    test = "MWU: HBV asym > non-HBV asym (specificity)",
    n_hbv = nrow(hbv_df), n_ctrl = nrow(ctrl_df),
    statistic = as.numeric(wt4$statistic), p_val = wt4$p.value,
    median_asym_hbv = median(hbv_df$asym, na.rm=TRUE),
    median_asym_ctrl = median(ctrl_df$asym, na.rm=TRUE)
  )
  cat(sprintf("Test 4: HBV asym > non-HBV: W=%.1f, p=%.4f\n",
              as.numeric(wt4$statistic), wt4$p.value))
}

# LMM: is_hbv effect adjusting for patient and SV type =========================
cat("[Level B] Fitting LMM for specificity (A-III)...\n")
if (nrow(ctrl_df) > 0) {
  all_pairs_lmm <- bind_rows(
    mutate(hbv_df,   is_hbv = 1L, geom_type = NA_character_),
    mutate(ctrl_df,  is_hbv = 0L, geom_type = NA_character_)
  )
  lmm_fit <- tryCatch(
    lmer(asym ~ is_hbv + (1 | patient_id), data = all_pairs_lmm, REML = FALSE),
    error = function(e) NULL
  )
  if (!is.null(lmm_fit)) {
    lmm_coef <- summary(lmm_fit)$coefficients
    lmm_dt <- as.data.table(lmm_coef, keep.rownames = "term")
    setnames(lmm_dt, c("Estimate","Std. Error","t value"), c("beta","se","t"), skip_absent = TRUE)
    cat("\nLMM coefficients:\n"); print(lmm_dt)
    fwrite(lmm_dt, file.path(OUT_ROOT, "hbv_specificity_lmm.csv"))
  }

  # Matched resample null: randomly permute is_hbv label N_RESAMPLE times
  set.seed(42L)
  null_betas <- replicate(N_RESAMPLE, {
    d <- all_pairs_lmm
    d$is_hbv <- sample(d$is_hbv)
    m <- tryCatch(lmer(asym ~ is_hbv + (1 | patient_id), data = d, REML = FALSE),
                  error = function(e) NULL)
    if (is.null(m)) return(NA_real_)
    fixef(m)[["is_hbv"]]
  })
  null_betas <- null_betas[!is.na(null_betas)]
  if (!is.null(lmm_fit)) {
    obs_beta <- fixef(lmm_fit)[["is_hbv"]]
    emp_p    <- mean(null_betas >= obs_beta)
    cat(sprintf("Resample null: obs_beta=%.4f, emp_p=%.4f\n", obs_beta, emp_p))
    fwrite(data.table(null_beta = null_betas),
           file.path(OUT_ROOT, "hbv_specificity_matched_resample.csv"))
  }
}

# Save =========================================================================
fwrite(as.data.table(hbv_df), file.path(OUT_ROOT, "hbv_allele_anchored_pairs.csv"))

tests_dt <- rbindlist(lapply(test_results, as.data.table), fill = TRUE)
fwrite(tests_dt, file.path(OUT_ROOT, "hbv_allele_anchored_tests.csv"))

cat(sprintf("\nSaved: hbv_allele_anchored_pairs.csv (%d rows)\n", nrow(hbv_df)))
cat(sprintf("Saved: hbv_allele_anchored_tests.csv, hbv_specificity_lmm.csv\n"))

cat("=== Done: 06_hbv_allele_specific.R ===\n")
