#!/usr/bin/env Rscript
# tcga_scna_meth_comparison.R
# For each Gold-tier aDMR locus, compare TCGA-LIHC 450k promoter β
# between tumor samples WITH vs WITHOUT a somatic copy-number aberration (SCNA)
# in the same ±50 kb window.
#
# Rationale: If SV-associated allele-specific methylation at Gold loci reflects
# a cis-regulatory mechanism, TCGA samples with focal SCNA at those loci should
# show a similar methylation shift compared to SCNA-free tumors.
#
# Usage: mamba run -n renv Rscript post_processing/tcga_scna_meth_comparison.R \
#          2>&1 | tee logs/tcga_scna_meth.log

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUTDIR    <- "/node200data/kachungk/hcc_data/DMR_SVs"
CACHE_DIR <- file.path(OUTDIR, "tcga_cache")
METH_RDS  <- file.path(CACHE_DIR, "tcga_lihc_meth450k_se.rds")
CNV_RDS   <- file.path(CACHE_DIR, "tcga_lihc_cnv_segment.rds")
GOLD_CSV  <- file.path(OUTDIR, "04.final_candidate/gold_tier_final.csv")
SILV_CSV  <- file.path(OUTDIR, "04.final_candidate/silver_tier.csv")
OUT_CSV   <- file.path(OUTDIR, "result/tcga_scna_vs_meth.csv")
FIG_PNG   <- file.path(OUTDIR, "figs/png/fig_tcga_scna_meth.png")
FIG_PDF   <- file.path(OUTDIR, "figs/panels/fig_tcga_scna_meth.pdf")
DEC_LOG   <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"

WINDOW_BP <- 50000L   # ±50 kb window for SCNA matching
MIN_N_PER_GROUP <- 5L # min samples per group for per-locus test

dir.create(file.path(OUTDIR, "figs/png"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "figs/panels"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "result"),      showWarnings = FALSE, recursive = TRUE)

# ── Module 1: Load Gold-tier aDMR loci ────────────────────────────────────────
message("=== Module 1: Gold-tier aDMR loci ===")
gold_dt <- fread(GOLD_CSV)
silv_dt <- fread(SILV_CSV)

# Use only loci with position info
loci <- rbind(
  gold_dt[!is.na(admr_chr) & !is.na(admr_start) & !is.na(admr_end),
          .(admr_chr, admr_start, admr_end, tier_class, sv_minus_wt,
            patient_code, bp_dist)],
  silv_dt[!is.na(admr_chr) & !is.na(admr_start) & !is.na(admr_end),
          .(admr_chr, admr_start, admr_end, tier_class, sv_minus_wt,
            patient_code, bp_dist)],
  fill = TRUE
)
setnames(loci, c("admr_chr", "admr_start", "admr_end"),
                c("seqnames", "start", "end"))

# Collapse to unique loci (same position across patients = same aDMR)
loci_uniq <- unique(loci[, .(seqnames, start, end)])
loci_uniq[, locus_id := sprintf("%s:%d-%d", seqnames, start, end)]
message(sprintf("  Unique Gold+Silver aDMR loci: %d", nrow(loci_uniq)))

# Per-locus HCC methylation direction (median sv_minus_wt across patients)
loci_dir <- loci[, .(
  hcc_median_svwt = median(sv_minus_wt, na.rm = TRUE),
  hcc_dir         = ifelse(median(sv_minus_wt, na.rm = TRUE) > 0, "Hyper", "Hypo"),
  n_patients      = .N
), by = .(seqnames, start, end)]
loci_dir[, locus_id := sprintf("%s:%d-%d", seqnames, start, end)]

loci_gr <- GRanges(
  seqnames = loci_uniq$seqnames,
  ranges   = IRanges(start = loci_uniq$start, end = loci_uniq$end),
  locus_id = loci_uniq$locus_id
)

# ── Module 2: Load TCGA-LIHC 450k methylation (cached) ───────────────────────
message("\n=== Module 2: TCGA-LIHC 450k methylation ===")
if (!file.exists(METH_RDS)) stop("Methylation SE not found — run fig8_tcga_validation.R first.")
meth_se <- readRDS(METH_RDS)
message(sprintf("  Loaded: %d probes × %d samples", nrow(meth_se), ncol(meth_se)))

# Sample type labels
cd <- as.data.table(colData(meth_se))
cd[, patient_id := substr(barcode, 1, 12)]
if ("sample_type" %in% names(cd)) {
  cd[, sample_type2 := fcase(
    sample_type == "Primary Tumor",       "Tumor",
    sample_type == "Solid Tissue Normal", "Normal",
    default = "Other"
  )]
} else {
  cd[, sample_code := sub(".*-([0-9]{2})[A-Z].*", "\\1", barcode)]
  cd[, sample_type2 := fcase(
    sample_code == "01", "Tumor",
    sample_code == "11", "Normal",
    default = "Other"
  )]
}
tumor_barcodes <- cd[sample_type2 == "Tumor", barcode]
message(sprintf("  Tumor samples: %d", length(tumor_barcodes)))

# Find probes overlapping Gold+Silver loci (±500 bp padding)
probe_gr <- rowRanges(meth_se)
if (is.null(probe_gr) || length(probe_gr) == 0) {
  rd <- as.data.table(rowData(meth_se))
  probe_gr <- GRanges(
    seqnames = rd$Chromosome,
    ranges   = IRanges(start = rd$Start, width = 1L)
  )
  names(probe_gr) <- rownames(meth_se)
}

# Pad loci ±500 bp for probe lookup
loci_padded <- GenomicRanges::resize(loci_gr, width = width(loci_gr) + 1000L, fix = "center")
probe_hits   <- findOverlaps(loci_padded, probe_gr, ignore.strand = TRUE)

if (length(probe_hits) == 0) {
  stop("No 450k probes overlap Gold+Silver loci (±500bp). Check seqname format (chrX vs X).")
}

probe_locus_map <- data.table(
  locus_id = loci_gr$locus_id[queryHits(probe_hits)],
  probe_id  = names(probe_gr)[subjectHits(probe_hits)]
)
n_loci_with_probes <- uniqueN(probe_locus_map$locus_id)
message(sprintf("  Loci with ≥1 probe: %d / %d", n_loci_with_probes, length(loci_gr)))

# Extract beta matrix for matched probes in tumor samples only
probe_ids_keep <- unique(probe_locus_map$probe_id)
tumor_idx      <- which(colnames(meth_se) %in% tumor_barcodes)
beta_mat <- assay(meth_se)[probe_ids_keep, tumor_idx, drop = FALSE]

# Aggregate to locus-level mean β per tumor sample
beta_dt   <- as.data.table(beta_mat, keep.rownames = "probe_id")
beta_long <- melt(beta_dt, id.vars = "probe_id",
                  variable.name = "barcode", value.name = "beta")
beta_long <- merge(beta_long, probe_locus_map, by = "probe_id", allow.cartesian = TRUE)
locus_beta <- beta_long[!is.na(beta),
  .(mean_beta = mean(beta, na.rm = TRUE), n_probes = .N),
  by = .(locus_id, barcode)]
message(sprintf("  Locus-sample β values: %d", nrow(locus_beta)))
rm(beta_mat, beta_dt, beta_long); gc()

# ── Module 3: TCGA-LIHC Copy Number Segment ───────────────────────────────────
message("\n=== Module 3: TCGA-LIHC CN segment data ===")

load_cnv <- function() {
  if (file.exists(CNV_RDS)) {
    message("  [cache] Loading ", basename(CNV_RDS))
    return(readRDS(CNV_RDS))
  }
  message("  Querying GDC for TCGA-LIHC CNV segments ...")
  cnv_query <- tryCatch(
    GDCquery(
      project       = "TCGA-LIHC",
      data.category = "Copy Number Variation",
      data.type     = "Copy Number Segment"
    ),
    error = function(e) {
      message("  [WARN] GDCquery CNV failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(cnv_query)) return(NULL)

  GDCdownload(cnv_query, directory = CACHE_DIR)
  cnv_se <- GDCprepare(cnv_query, directory = CACHE_DIR)
  saveRDS(cnv_se, CNV_RDS)
  message("  Saved CN segment SE to ", basename(CNV_RDS))
  cnv_se
}

cnv_obj <- load_cnv()

if (is.null(cnv_obj)) {
  message("[WARN] CNV download failed — using methylation quartile stratification as fallback")
  USE_CNV <- FALSE
} else {
  USE_CNV <- TRUE
  message(sprintf("  CNV object class: %s", class(cnv_obj)[1]))

  # Parse CN segment into long data.table
  if (is.data.frame(cnv_obj) || is.data.table(cnv_obj)) {
    cnv_dt <- as.data.table(cnv_obj)
  } else if (inherits(cnv_obj, "SummarizedExperiment")) {
    cnv_dt <- as.data.table(assay(cnv_obj), keep.rownames = "probe_id")
    cnv_dt <- melt(cnv_dt, id.vars = "probe_id", variable.name = "barcode",
                   value.name = "segment_mean")
  } else {
    # GDCprepare may return a list of data frames (one per sample)
    cnv_dt <- rbindlist(lapply(names(cnv_obj), function(s) {
      df <- as.data.table(cnv_obj[[s]])
      df[, barcode := s]
      df
    }), fill = TRUE)
  }

  # Standardise column names (GDC uses Chromosome/Start/End/Segment_Mean)
  setnames(cnv_dt,
    intersect(c("Chromosome","Start","End","Segment_Mean","Sample"),
              names(cnv_dt)),
    intersect(c("chr","seg_start","seg_end","log2_ratio","barcode"),
              c("chr","seg_start","seg_end","log2_ratio","barcode")),
    skip_absent = TRUE
  )
  # If still not renamed, try lower-case aliases
  name_map <- c(
    chromosome = "chr", start = "seg_start", end = "seg_end",
    segment_mean = "log2_ratio", sample = "barcode",
    num_probes = "n_probes_seg"
  )
  for (old in names(name_map)) {
    if (old %in% tolower(names(cnv_dt)) && !name_map[[old]] %in% names(cnv_dt)) {
      setnames(cnv_dt, names(cnv_dt)[tolower(names(cnv_dt)) == old], name_map[[old]])
    }
  }

  needed <- c("chr","seg_start","seg_end","log2_ratio","barcode")
  missing_cols <- setdiff(needed, names(cnv_dt))
  if (length(missing_cols) > 0) {
    message("[WARN] Missing CN segment columns: ", paste(missing_cols, collapse=", "),
            " — falling back to methylation quartile stratification")
    USE_CNV <- FALSE
    message("[WARN] Falling back to methylation quartile stratification — set scna_method='meth-quartile' in output")
  } else {
    cnv_dt[, patient_id := substr(barcode, 1, 12)]
    # Keep only tumor samples
    cnv_dt <- cnv_dt[barcode %in% tumor_barcodes | patient_id %in% cd[sample_type2=="Tumor", patient_id]]
    message(sprintf("  CNV segments: %d rows, %d tumor barcodes",
                    nrow(cnv_dt), uniqueN(cnv_dt$barcode)))
  }
}

# ── Module 4: Classify TCGA samples by SCNA at each Gold locus ───────────────
message("\n=== Module 4: SCNA classification per locus ===")

if (USE_CNV) {
  # Build SCNA GRanges
  cnv_dt[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  cnv_gr <- GRanges(
    seqnames = cnv_dt$chr,
    ranges   = IRanges(start = cnv_dt$seg_start, end = cnv_dt$seg_end),
    log2_ratio = cnv_dt$log2_ratio,
    barcode    = cnv_dt$barcode
  )

  # Pad loci ±WINDOW_BP for SCNA matching
  loci_win <- GenomicRanges::resize(loci_gr, width = width(loci_gr) + 2L * WINDOW_BP, fix = "center")
  scna_hits <- findOverlaps(loci_win, cnv_gr, ignore.strand = TRUE)

  scna_dt <- data.table(
    locus_id   = loci_gr$locus_id[queryHits(scna_hits)],
    barcode    = cnv_gr$barcode[subjectHits(scna_hits)],
    log2_ratio = cnv_gr$log2_ratio[subjectHits(scna_hits)]
  )
  # Use 12-char patient_id for joining (meth and CNV barcodes have different suffixes)
  scna_dt[, patient_id := substr(barcode, 1, 12)]
  # A locus is "SCNA+" if |log2_ratio| > 0.3 (focal aberration threshold)
  SCNA_THR <- 0.3
  scna_max <- scna_dt[, .(max_abs_log2 = max(abs(log2_ratio), na.rm = TRUE)),
                      by = .(locus_id, patient_id)]
  scna_max[, scna_status := ifelse(max_abs_log2 >= SCNA_THR, "SCNA+", "SCNA-")]
  message(sprintf("  SCNA+ locus-patient pairs: %d", sum(scna_max$scna_status == "SCNA+")))
} else {
  # Fallback: classify by methylation quartile (top/bottom 25% = "extreme")
  message("  [Fallback] Using methylation quartile stratification (no CNV data)")
  scna_max <- locus_beta[, {
    q_hi <- quantile(mean_beta, 0.75, na.rm = TRUE)
    q_lo <- quantile(mean_beta, 0.25, na.rm = TRUE)
    .(barcode    = barcode,
      max_abs_log2 = NA_real_,
      scna_status = fcase(
        mean_beta >= q_hi, "Extreme-Hi",
        mean_beta <= q_lo, "Extreme-Lo",
        default           = "Mid"
      ))
  }, by = locus_id]
}

# ── Module 5: Wilcoxon per locus (SCNA+ vs SCNA-) ────────────────────────────
message("\n=== Module 5: Per-locus Wilcoxon test ===")

# Add 12-char patient_id to locus_beta for joining with SCNA data
locus_beta[, patient_id := substr(barcode, 1, 12)]
merged_dt <- merge(locus_beta, scna_max[scna_status %in% c("SCNA+", "SCNA-")],
                   by = c("locus_id", "patient_id"))

locus_test <- merged_dt[, {
  scna_pos <- mean_beta[scna_status == "SCNA+"]
  scna_neg <- mean_beta[scna_status == "SCNA-"]
  n_pos    <- length(scna_pos)
  n_neg    <- length(scna_neg)
  if (n_pos < MIN_N_PER_GROUP || n_neg < MIN_N_PER_GROUP) {
    .(n_scna_pos   = n_pos,
      n_scna_neg   = n_neg,
      mean_beta_pos = NA_real_,
      mean_beta_neg = NA_real_,
      delta_beta    = NA_real_,
      wilcox_p      = NA_real_,
      tested        = FALSE)
  } else {
    wt <- suppressWarnings(wilcox.test(scna_pos, scna_neg, exact = FALSE))
    .(n_scna_pos    = n_pos,
      n_scna_neg    = n_neg,
      mean_beta_pos = mean(scna_pos, na.rm = TRUE),
      mean_beta_neg = mean(scna_neg, na.rm = TRUE),
      delta_beta    = mean(scna_pos, na.rm = TRUE) - mean(scna_neg, na.rm = TRUE),
      wilcox_p      = wt$p.value,
      tested        = TRUE)
  }
}, by = locus_id]

locus_test_tested <- locus_test[tested == TRUE]
locus_test_tested[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]
message(sprintf("  Loci tested: %d / %d", nrow(locus_test_tested), nrow(locus_test)))
message(sprintf("  Loci FDR < 0.20: %d", sum(locus_test_tested$wilcox_fdr < 0.20, na.rm = TRUE)))
message(sprintf("  Loci FDR < 0.05: %d", sum(locus_test_tested$wilcox_fdr < 0.05, na.rm = TRUE)))

# Add HCC direction from Gold loci
locus_test_tested <- merge(locus_test_tested,
  loci_dir[, .(locus_id, hcc_dir, hcc_median_svwt, n_patients)],
  by = "locus_id", all.x = TRUE)

# Concordance: SCNA-associated β direction matches HCC allelic direction
locus_test_tested[, tcga_dir := fcase(
  delta_beta >  0, "Hyper",
  delta_beta <= 0, "Hypo",
  default = NA_character_
)]
locus_test_tested[, concordant := !is.na(tcga_dir) & !is.na(hcc_dir) &
                                    (tcga_dir == hcc_dir)]
n_conc <- sum(locus_test_tested$concordant, na.rm = TRUE)
n_tot  <- sum(!is.na(locus_test_tested$tcga_dir) & !is.na(locus_test_tested$hcc_dir))
message(sprintf("  Directional concordance (TCGA SCNA vs HCC allelic): %d / %d",
                n_conc, n_tot))

# ── Save ──────────────────────────────────────────────────────────────────────
full_res <- merge(locus_test, loci_dir[, .(locus_id, hcc_dir, hcc_median_svwt, n_patients)],
                  by = "locus_id", all.x = TRUE)
full_res[, scna_method := ifelse(USE_CNV, "CNV-based", "meth-quartile")]
fwrite(full_res, OUT_CSV)
message(sprintf("\nResults saved: %s", OUT_CSV))

# ── Figures ───────────────────────────────────────────────────────────────────
theme_hcc <- theme_classic(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_blank(),
    legend.background = element_blank(),
    plot.title  = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey50", size = 9)
  )

if (nrow(locus_test_tested) < 3) {
  message("[WARN] Too few tested loci for figures — skipping visualization")
} else {
  # Panel A: Δβ (SCNA+ − SCNA-) forest plot, ordered by delta
  dat_a <- locus_test_tested[!is.na(delta_beta)][order(delta_beta)]
  dat_a[, locus_short := sub("chr(.*):.*", "chr\\1", locus_id)]
  dat_a[, sig_label := fcase(
    wilcox_fdr < 0.05,  "*",
    wilcox_fdr < 0.20,  ".",
    default             = ""
  )]
  dat_a[, concordant_label := ifelse(concordant, "concordant", "discordant")]

  p_forest <- ggplot(dat_a, aes(x = delta_beta,
                                y = reorder(locus_short, delta_beta),
                                colour = hcc_dir)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_segment(aes(xend = 0,
                     yend = reorder(locus_short, delta_beta)),
                 linewidth = 0.5, alpha = 0.5) +
    geom_point(size = 2.5) +
    geom_text(aes(label = sig_label),
              nudge_x = sign(dat_a$delta_beta) * 0.003,
              size = 4, vjust = 0.5) +
    scale_colour_manual(
      values   = c(Hyper = "#D6604D", Hypo = "#4393C3"),
      na.value = "grey50",
      name     = "HCC allelic\ndirection"
    ) +
    labs(
      title    = "A. SCNA+ vs SCNA- methylation at Gold+Silver aDMRs",
      subtitle = sprintf(
        "Δβ = mean(SCNA+) − mean(SCNA-) | concordance = %d/%d | *FDR<0.05  .FDR<0.20",
        n_conc, n_tot),
      x = "Δβ (SCNA+ − SCNA-)", y = NULL
    ) +
    theme_hcc +
    theme(axis.text.y = element_text(size = 7))

  # Panel B: volcano-style (Δβ vs -log10 FDR)
  p_volcano <- ggplot(locus_test_tested[!is.na(delta_beta) & !is.na(wilcox_fdr)],
                      aes(x = delta_beta, y = -log10(wilcox_fdr + 1e-10),
                          colour = hcc_dir, shape = concordant)) +
    geom_hline(yintercept = -log10(0.20), linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_manual(
      values   = c(Hyper = "#D6604D", Hypo = "#4393C3"),
      na.value = "grey50", name = "HCC direction"
    ) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                       name = "Concordant") +
    labs(
      title    = "B. Volcano: SCNA-associated methylation shift",
      subtitle = "Dashed line: FDR = 0.20",
      x = "Δβ (SCNA+ − SCNA-)", y = "-log10(FDR)"
    ) +
    theme_hcc

  # Panel C: violin by SCNA status (pooled across all tested loci)
  if (nrow(merged_dt) > 0) {
    p_violin <- ggplot(merged_dt, aes(x = scna_status, y = mean_beta, fill = scna_status)) +
      geom_violin(alpha = 0.7, colour = "grey40", trim = TRUE, scale = "width") +
      geom_boxplot(width = 0.08, fill = "white", colour = "grey30", alpha = 0.9,
                   outlier.size = 0.3) +
      scale_fill_manual(
        values = c("SCNA+" = "#E24B4A", "SCNA-" = "#3B8BD4"),
        guide = "none"
      ) +
      labs(
        title    = "C. 450k β at Gold+Silver loci (pooled)",
        subtitle = sprintf("All tested locus-sample pairs (n=%d)", nrow(merged_dt)),
        x = NULL, y = "Mean promoter β"
      ) +
      theme_hcc
  } else {
    p_violin <- ggplot() + labs(title = "C. No merged data available") + theme_hcc
  }

  fig_main <- (p_forest | (p_volcano / p_violin)) +
    plot_annotation(
      title   = "TCGA-LIHC: SCNA at Gold-tier aDMR loci associates with methylation shift",
      caption = sprintf(
        "n=%d TCGA Tumor samples | window=±%dkb | FDR BH | concordance (HCC allelic vs TCGA SCNA dir): %d/%d",
        length(tumor_barcodes), WINDOW_BP / 1000L, n_conc, n_tot
      ),
      tag_levels = "A",
      theme      = theme(plot.title = element_text(face = "bold", size = 13))
    ) +
    plot_layout(widths = c(2, 1))

  ggsave(FIG_PNG, fig_main, width = 14, height = 8, dpi = 200)
  ggsave(FIG_PDF, fig_main, width = 14, height = 8)
  message(sprintf("Figure saved: %s", FIG_PNG))
}

# ── Decision log ──────────────────────────────────────────────────────────────
cat(sprintf(
  "[%s] TCGA SCNA-meth: %d Gold+Silver loci; %d tested; FDR<0.05 n=%d; concordance %d/%d; USE_CNV=%s\n",
  format(Sys.Date()),
  nrow(loci_uniq), nrow(locus_test_tested),
  if (nrow(locus_test_tested) > 0) sum(locus_test_tested$wilcox_fdr < 0.05, na.rm=TRUE) else 0L,
  n_conc, n_tot, USE_CNV
), file = DEC_LOG, append = TRUE)

message("\nDone.")
