#!/usr/bin/env Rscript
# fig8_tcga_validation.R
# TCGA-LIHC validation of Gold gene promoter methylation:
#   - Tumor vs Normal promoter β (Wilcoxon)
#   - Methylation vs expression correlation (Spearman)
#   - Survival analysis by promoter methylation level (KM + logrank)
#
# Usage:
#   mamba run -n renv Rscript viz/fig8_tcga_validation.R 2>&1 | tee figs/logs/fig8_tcga.log
#
# Downloads are cached; re-runs skip GDCdownload if RDS exists.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(survival)
  library(broom)
  library(limma)
  library(edgeR)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

# Paths ========================================================================
OUTDIR    <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")
FIGDIR    <- file.path(OUTDIR, "figs")
PROM_RDS  <- file.path(OUTDIR, "canonical_promoters_hg38.gencode_v49.rds")
GENE_LIST <- file.path(OUTDIR, "result", "genes_gold_silver_admr.csv")
GOLD_CSV   <- file.path(OUTDIR, "result", "gold_tier_final.csv")
SILV_CSV   <- file.path(OUTDIR, "result", "silver_tier.csv")
GOLD_FINAL <- file.path(OUTDIR, "04.final_candidate", "gold_tier_final.csv")
SILV_FINAL <- file.path(OUTDIR, "04.final_candidate", "silver_tier.csv")
CACHE_DIR <- file.path(OUTDIR, "tcga_cache")

for (d in c(CACHE_DIR,
            file.path(FIGDIR, "panels"),
            file.path(FIGDIR, "rds"),
            file.path(FIGDIR, "png"),
            file.path(FIGDIR, "logs"),
            file.path(OUTDIR, "result"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Theme ========================================================================
theme_hcc <- theme_classic(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_blank(),
    legend.background = element_blank()
  )

# Gold gene list (promoter-overlap only) =======================================
gene_dt <- fread(GENE_LIST)
gold_genes_prom <- unique(gene_dt[source == "promoter", gene_name])
gold_genes_prom <- gold_genes_prom[!grepl("^ENSG", gold_genes_prom)]  # drop ENSG IDs
message(sprintf("Gold+Silver promoter genes for validation: %d", length(gold_genes_prom)))
message("  ", paste(gold_genes_prom, collapse = ", "))

# Load our HCC Δβ direction per gene (for concordance panel) ===================
load_hcc_direction <- function() {
  gold_df <- if (file.exists(GOLD_FINAL)) fread(GOLD_FINAL) else data.table()
  silv_df <- if (file.exists(SILV_FINAL)) fread(SILV_FINAL) else data.table()
  all_df <- rbind(gold_df, silv_df, fill = TRUE)
  if (nrow(all_df) == 0 || !"sv_minus_wt" %in% names(all_df)) return(NULL)
  all_df <- all_df[!is.na(admr_chr) & !is.na(sv_minus_wt)]
  if (nrow(all_df) == 0 || !file.exists(PROM_RDS)) return(NULL)
  # Overlap gold promoters with aDMRs to get per-gene HCC allelic direction
  admr_gr <- GRanges(seqnames = all_df$admr_chr,
                     ranges   = IRanges(start = all_df$admr_start, end = all_df$admr_end))
  prom_gr   <- readRDS(PROM_RDS)
  gold_prom <- prom_gr[prom_gr$gene_name %in% gold_genes_prom]
  hits <- findOverlaps(gold_prom, admr_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(NULL)
  hit_dt <- data.table(
    gene_name   = gold_prom$gene_name[queryHits(hits)],
    sv_minus_wt = all_df$sv_minus_wt[subjectHits(hits)]
  )
  hit_dt[, .(hcc_dir   = ifelse(median(sv_minus_wt, na.rm = TRUE) > 0, "Hyper", "Hypo"),
             hcc_delta = median(sv_minus_wt, na.rm = TRUE)),
         by = gene_name]
}

# Module 1: Download TCGA-LIHC data (cached) ===================================
message("\n=== Module 1: Download TCGA-LIHC data ===")

METH_RDS <- file.path(CACHE_DIR, "tcga_lihc_meth450k_se.rds")
RNA_RDS  <- file.path(CACHE_DIR, "tcga_lihc_rnaseq_se.rds")
CLIN_RDS <- file.path(CACHE_DIR, "tcga_lihc_clinical.rds")

download_and_prepare <- function(query, rds_path, label) {
  if (file.exists(rds_path)) {
    message(sprintf("  [cache] Loading %s from %s", label, basename(rds_path)))
    return(readRDS(rds_path))
  }
  message(sprintf("  Downloading %s ...", label))
  GDCdownload(query, directory = CACHE_DIR)
  se <- GDCprepare(query, directory = CACHE_DIR)
  saveRDS(se, rds_path)
  message(sprintf("  Saved %s → %s", label, rds_path))
  se
}

# Methylation query: Primary Tumor + Solid Tissue Normal
meth_query <- GDCquery(
  project       = "TCGA-LIHC",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)
meth_res <- getResults(meth_query)
message(sprintf("  Methylation: %d Primary Tumor, %d Solid Tissue Normal",
                sum(meth_res$sample_type == "Primary Tumor"),
                sum(meth_res$sample_type == "Solid Tissue Normal")))

# RNA-seq query
rna_query <- GDCquery(
  project             = "TCGA-LIHC",
  data.category       = "Transcriptome Profiling",
  data.type           = "Gene Expression Quantification",
  workflow.type       = "STAR - Counts",
  sample.type         = c("Primary Tumor", "Solid Tissue Normal")
)
rna_res <- getResults(rna_query)
message(sprintf("  RNA-seq: %d Primary Tumor, %d Solid Tissue Normal",
                sum(rna_res$sample_type == "Primary Tumor"),
                sum(rna_res$sample_type == "Solid Tissue Normal")))

meth_se <- download_and_prepare(meth_query, METH_RDS, "450k methylation")
rna_se  <- download_and_prepare(rna_query,  RNA_RDS,  "RNA-seq")

# Clinical data — extracted from RNA-seq colData (already enriched by GDCprepare)
# GDCquery_clinic has a data.table incompatibility bug; colData(rna_se) is equivalent
get_clinical <- function() {
  if (file.exists(CLIN_RDS)) {
    message("  [cache] Loading clinical from ", basename(CLIN_RDS))
    return(readRDS(CLIN_RDS))
  }
  # Try GDCquery_clinic first; fall back to rna_se colData
  clin <- tryCatch({
    suppressWarnings(as.data.frame(GDCquery_clinic("TCGA-LIHC", type = "clinical")))
  }, error = function(e) {
    message("  [WARN] GDCquery_clinic failed (", conditionMessage(e), "); using rna_se colData")
    as.data.frame(colData(rna_se))
  })
  saveRDS(clin, CLIN_RDS)
  clin
}
clin <- get_clinical()

# Module 2: Methylation processing =============================================
message("\n=== Module 2: Methylation processing ===")

METH_BETA_OUT <- file.path(OUTDIR, "result", "gold_meth_beta.csv.gz")

process_methylation <- function() {
  if (file.exists(METH_BETA_OUT)) {
    message("  [cache] Loading ", basename(METH_BETA_OUT))
    return(fread(METH_BETA_OUT))
  }

  message("  Building probe GRanges from rowRanges ...")
  probe_gr <- rowRanges(meth_se)
  if (is.null(probe_gr) || length(probe_gr) == 0) {
    # Fallback: rowData has coordinates
    rd <- as.data.table(rowData(meth_se))
    probe_gr <- GRanges(
      seqnames = rd$Chromosome,
      ranges   = IRanges(start = rd$Start, width = 1L)
    )
    names(probe_gr) <- rownames(meth_se)
  }
  message(sprintf("  Total probes: %d", length(probe_gr)))

  # Load Gold gene promoters
  prom_gr <- readRDS(PROM_RDS)
  gold_prom <- prom_gr[prom_gr$gene_name %in% gold_genes_prom]
  message(sprintf("  Gold gene promoters: %d ranges for %d genes",
                  length(gold_prom), uniqueN(gold_prom$gene_name)))

  # Overlap probes with Gold promoters
  prom_hits <- findOverlaps(gold_prom, probe_gr, ignore.strand = TRUE)
  probe_gene_map <- data.table(
    probe_id  = names(probe_gr)[subjectHits(prom_hits)],
    gene_name = gold_prom$gene_name[queryHits(prom_hits)]
  )
  message(sprintf("  Probes in Gold promoters: %d (unique probes: %d)",
                  nrow(probe_gene_map), uniqueN(probe_gene_map$probe_id)))

  if (nrow(probe_gene_map) == 0) {
    warning("No probes overlapping Gold gene promoters — check seqname format")
    return(NULL)
  }

  # Extract beta matrix for those probes
  probe_ids <- unique(probe_gene_map$probe_id)
  beta_mat  <- assay(meth_se)[probe_ids, , drop = FALSE]

  # Sample metadata: barcode → sample type, patient ID
  cd <- as.data.table(colData(meth_se))
  cd[, patient_id := substr(barcode, 1, 12)]
  # Prefer the explicit sample_type column; fall back to barcode position 14-15
  if ("sample_type" %in% names(cd)) {
    cd[, sample_type2 := fcase(
      sample_type == "Primary Tumor",        "Tumor",
      sample_type == "Solid Tissue Normal",  "Normal",
      default = "Other"
    )]
  } else {
    # 4th field of barcode is the sample code: 01x = Tumor, 11x = Normal
    cd[, sample_code := sub(".*-([0-9]{2})[A-Z].*", "\\1", barcode)]
    cd[, sample_type2 := fcase(
      sample_code == "01", "Tumor",
      sample_code == "11", "Normal",
      default = "Other"
    )]
  }
  message(sprintf("  Sample types: %s",
    paste(cd[, .N, by = sample_type2][order(-N), sprintf("%s=%d", sample_type2, N)],
          collapse = ", ")))

  # Aggregate: mean β per gene per sample
  beta_dt <- as.data.table(beta_mat, keep.rownames = "probe_id")
  beta_long <- melt(beta_dt, id.vars = "probe_id", variable.name = "barcode", value.name = "beta")
  beta_long <- merge(beta_long, probe_gene_map, by = "probe_id")
  gene_beta <- beta_long[, .(mean_beta = mean(beta, na.rm = TRUE),
                             n_probes  = .N),
                         by = .(gene_name, barcode)]
  gene_beta <- merge(gene_beta, cd[, .(barcode, patient_id, sample_type2)],
                     by = "barcode")

  fwrite(gene_beta, METH_BETA_OUT, compress = "gzip")
  message(sprintf("  Saved: %d gene-sample pairs → %s", nrow(gene_beta), METH_BETA_OUT))
  gene_beta
}

gene_beta <- process_methylation()

if (is.null(gene_beta) || nrow(gene_beta) == 0) {
  stop("Methylation processing failed — no gene-sample beta values. Check probe/promoter seqname format.")
}
message(sprintf("  Genes with ≥1 probe: %d; samples: %d",
                uniqueN(gene_beta$gene_name), uniqueN(gene_beta$barcode)))

# Module 3: Expression processing ==============================================
message("\n=== Module 3: Expression processing ===")

EXPR_OUT <- file.path(OUTDIR, "result", "gold_expr_tpm.csv.gz")

process_expression <- function() {
  if (file.exists(EXPR_OUT)) {
    message("  [cache] Loading ", basename(EXPR_OUT))
    return(fread(EXPR_OUT))
  }

  message("  Extracting TPM for Gold genes ...")
  # STAR - Counts SE: assay names include tpm_unstrand, unstranded, etc.
  anames <- assayNames(rna_se)
  message("  Available assays: ", paste(anames, collapse = ", "))
  tpm_assay <- if ("tpm_unstrand" %in% anames) "tpm_unstrand" else anames[1]

  # rowData has gene_id and gene_name
  rd_rna <- as.data.table(rowData(rna_se))
  target_rows <- which(rd_rna$gene_name %in% gold_genes_prom)
  message(sprintf("  Gold gene rows found: %d", length(target_rows)))

  if (length(target_rows) == 0) {
    # Try case-insensitive
    target_rows <- which(toupper(rd_rna$gene_name) %in% toupper(gold_genes_prom))
    message(sprintf("  (case-insensitive) Gold gene rows found: %d", length(target_rows)))
  }

  tpm_mat <- assay(rna_se, tpm_assay)[target_rows, , drop = FALSE]
  rownames(tpm_mat) <- rd_rna$gene_name[target_rows]

  # Sample metadata
  cd_rna <- as.data.table(colData(rna_se))
  cd_rna[, patient_id := substr(barcode, 1, 12)]
  if ("sample_type" %in% names(cd_rna)) {
    cd_rna[, sample_type2 := fcase(
      sample_type == "Primary Tumor",       "Tumor",
      sample_type == "Solid Tissue Normal", "Normal",
      default = "Other"
    )]
  } else {
    cd_rna[, sample_code := sub(".*-([0-9]{2})[A-Z].*", "\\1", barcode)]
    cd_rna[, sample_type2 := fcase(
      sample_code == "01", "Tumor",
      sample_code == "11", "Normal",
      default = "Other"
    )]
  }

  expr_dt <- as.data.table(log2(tpm_mat + 1), keep.rownames = "gene_name")
  expr_long <- melt(expr_dt, id.vars = "gene_name", variable.name = "barcode",
                    value.name = "log2tpm")
  expr_long <- merge(expr_long, cd_rna[, .(barcode, patient_id, sample_type2)],
                     by = "barcode")

  fwrite(expr_long, EXPR_OUT, compress = "gzip")
  message(sprintf("  Saved: %d gene-sample pairs → %s", nrow(expr_long), EXPR_OUT))
  expr_long
}

expr_long <- process_expression()

# Module 4: Statistical analyses ===============================================
message("\n=== Module 4: Statistical analyses ===")

# 4a. Tumor vs Normal Δβ per gene (Wilcoxon)
message("  4a. Tumor vs Normal Δβ ...")
tn_beta <- gene_beta[sample_type2 %in% c("Tumor", "Normal")]

tn_stats <- tn_beta[, {
  tumor_vals  <- mean_beta[sample_type2 == "Tumor"  & !is.na(mean_beta)]
  normal_vals <- mean_beta[sample_type2 == "Normal" & !is.na(mean_beta)]
  if (length(tumor_vals) < 3 || length(normal_vals) < 2) {
    .(n_tumor = length(tumor_vals), n_normal = length(normal_vals),
      mean_tumor = NA_real_, mean_normal = NA_real_,
      delta_beta = NA_real_, wilcox_p = NA_real_)
  } else {
    wt <- suppressWarnings(wilcox.test(tumor_vals, normal_vals, exact = FALSE))
    .(n_tumor   = length(tumor_vals),
      n_normal  = length(normal_vals),
      mean_tumor  = mean(tumor_vals,  na.rm = TRUE),
      mean_normal = mean(normal_vals, na.rm = TRUE),
      delta_beta  = mean(tumor_vals, na.rm = TRUE) - mean(normal_vals, na.rm = TRUE),
      wilcox_p    = wt$p.value)
  }
}, by = gene_name]
tn_stats[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]
tn_stats[, tcga_dir := ifelse(delta_beta > 0, "Hyper", "Hypo")]
setorder(tn_stats, wilcox_fdr)
message("  Genes with |Δβ| ≥ 0.05: ", sum(abs(tn_stats$delta_beta) >= 0.05, na.rm = TRUE))

fwrite(tn_stats, file.path(OUTDIR, "result", "tcga_tn_delta_beta.csv"))

# 4b. Methylation–expression correlation (tumor only)
message("  4b. Methylation–expression correlation (Spearman) ...")
tumor_meth <- gene_beta[sample_type2 == "Tumor",
                        .(gene_name, patient_id, mean_beta)]
tumor_expr <- expr_long[sample_type2 == "Tumor",
                        .(gene_name, patient_id, log2tpm)]

# Match on patient_id (12-char barcode prefix)
meth_expr <- merge(tumor_meth, tumor_expr, by = c("gene_name", "patient_id"))

cor_stats <- meth_expr[, {
  ok <- !is.na(mean_beta) & !is.na(log2tpm) & is.finite(mean_beta) & is.finite(log2tpm)
  n_ok <- sum(ok)
  if (n_ok < 10) {
    .(n = n_ok, spearman_rho = NA_real_, spearman_p = NA_real_)
  } else {
    ct <- cor.test(mean_beta[ok], log2tpm[ok], method = "spearman", exact = FALSE)
    .(n = n_ok, spearman_rho = ct$estimate, spearman_p = ct$p.value)
  }
}, by = gene_name]
cor_stats[, spearman_fdr := p.adjust(spearman_p, method = "BH")]
setorder(cor_stats, spearman_rho)
message("  Genes with ρ < -0.2: ", sum(cor_stats$spearman_rho < -0.2, na.rm = TRUE))

fwrite(cor_stats, file.path(OUTDIR, "result", "tcga_meth_expr_cor.csv"))

# 4c. Survival analysis: KM by median methylation split
message("  4c. Survival analysis ...")

# Prepare clinical
clin_dt <- as.data.table(clin)

# Normalise patient ID column — works for both GDCquery_clinic and rna_se colData
id_candidates <- c("bcr_patient_barcode", "submitter_id", "patient_id",
                   "case_submitter_id", "cases.submitter_id", "barcode")
id_col <- intersect(id_candidates, names(clin_dt))[1]
if (!is.na(id_col) && id_col != "bcr_patient_barcode")
  clin_dt[, bcr_patient_barcode := get(id_col)]
# If from colData, barcode is 28-char; trim to 12-char patient ID
if ("bcr_patient_barcode" %in% names(clin_dt))
  clin_dt[, bcr_patient_barcode := substr(bcr_patient_barcode, 1, 12)]

# OS columns — GDCquery_clinic uses days_to_death / days_to_last_follow_up
# colData from rna_se may use same or slightly different names
time_col   <- intersect(c("days_to_death", "days_to_last_follow_up",
                          "days_to_last_followup"), names(clin_dt))
status_col <- intersect(c("vital_status"), names(clin_dt))

if (length(time_col) == 0 || length(status_col) == 0) {
  message("  [WARN] Clinical OS columns not found; skipping survival analysis")
  surv_res <- NULL
} else {
  time_col <- time_col[1]
  clin_surv <- clin_dt[, .(
    patient_id  = bcr_patient_barcode,
    vital_status = get(status_col[1]),
    os_days = suppressWarnings(as.numeric(days_to_death))
  )]
  # Use days_to_last_followup for censored
  lf_col <- intersect(c("days_to_last_follow_up", "days_to_last_followup"), names(clin_dt))
  if (length(lf_col) > 0) {
    clin_surv[, lf_days := suppressWarnings(as.numeric(clin_dt[[lf_col[1]]]))]
    clin_surv[is.na(os_days), os_days := lf_days]
  }
  clin_surv[, os_event := as.integer(vital_status == "Dead" | vital_status == "dead")]
  clin_surv <- clin_surv[!is.na(os_days) & os_days >= 0]
  message(sprintf("  Clinical: %d patients with OS data", nrow(clin_surv)))

  # Per-gene KM
  surv_res <- lapply(gold_genes_prom, function(gn) {
    gb <- tumor_meth[gene_name == gn]
    if (nrow(gb) < 20) return(NULL)
    gb <- merge(gb, clin_surv, by = "patient_id")
    gb <- gb[!is.na(os_days) & !is.na(os_event) & is.finite(os_days) & !is.na(mean_beta)]
    if (nrow(gb) < 20) return(NULL)
    med_b <- median(gb$mean_beta, na.rm = TRUE)
    gb[, meth_group := ifelse(mean_beta >= med_b, "High", "Low")]
    fit <- survfit(Surv(os_days, os_event) ~ meth_group, data = gb)
    lr  <- survdiff(Surv(os_days, os_event) ~ meth_group, data = gb)
    lr_p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
    list(gene = gn, fit = fit, data = gb, logrank_p = lr_p, n = nrow(gb))
  })
  surv_res <- Filter(Negate(is.null), surv_res)

  if (length(surv_res) == 0) {
    message("  [WARN] No genes passed survival filters (n≥20 matched samples)")
    surv_summary <- data.table(gene_name = character(), logrank_p = numeric(),
                               n = integer(), logrank_fdr = numeric())
  } else {
    surv_summary <- rbindlist(lapply(surv_res, function(x) {
      data.table(gene_name = x$gene, logrank_p = x$logrank_p, n = x$n)
    }))
    surv_summary[, logrank_fdr := p.adjust(logrank_p, method = "BH")]
    setorder(surv_summary, logrank_p)
    message("  Genes with logrank p < 0.05: ",
            sum(surv_summary$logrank_p < 0.05, na.rm = TRUE))
  }
  fwrite(surv_summary, file.path(OUTDIR, "result", "tcga_survival_summary.csv"))
} # end survival block

# Module 5: Figures ============================================================
message("\n=== Module 5: Figures ===")

save_panel <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(FIGDIR, "panels", paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIGDIR, "png",    paste0(name, ".png")), p, width = w, height = h, dpi = 150)
}

# Panel A: Tumor vs Normal promoter β ==========================================
panel_A <- function() {
  dat <- gene_beta[sample_type2 %in% c("Tumor", "Normal")]
  dat[, gene_label := gene_name]

  # Order genes by Δβ
  ord <- tn_stats[order(delta_beta), gene_name]
  dat[, gene_label := factor(gene_name, levels = intersect(ord, unique(gene_name)))]
  dat[, sample_label := factor(sample_type2, levels = c("Normal", "Tumor"))]

  ggplot(dat, aes(x = sample_label, y = mean_beta, fill = sample_label)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8, width = 0.5) +
    facet_wrap(~ gene_label, nrow = 3, scales = "free_y") +
    scale_fill_manual(values = c(Normal = "#4393C3", Tumor = "#D6604D")) +
    labs(title = "TCGA-LIHC: Promoter methylation (Tumor vs Normal)",
         subtitle = sprintf("n = %d Tumor, %d Normal; 450k; Gold+Silver genes",
                            uniqueN(dat[sample_type2 == "Tumor", barcode]),
                            uniqueN(dat[sample_type2 == "Normal", barcode])),
         x = NULL, y = "Mean promoter β", fill = NULL) +
    theme_hcc +
    theme(strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# Panel B: Forest plot — Spearman ρ ============================================
panel_B <- function() {
  dat <- cor_stats[!is.na(spearman_rho)]
  dat <- merge(dat, tn_stats[, .(gene_name, delta_beta, tcga_dir)],
               by = "gene_name", all.x = TRUE)
  dat[, gene_name := factor(gene_name, levels = dat[order(spearman_rho), gene_name])]

  ggplot(dat, aes(x = spearman_rho, y = gene_name, colour = tcga_dir)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(aes(size = -log10(spearman_p + 1e-10))) +
    geom_errorbarh(aes(xmin = spearman_rho - 0.1, xmax = spearman_rho + 0.1), height = 0) +
    scale_colour_manual(values = c(Hyper = "#D6604D", Hypo = "#4393C3"),
                        na.value = "grey50",
                        name = "TCGA T vs N") +
    scale_size_continuous(name = "-log10(p)", range = c(2, 6)) +
    labs(title = "Methylation–expression correlation (Spearman ρ)",
         subtitle = "Tumor samples; promoter β vs log2(TPM+1)",
         x = "Spearman ρ", y = NULL) +
    theme_hcc +
    theme(legend.position = "right")
}

# Panel C: KM survival curves (top 3 by logrank p) =============================
km_ggplot <- function(sv_data, gene, logrank_p) {
  fit <- survfit(Surv(os_days, os_event) ~ meth_group, data = sv_data)
  km_df <- setDT(broom::tidy(fit))
  km_df[, group := sub("meth_group=", "", strata)]
  n_high <- sum(sv_data$meth_group == "High")
  n_low  <- sum(sv_data$meth_group == "Low")

  ggplot(km_df, aes(x = time / 30.4, y = estimate, colour = group, fill = group)) +
    geom_step(linewidth = 0.8) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    scale_colour_manual(values = c(High = "#D6604D", Low = "#4393C3"),
                        labels = c(High = sprintf("High (n=%d)", n_high),
                                   Low  = sprintf("Low  (n=%d)", n_low))) +
    scale_fill_manual(values = c(High = "#D6604D", Low = "#4393C3"), guide = "none") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("Log-rank p = %s",
                             ifelse(logrank_p < 0.001, "<0.001",
                                    sprintf("%.3f", logrank_p))),
             size = 3.5) +
    labs(title = gene,
         x = "Time (months)", y = "Overall survival", colour = NULL) +
    theme_hcc
}

panel_C_list <- NULL
if (!is.null(surv_res) && length(surv_res) > 0) {
  surv_summary_local <- rbindlist(lapply(surv_res, function(x) {
    data.table(gene_name = x$gene, logrank_p = x$logrank_p, n = x$n,
               idx = which(sapply(surv_res, `[[`, "gene") == x$gene))
  }))
  setorder(surv_summary_local, logrank_p)
  top3 <- head(surv_summary_local, 3)
  panel_C_list <- lapply(seq_len(nrow(top3)), function(i) {
    sr <- surv_res[[top3$idx[i]]]
    km_ggplot(sr$data, sr$gene, sr$logrank_p)
  })
}

# Panel D: Concordance tile — HCC direction vs TCGA Δβ =========================
panel_D <- function() {
  hcc_dir <- load_hcc_direction()
  if (is.null(hcc_dir) || nrow(hcc_dir) == 0) {
    # fallback: use Gold tier with known directionality from fig7/fig5 results
    message("  [WARN] HCC direction not available from gold/silver CSVs; using TCGA-only panel D")
    dat <- tn_stats[!is.na(delta_beta), .(gene_name, delta_beta, tcga_dir,
                                          sig = wilcox_fdr < 0.20)]
    ggplot(dat, aes(x = 1, y = gene_name, fill = delta_beta)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("Δβ=%.3f", delta_beta)), size = 3) +
      scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                           midpoint = 0, name = "TCGA Δβ\n(T-N)") +
      labs(title = "TCGA-LIHC promoter Δβ (Tumor − Normal)",
           subtitle = "Gold+Silver aDMR genes", x = NULL, y = NULL) +
      theme_hcc +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    dat <- merge(tn_stats[, .(gene_name, delta_beta, tcga_dir)],
                 hcc_dir, by = "gene_name", all.x = TRUE)
    dat[, concordant := (tcga_dir == hcc_dir)]
    dat[, gene_name := factor(gene_name, levels = dat[order(delta_beta), gene_name])]
    tile_dat <- melt(dat[, .(gene_name, `TCGA Δβ` = delta_beta,
                              `HCC Δβ` = hcc_delta)],
                     id.vars = "gene_name", variable.name = "Cohort",
                     value.name = "delta_beta")
    ggplot(tile_dat, aes(x = Cohort, y = gene_name, fill = delta_beta)) +
      geom_tile(colour = "white") +
      scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                           midpoint = 0, name = "Δβ") +
      labs(title = "HCC vs TCGA: promoter methylation direction concordance",
           subtitle = "Gold+Silver aDMR genes; left = our HCC, right = TCGA-LIHC",
           x = NULL, y = NULL) +
      theme_hcc
  }
}

# Panel E: HCC allelic Δβ at promoter aDMRs vs TCGA Δβ =========================
panel_E <- function() {
  gold_df <- if (file.exists(GOLD_FINAL)) fread(GOLD_FINAL) else data.table()
  silv_df <- if (file.exists(SILV_FINAL)) fread(SILV_FINAL) else data.table()
  all_df <- rbind(
    gold_df[, .(admr_chr, admr_start, admr_end, sv_minus_wt,
                hp1_beta, hp2_beta, sv_hp, patient_code, tier_class)],
    silv_df[, .(admr_chr, admr_start, admr_end, sv_minus_wt,
                hp1_beta, hp2_beta, sv_hp, patient_code, tier_class)],
    fill = TRUE
  )
  all_df <- all_df[!is.na(admr_chr) & !is.na(sv_minus_wt)]
  if (nrow(all_df) == 0 || !file.exists(PROM_RDS))
    return(ggplot() + labs(title = "Panel E: No allelic data") + theme_hcc)

  admr_gr <- GRanges(seqnames = all_df$admr_chr,
                     ranges   = IRanges(start = all_df$admr_start, end = all_df$admr_end))
  prom_gr   <- readRDS(PROM_RDS)
  gold_prom <- prom_gr[prom_gr$gene_name %in% gold_genes_prom]
  hits <- findOverlaps(gold_prom, admr_gr, ignore.strand = TRUE)
  if (length(hits) == 0)
    return(ggplot() + labs(title = "Panel E: No promoter-aDMR overlaps") + theme_hcc)

  hit_dt <- data.table(
    gene_name    = gold_prom$gene_name[queryHits(hits)],
    sv_minus_wt  = all_df$sv_minus_wt[subjectHits(hits)],
    hp1_beta     = all_df$hp1_beta[subjectHits(hits)],
    hp2_beta     = all_df$hp2_beta[subjectHits(hits)],
    sv_hp        = all_df$sv_hp[subjectHits(hits)],
    patient_code = all_df$patient_code[subjectHits(hits)],
    tier_class   = all_df$tier_class[subjectHits(hits)]
  )
  hit_dt <- unique(hit_dt)
  message(sprintf("  Panel E: %d gene-patient aDMR records across %d genes",
                  nrow(hit_dt), uniqueN(hit_dt$gene_name)))

  # Per-gene summary
  gene_sum <- hit_dt[, .(
    hcc_median = median(sv_minus_wt, na.rm = TRUE),
    hcc_dir    = ifelse(median(sv_minus_wt, na.rm = TRUE) > 0, "Hyper", "Hypo"),
    n_pts      = .N
  ), by = gene_name]

  gene_sum <- merge(gene_sum,
                    tn_stats[, .(gene_name, tcga_delta = delta_beta, tcga_dir, wilcox_fdr)],
                    by = "gene_name", all.x = TRUE)
  gene_sum[, concordant := !is.na(tcga_dir) & (hcc_dir == tcga_dir)]

  n_conc <- sum(gene_sum$concordant, na.rm = TRUE)
  n_tcga <- sum(!is.na(gene_sum$tcga_delta))
  message(sprintf("  Panel E concordance: %d / %d genes (HCC allelic dir == TCGA dir)",
                  n_conc, n_tcga))

  gene_ord <- gene_sum[order(hcc_median), gene_name]
  hit_dt[, gene_name   := factor(gene_name, levels = gene_ord)]
  gene_sum[, gene_name := factor(gene_name, levels = gene_ord)]
  tcga_layer <- gene_sum[!is.na(tcga_delta)]
  conc_layer <- gene_sum[concordant == TRUE]

  ggplot(hit_dt, aes(x = gene_name, y = sv_minus_wt)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_jitter(aes(colour = sv_minus_wt > 0),
                width = 0.18, size = 1.8, alpha = 0.65, stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.42,
                 colour = "grey20", linewidth = 0.7, fatten = 1) +
    geom_point(data = tcga_layer,
               aes(x = gene_name, y = tcga_delta, fill = tcga_dir),
               shape = 23, size = 3.2, colour = "black", stroke = 0.8) +
    geom_text(data = conc_layer,
              aes(x = gene_name, y = Inf, label = "✓"),
              colour = "#2E7D32", size = 4.5, vjust = 1.4) +
    scale_colour_manual(
      values = c("TRUE" = "#D6604D", "FALSE" = "#4393C3"),
      labels = c("TRUE" = "Hyper (SV HP > WT)", "FALSE" = "Hypo (SV HP < WT)"),
      name   = "HCC allelic Δβ"
    ) +
    scale_fill_manual(
      values   = c("Hyper" = "#D6604D", "Hypo" = "#4393C3"),
      na.value = "grey70",
      name     = "TCGA direction (◇)"
    ) +
    labs(
      title    = "Allelic Δβ at promoter aDMRs vs TCGA Δβ",
      subtitle = sprintf(
        "Points = per-patient (SV-HP − WT-HP); crossbar = median; ◇ = TCGA T−N Δβ; ✓ = concordant (%d/%d genes)",
        n_conc, n_tcga),
      x = NULL,
      y = "Allelic Δβ (SV haplotype − WT)"
    ) +
    theme_hcc +
    theme(axis.text.x  = element_text(angle = 40, hjust = 1),
          legend.position = "bottom",
          legend.box      = "horizontal")
}

# Assemble figure ==============================================================
message("  Assembling panels ...")
library(patchwork)

pA <- tryCatch(panel_A(), error = function(e) {
  message("  [WARN] Panel A error: ", conditionMessage(e))
  ggplot() + labs(title = paste("Panel A error:", conditionMessage(e))) + theme_hcc
})

pB <- tryCatch(panel_B(), error = function(e) {
  message("  [WARN] Panel B error: ", conditionMessage(e))
  ggplot() + labs(title = paste("Panel B error:", conditionMessage(e))) + theme_hcc
})

pD <- tryCatch(panel_D(), error = function(e) {
  message("  [WARN] Panel D error: ", conditionMessage(e))
  ggplot() + labs(title = paste("Panel D error:", conditionMessage(e))) + theme_hcc
})

pE <- tryCatch(panel_E(), error = function(e) {
  message("  [WARN] Panel E error: ", conditionMessage(e))
  ggplot() + labs(title = paste("Panel E error:", conditionMessage(e))) + theme_hcc
})

save_panel(pA, "fig8A_tn_methylation", w = 14, h = 8)
save_panel(pB, "fig8B_meth_expr_cor", w = 7, h = 6)
save_panel(pD, "fig8D_concordance", w = 6, h = 6)
save_panel(pE, "fig8E_allelic_delta", w = 10, h = 6)

# Panel C: KM curves (up to 3)
if (!is.null(panel_C_list) && length(panel_C_list) > 0) {
  pC_combined <- wrap_plots(panel_C_list, nrow = 1)
  save_panel(pC_combined, "fig8C_survival_km",
             w = 4 * length(panel_C_list), h = 5)
}

# Combined figure: A+B / C+D / E
if (!is.null(panel_C_list) && length(panel_C_list) > 0) {
  pC_combined <- wrap_plots(panel_C_list, nrow = 1)
  fig8 <- (pA | pB) / (pC_combined | pD) / pE +
    plot_annotation(
      title = "Fig 8. TCGA-LIHC validation of Gold-tier SV-associated methylation",
      tag_levels = "A"
    ) +
    plot_layout(heights = c(2, 2, 1.5))
} else {
  fig8 <- (pA | pB) / (pD | pE) +
    plot_annotation(
      title = "Fig 8. TCGA-LIHC validation of Gold-tier SV-associated methylation",
      tag_levels = "A"
    ) +
    plot_layout(heights = c(2, 1.5))
}

ggsave(file.path(FIGDIR, "png", "fig8_combined.png"),
       fig8, width = 18, height = 18, dpi = 150)
ggsave(file.path(FIGDIR, "panels", "fig8_combined.pdf"),
       fig8, width = 18, height = 18)
message("  Saved: fig8_combined.png/pdf")

# Summary ======================================================================
message("\n=== Summary ===")
message(sprintf("Genes validated      : %d", uniqueN(gene_beta$gene_name)))
message(sprintf("Tumor samples (meth) : %d",
                uniqueN(gene_beta[sample_type2 == "Tumor", barcode])))
message(sprintf("Normal samples (meth): %d",
                uniqueN(gene_beta[sample_type2 == "Normal", barcode])))
message(sprintf("Genes |Δβ| ≥ 0.05   : %d", sum(abs(tn_stats$delta_beta) >= 0.05, na.rm = TRUE)))
message(sprintf("Genes with ρ < -0.20 : %d",
                sum(cor_stats$spearman_rho < -0.20, na.rm = TRUE)))
n_surv_sig <- if (!is.null(surv_res)) {
  sum(sapply(surv_res, `[[`, "logrank_p") < 0.05)
} else 0L
message(sprintf("Genes logrank p<0.05 : %d", n_surv_sig))

message("\nDone.")
