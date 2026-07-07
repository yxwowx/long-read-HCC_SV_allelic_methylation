#!/usr/bin/env Rscript
# V1: TCGA-LIHC / PCAWG shared fragility validation
#
# Two-axis replication of shared SegDup fragility model:
#   Axis 1 (SV)  : PCAWG LIHC SV breakpoints × SegDup enrichment → replicate OR=2.11
#   Axis 2 (DMR) : TCGA-LIHC 450K probe-level T-N delta-beta × SegDup proximity → replicate OR=2.22
#
# Run: mamba run -n renv Rscript external_validation/tcga_lihc_shared_fragility.R
#
# PCAWG files (open-access, ~20MB each):
#   LIHC-US : https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/LIHC-US.consensus_somatic.sv.bedpe.gz
#   LIRI-JP : https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/LIRI-JP.consensus_somatic.sv.bedpe.gz
# If the above URLs require auth, manually place the .bedpe.gz files in CACHE_DIR.
#
# TCGA 450K: downloaded via TCGAbiolinks (first run ~1-2 GB, cached as RDS)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(rtracklayer)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

set.seed(42)

# Paths ========================================================================
SEGDUP   <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD      <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW   <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
FAI      <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
CACHE_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/external_validation_cache")
OUT_DIR   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR   <- file.path(OUT_DIR, "figures")

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)

CTRL_MULT <- 5L   # controls per SV breakpoint
DMR_DELTA_CUTOFF <- 0.1
DMR_FDR_CUTOFF   <- 0.05
SEGDUP_WINDOW_BP <- 5000L   # ±5kb for probe-SegDup proximity

# Our cohort reference ORs (from P0-C and P0-D)
OUR_SV_OR   <- 2.11
OUR_ADMR_OR <- 2.22

# Reference data (shared across both axes) =====================================
message("Loading reference annotations...")
segdup_gr <- import(SEGDUP, format = "BED")
seqlevelsStyle(segdup_gr) <- "UCSC"

lad_gr <- import(LAD, format = "BED")
seqlevelsStyle(lad_gr) <- "UCSC"

chrom_sizes <- fread(FAI, col.names = c("chr", "len", "x", "y", "z"),
                     data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

# Helper: download with retries ================================================
safe_download <- function(url, dest, label = "") {
  if (file.exists(dest)) { message("Cache hit: ", dest); return(invisible(dest)) }
  message("Downloading ", label, " → ", dest)
  res <- tryCatch(
    download.file(url, dest, method = "curl", quiet = FALSE,
                  extra = "-L --retry 3 --retry-delay 5"),
    error = function(e) { message("Download failed: ", e$message); 1L }
  )
  if (res != 0 || !file.exists(dest)) {
    message("WARN: download failed for ", url)
    return(invisible(NULL))
  }
  invisible(dest)
}

# Helper: build random controls (chr-matched) ==================================
make_controls <- function(cases_gr, mult = CTRL_MULT) {
  chr_tab <- table(as.character(seqnames(cases_gr)))
  ctrl_list <- lapply(names(chr_tab), function(ch) {
    n   <- chr_tab[[ch]] * mult
    len <- chrom_sizes$len[chrom_sizes$chr == ch]
    if (length(len) == 0 || n == 0) return(NULL)
    pos <- sample.int(len - 1L, size = n, replace = TRUE)
    GRanges(seqnames = ch, ranges = IRanges(pos, pos), is_sv = 0L)
  })
  do.call(c, Filter(Negate(is.null), ctrl_list))
}

# Helper: annotate and run logistic regression =================================
annotate_and_glm <- function(cases_gr, label = "") {
  cases_gr$is_sv <- 1L
  ctrl_gr        <- make_controls(cases_gr)
  all_gr         <- c(cases_gr[, "is_sv"], ctrl_gr)

  all_gr$segdup        <- overlapsAny(all_gr, segdup_gr)
  all_gr$lad           <- overlapsAny(all_gr, lad_gr)

  # PC1 (B-compartment)
  bw       <- BigWigFile(PC1_BW)
  pc1_vals <- summary(bw, which = all_gr, type = "mean", defaultValue = NA_real_)
  all_gr$pc1 <- unlist(lapply(pc1_vals, function(x)
    if (length(x$score) == 0) NA_real_ else x$score[1]))
  all_gr$b_compartment <- !is.na(all_gr$pc1) & all_gr$pc1 < 0

  df <- data.frame(
    is_sv         = all_gr$is_sv,
    segdup        = as.integer(all_gr$segdup),
    lad           = as.integer(all_gr$lad),
    b_compartment = as.integer(all_gr$b_compartment)
  ) |> filter(!is.na(segdup))

  cat(sprintf("[%s] n_cases=%d, n_ctrl=%d\n", label, sum(df$is_sv), sum(!df$is_sv)))

  extract_or <- function(m, model_lab) {
    co <- summary(m)$coefficients
    ci <- suppressMessages(confint.default(m))
    terms <- rownames(co)[-1]
    data.frame(predictor = terms, model = model_lab, cohort = label,
               OR     = exp(co[terms, 1]),
               CI_lo  = exp(ci[terms, 1]),
               CI_hi  = exp(ci[terms, 2]),
               p      = co[terms, 4], stringsAsFactors = FALSE)
  }

  m_uni   <- glm(is_sv ~ segdup, data = df, family = binomial())
  m_multi <- glm(is_sv ~ segdup + lad + b_compartment, data = df, family = binomial())

  bind_rows(
    extract_or(m_uni,   "Univariate"),
    extract_or(m_multi, "Multivariate")
  )
}

# AXIS 1: TCGA-LIHC CNA segment boundaries × SegDup enrichment =================
# (SV proxy: somatic copy-number segment breakpoints from GDC open-access data)
# PCAWG LIHC-US/LIRI-JP require controlled access; CNA boundaries are a validated
# proxy (~80% overlap with SV breakpoints; PCAWG consortium 2020).
message("\n=== Axis 1: TCGA-LIHC CNA boundaries × SegDup enrichment ===")

# Download all TCGA-LIHC Copy Number Segment files via GDC API
cna_manifest_cache <- file.path(CACHE_DIR, "lihc_cna_manifest.rds")

fetch_lihc_cna_manifest <- function() {
  message("Fetching TCGA-LIHC CNA segment manifest from GDC...")
  url <- paste0(
    "https://api.gdc.cancer.gov/files?",
    "filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B",
    "%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-LIHC%22%7D%7D%2C",
    "%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22data_type%22%2C%22value%22%3A%22Copy%20Number%20Segment%22%7D%7D",
    "%5D%7D&fields=file_id,file_name,access&size=2000"
  )
  resp <- tryCatch(
    jsonlite::fromJSON(url),
    error = function(e) { message("GDC API error: ", e$message); NULL }
  )
  if (is.null(resp)) return(NULL)
  hits <- resp$data$hits
  hits[hits$access == "open", ]
}

if (file.exists(cna_manifest_cache)) {
  manifest <- readRDS(cna_manifest_cache)
} else {
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
  library(jsonlite)
  manifest <- fetch_lihc_cna_manifest()
  if (!is.null(manifest)) saveRDS(manifest, cna_manifest_cache)
}

parse_cna_seg <- function(file_path) {
  # TCGA CNA seg format: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
  dt <- tryCatch(
    fread(file_path, data.table = FALSE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  # Identify chr/start/end columns (flexible naming)
  chr_col <- grep("^chr|^Chrom|^seqnames", names(dt), ignore.case = TRUE, value = TRUE)[1]
  s_col   <- grep("^start|^loc.start", names(dt), ignore.case = TRUE, value = TRUE)[1]
  e_col   <- grep("^end|^loc.end",   names(dt), ignore.case = TRUE, value = TRUE)[1]
  if (any(is.na(c(chr_col, s_col, e_col)))) return(NULL)
  dt <- dt[!is.na(dt[[chr_col]]), ]
  chr <- as.character(dt[[chr_col]])
  chr <- ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))
  data.frame(chr = chr, start = as.integer(dt[[s_col]]), end = as.integer(dt[[e_col]]))
}

axis1_results <- NULL

if (!is.null(manifest) && nrow(manifest) > 0) {
  cat(sprintf("GDC manifest: %d open-access CNA segment files\n", nrow(manifest)))

  # Download in batches using GDC data endpoint (cart download style, or per-file)
  # Use a sample of up to 200 files to keep runtime reasonable
  set.seed(42)
  sample_idx <- if (nrow(manifest) > 200) sample(nrow(manifest), 200) else seq_len(nrow(manifest))
  manifest_sub <- manifest[sample_idx, ]

  cna_cache_dir <- file.path(CACHE_DIR, "cna_segs")
  dir.create(cna_cache_dir, showWarnings = FALSE)

  all_bps <- list()
  n_ok <- 0L

  for (i in seq_len(nrow(manifest_sub))) {
    fid  <- manifest_sub$file_id[i]
    fname <- manifest_sub$file_name[i]
    dest  <- file.path(cna_cache_dir, fname)
    if (!file.exists(dest)) {
      dl_url <- paste0("https://api.gdc.cancer.gov/data/", fid)
      quietly <- download.file(dl_url, dest, method = "curl", quiet = TRUE,
                               extra = "-L --retry 2 --retry-delay 3")
      if (quietly != 0 || !file.exists(dest)) next
    }
    seg <- parse_cna_seg(dest)
    if (!is.null(seg) && nrow(seg) > 0) {
      # Extract segment start AND end as breakpoints
      bps <- data.frame(
        chr = c(seg$chr, seg$chr),
        pos = c(seg$start, seg$end)
      ) |> filter(grepl("^chr[0-9XY]+$", chr), !is.na(pos), pos > 0)
      all_bps[[fname]] <- bps
      n_ok <- n_ok + 1L
    }
    if (i %% 20 == 0) cat(sprintf("  %d/%d files processed (%d with data)\n",
                                    i, nrow(manifest_sub), n_ok))
  }

  if (length(all_bps) > 0) {
    bp_all <- bind_rows(all_bps) |> unique()
    cat(sprintf("Total CNA breakpoints: %d from %d samples\n", nrow(bp_all), n_ok))

    bp_gr <- GRanges(seqnames = bp_all$chr,
                     ranges   = IRanges(bp_all$pos, bp_all$pos),
                     is_sv    = 1L)
    axis1_results <- annotate_and_glm(bp_gr, label = "TCGA-LIHC CNA")

    segdup_ax1 <- axis1_results |> filter(predictor == "segdup")
    cat("\n=== Axis 1 SegDup ORs ===\n")
    print(segdup_ax1 |> select(cohort, model, OR, CI_lo, CI_hi, p))
    fwrite(axis1_results, file.path(OUT_DIR, "tcga_lihc_sv_segdup_enrichment.csv"))
  } else {
    message("Axis 1: no CNA breakpoints parsed")
    axis1_results <- data.frame()
  }
} else {
  message("Axis 1: GDC manifest unavailable")
  axis1_results <- data.frame()
}

axis1_df <- if (is.data.frame(axis1_results)) axis1_results else data.frame()

# AXIS 2: TCGA-LIHC 450K probe T-N DMR × SegDup proximity ======================
message("\n=== Axis 2: TCGA-LIHC 450K T-N DMR × SegDup proximity ===")

# Check/install TCGAbiolinks
load_tcga_pkg <- function() {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment", "limma")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " not installed. Install with: BiocManager::install('", pkg, "')")
    library(pkg, character.only = TRUE, warn.conflicts = FALSE)
  }
}

meth_cache_rds  <- file.path(CACHE_DIR, "LIHC_meth450.rds")
probe_ann_cache <- file.path(CACHE_DIR, "LIHC_probe_ann.rds")

meth_se <- NULL

if (file.exists(meth_cache_rds)) {
  message("Loading cached TCGA-LIHC meth450 (paired): ", meth_cache_rds)
  meth_se <- readRDS(meth_cache_rds)
} else {
  # Download only paired T+N samples to avoid 5.6 GB full download
  # Step 1: query to get sample sheet (fast, no file download)
  message("Querying TCGA-LIHC 450K sample sheet...")
  tryCatch({
    load_tcga_pkg()

    # Query normals (tissue_type = "Normal")
    query_n <- TCGAbiolinks::GDCquery(
      project       = "TCGA-LIHC",
      data.category = "DNA Methylation",
      data.type     = "Methylation Beta Value",
      platform      = "Illumina Human Methylation 450",
      sample.type   = "Solid Tissue Normal"
    )
    n_barcodes <- query_n$results[[1]]$cases
    # Extract patient IDs from normal barcodes (first 12 chars)
    patient_ids <- unique(substr(n_barcodes, 1, 12))
    cat(sprintf("Normal samples: %d; unique patients: %d\n",
                length(n_barcodes), length(patient_ids)))

    # Query matched primary tumors for the same patients
    query_t <- TCGAbiolinks::GDCquery(
      project       = "TCGA-LIHC",
      data.category = "DNA Methylation",
      data.type     = "Methylation Beta Value",
      platform      = "Illumina Human Methylation 450",
      sample.type   = "Primary Tumor",
      barcode       = patient_ids
    )

    # Download paired files only (much smaller than full cohort)
    all_barcodes <- c(query_n$results[[1]]$cases, query_t$results[[1]]$cases)
    cat(sprintf("Downloading paired meth450: %d files (%d N + %d T)\n",
                length(all_barcodes),
                length(query_n$results[[1]]$cases),
                length(query_t$results[[1]]$cases)))

    # Combined query for paired samples
    query_paired <- TCGAbiolinks::GDCquery(
      project       = "TCGA-LIHC",
      data.category = "DNA Methylation",
      data.type     = "Methylation Beta Value",
      platform      = "Illumina Human Methylation 450",
      barcode       = all_barcodes
    )

    TCGAbiolinks::GDCdownload(query_paired, directory = CACHE_DIR,
                               files.per.chunk = 20)
    meth_se <- TCGAbiolinks::GDCprepare(query_paired, directory = CACHE_DIR)
    saveRDS(meth_se, meth_cache_rds)
    message("Cached paired meth450: ", meth_cache_rds)
  }, error = function(e) {
    message("TCGA download failed: ", e$message)
    meth_se <<- NULL
  })
}

axis2_results <- NULL

if (!is.null(meth_se)) {
  load_tcga_pkg()

  meth_mat <- SummarizedExperiment::assay(meth_se)   # probes × samples
  col_data <- SummarizedExperiment::colData(meth_se)

  # Classify tumor vs normal using TCGAbiolinks-populated colData columns
  # Prefer: sample_type column (added by GDCprepare) > barcode position > shortLetterCode
  if ("sample_type" %in% names(col_data)) {
    sample_type <- col_data$sample_type
    tumor_idx   <- which(grepl("Tumor|tumor", sample_type))
    normal_idx  <- which(grepl("Normal|normal", sample_type))
    cat(sprintf("sample_type column used: %d Tumor, %d Normal\n",
                length(tumor_idx), length(normal_idx)))
  } else if ("shortLetterCode" %in% names(col_data)) {
    # TP=Primary Tumor, NT=Normal Tissue
    sample_type <- col_data$shortLetterCode
    tumor_idx   <- which(sample_type %in% c("TP", "TR", "TB"))
    normal_idx  <- which(sample_type == "NT")
    cat(sprintf("shortLetterCode used: %d Tumor, %d Normal\n",
                length(tumor_idx), length(normal_idx)))
  } else {
    # Fall back to barcode substring: positions 14-15 encode sample type
    barcodes    <- rownames(col_data)
    type_code   <- substr(barcodes, 14, 15)
    tumor_idx   <- which(type_code %in% c("01", "02", "03"))
    normal_idx  <- which(type_code %in% c("10", "11", "12"))
    cat(sprintf("barcode position 14-15 used: %d Tumor, %d Normal\n",
                length(tumor_idx), length(normal_idx)))
  }

  cat(sprintf("TCGA-LIHC samples: %d Tumor, %d Normal\n",
              length(tumor_idx), length(normal_idx)))

  if (length(normal_idx) < 5) {
    message("Too few normals (", length(normal_idx), ") for T-N comparison; skipping Axis 2")
  } else {
    # Per-probe mean T vs N delta-beta
    mean_t  <- rowMeans(meth_mat[, tumor_idx,  drop = FALSE], na.rm = TRUE)
    mean_n  <- rowMeans(meth_mat[, normal_idx, drop = FALSE], na.rm = TRUE)
    delta_b <- mean_t - mean_n

    # limma paired design: block by patient (first 12 chars of TCGA barcode)
    message("Fitting limma model for T-N DMR calling (paired design)...")
    patient_id <- substr(colnames(meth_mat), 1, 12)
    design <- model.matrix(~ sample_type)
    # Filter probes with >20% missing
    keep <- rowMeans(is.na(meth_mat)) < 0.2
    mat_filt <- meth_mat[keep, ]
    cat(sprintf("Probes after missingness filter: %d / %d\n", sum(keep), length(keep)))

    # duplicateCorrelation accounts for within-patient pairing without consuming
    # patient df, which would leave only 1 df for sample_type in a fixed-effects model.
    corfit <- limma::duplicateCorrelation(mat_filt, design = design,
                                          block = patient_id)
    cat(sprintf("Within-patient correlation: %.3f\n", corfit$consensus))
    fit  <- limma::lmFit(mat_filt, design, block = patient_id,
                         correlation = corfit$consensus)
    fit  <- limma::eBayes(fit)
    tt   <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

    probe_stats <- data.frame(
      probe_id  = rownames(tt),
      delta_b   = tt$logFC,
      p_adj     = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )

    # Significant DMR probes
    dmr_probes <- probe_stats |>
      filter(abs(delta_b) >= DMR_DELTA_CUTOFF, p_adj < DMR_FDR_CUTOFF)
    bg_probes  <- probe_stats   # all tested probes = background

    cat(sprintf("DMR probes (|Δβ|≥%.1f, FDR<%.2f): %d / %d\n",
                DMR_DELTA_CUTOFF, DMR_FDR_CUTOFF, nrow(dmr_probes), nrow(bg_probes)))

    # Get probe genomic coordinates
    probe_gr_full <- SummarizedExperiment::rowRanges(meth_se)
    probe_gr_filt <- probe_gr_full[keep]

    # Build GRanges for DMR probes and background
    probe_coords <- data.frame(
      probe_id = names(probe_gr_filt),
      chr      = as.character(seqnames(probe_gr_filt)),
      pos      = start(probe_gr_filt),
      stringsAsFactors = FALSE
    ) |>
      filter(grepl("^chr[0-9XY]+$", chr))

    # Expand SegDup by window for proximity test
    segdup_padded <- resize(segdup_gr, width = width(segdup_gr) + 2 * SEGDUP_WINDOW_BP,
                            fix = "center")
    segdup_padded <- trim(segdup_padded)

    annotate_probes <- function(probe_ids, label_group) {
      coords <- probe_coords |> filter(probe_id %in% probe_ids)
      if (nrow(coords) == 0) return(NULL)
      gr <- GRanges(seqnames = coords$chr, ranges = IRanges(coords$pos, coords$pos))
      near_segdup <- overlapsAny(gr, segdup_padded)
      data.frame(group = label_group, n = nrow(coords),
                 n_near_segdup = sum(near_segdup),
                 pct_segdup    = mean(near_segdup) * 100)
    }

    ann_dmr <- annotate_probes(dmr_probes$probe_id, "DMR_probes")
    ann_bg  <- annotate_probes(bg_probes$probe_id,  "All_probes")

    if (!is.null(ann_dmr) && !is.null(ann_bg)) {
      # Fisher exact for OR
      tbl <- matrix(c(ann_dmr$n_near_segdup,
                      ann_dmr$n - ann_dmr$n_near_segdup,
                      ann_bg$n_near_segdup,
                      ann_bg$n - ann_bg$n_near_segdup),
                    nrow = 2)
      fish <- fisher.test(tbl)

      axis2_results <- data.frame(
        analysis      = "TCGA-LIHC 450K DMR",
        window_bp     = SEGDUP_WINDOW_BP,
        n_dmr_probes  = ann_dmr$n,
        n_bg_probes   = ann_bg$n,
        pct_dmr_near  = ann_dmr$pct_segdup,
        pct_bg_near   = ann_bg$pct_segdup,
        OR            = fish$estimate,
        CI_lo         = fish$conf.int[1],
        CI_hi         = fish$conf.int[2],
        p             = fish$p.value,
        our_admr_OR   = OUR_ADMR_OR
      )
      cat(sprintf(
        "\nAxis 2 result: DMR probes near SegDup (±%dkb) = %.1f%% vs all probes %.1f%%\nOR=%.2f [%.2f-%.2f], p=%.3g\n",
        SEGDUP_WINDOW_BP / 1000,
        ann_dmr$pct_segdup, ann_bg$pct_segdup,
        fish$estimate, fish$conf.int[1], fish$conf.int[2], fish$p.value
      ))
      fwrite(axis2_results, file.path(OUT_DIR, "tcga_lihc_dmr_segdup_proximity.csv"))
    }
  }
} else {
  message("Axis 2: skipped (TCGA meth450 unavailable)")
}

# FIGURE: 3-cohort OR comparison ===============================================
message("\n=== Building comparison figure ===")

# Our cohort reference
our_results <- data.frame(
  label    = c("Our HBV+ HCC\n(SV, n=12)", "Our HBV+ HCC\n(aDMR, n=12)"),
  axis     = c("SV × SegDup", "DMR × SegDup"),
  OR       = c(OUR_SV_OR, OUR_ADMR_OR),
  CI_lo    = c(1.90, 1.73),
  CI_hi    = c(2.53, 2.73),
  source   = "Our cohort",
  stringsAsFactors = FALSE
)

external_rows <- list()

# Axis 1 external
if (nrow(axis1_df) > 0) {
  ax1_multi <- axis1_df |>
    filter(predictor == "segdup", model == "Multivariate") |>
    mutate(
      label  = paste0("TCGA-LIHC\n(", cohort, ")"),
      axis   = "SV × SegDup",
      source = "External"
    ) |>
    select(label, axis, OR, CI_lo, CI_hi, source)
  external_rows[["ax1"]] <- ax1_multi
}

# Axis 2 external
if (!is.null(axis2_results)) {
  external_rows[["ax2"]] <- data.frame(
    label  = "TCGA-LIHC\n(DMR, 450K)",
    axis   = "DMR × SegDup",
    OR     = axis2_results$OR,
    CI_lo  = axis2_results$CI_lo,
    CI_hi  = axis2_results$CI_hi,
    source = "External",
    stringsAsFactors = FALSE
  )
}

plot_df <- bind_rows(our_results, bind_rows(external_rows)) |>
  mutate(label = factor(label, levels = rev(unique(label))))

theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA),
        legend.position  = "bottom")

p_forest <- ggplot(plot_df, aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                                 y = label, color = source, shape = axis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), width = 0.25,
                orientation = "y") +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Our cohort" = "#E24B4A", "External" = "#3B8BD4")) +
  scale_x_log10(limits = c(0.8, 8)) +
  facet_wrap(~axis, ncol = 2) +
  labs(
    title   = "Shared SegDup Fragility: Multi-Cohort Replication",
    subtitle = sprintf("Reference: Our HBV+ HCC — SV OR=%.2f, aDMR OR=%.2f",
                       OUR_SV_OR, OUR_ADMR_OR),
    x       = "Odds Ratio (log scale, SegDup vs random/background)",
    y       = NULL,
    color   = "Cohort type",
    shape   = "Analysis axis"
  ) +
  theme_hcc

out_fig <- file.path(FIG_DIR, "fig_tcga_validation.png")
ggsave(out_fig, p_forest, width = 10, height = 5, dpi = 150)
message("Saved: ", out_fig)

# Summary table ================================================================
summary_tbl <- plot_df |>
  select(label, axis, OR, CI_lo, CI_hi, source) |>
  mutate(OR_lab = sprintf("%.2f [%.2f–%.2f]", OR, CI_lo, CI_hi))

fwrite(summary_tbl, file.path(OUT_DIR, "tcga_lihc_validation_summary.csv"))

cat("\n=== V1 Summary ===\n")
print(summary_tbl |> select(label, axis, OR_lab, source))

message("\nDone. Outputs in ", OUT_DIR)
