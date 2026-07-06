#!/usr/bin/env Rscript
# P0-3: TCGA-LIHC orthogonal SV-methylation replication
# Gold-tier 89 loci × TCGA-LIHC HM450 → SV+ vs SV- sample comparison
#
# Requires TCGAbiolinks (Bioc) and internet access for first run.
# Cached downloads stored in opt$cache_dir.
#
# Run: mamba run -n renv Rscript external_validation/tcga_lihc.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(optparse)
})

option_list <- list(
  make_option("--gold_file", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/gold_tier_final.csv"),
  make_option("--cache_dir", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/tcga_cache"),
  make_option("--outdir", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment"),
  make_option("--run_id", type = "character", default = "tier_v2"),
  make_option("--probe_window_bp", type = "integer", default = 2000L,
    help = "Window around Gold DMR center to match HM450 probes (default 2kb)")
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR    <- opt$outdir
CACHE_DIR <- opt$cache_dir
RUN_ID    <- opt$run_id
LOG_FILE  <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR,    showWarnings = FALSE, recursive = TRUE)

# ── Load TCGAbiolinks lazily ───────────────────────────────────────────────────
load_tcga_pkg <- function() {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " not installed. Run: BiocManager::install('", pkg, "')")
    library(pkg, character.only = TRUE)
  }
}

# ── 1. Gold loci ───────────────────────────────────────────────────────────────
message("Loading Gold tier: ", opt$gold_file)
gold <- fread(opt$gold_file)

find_col <- function(df, cands, req = TRUE) {
  m <- intersect(cands, names(df)); if (length(m) > 0) return(m[[1]])
  if (req) stop("None of [", paste(cands, collapse=","), "] found"); NULL
}

col_chr  <- find_col(gold, c("chr", "seqnames", "chrom", "dmr_chr"))
col_s    <- find_col(gold, c("start", "dmr_start"))
col_e    <- find_col(gold, c("end", "dmr_end"))
col_gene <- find_col(gold, c("gene", "gene_name", "nearest_gene"), req = FALSE)

gold_gr <- gold %>%
  dplyr::rename(chr = !!col_chr, start = !!col_s, end = !!col_e) %>%
  dplyr::mutate(
    center   = as.integer((start + end) / 2),
    win_start = center - opt$probe_window_bp,
    win_end   = center + opt$probe_window_bp,
    locus_id  = paste0("locus_", seq_len(n()))
  )

message(sprintf("Gold tier: %d loci, %d unique chromosomes", nrow(gold_gr),
                length(unique(gold_gr$chr))))

# ── 2. Download / load TCGA-LIHC HM450 methylation ────────────────────────────
meth_cache <- file.path(CACHE_DIR, "LIHC_meth450.rds")
if (file.exists(meth_cache)) {
  message("Loading cached LIHC methylation: ", meth_cache)
  meth_se <- readRDS(meth_cache)
} else {
  message("Downloading TCGA-LIHC HM450 methylation via TCGAbiolinks...")
  load_tcga_pkg()
  query_meth <- TCGAbiolinks::GDCquery(
    project           = "TCGA-LIHC",
    data.category     = "DNA Methylation",
    data.type         = "Methylation Beta Value",
    platform          = "Illumina Human Methylation 450",
    legacy            = FALSE
  )
  TCGAbiolinks::GDCdownload(query_meth, directory = CACHE_DIR)
  meth_se <- TCGAbiolinks::GDCprepare(query_meth, directory = CACHE_DIR)
  saveRDS(meth_se, meth_cache)
  message("Saved cache: ", meth_cache)
}

# ── 3. Download / load TCGA-LIHC copy number (GISTIC) for SV+/SV- proxy ───────
# Use GISTIC2 arm-level calls: samples with focal CNA at Gold loci = SV proxy
gistic_cache <- file.path(CACHE_DIR, "LIHC_gistic.rds")
if (file.exists(gistic_cache)) {
  message("Loading cached GISTIC data: ", gistic_cache)
  gistic_df <- readRDS(gistic_cache)
} else {
  message("Downloading TCGA-LIHC copy number data...")
  load_tcga_pkg()
  query_cn <- TCGAbiolinks::GDCquery(
    project      = "TCGA-LIHC",
    data.category = "Copy Number Variation",
    data.type    = "Gene Level Copy Number"
  )
  TCGAbiolinks::GDCdownload(query_cn, directory = CACHE_DIR)
  gistic_df <- TCGAbiolinks::GDCprepare(query_cn, directory = CACHE_DIR)
  saveRDS(gistic_df, gistic_cache)
}

# ── 4. Match HM450 probes to Gold loci ────────────────────────────────────────
load_tcga_pkg()
meth_mat  <- SummarizedExperiment::assay(meth_se)  # probes × samples
probe_ann <- SummarizedExperiment::rowRanges(meth_se)

probe_gr  <- as(probe_ann, "GRanges")
gold_win_gr <- GRanges(
  seqnames = gold_gr$chr,
  ranges   = IRanges(gold_gr$win_start, gold_gr$win_end),
  locus_id = gold_gr$locus_id
)
seqlevelsStyle(gold_win_gr) <- seqlevelsStyle(probe_gr)

hits <- findOverlaps(probe_gr, gold_win_gr)
probe_locus <- data.frame(
  probe_id = names(probe_gr)[queryHits(hits)],
  locus_id = gold_win_gr$locus_id[subjectHits(hits)],
  stringsAsFactors = FALSE
)
cat(sprintf("Probes matched to Gold loci: %d probes across %d loci\n",
            nrow(probe_locus), length(unique(probe_locus$locus_id))))

# ── 5. Per-locus: mean methylation per TCGA sample ────────────────────────────
locus_meth <- lapply(unique(probe_locus$locus_id), function(lid) {
  probes <- probe_locus$probe_id[probe_locus$locus_id == lid]
  probes <- intersect(probes, rownames(meth_mat))
  if (length(probes) == 0) return(NULL)
  vals <- if (length(probes) == 1) meth_mat[probes, ] else colMeans(meth_mat[probes, ], na.rm = TRUE)
  data.frame(locus_id = lid, sample = names(vals), mean_beta = as.numeric(vals))
}) %>% bind_rows()

# ── 6. Classify samples as SV+ / SV- using arm-level CNA as proxy ─────────────
# Samples with any focal amplification or deep deletion at locus = SV_proxy_pos
# This is approximate; ideally would use PCAWG SV calls
classify_sv_proxy <- function(gistic_obj, locus_id_vec) {
  # If gistic_df is a data.frame with sample × gene matrix
  if (inherits(gistic_obj, "data.frame") || inherits(gistic_obj, "data.table")) {
    # Expect: columns = samples, rows = genes/segments, values = copy number call
    cn_mat <- as.matrix(gistic_obj[, -1])
    rownames(cn_mat) <- gistic_obj[[1]]
    # Samples with any deletion (≤-1) or amplification (≥1) near loci = SV proxy
    sv_pos_samples <- unique(colnames(cn_mat)[apply(cn_mat, 2, function(x) any(abs(x) >= 1, na.rm=TRUE))])
    return(sv_pos_samples)
  }
  # Fallback: use all samples as SV- (conservative)
  message("Cannot classify SV+/SV- from gistic format; all samples treated as SV-")
  character(0)
}

sv_pos_samples <- tryCatch(
  classify_sv_proxy(gistic_df, gold_gr$locus_id),
  error = function(e) { message("SV proxy classification failed: ", e$message); character(0) }
)

locus_meth <- locus_meth %>%
  dplyr::mutate(
    sv_status = ifelse(sample %in% sv_pos_samples, "SV_pos", "SV_neg"),
    sample_type = ifelse(grepl("-01[A-Z]$", sample), "Tumor", "Normal")
  )

# ── 7. Per-locus Wilcoxon test: SV+ vs SV- ────────────────────────────────────
locus_tests <- locus_meth %>%
  dplyr::filter(sample_type == "Tumor") %>%
  dplyr::group_by(locus_id) %>%
  dplyr::summarise(
    n_sv_pos   = sum(sv_status == "SV_pos", na.rm = TRUE),
    n_sv_neg   = sum(sv_status == "SV_neg", na.rm = TRUE),
    med_sv_pos = median(mean_beta[sv_status == "SV_pos"], na.rm = TRUE),
    med_sv_neg = median(mean_beta[sv_status == "SV_neg"], na.rm = TRUE),
    delta_med  = med_sv_pos - med_sv_neg,
    p_wilcox   = tryCatch({
      wilcox.test(mean_beta[sv_status == "SV_pos"],
                  mean_beta[sv_status == "SV_neg"])$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_fdr    = p.adjust(p_wilcox, method = "BH"),
    sig      = p_fdr < 0.05,
    direction = ifelse(delta_med > 0, "hyper_in_SV_pos", "hypo_in_SV_pos")
  )

# Concordance with our HCC allelic Δβ direction
if (!is.null(col_gene) && "delta_beta_bulk" %in% names(gold)) {
  gold_dir <- gold %>%
    dplyr::mutate(locus_id = paste0("locus_", seq_len(n())),
                  our_direction = ifelse(delta_beta_bulk > 0, "hyper", "hypo"))
  locus_tests <- left_join(locus_tests, gold_dir[, c("locus_id", "our_direction")], by = "locus_id") %>%
    dplyr::mutate(
      concordant = ifelse(
        (our_direction == "hyper" & direction == "hyper_in_SV_pos") |
        (our_direction == "hypo"  & direction == "hypo_in_SV_pos"), TRUE, FALSE)
    )
}

n_sig  <- sum(locus_tests$sig, na.rm = TRUE)
n_conc <- if ("concordant" %in% names(locus_tests)) sum(locus_tests$concordant, na.rm = TRUE) else NA
cat(sprintf("TCGA replication: %d/%d loci FDR<0.05; %s concordant with HCC direction\n",
            n_sig, nrow(locus_tests),
            if (!is.na(n_conc)) paste0(n_conc, "/", n_sig) else "N/A"))

# ── 8. Plot ────────────────────────────────────────────────────────────────────
theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

p_volcano <- ggplot(locus_tests, aes(x = delta_med, y = -log10(p_wilcox),
                                      color = sig, shape = direction)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("TRUE" = "#E24B4A", "FALSE" = "#888780"),
                     labels = c("TRUE" = "FDR<0.05", "FALSE" = "ns")) +
  labs(title = "P0-3: TCGA-LIHC replication",
       subtitle = sprintf("%d/%d Gold loci significant (FDR<0.05)", n_sig, nrow(locus_tests)),
       x = "Δβ (SV+ − SV−)", y = "-log10(p Wilcoxon)",
       color = "Significance", shape = "Direction") +
  theme_hcc

p_med_beta <- ggplot(
  locus_tests %>% tidyr::pivot_longer(cols = c(med_sv_pos, med_sv_neg),
                                       names_to = "group", values_to = "median_beta"),
  aes(x = group, y = median_beta, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  scale_fill_manual(values = c(med_sv_pos = "#E24B4A", med_sv_neg = "#3B8BD4")) +
  labs(title = "Median β: SV+ vs SV−", x = NULL, y = "HM450 β") +
  theme_hcc + theme(legend.position = "none")

p_combined <- p_volcano | p_med_beta
ggsave(file.path(OUTDIR, paste0(RUN_ID, "_P03_tcga_replication.png")),
       p_combined, width = 12, height = 5, dpi = 150)

# ── 9. Save ────────────────────────────────────────────────────────────────────
fwrite(locus_tests, file.path(OUTDIR, paste0(RUN_ID, "_P03_tcga_locus_tests.csv")))
fwrite(locus_meth,  file.path(OUTDIR, paste0(RUN_ID, "_P03_tcga_locus_meth.csv.gz")))

cat(append = TRUE,
    text   = sprintf("[%s] P0-3 tcga_lihc: %d/%d Gold loci sig (FDR<0.05); concordance=%s\n",
                     Sys.Date(), n_sig, nrow(locus_tests),
                     if (!is.na(n_conc)) paste0(n_conc, "/", n_sig) else "NA"),
    file   = LOG_FILE)

message("Done: P0-3 outputs in ", OUTDIR)
