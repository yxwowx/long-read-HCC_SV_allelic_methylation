#!/usr/bin/env Rscript
# P1-5: Replication timing / fragility annotation
# Tests whether non-boundary SVs are enriched at late-replicating / fragile loci
# (fragility-driven vs selection-driven model for non-boundary SV dominance).
#
# Reference data:
#   Repli-seq: ENCODE HepG2 (ENCFF001LVM or ENCFF001LVK) — 6-fraction BED
#   Common fragile sites (CFS): doi:10.1038/ng.1079 (Fungtammasan 2012) → hg38-lifted BED
#
# Run: mamba run -n renv Rscript pipeline/11_replication_timing.R --fetch_encode

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(rtracklayer)
  library(optparse)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

option_list <- list(
  make_option("--sv_file", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz"),
  make_option("--repliseq_bed", type = "character",
    default = "/node200data/kachungk/reference/GRCh38/encode/HepG2_repliseq_6frac.bed.gz",
    help = "ENCODE HepG2 Repli-seq 6-fraction BED (download if absent)"),
  make_option("--cfs_bed", type = "character",
    default = "/node200data/kachungk/reference/GRCh38/genomic_element/common_fragile_sites_hg38.bed",
    help = "Common fragile sites BED (hg38)"),
  make_option("--outdir", type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment"),
  make_option("--run_id", type = "character", default = "tier_v2"),
  make_option("--fetch_encode", action = "store_true", default = FALSE,
    help = "Download ENCODE HepG2 Repli-seq if --repliseq_bed does not exist"),
  make_option("--n_perm", type = "integer", default = 1000L,
    help = "Permutations for CFS enrichment test")
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR   <- opt$outdir
RUN_ID   <- opt$run_id
LOG_FILE <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ENCODE_REPLISEQ_URL <- "https://www.encodeproject.org/files/ENCFF001LVM/@@download/ENCFF001LVM.bed.gz"

find_col <- function(df, cands, req = TRUE) {
  m <- intersect(cands, names(df)); if (length(m) > 0) return(m[[1]])
  if (req) stop("None of [", paste(cands, collapse=","), "] found"); NULL
}

# ── 1. Load SV data ───────────────────────────────────────────────────────────
message("Reading SV file: ", opt$sv_file)
sv <- fread(opt$sv_file)

col_chr    <- find_col(sv, c("seqnames", "chr", "chrom", "CHROM", "bp_chr"))
col_pos    <- find_col(sv, c("pos", "start", "POS", "bp_start"))
col_tier   <- find_col(sv, c("stratification", "sv_tier", "tier"))
col_type   <- find_col(sv, c("sv_type", "geom_type", "SVTYPE"))
col_pt     <- find_col(sv, c("sample", "patient_code"))

sv <- sv %>%
  dplyr::rename(chrom = !!col_chr, pos = !!col_pos,
                sv_tier = !!col_tier, sv_type = !!col_type, sample = !!col_pt) %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X")))

sv_gr <- GRanges(seqnames = sv$chrom,
                 ranges   = IRanges(sv$pos, sv$pos),
                 sv_tier  = sv$sv_tier,
                 sv_type  = sv$sv_type,
                 sample   = sv$sample)

sv_gr$sv_arch <- ifelse(
  sv_gr$sv_tier %in% c("TAD_CTCF", "TAD_only", "CTCF_only",
                         "TAD+CTCF disrupting", "TAD-only", "CTCF-only"),
  "boundary", "non_boundary"
)

# ── 2. Fetch / load Repli-seq ─────────────────────────────────────────────────
if (!file.exists(opt$repliseq_bed)) {
  if (opt$fetch_encode) {
    message("Downloading ENCODE HepG2 Repli-seq: ", ENCODE_REPLISEQ_URL)
    dir.create(dirname(opt$repliseq_bed), showWarnings = FALSE, recursive = TRUE)
    download.file(ENCODE_REPLISEQ_URL, opt$repliseq_bed, method = "wget", quiet = FALSE)
  } else {
    stop("Repli-seq file not found: ", opt$repliseq_bed,
         "\nRun with --fetch_encode to download from ENCODE, or provide path via --repliseq_bed")
  }
}

message("Loading Repli-seq: ", opt$repliseq_bed)
# Expected format: BED4+ with 4th column = fraction (G1b, S1, S2, S3, S4, G2)
# or BED with score column for RT value (higher = earlier replication)
repli_raw <- fread(opt$repliseq_bed, header = FALSE)
if (ncol(repli_raw) >= 5) {
  names(repli_raw)[1:5] <- c("chrom", "start", "end", "name", "score")
} else if (ncol(repli_raw) >= 4) {
  names(repli_raw)[1:4] <- c("chrom", "start", "end", "score")
  repli_raw$name <- "RT"
} else {
  names(repli_raw)[1:3] <- c("chrom", "start", "end")
  repli_raw$score <- NA_real_
}

# If multi-fraction BED (ENCODE 6-fraction): use fraction label to assign RT score
# G1b=1 (earliest), S1=2, S2=3, S3=4, S4=5, G2=6 (latest)
frac_order <- c(G1b=1, S1=2, S2=3, S3=4, S4=5, G2=6)
if ("name" %in% names(repli_raw) && any(repli_raw$name %in% names(frac_order))) {
  repli_raw <- repli_raw %>%
    dplyr::mutate(rt_score = frac_order[name])
  message("6-fraction Repli-seq: fraction label → RT score (1=early, 6=late)")
} else {
  repli_raw <- repli_raw %>% dplyr::mutate(rt_score = as.numeric(score))
  message("Continuous Repli-seq: using score column as RT value")
}

repli_gr <- GRanges(seqnames = repli_raw$chrom,
                    ranges   = IRanges(repli_raw$start + 1L, repli_raw$end),
                    rt_score = repli_raw$rt_score)
seqlevelsStyle(repli_gr) <- seqlevelsStyle(sv_gr)

# ── 3. Annotate SVs with replication timing ───────────────────────────────────
hits_rt <- findOverlaps(sv_gr, repli_gr, select = "first")
sv_gr$rt_score <- repli_gr$rt_score[hits_rt]

sv_rt <- as.data.frame(sv_gr) %>%
  dplyr::filter(!is.na(rt_score)) %>%
  dplyr::mutate(
    rt_class = cut(rt_score,
      breaks = quantile(rt_score, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
      labels = c("Early (Q1)", "EarlyMid (Q2)", "LateMid (Q3)", "Late (Q4)"),
      include.lowest = TRUE)
  )

cat("\n=== RT score distribution by SV architecture ===\n")
print(sv_rt %>%
  dplyr::group_by(sv_arch, rt_class) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = rt_class, values_from = n, values_fill = 0))

# KW: is rt_score different between boundary vs non_boundary SVs?
kw_rt <- kruskal.test(rt_score ~ sv_arch, data = sv_rt)
cat(sprintf("\nKW (rt_score ~ sv_arch): p = %.4f\n", kw_rt$p.value))

# Wilcoxon: non_boundary > boundary in late replication?
wt_rt <- wilcox.test(
  sv_rt$rt_score[sv_rt$sv_arch == "non_boundary"],
  sv_rt$rt_score[sv_rt$sv_arch == "boundary"],
  alternative = "greater", exact = FALSE
)
cat(sprintf("Wilcoxon (non_boundary RT > boundary RT): p = %.4f\n", wt_rt$p.value))

# ── 4. Common fragile sites ────────────────────────────────────────────────────
if (file.exists(opt$cfs_bed)) {
  message("Loading CFS: ", opt$cfs_bed)
  cfs_gr <- tryCatch(import(opt$cfs_bed), error = function(e) { message(e$message); NULL })
  if (!is.null(cfs_gr)) {
    seqlevelsStyle(cfs_gr) <- seqlevelsStyle(sv_gr)
    sv_gr$in_cfs <- overlapsAny(sv_gr, cfs_gr)
    sv_cfs <- as.data.frame(sv_gr)

    cfs_tab <- sv_cfs %>%
      dplyr::group_by(sv_arch) %>%
      dplyr::summarise(
        n_total  = n(),
        n_in_cfs = sum(in_cfs, na.rm = TRUE),
        pct_cfs  = mean(in_cfs, na.rm = TRUE) * 100,
        .groups  = "drop"
      )
    cat("\n=== CFS overlap by SV architecture ===\n")
    print(cfs_tab)

    # Permutation test: is non_boundary enriched at CFS relative to boundary?
    n_nb <- sum(sv_cfs$sv_arch == "non_boundary")
    n_b  <- sum(sv_cfs$sv_arch == "boundary")
    obs_diff <- mean(sv_cfs$in_cfs[sv_cfs$sv_arch == "non_boundary"], na.rm = TRUE) -
                mean(sv_cfs$in_cfs[sv_cfs$sv_arch == "boundary"],     na.rm = TRUE)
    perm_diffs <- replicate(opt$n_perm, {
      idx <- sample(nrow(sv_cfs), n_nb)
      mean(sv_cfs$in_cfs[idx], na.rm = TRUE) - mean(sv_cfs$in_cfs[-idx], na.rm = TRUE)
    })
    p_perm_cfs <- mean(perm_diffs >= obs_diff)
    cat(sprintf("CFS enrichment permutation: obs_diff=%.4f  p_perm=%.4f\n", obs_diff, p_perm_cfs))
  }
} else {
  message("CFS BED not found at: ", opt$cfs_bed,
          "\nDownload hg38-lifted CFS from supplementary of Fungtammasan 2012 / Wenger 2019.")
  p_perm_cfs <- NA; cfs_tab <- NULL
}

# ── 5. Plots ──────────────────────────────────────────────────────────────────
p_rt <- ggplot(sv_rt, aes(x = sv_arch, y = rt_score, fill = sv_arch)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c(boundary = "#7F77DD", non_boundary = "#888780")) +
  labs(title = "P1-5: Replication timing × SV architecture",
       subtitle = sprintf("KW p=%.4f  |  Wilcoxon (non_boundary>boundary) p=%.4f",
                          kw_rt$p.value, wt_rt$p.value),
       x = "SV architecture", y = "Replication timing score (higher = later)") +
  theme_hcc + theme(legend.position = "none")

p_tier_rt <- ggplot(sv_rt, aes(x = sv_tier, y = rt_score, fill = sv_arch)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1) +
  scale_fill_manual(values = c(boundary = "#7F77DD", non_boundary = "#888780")) +
  labs(title = "RT score by SV tier", x = NULL, y = "RT score") +
  coord_flip() + theme_hcc

p_combined <- p_rt | p_tier_rt
ggsave(file.path(OUTDIR, paste0(RUN_ID, "_P15_replication_timing.png")),
       p_combined, width = 12, height = 5, dpi = 150)

# ── 6. Save ───────────────────────────────────────────────────────────────────
rt_summary <- sv_rt %>%
  dplyr::group_by(sv_arch, sv_tier) %>%
  dplyr::summarise(n=n(), med_rt=median(rt_score,na.rm=T), .groups="drop")
fwrite(rt_summary, file.path(OUTDIR, paste0(RUN_ID, "_P15_rt_summary.csv")))
if (!is.null(cfs_tab)) fwrite(cfs_tab, file.path(OUTDIR, paste0(RUN_ID, "_P15_cfs_overlap.csv")))

cat(append = TRUE,
    text   = sprintf("[%s] P1-5 replication_timing: KW p=%.4f; Wilcoxon non_boundary>boundary p=%.4f; CFS_perm_p=%s\n",
                     Sys.Date(), kw_rt$p.value, wt_rt$p.value,
                     ifelse(is.na(p_perm_cfs), "NA", round(p_perm_cfs, 4))),
    file   = LOG_FILE)

message("Done: P1-5 outputs in ", OUTDIR)
