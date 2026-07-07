#!/usr/bin/env Rscript
# A3: SV-exclusion sensitivity test for C14
#
# If aDMRs are truly enriched at SegDup independently of SV proximity,
# the SegDup OR should persist after excluding aDMRs near SVs.
# Tests three subsets: all, SV-far >10kb, SV-far >50kb.
#
# Output:
#   result/a3_sv_exclusion_segdup.csv      OR by distance threshold
#   result/figures/fig_a3_sv_exclusion.png  forest + bar plot

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

set.seed(42)

GOLD_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILV_FILE <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/silver_tier.csv")
SV_FILE   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/sv_tad_ctcf_annotation.csv.gz")
SEGDUP    <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD       <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
PC1_BW    <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
FAI       <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_CTRL_MULT <- 10

# 1. Load aDMR loci ============================================================
message("Loading aDMR loci...")
admr_cols <- c("tier_class","admr_chr","admr_start","admr_end","nCG")
gold <- fread(GOLD_FILE, data.table=FALSE) |> select(all_of(admr_cols))
silv <- fread(SILV_FILE, data.table=FALSE) |> select(all_of(admr_cols))

admr <- bind_rows(gold, silv) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr)) |>
  distinct(admr_chr, admr_start, admr_end, .keep_all=TRUE)

cat(sprintf("Total unique aDMR loci: %d\n", nrow(admr)))

admr_gr <- GRanges(
  seqnames = admr$admr_chr,
  ranges   = IRanges(admr$admr_start, admr$admr_end),
  tier     = admr$tier_class,
  locus_id = seq_len(nrow(admr))
)

# 2. Load SV breakpoints =======================================================
message("Loading SV breakpoints...")
sv <- fread(SV_FILE, data.table=FALSE) |>
  filter(grepl("^chr[0-9XY]+$", seqnames))

sv_gr <- GRanges(
  seqnames = sv$seqnames,
  ranges   = IRanges(sv$start, sv$start)  # breakpoints are single positions
)
cat(sprintf("SV breakpoints: %d\n", length(sv_gr)))

# 3. Compute minimum SV distance for each aDMR =================================
message("Computing min SV distance per aDMR...")
dist_mat <- distanceToNearest(admr_gr, sv_gr, ignore.strand=TRUE)
admr$min_sv_dist <- Inf
admr$min_sv_dist[queryHits(dist_mat)] <- mcols(dist_mat)$distance

cat(sprintf("aDMRs with SV within 10kb : %d\n", sum(admr$min_sv_dist <= 10000)))
cat(sprintf("aDMRs with SV within 50kb : %d\n", sum(admr$min_sv_dist <= 50000)))
cat(sprintf("aDMRs SV-far >10kb        : %d\n", sum(admr$min_sv_dist >  10000)))
cat(sprintf("aDMRs SV-far >50kb        : %d\n", sum(admr$min_sv_dist >  50000)))

# 4. Load fragility feature annotations ========================================
message("Loading fragility features...")
segdup_gr <- import(SEGDUP, format="BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr    <- import(LAD,    format="BED"); seqlevelsStyle(lad_gr)    <- "UCSC"
chrom_sizes <- fread(FAI, col.names=c("chr","len","x","y","z"), data.table=FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |> select(chr, len)

# Helper: build dataset (aDMR subset + matched controls) and annotate
build_df <- function(admr_sub) {
  if (nrow(admr_sub) < 5) return(NULL)
  gr_a <- GRanges(seqnames=admr_sub$admr_chr,
                  ranges=IRanges(admr_sub$admr_start, admr_sub$admr_end),
                  is_admr=1L)
  # generate controls
  ctrl_list <- lapply(seq_len(nrow(admr_sub)), function(i) {
    chr <- admr_sub$admr_chr[i]
    w   <- admr_sub$admr_end[i] - admr_sub$admr_start[i] + 1L
    len <- chrom_sizes$len[chrom_sizes$chr==chr]
    if (length(len)==0 || len - w < 1) return(NULL)
    starts <- sample.int(len - w, size=N_CTRL_MULT, replace=FALSE)
    GRanges(seqnames=chr, ranges=IRanges(starts, starts+w-1L), is_admr=0L)
  })
  gr_c <- do.call(c, Filter(Negate(is.null), ctrl_list))
  all  <- c(gr_a, gr_c)
  all$segdup        <- overlapsAny(all, segdup_gr)
  all$lad           <- overlapsAny(all, lad_gr)
  bw <- BigWigFile(PC1_BW)
  pc1v <- summary(bw, which=all, type="mean", defaultValue=NA_real_)
  all$pc1 <- unlist(lapply(pc1v, function(x) if(length(x$score)==0) NA_real_ else x$score[1]))
  all$b_compartment <- !is.na(all$pc1) & all$pc1 < 0
  data.frame(
    is_admr       = all$is_admr,
    segdup        = as.integer(all$segdup),
    lad           = as.integer(all$lad),
    b_compartment = as.integer(all$b_compartment),
    obs_weight    = ifelse(all$is_admr == 1L, 1, 1 / N_CTRL_MULT)
  )
}

# Helper: run Fisher + GLM on a df
run_enrichment <- function(df, label) {
  if (is.null(df) || nrow(df) < 10) return(NULL)
  # Fisher for SegDup
  tab <- table(is_admr=df$is_admr, segdup=df$segdup)
  if (nrow(tab)<2 || ncol(tab)<2) return(NULL)
  ft  <- fisher.test(tab, simulate.p.value=TRUE, B=5000)
  # GLM
  m   <- glm(is_admr ~ segdup + lad + b_compartment, data=df, weights=obs_weight, family=binomial())
  co  <- summary(m)$coefficients
  ci  <- confint.default(m)
  glm_or    <- exp(co["segdup",1])
  glm_ci_lo <- exp(ci["segdup",1])
  glm_ci_hi <- exp(ci["segdup",2])
  glm_p     <- co["segdup",4]
  data.frame(
    subset     = label,
    n_admr     = sum(df$is_admr),
    pct_admr_segdup = 100 * mean(df$segdup[df$is_admr==1]),
    pct_ctrl_segdup = 100 * mean(df$segdup[df$is_admr==0]),
    fisher_OR  = ft$estimate,
    fisher_p   = ft$p.value,
    glm_OR     = glm_or,
    glm_CI_lo  = glm_ci_lo,
    glm_CI_hi  = glm_ci_hi,
    glm_p      = glm_p,
    stringsAsFactors=FALSE
  )
}

# 5. Run enrichment for each subset ============================================
thresholds <- list(
  "All aDMRs"      = admr,
  "SV-far > 10kb"  = admr[admr$min_sv_dist > 10000, ],
  "SV-far > 50kb"  = admr[admr$min_sv_dist > 50000, ]
)

results <- lapply(names(thresholds), function(nm) {
  message(sprintf("Processing subset: %s (n=%d)...", nm, nrow(thresholds[[nm]])))
  df <- build_df(thresholds[[nm]])
  run_enrichment(df, nm)
})

res_df <- bind_rows(Filter(Negate(is.null), results)) |>
  mutate(
    sig_fisher = cut(fisher_p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")),
    sig_glm    = cut(glm_p,    c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns"))
  )

cat("\n=== A3 SV-exclusion results ===\n")
print(res_df |> select(subset, n_admr, glm_OR, glm_CI_lo, glm_CI_hi, glm_p, sig_glm))

fwrite(res_df, file.path(OUT_DIR, "a3_sv_exclusion_segdup.csv"))
message("Wrote: a3_sv_exclusion_segdup.csv")

# 6. Figure ====================================================================
res_df$subset <- factor(res_df$subset,
  levels=c("All aDMRs","SV-far > 10kb","SV-far > 50kb"))

pA <- ggplot(res_df, aes(x=glm_OR, xmin=glm_CI_lo, xmax=glm_CI_hi, y=subset)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=glm_CI_lo, xmax=glm_CI_hi), height=0.25, color="#d73027") +
  geom_point(size=4, color="#d73027") +
  geom_text(aes(x=glm_CI_hi*1.12,
                label=sprintf("OR=%.2f\nn=%d\n%s", glm_OR, n_admr, sig_glm)),
            hjust=0, size=3) +
  scale_x_log10() +
  expand_limits(x=max(res_df$glm_CI_hi)*1.8) +
  labs(
    title    = "A3: aDMR-SegDup enrichment after SV exclusion",
    subtitle = "SegDup OR (multivariate GLM) across SV-distance thresholds",
    x = "Odds Ratio (log scale)", y = NULL,
    caption  = "Model: is_aDMR ~ segdup + lad + b_compartment (matched random controls)"
  ) +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"))

pct_df <- tidyr::pivot_longer(
  res_df |> select(subset, pct_admr_segdup, pct_ctrl_segdup),
  cols=c(pct_admr_segdup, pct_ctrl_segdup),
  names_to="group", values_to="pct"
) |> mutate(group=recode(group,
  pct_admr_segdup="aDMR", pct_ctrl_segdup="Random control"))

pB <- ggplot(pct_df, aes(x=subset, y=pct, fill=group)) +
  geom_col(position=position_dodge(0.7), width=0.65, alpha=0.85) +
  scale_fill_manual(values=c("aDMR"="#d73027","Random control"="#6baed6")) +
  labs(title="% SegDup overlap by SV-exclusion subset",
       x=NULL, y="% overlapping SegDup", fill=NULL) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=15, hjust=1))

combined <- pA / pB + plot_layout(heights=c(2,1))
ggsave(file.path(FIG_DIR, "fig_a3_sv_exclusion.png"), combined,
       width=9, height=8, dpi=150)
message("Saved: fig_a3_sv_exclusion.png")

cat(sprintf("\n=== A3 SUMMARY ===\n"))
for (i in seq_len(nrow(res_df))) {
  cat(sprintf("%-20s  n=%3d  GLM OR=%.2f [%.2f-%.2f]  p=%.3g  %s\n",
    as.character(res_df$subset[i]), res_df$n_admr[i],
    res_df$glm_OR[i], res_df$glm_CI_lo[i], res_df$glm_CI_hi[i],
    res_df$glm_p[i], as.character(res_df$sig_glm[i])))
}
cat("Done.\n")
