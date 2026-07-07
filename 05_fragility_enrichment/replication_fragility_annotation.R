#!/usr/bin/env Rscript
# Replication timing / fragile site annotation of SV breakpoints
# Proxies: LaminB1-LAD (late replication), HepG2 MicroC PC1 (B-compartment),
#          Segmental duplications, RepeatMasker density
# Usage: mamba run -n renv Rscript replication_fragility_annotation.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(rtracklayer)
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
REF      <- Sys.getenv("REFERENCE_DIR")
LAD_BED  <- file.path(REF, "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
SEGDUP   <- file.path(REF, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
RMSK     <- file.path(REF, "LOLACore_180423/hg38/ucsc_features/regions/rmsk.bed")
PC1_BW   <- file.path(REF, "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
SV_CSV   <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/sv_tad_ctcf_annotation.v2.csv.gz")
GOLD_CSV <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILV_CSV <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/04.final_candidate/silver_tier.csv")
OUT_DIR  <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")

# Load annotation tracks =======================================================
cat("Loading annotation tracks...\n")

lad_gr <- tryCatch({
  dt <- fread(LAD_BED, header=FALSE, col.names=c("chr","start","end"))
  makeGRangesFromDataFrame(dt)
}, error=function(e) { cat("LAD load error:", conditionMessage(e), "\n"); NULL })

segdup_gr <- tryCatch({
  dt <- fread(SEGDUP, header=FALSE, col.names=c("chr","start","end"))
  makeGRangesFromDataFrame(dt)
}, error=function(e) NULL)

rmsk_gr <- tryCatch({
  dt <- fread(RMSK, header=FALSE, col.names=c("chr","start","end"))
  makeGRangesFromDataFrame(dt)
}, error=function(e) NULL)

cat(sprintf("LAD regions: %d | SegDup: %d | rmsk: %d\n",
            length(lad_gr), length(segdup_gr), length(rmsk_gr)))

# Load SV breakpoints ==========================================================
cat("Loading SV breakpoints...\n")
sv <- fread(SV_CSV)
sv <- sv[!is.na(start) & seqnames %in% paste0("chr", c(1:22, "X","Y"))]

# Deduplicate: one row per unique breakpoint (bp_id)
sv_bp <- unique(sv[, .(bp_id, seqnames, start, sv_tier, stratification,
                        cnv_class, svtype, sample, patient_name)])
cat(sprintf("SV breakpoints: %d unique (from %d rows)\n", nrow(sv_bp), nrow(sv)))

sv_gr <- GRanges(sv_bp$seqnames, IRanges(sv_bp$start, sv_bp$start))
mcols(sv_gr) <- sv_bp[, .(bp_id, sv_tier, stratification, cnv_class, svtype, sample)]

# Annotate: LAD overlap ========================================================
cat("Annotating LAD overlap...\n")
sv_bp[, lad_overlap := countOverlaps(sv_gr, lad_gr) > 0]

# Annotate: SegDup overlap =====================================================
cat("Annotating SegDup overlap...\n")
sv_bp[, segdup_overlap := countOverlaps(sv_gr, segdup_gr) > 0]

# Annotate: repeat density in +/-50kb ==========================================
cat("Annotating repeat density...\n")
sv_win <- GRanges(sv_bp$seqnames, IRanges(pmax(1, sv_bp$start - 50000),
                                           sv_bp$start + 50000))
sv_bp[, repeat_density := countOverlaps(sv_win, rmsk_gr)]

# Annotate: HepG2 PC1 (B-compartment = negative) ===============================
cat("Annotating HepG2 PC1...\n")
pc1_bw <- import(PC1_BW, as = "GRanges")
pc1_hits <- findOverlaps(sv_gr, pc1_bw)
# Take the PC1 score of the overlapping 50kb bin
sv_bp[, pc1_score := NA_real_]
sv_bp[queryHits(pc1_hits), pc1_score := pc1_bw$score[subjectHits(pc1_hits)]]
sv_bp[, b_compartment := pc1_score < 0]   # TRUE = B compartment = late replicating

cat(sprintf("PC1 annotated: %d/%d breakpoints\n",
            sv_bp[!is.na(pc1_score), .N], nrow(sv_bp)))

# Group tiers: boundary (2-5) vs non_boundary (6) ==============================
# sv_tier: 2=TAD+CTCF, 3=TAD-only, 4=CTCF-only, 5=near-boundary, 6=non-boundary
sv_bp[, tier_group := ifelse(sv_tier == 6, "Non-boundary", "Boundary")]
sv_bp[, sv_tier_label := factor(paste0("Tier ", sv_tier),
                                 levels = paste0("Tier ", 2:6))]

# Statistical tests ============================================================
cat("\n=== Statistical Tests ===\n")

# LAD enrichment: Fisher test by tier group
lad_tab <- sv_bp[!is.na(lad_overlap),
                  table(tier_group, lad_overlap)]
ft_lad  <- fisher.test(lad_tab)
cat(sprintf("LAD overlap: Non-boundary %.1f%% vs Boundary %.1f%% | Fisher p=%.4f OR=%.2f\n",
            100*mean(sv_bp[tier_group=="Non-boundary", lad_overlap]),
            100*mean(sv_bp[tier_group=="Boundary",     lad_overlap]),
            ft_lad$p.value, ft_lad$estimate))

# B-compartment: Fisher test
bc_tab <- sv_bp[!is.na(b_compartment),
                 table(tier_group, b_compartment)]
ft_bc  <- fisher.test(bc_tab)
cat(sprintf("B-compartment: Non-boundary %.1f%% vs Boundary %.1f%% | Fisher p=%.4f OR=%.2f\n",
            100*mean(sv_bp[tier_group=="Non-boundary" & !is.na(b_compartment), b_compartment]),
            100*mean(sv_bp[tier_group=="Boundary"     & !is.na(b_compartment), b_compartment]),
            ft_bc$p.value, ft_bc$estimate))

# PC1 score: Wilcoxon
wt_pc1 <- wilcox.test(pc1_score ~ tier_group, data = sv_bp[!is.na(pc1_score)])
cat(sprintf("PC1 score: Non-boundary median=%.3f vs Boundary median=%.3f | Wilcoxon p=%.4f\n",
            sv_bp[tier_group=="Non-boundary" & !is.na(pc1_score), median(pc1_score)],
            sv_bp[tier_group=="Boundary"     & !is.na(pc1_score), median(pc1_score)],
            wt_pc1$p.value))

# SegDup: Fisher test
sd_tab <- sv_bp[!is.na(segdup_overlap), table(tier_group, segdup_overlap)]
ft_sd  <- fisher.test(sd_tab)
cat(sprintf("SegDup overlap: Non-boundary %.1f%% vs Boundary %.1f%% | Fisher p=%.4f OR=%.2f\n",
            100*mean(sv_bp[tier_group=="Non-boundary", segdup_overlap]),
            100*mean(sv_bp[tier_group=="Boundary",     segdup_overlap]),
            ft_sd$p.value, ft_sd$estimate))

# Repeat density: Wilcoxon
wt_rep <- wilcox.test(repeat_density ~ tier_group, data = sv_bp)
cat(sprintf("Repeat density: Non-boundary median=%.0f vs Boundary median=%.0f | Wilcoxon p=%.4f\n",
            sv_bp[tier_group=="Non-boundary", median(repeat_density)],
            sv_bp[tier_group=="Boundary",     median(repeat_density)],
            wt_rep$p.value))

# Annotate recurrent Gold/Silver loci ==========================================
cat("\n=== Recurrent Loci Annotation ===\n")
gold   <- fread(GOLD_CSV)
silver <- fread(SILV_CSV)
both   <- rbind(gold, silver, fill=TRUE)

# Define the 6 zoom loci
zoom_loci <- data.table(
  label = c("chr2:87.35Mb(n=7)","chrY:20.1Mb(n=7)","chr3:195.5Mb(n=3,HYPER)",
            "chr2:91.5Mb(n=3)","chrX:PAR1(n=3)","chr5:tel(n=4)"),
  chr   = c("chr2","chrY","chr3","chr2","chrX","chr5"),
  pos   = c(87390000, 20200000, 195480000, 91480000, 500000, 720000)
)

zoom_gr <- GRanges(zoom_loci$chr, IRanges(zoom_loci$pos, zoom_loci$pos))

# LAD
zoom_loci[, lad    := countOverlaps(zoom_gr, lad_gr) > 0]
# B compartment
zoom_hits <- findOverlaps(zoom_gr, pc1_bw)
zoom_loci[, pc1   := NA_real_]
zoom_loci[queryHits(zoom_hits), pc1 := pc1_bw$score[subjectHits(zoom_hits)]]
zoom_loci[, b_comp := pc1 < 0]
# SegDup
zoom_loci[, segdup := countOverlaps(zoom_gr, segdup_gr) > 0]
# Repeat density ±100kb
zoom_win <- GRanges(zoom_loci$chr, IRanges(pmax(1, zoom_loci$pos - 100000),
                                            zoom_loci$pos + 100000))
zoom_loci[, repeat_dens := countOverlaps(zoom_win, rmsk_gr)]

cat("Recurrent locus annotations:\n")
print(zoom_loci[, .(label, lad, b_comp, pc1, segdup, repeat_dens)])

# Save results table ===========================================================
fwrite(sv_bp[, .(bp_id, seqnames, start, sv_tier, tier_group, stratification,
                  cnv_class, svtype, sample, lad_overlap, b_compartment,
                  pc1_score, segdup_overlap, repeat_density)],
       file.path(OUT_DIR, "result/sv_fragility_annotation.csv"))

fwrite(zoom_loci, file.path(OUT_DIR, "result/recurrent_loci_fragility.csv"))

# Figures ======================================================================
cat("\nBuilding figures...\n")

# Panel A: PC1 violin by tier group ============================================
sv_plot <- sv_bp[!is.na(pc1_score)]
pa <- ggplot(sv_plot, aes(x = tier_group, y = pc1_score, fill = tier_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Non-boundary" = "#d73027", "Boundary" = "#4393c3")) +
  annotate("text", x = 1.5, y = max(sv_plot$pc1_score, na.rm=TRUE) * 0.9,
           label = sprintf("Wilcoxon p=%.4f", wt_pc1$p.value), size = 3) +
  labs(x = NULL, y = "HepG2 MicroC PC1\n(negative = B compartment = late replicating)",
       title = "A. Chromatin compartment at SV breakpoints",
       fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

# Panel B: Stacked bar — LAD + B-comp + SegDup rate by tier group ==============
sv_bp2 <- sv_bp[!is.na(lad_overlap) & !is.na(b_compartment)]
bar_dt <- rbind(
  sv_bp2[, .(tier_group, feature = "LAD overlap",     value = lad_overlap)],
  sv_bp2[, .(tier_group, feature = "B-compartment",   value = b_compartment)],
  sv_bp2[, .(tier_group, feature = "SegDup overlap",  value = segdup_overlap)]
)
bar_sum <- bar_dt[, .(pct = mean(value, na.rm=TRUE) * 100, n = .N),
                   by = .(tier_group, feature)]
bar_sum[, feature := factor(feature, levels = c("LAD overlap","B-compartment","SegDup overlap"))]

pb <- ggplot(bar_sum, aes(x = feature, y = pct, fill = tier_group)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Non-boundary" = "#d73027", "Boundary" = "#4393c3")) +
  labs(x = NULL, y = "% breakpoints overlapping",
       title = "B. Late-replication / fragility markers",
       fill = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Panel C: Repeat density boxplot by tier group ================================
pc <- ggplot(sv_bp, aes(x = tier_group, y = log10(repeat_density + 1), fill = tier_group)) +
  geom_violin(alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Non-boundary" = "#d73027", "Boundary" = "#4393c3")) +
  annotate("text", x = 1.5, y = max(log10(sv_bp$repeat_density+1)) * 0.95,
           label = sprintf("Wilcoxon p=%.4f", wt_rep$p.value), size = 3) +
  labs(x = NULL, y = "log10(repeat elements in ±50kb + 1)",
       title = "C. Repeat element density at SV breakpoints",
       fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

# Panel D: Per-tier breakdown of B-compartment + LAD rate ======================
tier_sum <- sv_bp[!is.na(b_compartment),
                   .(pct_b   = mean(b_compartment, na.rm=TRUE) * 100,
                     pct_lad = mean(lad_overlap, na.rm=TRUE)   * 100,
                     n = .N),
                   by = sv_tier_label]
tier_long <- melt(tier_sum, id.vars = c("sv_tier_label","n"),
                  measure.vars = c("pct_b","pct_lad"),
                  variable.name = "feature", value.name = "pct")
tier_long[, feature := ifelse(feature == "pct_b", "B-compartment", "LAD")]

pd <- ggplot(tier_long, aes(x = sv_tier_label, y = pct, color = feature, group = feature)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text(data = tier_sum, aes(x = sv_tier_label, y = -5, label = paste0("n=",n)),
            inherit.aes = FALSE, size = 2.5) +
  scale_color_manual(values = c("B-compartment" = "#7b3294", "LAD" = "#1a9641")) +
  labs(x = "SV tier (2=TAD+CTCF → 6=non-boundary)",
       y = "% breakpoints",
       title = "D. Fragility rate across tier gradient",
       color = NULL) +
  theme_bw(base_size = 10)

# Panel E: Recurrent loci dot plot =============================================
zoom_long <- melt(zoom_loci[, .(label, lad, b_comp, segdup)],
                  id.vars = "label", variable.name = "feature", value.name = "present")
zoom_long[, feature := factor(feature,
  levels = c("lad","b_comp","segdup"),
  labels = c("LAD overlap","B-compartment","SegDup overlap"))]
zoom_long[, label := factor(label, levels = rev(zoom_loci$label))]

pe <- ggplot(zoom_long, aes(x = feature, y = label, fill = present)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(present, "✓", "–")), size = 4) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#f7f7f7")) +
  labs(x = NULL, y = NULL,
       title = "E. Recurrent Gold/Silver loci: fragility context",
       fill = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none")

# Combine and save =============================================================
top <- (pa | pb | pc) + plot_layout(widths = c(1,1,1))
bot <- (pd | pe) + plot_layout(widths = c(2, 1.2))
combined <- top / bot + plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title   = "Replication timing & fragile site annotation of SV breakpoints",
    subtitle = "Non-boundary SVs (tier 6) vs Boundary-disrupting SVs (tiers 2–5)",
    theme = theme(plot.title = element_text(size=13, face="bold"),
                  plot.subtitle = element_text(size=10))
  )

out_png <- file.path(OUT_DIR, "figs/png/fig_replication_fragility.png")
ggsave(out_png, combined, width = 14, height = 9, dpi = 150)
cat(sprintf("Saved: %s\n", out_png))

cat("Done.\n")
