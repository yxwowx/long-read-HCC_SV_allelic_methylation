#!/usr/bin/env Rscript
# sv_clonality.R — SV clonality analysis using SVclone CCF estimates
#
# Uses SVclone CCF (purity-corrected) and CCF-based timing categories.
#
# Prerequisites:
#   1. bash post_processing/run_svclone_ccf.sh
#   2. mamba run -n renv Rscript post_processing/mutationtimer_sv_timing.R
#
# Input:
#   result/svclone/sv_timing_per_patient.csv — SVclone CCF + timing per SV
#   result/sv_fragility_annotation.csv       — SegDup/fragility annotations (P-codes)
#   ~/patient_code_mapping.csv               — JJT → P1 sample name mapping
#
# Output:
#   result/sv_clonality_fragility.csv
#   result/figures/fig_sv_clonality.png/pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

BASE     <- "/node200data/kachungk/hcc_data/DMR_SVs"
SVCLONE  <- file.path(BASE, "result/svclone")
FRAG_F   <- file.path(BASE, "result/sv_fragility_annotation.csv")
TAD_F    <- file.path(BASE, "sv_tad_ctcf_annotation.csv.gz")
MAP_F    <- path.expand("~/patient_code_mapping.csv")
OUT_DIR  <- file.path(BASE, "result")
FIG_DIR  <- file.path(OUT_DIR, "figures")
LOG_FILE <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(FIG_DIR, showWarnings = FALSE)

TIMING_LEVELS <- c("early clonal", "clonal [NA]", "late clonal", "subclonal")

# 1. Load inputs =====================================================
timing_f <- file.path(SVCLONE, "sv_timing_per_patient.csv")
if (!file.exists(timing_f)) stop("sv_timing_per_patient.csv not found. Run prerequisite scripts.")

message("Loading SVclone timing ...")
sv_dt <- fread(timing_f)
sv_dt[, timing := factor(timing, levels = TIMING_LEVELS)]

# Patient code mapping
pmap <- fread(MAP_F)
pmap[, sample_short := sub("_HCC$", "", Samples_ID)]
code_map <- pmap[, .(sample = sample_short, patient_code)]
sv_dt <- merge(sv_dt, code_map, by = "sample", all.x = TRUE)

# 2. Load fragility annotation (uses patient_code = P1..P12) ============
message("Loading fragility annotation ...")
frag <- fread(FRAG_F)
frag[, start := as.integer(start)]

# 3. Positional join SVclone → fragility (within ±50 bp) ================
message("Positional join SVclone ↔ fragility ...")
sv_dt[, seqnames := as.character(chr1)]
sv_dt[, bp_pos   := as.integer(pos1)]

# Rename patient_code in frag for join
frag2 <- frag[, .(bp_id, sample, seqnames, start, svtype,
                  segdup_overlap, lad_overlap, b_compartment, pc1_score)]
setkey(frag2, sample, seqnames, start)
# Join on patient_code (P1/P2..) which equals frag$sample
setkey(sv_dt, patient_code, seqnames, bp_pos)

# Rolling join: for each SVclone record find nearest fragility bp within 50 bp
matched <- frag2[sv_dt, on = .(sample = patient_code, seqnames, start = bp_pos), roll = 50, nomatch = NA]
n_matched <- sum(!is.na(matched$bp_id))
cat(sprintf("Matched SVclone → bp_id: %d / %d SVclone records (%.0f%%)\n",
            n_matched, nrow(sv_dt), 100 * n_matched / nrow(sv_dt)))

# 4. Load VAF from TAD file and join ====================================
message("Loading VAF from TAD annotation ...")
tad_vaf <- fread(TAD_F, select = c("bp_id", "VAF", "hVAF_H0", "hVAF_H1", "hVAF_H2", "vaf_concordant"))
tad_vaf <- unique(tad_vaf, by = "bp_id")
matched <- merge(matched, tad_vaf, by = "bp_id", all.x = TRUE)

# 5. Build final table ===================================================
df <- matched[, .(
  bp_id, seqnames, start, svtype,
  sample = i.sample,          # actual patient name (JJT etc.)
  patient_code = sample,      # P-code (P1 etc.) from frag join
  segdup_overlap = as.logical(segdup_overlap),
  lad_overlap, b_compartment, pc1_score,
  VAF, hVAF_H0, hVAF_H1, hVAF_H2, vaf_concordant,
  CCF = ccf, cellular_prevalence, purity,
  timing, most_likely_assignment
)]
df[, segdup_lab := ifelse(segdup_overlap, "SegDup SV", "Non-SegDup SV")]

# 5. Statistical tests ===================================================
df_test <- df[!is.na(CCF) & !is.na(segdup_overlap)]

wt_ccf <- if (nrow(df_test) > 5)
  wilcox.test(CCF ~ segdup_overlap, data = df_test, exact = FALSE) else list(p.value = NA)

summary_stats <- df_test[, .(
  n          = .N,
  median_CCF = round(median(CCF, na.rm = TRUE), 3),
  pct_clonal = round(100 * mean(timing %in% c("early clonal","clonal [NA]"), na.rm = TRUE), 1)
), by = segdup_lab]

cat("\n=== SVclone CCF by SegDup ===\n"); print(summary_stats)
cat(sprintf("Wilcoxon CCF SegDup vs Non-SegDup: p=%.4g\n", wt_ccf$p.value))

timing_tab <- df[!is.na(timing), .N, by = .(segdup_lab, timing)][order(segdup_lab, timing)]
cat("\n=== Timing × SegDup ===\n")
print(dcast(timing_tab, segdup_lab ~ timing, value.var = "N", fill = 0L))

# 6. Save CSV ===================================================
fwrite(df, file.path(OUT_DIR, "sv_clonality_fragility.csv"))
message("Wrote: sv_clonality_fragility.csv")

# 7. Figures ====================================================
pal <- c("SegDup SV" = "#d73027", "Non-SegDup SV" = "#4393c3")
pal_timing <- c(
  "early clonal" = "#d73027", "clonal [NA]" = "#fc8d59",
  "late clonal"  = "#4575b4", "subclonal"   = "#91bfdb"
)

sig_ccf <- ifelse(is.na(wt_ccf$p.value), "nd",
           ifelse(wt_ccf$p.value < 0.001, "***",
           ifelse(wt_ccf$p.value < 0.01, "**",
           ifelse(wt_ccf$p.value < 0.05, "*", "ns"))))

# Panel A: CCF violin by SegDup
pA <- ggplot(df_test, aes(x = segdup_lab, y = CCF, fill = segdup_lab)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1.5, y = 1.05,
           label = sprintf("%s\np=%.3g", sig_ccf, wt_ccf$p.value), size = 3.5) +
  scale_fill_manual(values = pal) +
  labs(title = "SVclone CCF by SegDup Overlap",
       subtitle = sprintf("n=%d SVs matched to fragility annotation", nrow(df_test)),
       x = NULL, y = "SVclone CCF (purity-corrected)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# Panel B: timing per patient stacked bar
pct_timing <- df[!is.na(timing), .N, by = .(sample, timing)]
pct_timing[, total := sum(N), by = sample]
pct_timing[, pct   := 100 * N / total]

pB <- ggplot(pct_timing, aes(x = sample, y = pct, fill = timing)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = pal_timing, name = "Timing",
                    drop = FALSE, limits = rev(TIMING_LEVELS)) +
  labs(title = "SV Timing by Patient", x = "Patient", y = "% SVs") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

# Panel C: timing × SegDup
if (nrow(timing_tab) > 0) {
  timing_tab[, total := sum(N), by = segdup_lab]
  timing_tab[, pct   := 100 * N / total]
  timing_tab[, timing := factor(timing, levels = TIMING_LEVELS)]
  pC <- ggplot(timing_tab, aes(x = segdup_lab, y = pct, fill = timing)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = pal_timing) +
    labs(title = "Timing × SegDup", x = NULL, y = "% SVs", fill = "Timing") +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")
} else {
  pC <- ggplot() + theme_void() + labs(title = "No timing data")
}

# Panel D: ECDF of CCF
pD <- ggplot(df_test, aes(x = CCF, color = segdup_lab)) +
  stat_ecdf(linewidth = 1) +
  scale_color_manual(values = pal) +
  labs(title = "ECDF of SVclone CCF", x = "CCF", y = "Cumulative proportion", color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

combined <- (pA | pB) / (pC | pD)
ggsave(file.path(FIG_DIR, "fig_sv_clonality.png"), combined, width = 12, height = 9, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_sv_clonality.pdf"), combined, width = 12, height = 9)
message("Saved: fig_sv_clonality.png/pdf")

# 8. Log ====================================================
cat(sprintf("[%s] sv_clonality: n=%d matched SVs, median_CCF_segdup=%.3f, median_CCF_nonsegdup=%.3f, p_ccf=%.4g\n",
            Sys.Date(), nrow(df_test),
            summary_stats[segdup_lab == "SegDup SV",     median_CCF],
            summary_stats[segdup_lab == "Non-SegDup SV", median_CCF],
            wt_ccf$p.value),
    file = LOG_FILE, append = TRUE)
message("Done.")
