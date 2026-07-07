#!/usr/bin/env Rscript
# P0-C: Multivariate logistic regression — independent fragility effects on SV breakpoints
#
# Design:
#   Cases   = somatic SV breakpoints (from sv_fragility_annotation.csv)
#   Controls = 5× random genome positions (same chr distribution, shuffled)
#   Response = is_sv (1/0)
#   Predictors = segdup_overlap + lad_overlap + b_compartment + pc1_score + repeat_density
#   Goal: confirm SegDup has independent effect after adjusting for LAD/B-comp confounders
#
# Output:
#   result/fragility_multivariate_glm.csv     per-predictor OR, CI, p
#   result/figures/fig_fragility_forest.png   forest plot

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

set.seed(42)

FRAG    <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result/sv_fragility_annotation.csv")
SEGDUP  <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
LAD     <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed")
RMSK    <- file.path(Sys.getenv("REFERENCE_DIR"), "LOLACore_180423/hg38/ucsc_features/regions/rmsk.bed")
PC1_BW  <- file.path(Sys.getenv("REFERENCE_DIR"), "3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw")
FAI     <- file.path(Sys.getenv("REFERENCE_DIR"), "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

CTRL_MULT <- 5   # controls per case

# 1. Load SV breakpoints (cases) ===============================================
message("Loading SV breakpoints...")
sv <- fread(FRAG, data.table = FALSE) |>
  filter(!is.na(seqnames), grepl("^chr[0-9XY]+$", seqnames)) |>
  distinct(bp_id, .keep_all = TRUE)
cat(sprintf("SV breakpoints: %d\n", nrow(sv)))

sv_gr <- GRanges(seqnames = sv$seqnames,
                 ranges   = IRanges(sv$start, sv$start),
                 is_sv    = 1L)

# 2. Generate random control regions (same chr distribution) ===================
message("Generating random controls...")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z"), data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

# match chromosome distribution of SVs
sv_chr_tab <- table(seqnames(sv_gr))
n_ctrl <- nrow(sv) * CTRL_MULT

ctrl_list <- lapply(names(sv_chr_tab), function(chr) {
  n   <- sv_chr_tab[[chr]] * CTRL_MULT
  len <- chrom_sizes$len[chrom_sizes$chr == chr]
  if (length(len) == 0 || n == 0) return(NULL)
  pos <- sample.int(len - 1L, size = n, replace = TRUE)
  GRanges(seqnames = chr, ranges = IRanges(pos, pos), is_sv = 0L)
})
ctrl_gr <- do.call(c, Filter(Negate(is.null), ctrl_list))
cat(sprintf("Controls generated: %d\n", length(ctrl_gr)))

all_gr <- c(sv_gr, ctrl_gr)

# 3. Annotate with fragility features ==========================================
message("Annotating with fragility features...")

# SegDup
segdup_gr <- import(SEGDUP, format = "BED")
seqlevelsStyle(segdup_gr) <- "UCSC"
all_gr$segdup <- overlapsAny(all_gr, segdup_gr)

# LAD
lad_gr <- import(LAD, format = "BED")
seqlevelsStyle(lad_gr) <- "UCSC"
all_gr$lad <- overlapsAny(all_gr, lad_gr)

# rmsk repeat density (count overlapping repeats per 10kb window)
rmsk_gr <- import(RMSK, format = "BED")
seqlevelsStyle(rmsk_gr) <- "UCSC"
windows <- GRanges(seqnames(all_gr),
                   IRanges(pmax(1, start(all_gr) - 5000),
                                start(all_gr) + 5000))
repeat_hits <- countOverlaps(windows, rmsk_gr)
all_gr$repeat_density <- repeat_hits

# PC1 from BigWig
message("Extracting PC1 scores from BigWig...")
bw         <- BigWigFile(PC1_BW)
pc1_vals   <- summary(bw, which = all_gr, type = "mean", defaultValue = NA_real_)
all_gr$pc1 <- unlist(lapply(pc1_vals, function(x) if (length(x$score) == 0) NA_real_ else x$score[1]))

# B-compartment flag (PC1 < 0)
all_gr$b_compartment <- !is.na(all_gr$pc1) & all_gr$pc1 < 0

cat(sprintf("Annotated %d regions\n", length(all_gr)))

# 4. Build data frame for logistic regression ==================================
df <- data.frame(
  is_sv           = all_gr$is_sv,
  segdup          = as.integer(all_gr$segdup),
  lad             = as.integer(all_gr$lad),
  b_compartment   = as.integer(all_gr$b_compartment),
  pc1_score       = all_gr$pc1,
  repeat_density  = all_gr$repeat_density
) |> filter(!is.na(segdup))

# Inverse-frequency weights: each control counts as 1/CTRL_MULT to balance the 1:CTRL_MULT
# design without inflating precision by the sampling multiplier.
df$obs_weight <- ifelse(df$is_sv == 1L, 1, 1 / CTRL_MULT)

cat(sprintf("Complete cases for regression: %d (SV=%d, ctrl=%d; weighted 1:1 effective ratio)\n",
            nrow(df), sum(df$is_sv), sum(!df$is_sv)))

# 5. Univariate ORs ============================================================
uni_features <- c("segdup", "lad", "b_compartment", "repeat_density")
uni_res <- lapply(uni_features, function(f) {
  m  <- glm(as.formula(paste("is_sv ~", f)), data = df, weights = obs_weight, family = binomial())
  co <- summary(m)$coefficients
  or <- exp(co[2, 1])
  ci <- exp(confint.default(m)[2, ])
  data.frame(predictor = f, model = "Univariate",
             OR = or, CI_lo = ci[1], CI_hi = ci[2], p = co[2, 4])
}) |> bind_rows()

# 6. Multivariate logistic regression ==========================================
message("Fitting multivariate logistic regression...")

# Model 1: binary features only (weighted to correct 1:CTRL_MULT sampling)
m_bin  <- glm(is_sv ~ segdup + lad + b_compartment, data = df, weights = obs_weight, family = binomial())
# Model 2: add continuous predictors
df_cc  <- df |> filter(!is.na(pc1_score))
m_full <- glm(is_sv ~ segdup + lad + b_compartment + repeat_density,
              data = df_cc, weights = obs_weight, family = binomial())

extract_or <- function(model, model_label, df_used) {
  co <- summary(model)$coefficients
  ci <- confint.default(model)
  terms <- rownames(co)[-1]
  data.frame(
    predictor = terms,
    model     = model_label,
    OR        = exp(co[terms, 1]),
    CI_lo     = exp(ci[terms, 1]),
    CI_hi     = exp(ci[terms, 2]),
    p         = co[terms, 4]
  )
}

multi_bin  <- extract_or(m_bin,  "Multivariate (binary)",   df)
multi_full <- extract_or(m_full, "Multivariate (full)",     df_cc)

result_df <- bind_rows(uni_res, multi_bin, multi_full) |>
  mutate(
    sig     = cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf),
                  labels = c("***", "**", "*", "ns")),
    OR_lab  = sprintf("%.2f (%.2f–%.2f) %s", OR, CI_lo, CI_hi, sig)
  )

cat("\n=== Multivariate GLM Results ===\n")
print(result_df |> select(predictor, model, OR, CI_lo, CI_hi, p, sig))

fwrite(result_df, file.path(OUT_DIR, "fragility_multivariate_glm.csv"))
message("Wrote: fragility_multivariate_glm.csv")

# AIC comparison
cat(sprintf("\nAIC — binary model: %.1f | full model: %.1f\n",
            AIC(m_bin), AIC(m_full)))

# 7. Forest plot ===============================================================
feature_labels <- c(
  segdup         = "SegDup overlap",
  lad            = "LAD overlap",
  b_compartment  = "B-compartment (PC1<0)",
  repeat_density = "Repeat density (10kb)"
)
plot_df <- result_df |>
  filter(model != "Univariate") |>
  mutate(
    predictor_lab = recode(predictor, !!!feature_labels),
    predictor_lab = factor(predictor_lab, levels = rev(unique(predictor_lab)))
  )

p <- ggplot(plot_df, aes(x = OR, xmin = CI_lo, xmax = CI_hi,
                          y = predictor_lab, color = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_text(aes(x = CI_hi * 1.05, label = sig),
            position = position_dodge(0.5), hjust = 0, size = 3.5) +
  scale_color_manual(values = c("Multivariate (binary)" = "#d73027",
                                "Multivariate (full)"   = "#1a9850")) +
  scale_x_log10() +
  labs(
    title    = "Fragility Features: Independent Effects on SV Breakpoints",
    subtitle = sprintf("Cases = SV breakpoints (n=%d); Controls = random regions (n=%d)",
                       sum(df$is_sv), sum(!df$is_sv)),
    x        = "Odds Ratio (log scale, vs random genome)",
    y        = NULL,
    color    = "Model", shape = "Model",
    caption  = "Weighted multivariate logistic regression (controls weight=1/CTRL_MULT); error bars = 95% CI"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "fig_fragility_forest.png"), p,
       width = 8, height = 5, dpi = 150)
message("Saved: fig_fragility_forest.png")

cat("\nDone.\n")
