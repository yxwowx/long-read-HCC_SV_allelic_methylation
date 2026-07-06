#!/usr/bin/env Rscript
# A1: CpG-density adjusted C14 logistic GLM
#
# Tests whether aDMR-SegDup enrichment (OR=2.22, C14) persists after
# adjusting for CpG density — identified as potential confounder in P1-B.
#
# Model: is_admr ~ segdup + log10(nCG_density) + lad + b_compartment
# Controls: matched random regions (same chr, same width), CpGs counted via BSgenome
#
# Output:
#   result/a1_cpg_adjusted_glm.csv        OR table (unadjusted vs CpG-adjusted)
#   result/figures/fig_a1_cpg_glm.png     side-by-side forest plot

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

GOLD_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/gold_tier_final.csv"
SILV_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/silver_tier.csv"
SEGDUP    <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
LAD       <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/laminB1Lads.bed"
PC1_BW    <- "/node200data/kachungk/reference/GRCh38/3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_cis_pc1.bw"
FAI       <- "/node200data/kachungk/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
OUT_DIR   <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_CTRL_MULT <- 10

# ── 1. Load unique aDMR loci with nCG ────────────────────────────────────────
message("Loading aDMR coordinates + nCG...")
admr_cols <- c("tier_class", "admr_chr", "admr_start", "admr_end", "nCG")
gold <- fread(GOLD_FILE, data.table = FALSE) |> select(all_of(admr_cols))
silv <- fread(SILV_FILE, data.table = FALSE) |> select(all_of(admr_cols))

admr <- bind_rows(gold, silv) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr)) |>
  distinct(admr_chr, admr_start, admr_end, .keep_all = TRUE)

admr$width_kb  <- (admr$admr_end - admr$admr_start + 1) / 1000
admr$nCG_dens  <- admr$nCG / pmax(admr$width_kb, 0.1)  # CpGs per kb
admr$log10_cpg <- log10(admr$nCG_dens + 1)

cat(sprintf("Unique aDMR loci: %d (Gold=%d, Silver=%d)\n",
            nrow(admr), sum(admr$tier_class=="Gold"), sum(admr$tier_class=="Silver")))
cat(sprintf("aDMR nCG density: median=%.1f CpGs/kb\n", median(admr$nCG_dens)))

admr_gr <- GRanges(
  seqnames = admr$admr_chr,
  ranges   = IRanges(admr$admr_start, admr$admr_end),
  tier     = admr$tier_class,
  is_admr  = 1L,
  nCG      = admr$nCG,
  nCG_dens = admr$nCG_dens,
  log10_cpg = admr$log10_cpg
)

# ── 2. Generate matched random control regions ────────────────────────────────
message("Generating matched controls...")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z"), data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

bsgenome <- BSgenome.Hsapiens.UCSC.hg38

ctrl_list <- lapply(seq_len(nrow(admr)), function(i) {
  chr <- admr$admr_chr[i]
  w   <- admr$admr_end[i] - admr$admr_start[i] + 1L
  len <- chrom_sizes$len[chrom_sizes$chr == chr]
  if (length(len) == 0) return(NULL)
  max_start <- len - w
  if (max_start < 1) return(NULL)
  starts <- sample.int(max_start, size = N_CTRL_MULT, replace = FALSE)
  GRanges(seqnames = chr,
          ranges   = IRanges(starts, starts + w - 1L),
          tier     = "Control",
          is_admr  = 0L)
})
ctrl_gr <- do.call(c, Filter(Negate(is.null), ctrl_list))
cat(sprintf("Controls generated: %d\n", length(ctrl_gr)))

# ── 3. Count CpG dinucleotides in control regions via BSgenome ───────────────
message("Counting CpGs in control regions (BSgenome)...")
ctrl_seqs      <- getSeq(bsgenome, ctrl_gr)
ctrl_cpg_cnt   <- vcountPattern("CG", ctrl_seqs)
ctrl_gr$nCG    <- ctrl_cpg_cnt
ctrl_gr$nCG_dens  <- ctrl_cpg_cnt / (width(ctrl_gr) / 1000)
ctrl_gr$log10_cpg <- log10(ctrl_gr$nCG_dens + 1)

cat(sprintf("Control nCG density: median=%.1f CpGs/kb\n", median(ctrl_gr$nCG_dens)))

all_gr <- c(admr_gr, ctrl_gr)

# ── 4. Annotate fragility features ───────────────────────────────────────────
message("Annotating segdup / LAD / B-compartment...")
segdup_gr <- import(SEGDUP, format="BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
lad_gr    <- import(LAD,    format="BED"); seqlevelsStyle(lad_gr)    <- "UCSC"

all_gr$segdup <- overlapsAny(all_gr, segdup_gr)
all_gr$lad    <- overlapsAny(all_gr, lad_gr)

bw <- BigWigFile(PC1_BW)
pc1_vals <- summary(bw, which = all_gr, type = "mean", defaultValue = NA_real_)
all_gr$pc1 <- unlist(lapply(pc1_vals, function(x)
  if (length(x$score)==0) NA_real_ else x$score[1]))
all_gr$b_compartment <- !is.na(all_gr$pc1) & all_gr$pc1 < 0

df <- data.frame(
  is_admr       = all_gr$is_admr,
  tier          = all_gr$tier,
  segdup        = as.integer(all_gr$segdup),
  lad           = as.integer(all_gr$lad),
  b_compartment = as.integer(all_gr$b_compartment),
  log10_cpg     = all_gr$log10_cpg,
  nCG_dens      = all_gr$nCG_dens
)
df$obs_weight <- ifelse(df$is_admr == 1L, 1, 1 / N_CTRL_MULT)

# ── 5. Logistic GLMs: unadjusted vs CpG-adjusted ────────────────────────────
message("Fitting logistic GLMs...")

fit_glm <- function(formula_str, data, label) {
  m  <- glm(as.formula(formula_str), data = data, weights = obs_weight, family = binomial())
  co <- summary(m)$coefficients
  ci <- confint.default(m)
  terms <- rownames(co)[-1]
  data.frame(
    feature = terms,
    model   = label,
    OR      = exp(co[terms, 1]),
    CI_lo   = exp(ci[terms, 1]),
    CI_hi   = exp(ci[terms, 2]),
    p       = co[terms, 4],
    stringsAsFactors = FALSE
  )
}

df_cc <- df |> filter(!is.na(log10_cpg))

res_base <- fit_glm("is_admr ~ segdup + lad + b_compartment", df_cc,
                    "Unadjusted (C14 base)")
res_adj  <- fit_glm("is_admr ~ segdup + log10_cpg + lad + b_compartment", df_cc,
                    "CpG-density adjusted")

res_all <- bind_rows(res_base, res_adj) |>
  mutate(
    sig = cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","ns")),
    feature_lab = recode(feature,
      segdup        = "SegDup overlap",
      log10_cpg     = "log10(CpG density)",
      lad           = "LAD overlap",
      b_compartment = "B-compartment")
  )

cat("\n=== A1 GLM results ===\n")
print(res_all |> select(model, feature_lab, OR, CI_lo, CI_hi, p, sig))

# ── 6. Gold-only adjusted model ──────────────────────────────────────────────
df_gold_ctrl <- df_cc |> filter(is_admr == 0 | tier == "Gold")
res_gold <- fit_glm("is_admr ~ segdup + log10_cpg + lad + b_compartment",
                    df_gold_ctrl, "CpG-adjusted (Gold only)")  # obs_weight used via closure

# ── 7. Save results ───────────────────────────────────────────────────────────
out <- bind_rows(res_all, res_gold |>
  mutate(sig = cut(p, c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")),
         feature_lab = recode(feature,
           segdup="SegDup overlap", log10_cpg="log10(CpG density)",
           lad="LAD overlap", b_compartment="B-compartment")))

fwrite(out, file.path(OUT_DIR, "a1_cpg_adjusted_glm.csv"))
message("Wrote: a1_cpg_adjusted_glm.csv")

# ── 8. Figure: side-by-side forest ───────────────────────────────────────────
feat_order <- c("SegDup overlap","log10(CpG density)","LAD overlap","B-compartment")

plot_df <- out |>
  filter(model %in% c("Unadjusted (C14 base)", "CpG-density adjusted")) |>
  filter(feature_lab %in% c("SegDup overlap","LAD overlap","B-compartment")) |>
  mutate(feature_lab = factor(feature_lab, levels = rev(c("SegDup overlap","LAD overlap","B-compartment"))),
         model = factor(model, levels = c("Unadjusted (C14 base)","CpG-density adjusted")))

seg_or_base <- res_base$OR[res_base$feature == "segdup"]
seg_or_adj  <- res_adj$OR[res_adj$feature  == "segdup"]
sig_base    <- as.character(cut(res_base$p[res_base$feature=="segdup"],
                c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")))
sig_adj     <- as.character(cut(res_adj$p[res_adj$feature=="segdup"],
                c(-Inf,0.001,0.01,0.05,Inf), labels=c("***","**","*","ns")))

subtitle_txt <- sprintf(
  "SegDup OR: unadjusted=%.2f %s → CpG-adjusted=%.2f %s",
  seg_or_base, sig_base, seg_or_adj, sig_adj)

pA <- ggplot(plot_df, aes(x = OR, y = feature_lab, color = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3.5, position = position_dodge(0.5)) +
  geom_text(aes(x = CI_hi * 1.15, label = sig),
            position = position_dodge(0.5), hjust = 0, size = 3) +
  scale_color_manual(values = c("Unadjusted (C14 base)" = "#6baed6",
                                 "CpG-density adjusted"  = "#d73027")) +
  scale_x_log10() +
  expand_limits(x = max(plot_df$CI_hi, na.rm=TRUE) * 1.6) +
  labs(
    title    = "A1: C14 SegDup enrichment after CpG-density adjustment",
    subtitle = subtitle_txt,
    x = "Odds Ratio (log scale)", y = NULL,
    color = "Model", shape = "Model",
    caption = "Logistic GLM: is_aDMR ~ segdup + lad + b_compartment [± log10(CpG/kb)]"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10))

# CpG density distribution: aDMR vs control
cpg_dens_df <- data.frame(
  group    = ifelse(df_cc$is_admr == 1, "aDMR", "Control"),
  nCG_dens = df_cc$nCG_dens,
  segdup   = factor(df_cc$segdup, labels = c("non-SegDup","SegDup"))
)
pB <- ggplot(cpg_dens_df |> filter(nCG_dens < quantile(nCG_dens, 0.99)),
             aes(x = nCG_dens, fill = group)) +
  geom_density(alpha = 0.5, bw = "SJ") +
  scale_fill_manual(values = c("aDMR" = "#d73027", "Control" = "#6baed6")) +
  facet_wrap(~segdup) +
  labs(title = "CpG density distribution by group",
       x = "CpG density (CpGs/kb)", y = "Density", fill = NULL) +
  theme_classic(base_size = 11)

combined <- pA / pB + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG_DIR, "fig_a1_cpg_glm.png"), combined,
       width = 10, height = 9, dpi = 150)
message("Saved: fig_a1_cpg_glm.png")

# summary to stdout
cat(sprintf("\n=== A1 SUMMARY ===\n"))
cat(sprintf("SegDup OR unadjusted : %.2f [%.2f-%.2f] %s\n",
            seg_or_base,
            res_base$CI_lo[res_base$feature=="segdup"],
            res_base$CI_hi[res_base$feature=="segdup"],
            sig_base))
cat(sprintf("SegDup OR CpG-adjusted: %.2f [%.2f-%.2f] %s\n",
            seg_or_adj,
            res_adj$CI_lo[res_adj$feature=="segdup"],
            res_adj$CI_hi[res_adj$feature=="segdup"],
            sig_adj))
cat("Done.\n")
