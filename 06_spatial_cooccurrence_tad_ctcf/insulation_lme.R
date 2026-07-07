#!/usr/bin/env Rscript
# P2-7: Insulation score × HP|Δβ| LME
# Tests whether continuous Micro-C insulation score predicts HP|Δβ| better than
# the categorical tier hierarchy (C3 robustness: binning artifact vs true null).
#
# Model: abs_db ~ log2_ins + sv_arch + abs_dist_kb + (1|sample)
# Compared to: abs_db ~ sv_tier + abs_dist_kb + (1|sample)  [existing Layer 4]
#
# Run: mamba run -n renv Rscript pipeline/13_insulation_lme.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)  # p-values for lmer
  library(GenomicRanges)
  library(optparse)
})
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

DMR_SVS_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")

option_list <- list(
  make_option("--phased_ov", type = "character",
    default = file.path(DMR_SVS_DIR, "03.haplotype_sv_admr_analysis/all_hp_admr_tier.csv.gz")),
  make_option("--insulation", type = "character",
    default = file.path(DMR_SVS_DIR, "02.sv_dmr_enrichment/tad_ctcf_validation/insulation_8kb.tsv.gz")),
  make_option("--sv_file", type = "character",
    default = file.path(DMR_SVS_DIR, "sv_tad_ctcf_annotation.v2.csv.gz")),
  make_option("--outdir", type = "character",
    default = file.path(DMR_SVS_DIR, "03.haplotype_sv_admr_analysis")),
  make_option("--run_id", type = "character", default = "tier_v2"),
  make_option("--ins_window", type = "integer", default = 240000L,
    help = "Insulation window to use (80000/240000/480000, default 240000)"),
  make_option("--min_abs_db", type = "double", default = 0.05)
)
opt <- parse_args(OptionParser(option_list = option_list))

OUTDIR   <- opt$outdir
RUN_ID   <- opt$run_id
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

INS_COL <- paste0("log2_insulation_score_", opt$ins_window)

find_col <- function(df, cands, req = TRUE) {
  m <- intersect(cands, names(df)); if (length(m) > 0) return(m[[1]])
  if (req) stop("None of [", paste(cands, collapse=","), "] found"); NULL
}

# 1. Load insulation scores ====================================================
message("Reading insulation: ", opt$insulation)
ins <- fread(opt$insulation)
if (!INS_COL %in% names(ins))
  stop("Column '", INS_COL, "' not found. Available: ", paste(grep("log2", names(ins), value=T), collapse=", "))

ins_gr <- GRanges(seqnames = ins$chrom,
                  ranges   = IRanges(ins$start + 1L, ins$end),
                  log2_ins = ins[[INS_COL]])
ins_gr <- ins_gr[!is.nan(ins_gr$log2_ins) & !is.na(ins_gr$log2_ins)]

# 2. Load SV annotations =======================================================
message("Reading SV file: ", opt$sv_file)
sv <- fread(opt$sv_file)

col_chr_sv <- find_col(sv, c("seqnames", "chr", "chrom", "CHROM", "bp_chr"))
col_pos_sv <- find_col(sv, c("pos", "start", "POS", "bp_start"))
col_bpid   <- find_col(sv, c("bp_id", "sourceId", "ID"))
col_tier   <- find_col(sv, c("stratification", "sv_tier", "tier"))
col_pt_sv  <- find_col(sv, c("sample", "patient_code"))

sv_gr <- GRanges(seqnames = sv[[col_chr_sv]],
                 ranges   = IRanges(sv[[col_pos_sv]], sv[[col_pos_sv]]),
                 bp_id    = sv[[col_bpid]],
                 sv_tier  = sv[[col_tier]],
                 sample   = sv[[col_pt_sv]])
seqlevelsStyle(sv_gr) <- seqlevelsStyle(ins_gr)

# Annotate each SV with its local insulation score
hits_ins <- findOverlaps(sv_gr, ins_gr, select = "first")
sv_gr$log2_ins <- ins_gr$log2_ins[hits_ins]

sv_ann <- as.data.frame(sv_gr) %>%
  dplyr::filter(!is.na(log2_ins)) %>%
  dplyr::mutate(
    sv_arch = ifelse(sv_tier %in% c("TAD_CTCF", "TAD_only", "CTCF_only",
                                     "TAD+CTCF disrupting", "TAD-only", "CTCF-only"),
                     "boundary", "non_boundary")
  )

message(sprintf("SVs with insulation score: %d / %d", nrow(sv_ann), length(sv_gr)))

# 3. Load phased overlap =======================================================
message("Reading phased overlap: ", opt$phased_ov)
dat <- fread(opt$phased_ov)

col_abs_db <- find_col(dat, c("abs_hp_delta", "abs_delta_beta", "abs_db", "hp_abs_diff"))
col_dist   <- find_col(dat, c("dist", "distance", "abs_dist"))
col_bpid2  <- find_col(dat, c("bp_id", "sourceId", "sv_id"))
col_pt     <- find_col(dat, c("sample", "patient_code"))

dat <- dat %>%
  dplyr::rename(abs_db = !!col_abs_db, dist   = !!col_dist,
                bp_id  = !!col_bpid2,  sample = !!col_pt) %>%
  dplyr::filter(!is.na(abs_db), abs_db >= opt$min_abs_db) %>%
  dplyr::mutate(abs_dist_kb = abs(dist) / 1e3)

# Join insulation scores from SV annotations
dat <- dat %>%
  dplyr::left_join(
    sv_ann %>% dplyr::select(bp_id, log2_ins, sv_tier, sv_arch),
    by = "bp_id"
  ) %>%
  dplyr::filter(!is.na(log2_ins))

message(sprintf("Analysis rows after join: %d", nrow(dat)))

# 4. LME models ================================================================
# Model A: continuous insulation + arch + distance
# Model B: categorical tier + distance (existing Layer 4 equivalent)
# Model C: insulation + tier + distance (full)
# Model null: distance only

lme_safe <- function(formula, data) {
  tryCatch(lmer(formula, data = data, REML = FALSE,
                control = lmerControl(optimizer = "bobyqa")),
           error = function(e) { message("lmer failed: ", e$message); NULL })
}

m_null <- lme_safe(abs_db ~ abs_dist_kb + (1|sample), dat)
m_ins  <- lme_safe(abs_db ~ log2_ins + abs_dist_kb + (1|sample), dat)
m_tier <- lme_safe(abs_db ~ sv_tier  + abs_dist_kb + (1|sample), dat)
m_arch <- lme_safe(abs_db ~ sv_arch  + abs_dist_kb + (1|sample), dat)
m_full <- lme_safe(abs_db ~ log2_ins + sv_tier + abs_dist_kb + (1|sample), dat)

model_list <- list(null = m_null, ins_only = m_ins, tier = m_tier,
                   arch = m_arch, full = m_full)
model_list <- Filter(Negate(is.null), model_list)

# AIC comparison
aic_df <- data.frame(
  model = names(model_list),
  AIC   = sapply(model_list, AIC),
  BIC   = sapply(model_list, BIC)
) %>%
  dplyr::mutate(
    ΔAIC       = AIC - min(AIC),
    best_model = AIC == min(AIC)
  ) %>%
  dplyr::arrange(AIC)

cat("\n=== LME AIC comparison ===\n")
print(aic_df)

# Likelihood ratio test: ins_only vs null
lrt_ins <- if (!is.null(m_null) && !is.null(m_ins)) {
  tryCatch(anova(m_null, m_ins), error = function(e) NULL)
} else NULL
lrt_tier <- if (!is.null(m_null) && !is.null(m_tier)) {
  tryCatch(anova(m_null, m_tier), error = function(e) NULL)
} else NULL

cat("\n=== LRT: insulation vs null ===\n")
if (!is.null(lrt_ins)) print(lrt_ins)
cat("\n=== LRT: tier vs null ===\n")
if (!is.null(lrt_tier)) print(lrt_tier)

# Fixed-effect coefficients for best model
best_name <- aic_df$model[1]
best_mod  <- model_list[[best_name]]
cat(sprintf("\n=== Best model: %s (AIC=%.1f) ===\n", best_name, aic_df$AIC[1]))
if (!is.null(best_mod)) print(summary(best_mod)$coefficients)

# 5. Insulation scatter per tier ===============================================
ins_summary <- dat %>%
  dplyr::group_by(sv_tier, sv_arch) %>%
  dplyr::summarise(
    n          = n(),
    med_ins    = median(log2_ins, na.rm = TRUE),
    med_abs_db = median(abs_db,   na.rm = TRUE),
    cor_rho    = cor(log2_ins, abs_db, method = "spearman", use = "complete.obs"),
    .groups    = "drop"
  )
cat("\n=== Insulation × |Δβ| correlation per tier ===\n")
print(ins_summary)

# 6. Plots =====================================================================
p_scatter <- ggplot(dat %>% dplyr::sample_n(min(5000, nrow(.))),
                    aes(x = log2_ins, y = abs_db, color = sv_arch)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  facet_wrap(~sv_arch) +
  scale_color_manual(values = c(boundary = "#7F77DD", non_boundary = "#888780")) +
  labs(title = "P2-7: Insulation score × HP|Δβ|",
       subtitle = sprintf("Best LME: %s (ΔAIC vs null = %.1f)",
                          best_name, aic_df$ΔAIC[aic_df$model == "null"]),
       x = sprintf("log2 insulation score (%dkb window)", opt$ins_window / 1e3),
       y = "HP |Δβ|") +
  theme_hcc + theme(legend.position = "none")

p_aic <- ggplot(aic_df, aes(x = reorder(model, AIC), y = ΔAIC,
                              fill = best_model)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#E24B4A", "FALSE" = "#3B8BD4")) +
  coord_flip() +
  labs(title = "LME AIC comparison",
       x = "Model", y = "ΔAIC (lower = better)") +
  theme_hcc + theme(legend.position = "none")

p_tier_scatter <- ggplot(dat %>% dplyr::sample_n(min(5000, nrow(.))),
                          aes(x = log2_ins, y = abs_db)) +
  geom_point(alpha = 0.2, size = 0.7, color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "#E24B4A", linewidth = 1) +
  facet_wrap(~sv_tier, scales = "free_x") +
  labs(title = "By SV tier",
       x = sprintf("log2 insulation (%dkb)", opt$ins_window / 1e3), y = "HP |Δβ|") +
  theme_hcc

p_combined <- (p_scatter | p_aic) / p_tier_scatter +
  plot_layout(heights = c(1, 1.5))
ggsave(file.path(OUTDIR, paste0(RUN_ID, "_P27_insulation_lme.png")),
       p_combined, width = 14, height = 10, dpi = 150)

# 7. Save ======================================================================
fwrite(aic_df,       file.path(OUTDIR, paste0(RUN_ID, "_P27_lme_aic.csv")))
fwrite(ins_summary,  file.path(OUTDIR, paste0(RUN_ID, "_P27_ins_tier_cor.csv")))

message("Done: P2-7 outputs in ", OUTDIR)
