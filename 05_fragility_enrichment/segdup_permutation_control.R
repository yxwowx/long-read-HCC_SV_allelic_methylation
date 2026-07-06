#!/usr/bin/env Rscript
# A2: Permuted SegDup mask control for C13/C14
#
# Tests whether observed SV/aDMR-SegDup enrichment is specific to SegDup
# genomic positions, vs. any similarly-sized genomic feature.
# SegDup regions are shuffled 1000 times within chromosomes (preserving
# total coverage per chromosome). Observed OR is compared to null distribution.
#
# Output:
#   result/a2_segdup_permutation.csv         empirical OR + p-values
#   result/figures/fig_a2_segdup_permutation.png  null distribution plots

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

GOLD_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/gold_tier_final.csv"
SILV_FILE <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/silver_tier.csv"
SV_FRAG   <- "/node200data/kachungk/hcc_data/DMR_SVs/result/sv_fragility_annotation.csv"
SEGDUP    <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
FAI       <- "/node200data/kachungk/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
OUT_DIR   <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

N_PERM <- 1000

# ── 1. Load SVs (non-boundary), aDMRs, SegDup ───────────────────────────────
message("Loading SVs, aDMRs, SegDup...")

sv_frag <- fread(SV_FRAG, data.table=FALSE) |>
  filter(grepl("^chr[0-9XY]+$", seqnames))
sv_gr <- GRanges(seqnames=sv_frag$seqnames, ranges=IRanges(sv_frag$start, sv_frag$start))
cat(sprintf("Non-boundary SVs: %d\n", length(sv_gr)))

admr_cols <- c("admr_chr","admr_start","admr_end")
gold <- fread(GOLD_FILE, data.table=FALSE) |> select(all_of(admr_cols)) |> mutate(tier="Gold")
silv <- fread(SILV_FILE, data.table=FALSE) |> select(all_of(admr_cols)) |> mutate(tier="Silver")
admr <- bind_rows(gold, silv) |>
  filter(grepl("^chr[0-9XY]+$", admr_chr)) |>
  distinct(admr_chr, admr_start, admr_end, .keep_all=TRUE)
admr_gr <- GRanges(seqnames=admr$admr_chr, ranges=IRanges(admr$admr_start, admr$admr_end))
cat(sprintf("Gold+Silver aDMRs: %d\n", length(admr_gr)))

segdup_gr  <- import(SEGDUP, format="BED"); seqlevelsStyle(segdup_gr) <- "UCSC"
segdup_gr  <- keepStandardChromosomes(segdup_gr, pruning.mode="coarse")
cat(sprintf("SegDup regions: %d\n", length(segdup_gr)))

chrom_sizes <- fread(FAI, col.names=c("chr","len","x","y","z"), data.table=FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |> select(chr, len)
chrom_len <- setNames(chrom_sizes$len, chrom_sizes$chr)

# ── 2. Observed OR (SV and aDMR vs SegDup) ──────────────────────────────────
calc_or <- function(query_gr, feature_gr) {
  n_feat   <- sum(overlapsAny(query_gr, feature_gr))
  n_nofeat <- length(query_gr) - n_feat
  # background: fraction of genome covered
  genome_len <- sum(chrom_len[names(chrom_len) %in% as.character(unique(seqnames(query_gr)))])
  feat_len   <- sum(width(reduce(feature_gr[seqnames(feature_gr) %in%
    as.character(unique(seqnames(query_gr)))])))
  bg_rate <- feat_len / genome_len
  bg_nfeat   <- round(length(query_gr) * bg_rate)
  bg_nofeat  <- length(query_gr) - bg_nfeat
  if (bg_nfeat < 1 || bg_nofeat < 1) return(NA_real_)
  tab <- matrix(c(n_feat, n_nofeat, bg_nfeat, bg_nofeat), nrow=2,
                dimnames=list(c("observed","background"), c("feat","nofeat")))
  ft <- fisher.test(tab, simulate.p.value=TRUE, B=2000)
  ft$estimate
}

obs_sv_or   <- calc_or(sv_gr,   segdup_gr)
obs_admr_or <- calc_or(admr_gr, segdup_gr)
cat(sprintf("Observed SV OR vs SegDup  : %.3f\n", obs_sv_or))
cat(sprintf("Observed aDMR OR vs SegDup: %.3f\n", obs_admr_or))

# ── 3. Permutation: shuffle SegDup within chromosomes ───────────────────────
message(sprintf("Running %d permutations...", N_PERM))

shuffle_segdup <- function(segdup_gr, chrom_len) {
  # Per-chromosome circular shift: one uniform random offset per chromosome.
  # All intervals on a chromosome shift together, preserving relative distances
  # and total coverage per chromosome (no overlap creation, no coverage change).
  chrs   <- as.character(seqnames(segdup_gr))
  starts <- start(segdup_gr)
  widths <- width(segdup_gr)
  chr_uniq <- unique(chrs)
  shifts <- setNames(
    vapply(chr_uniq, function(ch) {
      cl <- chrom_len[ch]
      if (is.na(cl) || cl < 1L) return(0L)
      sample.int(cl, 1L)
    }, integer(1L)),
    chr_uniq
  )
  new_starts <- as.integer((starts + shifts[chrs]) %% pmax(chrom_len[chrs], 1L))
  new_starts <- pmax(1L, new_starts)
  GRanges(seqnames = chrs,
          ranges   = IRanges(new_starts, width = widths))
}

null_sv_ors   <- numeric(N_PERM)
null_admr_ors <- numeric(N_PERM)

pb_step <- max(1, N_PERM %/% 10)
for (i in seq_len(N_PERM)) {
  if (i %% pb_step == 0) message(sprintf("  Permutation %d/%d", i, N_PERM))
  perm_segdup      <- shuffle_segdup(segdup_gr, chrom_len)
  null_sv_ors[i]   <- tryCatch(calc_or(sv_gr,   perm_segdup), error=function(e) NA_real_)
  null_admr_ors[i] <- tryCatch(calc_or(admr_gr, perm_segdup), error=function(e) NA_real_)
}

null_sv_ors   <- null_sv_ors[!is.na(null_sv_ors)]
null_admr_ors <- null_admr_ors[!is.na(null_admr_ors)]

emp_p_sv   <- mean(null_sv_ors   >= obs_sv_or)
emp_p_admr <- mean(null_admr_ors >= obs_admr_or)

cat(sprintf("\nNull SV OR   : median=%.3f, 95%%ile=%.3f, 99%%ile=%.3f\n",
  median(null_sv_ors), quantile(null_sv_ors,0.95), quantile(null_sv_ors,0.99)))
cat(sprintf("Observed SV OR=%.3f → empirical p=%.4f\n", obs_sv_or, emp_p_sv))
cat(sprintf("\nNull aDMR OR : median=%.3f, 95%%ile=%.3f, 99%%ile=%.3f\n",
  median(null_admr_ors), quantile(null_admr_ors,0.95), quantile(null_admr_ors,0.99)))
cat(sprintf("Observed aDMR OR=%.3f → empirical p=%.4f\n", obs_admr_or, emp_p_admr))

# ── 4. Save results ───────────────────────────────────────────────────────────
res <- data.frame(
  analysis           = c("SV (C13)", "aDMR (C14)"),
  observed_OR        = c(obs_sv_or,   obs_admr_or),
  null_median_OR     = c(median(null_sv_ors),   median(null_admr_ors)),
  null_95pct_OR      = c(quantile(null_sv_ors,0.95), quantile(null_admr_ors,0.95)),
  null_99pct_OR      = c(quantile(null_sv_ors,0.99), quantile(null_admr_ors,0.99)),
  empirical_p        = c(emp_p_sv, emp_p_admr),
  n_permutations     = c(length(null_sv_ors), length(null_admr_ors)),
  obs_exceeds_99pct  = c(obs_sv_or   > quantile(null_sv_ors,0.99),
                          obs_admr_or > quantile(null_admr_ors,0.99))
)

fwrite(res, file.path(OUT_DIR, "a2_segdup_permutation.csv"))
message("Wrote: a2_segdup_permutation.csv")

# save null distributions for reproducibility
null_dist <- data.frame(
  permutation = c(seq_along(null_sv_ors), seq_along(null_admr_ors)),
  analysis    = c(rep("SV (C13)",   length(null_sv_ors)),
                  rep("aDMR (C14)", length(null_admr_ors))),
  null_OR     = c(null_sv_ors, null_admr_ors)
)
fwrite(null_dist, file.path(OUT_DIR, "a2_segdup_permutation_null_dist.csv"))

# ── 5. Figure ─────────────────────────────────────────────────────────────────
obs_lines <- data.frame(
  analysis   = c("SV (C13)", "aDMR (C14)"),
  observed   = c(obs_sv_or,   obs_admr_or),
  emp_p      = c(emp_p_sv,    emp_p_admr),
  label      = c(sprintf("Observed OR=%.2f\nemp. p=%.4f", obs_sv_or,   emp_p_sv),
                 sprintf("Observed OR=%.2f\nemp. p=%.4f", obs_admr_or, emp_p_admr))
)
null_dist$analysis <- factor(null_dist$analysis,
                              levels=c("SV (C13)","aDMR (C14)"))
obs_lines$analysis <- factor(obs_lines$analysis,
                              levels=c("SV (C13)","aDMR (C14)"))

pA <- ggplot(null_dist, aes(x=null_OR)) +
  geom_histogram(bins=50, fill="#6baed6", color="white", alpha=0.8) +
  geom_vline(data=obs_lines, aes(xintercept=observed), color="#d73027", linewidth=1.2) +
  geom_text(data=obs_lines, aes(x=observed, y=Inf, label=label),
            vjust=1.5, hjust=-0.1, color="#d73027", size=3.2) +
  facet_wrap(~analysis, scales="free") +
  labs(
    title    = "A2: Permuted SegDup null distribution",
    subtitle = sprintf("%d circular shift permutations (per-chromosome uniform offset)", N_PERM),
    x        = "Null OR (permuted SegDup vs observed SV/aDMR)",
    y        = "Count",
    caption  = "Red line = observed OR. emp. p = fraction of permutations ≥ observed OR."
  ) +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"),
        strip.background=element_blank(), strip.text=element_text(face="bold"))

# summary table panel
sum_df <- res |>
  mutate(sig = ifelse(empirical_p < 0.001, "***",
                ifelse(empirical_p < 0.01, "**",
                  ifelse(empirical_p < 0.05, "*", "ns"))),
         label = sprintf("OR=%.2f (emp.p=%.4f %s)", observed_OR, empirical_p, sig))

pB <- ggplot(sum_df, aes(x=observed_OR, y=analysis)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_vline(aes(xintercept=null_99pct_OR), linetype="dotted", color="steelblue") +
  geom_point(size=5, color="#d73027") +
  geom_text(aes(x=observed_OR+0.05, label=label), hjust=0, size=3.5) +
  scale_x_continuous(limits=c(0.8, max(sum_df$observed_OR)*1.5)) +
  labs(title="Observed OR vs null 99th percentile (dotted)",
       x="Odds Ratio", y=NULL,
       caption="Dotted = null 99th percentile") +
  theme_classic(base_size=12)

combined <- pA / pB + plot_layout(heights=c(3,1))
ggsave(file.path(FIG_DIR, "fig_a2_segdup_permutation.png"), combined,
       width=10, height=9, dpi=150)
message("Saved: fig_a2_segdup_permutation.png")

cat("\n=== A2 SUMMARY ===\n")
print(res)
cat("Done.\n")
