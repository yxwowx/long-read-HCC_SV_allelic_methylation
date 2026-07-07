#!/usr/bin/env Rscript
# 03_segdup_permutation.R
# Circular permutation null for SegDup enrichment (C24, §4.2 fix)
# Per-chromosome circular shift with fast findInterval overlap counting.
# Non-GRanges inner loop: ~0.1s/iter vs ~5s for GRanges %over%, so 1000 iters ≈ 2 min.
# Output: SV_aDMR/segdup_permutation_null.csv

suppressPackageStartupMessages({
  library(data.table)
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

DATA_ROOT    <- Sys.getenv("HCC_DATA_DIR")
DMR_SVS_ROOT <- file.path(DATA_ROOT, "DMR_SVs")
OUT_ROOT     <- file.path(DATA_ROOT, "SV_aDMR")

REF_ROOT   <- Sys.getenv("REFERENCE_DIR")
SEGDUP_BED <- file.path(REF_ROOT, "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")

N_PERM      <- 1000L
MAIN_CHROMS <- paste0("chr", c(1:22, "X"))
set.seed(2026)

# Load SegDup; reduce to non-overlapping intervals per chrom ===================
cat("[0] Loading SegDup...\n")
segdup_gr  <- import.bed(SEGDUP_BED) |> keepStandardChromosomes(pruning.mode = "coarse")
segdup_gr  <- segdup_gr[seqnames(segdup_gr) %in% MAIN_CHROMS]
segdup_red <- reduce(segdup_gr)  # non-overlapping, for fast findInterval overlap

# Per-chromosome sorted SegDup lookup (plain integer vectors, no GRanges)
seg_by_chr <- lapply(MAIN_CHROMS, function(chr) {
  g <- segdup_red[seqnames(segdup_red) == chr]
  if (length(g) == 0) return(list(start = integer(0), end = integer(0)))
  list(start = start(g), end = end(g))
})
names(seg_by_chr) <- MAIN_CHROMS

# Observed genome SegDup fraction (for approximate OR)
genome_bp <- sum(as.numeric(CHROM_LENS[MAIN_CHROMS]))
bg_frac   <- sum(as.numeric(width(segdup_red))) / genome_bp

# Fast interval-in-sorted-set overlap count ====================================
# Correct for non-overlapping (reduced) intervals:
# If the last SegDup with start ≤ query_end doesn't end ≥ query_start,
# no earlier SegDup can overlap either (proved by the non-overlap property).
count_in_segdup <- function(q_starts, q_ends, seg_starts, seg_ends) {
  if (length(q_starts) == 0 || length(seg_starts) == 0) return(0L)
  idx <- findInterval(q_ends, seg_starts)   # last seg with start <= q_end
  hits <- idx >= 1L & seg_ends[pmax(1L, idx)] >= q_starts
  sum(hits)
}

# Load case sets ===============================================================
cat("[1] Loading somatic aDMR...\n")
admr_dt  <- fread(file.path(OUT_ROOT, "somatic_admr_annotated.csv.gz"))
admr_pos <- admr_dt[seqnames %in% MAIN_CHROMS, .(seqnames, start, end)]
# Pre-split by chromosome (plain integer vectors)
admr_by_chr <- lapply(MAIN_CHROMS, function(chr) {
  dt <- admr_pos[seqnames == chr]
  list(start = dt$start, end = dt$end,
       width = dt$end - dt$start,
       chr_len = CHROM_LENS[[chr]])
})
names(admr_by_chr) <- MAIN_CHROMS
cat(sprintf("  %d somatic aDMR\n", nrow(admr_pos)))

cat("[2] Loading SV breakpoints...\n")
sv_dt  <- fread(file.path(DMR_SVS_ROOT, "sv_tad_ctcf_annotation.csv.gz"))
sv_pos <- sv_dt[seqnames %in% MAIN_CHROMS, .(seqnames, start, end = start + 1L)]
sv_by_chr <- lapply(MAIN_CHROMS, function(chr) {
  dt <- sv_pos[seqnames == chr]
  list(start = dt$start, end = dt$end,
       width = rep(1L, nrow(dt)),
       chr_len = CHROM_LENS[[chr]])
})
names(sv_by_chr) <- MAIN_CHROMS
cat(sprintf("  %d SV breakpoints\n", nrow(sv_pos)))

# Observed overlaps ============================================================
obs_admr <- sum(mapply(function(ac, sc) {
  count_in_segdup(ac$start, ac$end, sc$start, sc$end)
}, admr_by_chr, seg_by_chr))

obs_sv <- sum(mapply(function(sc_sv, sc) {
  count_in_segdup(sc_sv$start, sc_sv$end, sc$start, sc$end)
}, sv_by_chr, seg_by_chr))

cat(sprintf("Observed: aDMR=%d (%.2f%%), SV=%d (%.2f%%)\n",
            obs_admr, obs_admr/nrow(admr_pos)*100,
            obs_sv,   obs_sv/nrow(sv_pos)*100))

# Circular permutation (per-chrom uniform random shift) ========================
cat(sprintf("[3] Running %d circular permutations...\n", N_PERM))

one_perm <- function(pos_by_chr, seg_by_chr) {
  total <- 0L
  for (chr in MAIN_CHROMS) {
    ac <- pos_by_chr[[chr]]
    sc <- seg_by_chr[[chr]]
    if (length(ac$start) == 0 || length(sc$start) == 0) next
    shift_amt  <- sample.int(ac$chr_len - 1L, 1L)
    new_starts <- (ac$start - 1L + shift_amt) %% ac$chr_len + 1L
    new_ends   <- pmin(new_starts + ac$width, ac$chr_len)
    total      <- total + count_in_segdup(new_starts, new_ends, sc$start, sc$end)
  }
  total
}

cat("  Permuting somatic aDMR...\n")
t0 <- proc.time()
null_admr <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  if (i %% 100 == 0) cat(sprintf("  perm %d/%d (%.0fs elapsed)\n", i, N_PERM,
                                   (proc.time() - t0)[[3]]))
  null_admr[i] <- one_perm(admr_by_chr, seg_by_chr)
}

cat("  Permuting SV breakpoints...\n")
null_sv <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  null_sv[i] <- one_perm(sv_by_chr, seg_by_chr)
}

# Results ======================================================================
n_admr <- nrow(admr_pos)
n_sv   <- nrow(sv_pos)

emp_p_admr   <- mean(null_admr >= obs_admr)
emp_p_sv     <- mean(null_sv   >= obs_sv)

obs_or_admr  <- (obs_admr / n_admr)    / bg_frac
obs_or_sv    <- (obs_sv   / n_sv)      / bg_frac
null_or_admr <- (null_admr / n_admr)   / bg_frac
null_or_sv   <- (null_sv   / n_sv)     / bg_frac

cat(sprintf("\n=== Circular Permutation Results ===\n"))
cat(sprintf("  Somatic aDMR: obs_or=%.3f, null_mean_or=%.3f±%.3f, emp_p=%.4f\n",
            obs_or_admr, mean(null_or_admr), sd(null_or_admr), emp_p_admr))
cat(sprintf("  SV bp:        obs_or=%.3f, null_mean_or=%.3f±%.3f, emp_p=%.4f\n",
            obs_or_sv, mean(null_or_sv), sd(null_or_sv), emp_p_sv))

# Save =========================================================================
null_dt <- rbind(
  data.table(set = "somatic_admr", perm = seq_len(N_PERM),
             n_overlap = null_admr, or_approx = null_or_admr),
  data.table(set = "sv_breakpoint", perm = seq_len(N_PERM),
             n_overlap = null_sv,   or_approx = null_or_sv)
)

summ_dt <- data.table(
  set           = c("somatic_admr", "sv_breakpoint"),
  obs_overlap   = c(obs_admr, obs_sv),
  obs_or_approx = c(obs_or_admr, obs_or_sv),
  null_mean     = c(mean(null_admr), mean(null_sv)),
  null_sd       = c(sd(null_admr), sd(null_sv)),
  null_p99      = c(quantile(null_admr, 0.99), quantile(null_sv, 0.99)),
  null_or_mean  = c(mean(null_or_admr), mean(null_or_sv)),
  empirical_p   = c(emp_p_admr, emp_p_sv),
  n_perm        = c(N_PERM, N_PERM),
  n_case        = c(n_admr, n_sv)
)

fwrite(null_dt,  file.path(OUT_ROOT, "segdup_permutation_null.csv"))
fwrite(summ_dt,  file.path(OUT_ROOT, "segdup_permutation_summary.csv"))

cat("\nSaved: segdup_permutation_null.csv, segdup_permutation_summary.csv\n")
print(summ_dt)

cat("=== Done: 03_segdup_permutation.R ===\n")
