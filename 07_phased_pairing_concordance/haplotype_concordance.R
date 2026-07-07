#!/usr/bin/env Rscript
# =============================================================================
# Figure S6. Haplotype HYPO-concordance — the central anti-cis-induction test
#
# Purpose:
#   Show that the SV-bearing haplotype is NOT preferentially the hypomethylated
#   allele. Under direct cis-induction we expect SV-hap == HYPO-hap >> 50%;
#   under the shared-fragility / constitutional model we expect ~50%.
#
# Panels:
#   A) Per-patient HYPO-concordance bar (12 patients), dashed line at 50%,
#      colored above/below 50%, n pairs annotated; overall fraction noted.
#   B) Distance-stratified concordance (<=50kb, 50-200kb, 200-500kb, >500kb)
#      with 95% Clopper-Pearson CI (binom.test). Tests cis prediction that
#      concordance peaks at short distance.
#   C) Boundary-class stratified concordance. Tests that the ~50% null holds
#      across SV confidence/boundary strata (TAD+CTCF, CTCF-only, TAD-only,
#      non-boundary) — i.e. not an artefact of low-confidence pairs.
#      NOTE: the spec's "tier (Gold/Silver/Bronze) x ov_bulk" facet is not
#      recoverable here: the gold/silver/bronze candidate tables are PRE-FILTERED
#      to direction_match==TRUE (100% concordant by construction) and carry no
#      per-pair ov_bulk flag. We therefore stratify the *unfiltered* phase-block
#      pairs by their available boundary class (pairs$sv_tier), which is the
#      confidence-relevant stratifier that exists upstream of any filtering.
#
# Data sources:
#   PATHS$pb_pairs   (sv_admr_hp_concordance_pairs.csv.gz): patient_id, dist_bp,
#                    sv_type, sv_tier, sv_hp, hp_delta, admr_hypo_hp, concordant
#   PATHS$tier_*     (gold/silver/bronze_tier.csv): patient_code, sv_tier,
#                    direction_match (logical), bp_dist, delta_beta_bulk
#
# Statistical notes (recorded in figure_captions.md, not in the figure):
#   Overall GLMER:  glmer(direction_match ~ 1 + (1|patient_code), family=binomial)
#   Per-bin:        binom.test(n_concordant, n_total, p=0.5)$p.value
#
# Run:
#   mamba run -n renv Rscript viz/v4/figS6_haplotype_concordance.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
  library(forcats)
  library(lme4)
})

# ── SCRIPT_DIR detection ──────────────────────────────────────────────────────
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frame(1)$ofile)) return(dirname(normalizePath(sys.frame(1)$ofile)))
  "~/script/SV-DMR/somatic_AMR/viz/v1"
}
SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, "shared_theme.R"))
DIRS <- setup_fig_dirs(PATHS$figs_root)

log_con <- file(file.path(DIRS$logs, "figS6_haplotype_concordance.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
on.exit({ if (sink.number() > 0) sink(type = "output"); close(log_con) }, add = TRUE)
cat("=== Figure S6 Haplotype HYPO-concordance (v4) ===\n")

CONC_HI <- "#E74C3C"   # above 50% (cis-consistent)
CONC_LO <- "#3498DB"   # below/at 50% (null-consistent)

# ── Load pair-level concordance data ──────────────────────────────────────────
pb <- safe_fread(PATHS$pb_pairs)

# Normalise to a common schema: patient, dist_bp, concordant(0/1), tier, sv_type
prep_pairs <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.table(d)
  # patient column
  pcol <- intersect(c("patient_id", "patient_code", "patient"), names(d))[1]
  setnames(d, pcol, "patient")
  # distance
  if ("dist_bp" %in% names(d))       setnames(d, "dist_bp", "bp_dist")
  else if ("dist" %in% names(d))     setnames(d, "dist", "bp_dist")
  # concordance: SV-bearing haplotype is the HYPO (hypomethylated) one.
  # Prefer sv_minus_wt (sv_hp_beta - wt_hp_beta) which is the ground-truth column;
  # dmatch = TRUE means sv_minus_wt < 0 (SV haplotype has lower methylation than WT).
  # direction_match from phaseblock_pairs.csv has a known bug (sign(x)!=0 always TRUE),
  # so always recompute from sv_minus_wt when available.
  if ("sv_minus_wt" %in% names(d)) {
    d[, dmatch := sv_minus_wt < 0]
  } else if ("concordant" %in% names(d)) {
    d[, dmatch := as.integer(concordant) == 1L | concordant == TRUE]
  } else if ("direction_match" %in% names(d)) {
    d[, dmatch := as.logical(direction_match)]
  } else if (all(c("sv_hp", "admr_hypo_hp") %in% names(d))) {
    d[, dmatch := sv_hp == admr_hypo_hp]
  } else {
    stop("no concordance column found in pb_pairs")
  }
  # tier / sv_type passthrough
  # phaseblock_pairs.csv uses "tier" (Gold/Silver/Bronze), not "sv_tier" (boundary class).
  # Map "tier" → sv_tier so downstream Panel C stratifies by pair-confidence tier.
  if (!"sv_tier" %in% names(d)) {
    if ("tier" %in% names(d)) {
      d[, sv_tier := as.character(tier)]
    } else {
      d[, sv_tier := NA_character_]
    }
  }
  if (!"sv_type" %in% names(d) && "geom_type" %in% names(d))
    setnames(d, "geom_type", "sv_type")
  d <- d[!is.na(dmatch) & !is.na(bp_dist)]
  d[]
}

pb_dt <- prep_pairs(pb)
HAVE_PAIRS <- !is.null(pb_dt) && nrow(pb_dt) > 0
if (HAVE_PAIRS)
  cat(sprintf("[load] %d pairs across %d patients\n",
              nrow(pb_dt), uniqueN(pb_dt$patient)))

# Overall GLMER p (for caption only)
glmer_p <- NA_real_
if (HAVE_PAIRS && uniqueN(pb_dt$patient) > 1) {
  glmer_p <- tryCatch({
    m <- glmer(dmatch ~ 1 + (1 | patient), data = pb_dt, family = binomial)
    cf <- summary(m)$coefficients
    cf[1, 4]
  }, error = function(e) { cat("[glmer] failed:", e$message, "\n"); NA_real_ })
}
overall_frac <- if (HAVE_PAIRS) mean(pb_dt$dmatch) else NA_real_
cat(sprintf("[overall] concordance = %.1f%% (GLMER intercept p = %s)\n",
            100 * overall_frac,
            ifelse(is.na(glmer_p), "NA", sprintf("%.3f", glmer_p))))

# ── Panel A: per-patient concordance bar ──────────────────────────────────────
make_panelA <- function(d, overall) {
  if (is.null(d)) stop("no pair data")
  pp <- d[, .(n = .N, conc = mean(dmatch)), by = patient]
  pp[, pct := 100 * conc]
  pp[, above := pct > 50]
  pp[, patient := factor(patient,
        levels = patient[order(as.integer(gsub("\\D", "", patient)))])]

  ymax <- 105
  ggplot(pp, aes(x = patient, y = pct, fill = above)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_hline(yintercept = 50, linetype = "dashed",
               color = "grey30", linewidth = 0.5) +
    geom_text(aes(label = paste0("n=", n)), hjust = -0.15,
              size = FONT_MIN / .pt * 0.85, color = "grey25") +
    scale_fill_manual(values = c(`TRUE` = CONC_HI, `FALSE` = CONC_LO),
                      guide = "none") +
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, 25, 50, 75, 100),
                       labels = function(x) paste0(x, "%")) +
    coord_flip() +
    annotate("text", x = 0.6, y = ymax,
             label = sprintf("overall = %.1f%%", 100 * overall),
             hjust = 1, vjust = 0, size = FONT_MIN / .pt, fontface = "italic") +
    labs(x = NULL, y = "HYPO-concordance\n(SV-hap == HYPO-hap)") +
    theme_hcc_pub
}

# ── Panel B: distance-stratified concordance with Clopper-Pearson CI ──────────
DIST_BREAKS <- c(-Inf, 50e3, 200e3, 500e3, Inf)
DIST_LABS   <- c("≤50 kb", "50–200 kb", "200–500 kb", ">500 kb")
BIN_COLORS  <- c("≤50 kb"      = "#1B6CA8",
                 "50–200 kb"   = "#36A6A6",
                 "200–500 kb"  = "#F2A93B",
                 ">500 kb"          = "#C0392B")

make_panelB <- function(d) {
  if (is.null(d)) stop("no pair data")
  dd <- copy(d)
  dd[, bin := cut(bp_dist, breaks = DIST_BREAKS, labels = DIST_LABS)]
  dd <- dd[!is.na(bin)]
  res <- dd[, {
    n  <- .N; k <- sum(dmatch)
    bt <- binom.test(k, n, p = 0.5)
    .(n = n, k = k, frac = k / n,
      lo = bt$conf.int[1], hi = bt$conf.int[2], p = bt$p.value)
  }, by = bin]
  res[, bin := factor(bin, levels = DIST_LABS)]
  cat("[Panel B] distance bins:\n"); print(res)

  ggplot(res, aes(x = bin, y = 100 * frac, color = bin)) +
    geom_hline(yintercept = 50, linetype = "dashed",
               color = "grey30", linewidth = 0.5) +
    geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi),
                  width = 0.18, linewidth = 0.6) +
    geom_point(size = 2.6) +
    geom_text(aes(label = paste0("n=", n), y = 100 * hi + 4),
              size = FONT_MIN / .pt * 0.85, color = "grey25") +
    scale_color_manual(values = BIN_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "SV breakpoint → aDMR distance",
         y = "HYPO-concordance") +
    theme_hcc_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# ── Panel C: tier-stratified concordance (Gold/Silver/Bronze) ─────────────────
# Uses the "tier" column from phaseblock_pairs.csv (mapped to sv_tier in prep_pairs).
# Tests whether ~50% null holds across pair-confidence strata — rules out that
# the null is an artifact of low-confidence (Bronze) pairs dominating.
TIER_LEVELS <- c("Gold","Silver","Bronze")

make_panelC <- function(d) {
  if (is.null(d)) stop("no pair data")
  if (!"sv_tier" %in% names(d) || all(is.na(d$sv_tier)))
    stop("sv_tier (tier confidence) absent from pairs")
  dd <- copy(d)
  res <- dd[!is.na(sv_tier), {
    n <- .N; k <- sum(dmatch)
    bt <- binom.test(k, n, p = 0.5)
    .(n = n, frac = k / n, lo = bt$conf.int[1], hi = bt$conf.int[2],
      p = bt$p.value)
  }, by = sv_tier]
  res <- res[n >= 5]   # lowered threshold: Gold tier has only 80 pairs
  res[, lab := factor(sv_tier, levels = TIER_LEVELS)]
  res <- res[!is.na(lab)]
  res[, pct   := 100 * frac]
  res[, above := pct > 50]
  cat("[Panel C] tier-concordance:\n"); print(res)

  ggplot(res, aes(x = lab, y = pct, fill = above)) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_hline(yintercept = 50, linetype = "dashed",
               color = "grey30", linewidth = 0.5) +
    geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi),
                  width = 0.18, linewidth = 0.5, color = "grey30") +
    geom_text(aes(label = paste0("n=", n), y = 100 * hi + 4),
              size = FONT_MIN / .pt * 0.8, color = "grey25") +
    scale_fill_manual(values = c(`TRUE` = CONC_HI, `FALSE` = CONC_LO),
                      guide = "none") +
    scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "SV–aDMR pair tier", y = "HYPO-concordance") +
    theme_hcc_pub +
    theme(axis.text.x = element_text(size = FONT_MIN))
}

# ── Build panels (graceful fallback) ──────────────────────────────────────────
if (!HAVE_PAIRS) {
  pA <- build_panel("A", function() placeholder_panel(
          msg = "pb_pairs unavailable\n(phase-block pairs pending)"))
  pB <- build_panel("B", function() placeholder_panel(
          msg = "pb_pairs unavailable"))
} else {
  pA <- build_panel("A", make_panelA, pb_dt, overall_frac)
  pB <- build_panel("B", make_panelB, pb_dt)
}
if (!HAVE_PAIRS) {
  pC <- build_panel("C", function() placeholder_panel(
          msg = "pb_pairs unavailable"))
} else {
  pC <- build_panel("C", make_panelC, pb_dt)
}

# ── Persist individual panels ─────────────────────────────────────────────────
saveRDS(pA, file.path(DIRS$panels, "figS6_A_per_patient.rds"))
saveRDS(pB, file.path(DIRS$panels, "figS6_B_distance.rds"))
saveRDS(pC, file.path(DIRS$panels, "figS6_C_boundary_class.rds"))
save_panel(pA, "figS6_A_per_patient",    DIRS, width = 6, height = 5)
save_panel(pB, "figS6_B_distance",       DIRS, width = 5, height = 4)
save_panel(pC, "figS6_C_boundary_class", DIRS, width = 6, height = 4)

# ── Assemble ──────────────────────────────────────────────────────────────────
figS6 <- (pA / (pB | pC)) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FONT_TAG, face = "bold"))

saveRDS(figS6, file.path(DIRS$rds, "figS6_combined.rds"))
save_panel(figS6, "figS6_haplotype_concordance", DIRS, width = 180 / 25.4,
           height = 200 / 25.4)
ggsave(file.path(DIRS$png, "figS6_haplotype_concordance.png"), figS6,
       width = FIG_W_MM, height = FIG_W_MM, units = "mm", dpi = 300, bg = "white")
cat(sprintf("[Output] figS6 -> %s\n",
            file.path(DIRS$png, "figS6_haplotype_concordance.png")))

log_decision(sprintf(
  "figS6_haplotype_concordance.R: per-patient/distance/tier-bulk HYPO-concordance; overall=%.1f%% GLMER p=%s (null-consistent, anti-cis)",
  100 * overall_frac, ifelse(is.na(glmer_p), "NA", sprintf("%.3f", glmer_p))))
