#!/usr/bin/env Rscript
# Figure S6. RepeatMasker Stratification (Informative Negative Control)
#
# Addresses the alternative hypothesis: SegDup aDMR co-enrichment reflects
# repeat-element hypomethylation rather than structural fragility.
#
# Panels:
#   A) OR by repeat-density quartile + Cochran-Q test (p=0.278, I²=22%)
#   B) OR by dominant repeat class (>=2 aDMR loci)
#   C) LINE proportion at SegDup aDMRs vs background (depleted, p<0.0001)
#
# NOTE: this script only re-exports a cached combined-figure RDS (the fast
# path below). The from-scratch computation that produces panels A-C lived in
# a separate legacy visualization script that is not part of this repo; if no
# cached RDS is found, a placeholder figure is written instead of guessing at
# the analysis.
#
# Run:
#   mamba run -n renv Rscript 05_fragility_enrichment/repeatmask.R

suppressPackageStartupMessages({
  library(data.table)
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

# Minimal local helper: this script's figure-export plumbing originally came
# from a private viz/shared_theme.R framework that is not included in this repo.
setup_fig_dirs <- function(figs_root) {
  dirs <- file.path(figs_root, c("panels", "rds", "png", "logs"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(panels = file.path(figs_root, "panels"),
       rds    = file.path(figs_root, "rds"),
       png    = file.path(figs_root, "png"),
       logs   = file.path(figs_root, "logs"))
}

DMR_SVS_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs")
PATHS <- list(figs_root = file.path(DMR_SVS_DIR, "figs/v4"))

DIRS   <- setup_fig_dirs(PATHS$figs_root)
V2_RDS <- file.path(DMR_SVS_DIR, "figs/v2/rds/figS11_repeatmask_stratification.rds")

log_con <- file(file.path(DIRS$logs, "figS6_repeatmask.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
on.exit({ if (sink.number() > 0) sink(type = "output"); close(log_con) }, add = TRUE)
cat("=== Figure S6 RepeatMasker Stratification (v4) ===\n")

# Fast path: use cached v2 combined figure =====================================
if (file.exists(V2_RDS)) {
  cat("[figS6] Loading from v2 cached RDS:", V2_RDS, "\n")
  raw_obj <- readRDS(V2_RDS)
  # v2 RDS may be a list with $fig element, or a direct patchwork/ggplot
  fig_obj <- if (is.list(raw_obj) && !inherits(raw_obj, "gg") && "fig" %in% names(raw_obj))
               raw_obj$fig else raw_obj
  if (!inherits(fig_obj, c("gg","patchwork"))) {
    cat("[figS6] WARNING: cached object is not gg/patchwork — writing placeholder\n")
    fig_obj <- ggplot() + annotate("text", x=0.5, y=0.5,
      label="figS6: v2 cache loaded but not a ggplot/patchwork object", size=4) + theme_void()
  }
  # Re-annotate title
  fig_obj <- fig_obj + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
  saveRDS(fig_obj, file.path(DIRS$rds, "figS6_combined.rds"))
  ggsave(file.path(DIRS$png, "figS6.png"), fig_obj,
         width = 13, height = 9, dpi = 300, bg = "white")
  cat(sprintf("[Output] figS6.png (from v2 cache) → %s\n", file.path(DIRS$png, "figS6.png")))
  quit(save = "no")
}

# No cache found: the from-scratch RepeatMasker stratification computation
# (OR by quartile, repeat-class breakdown, LINE depletion test) is not
# included in this repo — write a placeholder instead of guessing at it. =======
cat("[figS6] v2 RDS not found; from-scratch computation is not part of this repo\n")
placeholder <- ggplot() +
  annotate("text", x=0.5, y=0.5,
           label="figS6: cached RDS not found; from-scratch RepeatMasker\nstratification computation is not included in this repo",
           size=4, color="grey40") + theme_void()
ggsave(file.path(DIRS$png, "figS6.png"), placeholder,
       width=13, height=9, dpi=300, bg="white")
cat(sprintf("[Output] figS6.png (placeholder) → %s\n", file.path(DIRS$png, "figS6.png")))
