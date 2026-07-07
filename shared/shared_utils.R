# shared/shared_utils.R ========================================================
# Shared constants, palettes, and utility functions.
# Source at the top of each pipeline script (after the library block), finding
# this repo's shared/ directory regardless of the calling script's own depth:
#   REPO_ROOT <- local({
#     d <- dirname(normalizePath(sub("--file=", "",
#       grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
#     while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
#     d
#   })
#   source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

# Paths ========================================================================
# Locate this file's own path (works when source()'d, independent of caller)
# and load repo-root .env so all path variables below come from it.
.this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
.repo_root <- if (!is.na(.this_file)) dirname(dirname(.this_file)) else normalizePath(".")
dotenv::load_dot_env(file.path(.repo_root, ".env"))

PATIENT_MAP_PATH <- Sys.getenv("PATIENT_CODE_MAP")
BLACKLIST_PATH   <- file.path(Sys.getenv("REFERENCE_DIR"), "ensembl/ENCFF356LFX_ENCODE-blacklist.bed.gz")

# CNV quality thresholds =======================================================
CN_NORMAL_RANGE     <- c(1.7, 2.3)
MIN_DEPTH_WINDOW    <- 5L
MIN_BAF_COUNT       <- 5L
MIN_SVLEN_FOR_CN    <- 1000L
VAF_CONCORDANCE_THR <- 0.3

# GRCh38 chromosome lengths ====================================================
CHROM_LENS <- c(
  chr1 = 248956422L, chr2 = 242193529L, chr3 = 198295559L, chr4 = 190214555L,
  chr5 = 181538259L, chr6 = 170805979L, chr7 = 159345973L, chr8 = 145138636L,
  chr9 = 138394717L, chr10 = 133797422L, chr11 = 135086622L, chr12 = 133275309L,
  chr13 = 114364328L, chr14 = 107043718L, chr15 = 101991189L, chr16 = 90338345L,
  chr17 = 83257441L,  chr18 = 80373285L,  chr19 = 58617616L,  chr20 = 64444167L,
  chr21 = 46709983L,  chr22 = 50818468L,  chrX = 156040895L
)

# Common ggplot2 theme =========================================================
theme_hcc <- ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title         = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle      = ggplot2::element_text(color = "grey50", size = 10),
    axis.title         = ggplot2::element_text(size = 10),
    axis.text          = ggplot2::element_text(size = 9),
    legend.position    = "bottom",
    legend.text        = ggplot2::element_text(size = 9),
    legend.title       = ggplot2::element_text(size = 9, face = "bold"),
    strip.background   = ggplot2::element_rect(fill = "grey95", color = NA),
    strip.text         = ggplot2::element_text(face = "bold", size = 10),
    panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3)
  )

# Color palettes ===============================================================

# SV CNV class (matches $cnv_class column values: copy_neutral/gaining/losing/insertion/COM)
SV_COLORS <- c(
  copy_neutral = "#3B8BD4",
  copy_gaining = "#E24B4A",
  copy_losing  = "#BA7517",
  insertion    = "#7F77DD",
  COM          = "#888780"
)

# Per geom_type (DEL/DUP/INV/TRA/INS/COM)
SV_TYPE_LEVELS <- c("DEL", "DUP", "INV", "TRA", "INS", "COM")
SV_TYPE_COLORS <- c(
  DEL = "#BA7517", DUP = "#E24B4A",
  INV = "#3B8BD4", TRA = "#7F77DD", INS = "#1D9E75", COM = "#888780"
)

# TAD/CTCF stratification (matches $stratification column values from 02_sv_annotation_stratify.R)
STRAT_LEVELS <- c("TAD+CTCF disrupting", "CTCF-only", "TAD-only",
                  "Copy-neutral", "Non-boundary")
STRAT_COLORS <- c(
  "TAD+CTCF disrupting" = "#C0392B",
  "CTCF-only"           = "#E67E22",
  "TAD-only"            = "#F1C40F",
  "Copy-neutral"        = "#3B8BD4",
  "Non-boundary"        = "#95A5A6"
)

# Recoded (underscore) tier labels — for use after applying TIER_RECODE
SV_TIER_LEVELS <- c("HBV_associated", "TAD_CTCF", "TAD_only",
                    "CTCF_only", "copy_neutral", "non_boundary")
SV_TIER_COLORS <- c(
  TAD_CTCF     = "#7F77DD",
  TAD_only     = "#3B8BD4",
  CTCF_only    = "#1D9E75",
  copy_neutral = "#BA7517",
  non_boundary = "#888780"
)

# Recode mapping: stratification labels → SV_TIER_COLORS names
TIER_RECODE <- c(
  "TAD+CTCF disrupting" = "TAD_CTCF",
  "TAD-only"            = "TAD_only",
  "CTCF-only"           = "CTCF_only",
  "Copy-neutral"        = "copy_neutral",
  "Non-boundary"        = "non_boundary",
  "HBV-associated"      = "HBV_associated"
)

# Haplotype and methylation direction
HP_COLORS   <- c(HP1 = "#3B8BD4", HP2 = "#E24B4A")
METH_COLORS <- c(HYPER = "#C0392B", HYPO = "#2980B9")

# Utility functions ============================================================

anonym_sample <- function(df, sample_col = "sample") {
  if (!exists("patient_code_rule", envir = .GlobalEnv)) {
    patient_code_rule <- data.table::fread(PATIENT_MAP_PATH)
    assign("patient_code_rule", patient_code_rule, envir = .GlobalEnv)
  }
  df |> dplyr::left_join(patient_code_rule, by = setNames("Samples_ID", sample_col))
}

filter_dmr_region <- function(dmr_gr, meth_diff_cutoff = 0.2, ncg_cutoff = 5) {
  if (!is.null(mcols(dmr_gr)$sample))
    message("Processing Sample:", mcols(dmr_gr)$sample[1])
  message("Raw DMR: ", length(unique(dmr_gr)), appendLF = FALSE)
  if (!exists("blacklist", envir = .GlobalEnv)) {
    bl_df <- data.table::fread(BLACKLIST_PATH, header = FALSE)
    blacklist <- GenomicRanges::makeGRangesFromDataFrame(
      bl_df, seqnames.field = "V1", start.field = "V2", end.field = "V3"
    )
    rm(bl_df)
    assign("blacklist", blacklist, envir = .GlobalEnv)
  }
  dmr_gr <- dmr_gr[!(dmr_gr %over% blacklist)]
  message("  DMR number after blacklist filtering: ", length(unique(dmr_gr)), appendLF = FALSE)
  dmr_gr <- dmr_gr[
    abs(mcols(dmr_gr)$diff.Methy) >= meth_diff_cutoff &
      mcols(dmr_gr)$nCG >= ncg_cutoff
  ]
  message("  DMR after QC: ", length(unique(dmr_gr)))
  dmr_gr
}

reciprocal_overlap_hits <- function(gr1, gr2, min_pct = 0.3, min_bp = 1L) {
  hits <- GenomicRanges::findOverlaps(gr1, gr2, minoverlap = min_bp)
  if (length(hits) == 0) return(hits)
  ol   <- GenomicRanges::pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
  pct1 <- width(ol) / width(gr1[queryHits(hits)])
  pct2 <- width(ol) / width(gr2[subjectHits(hits)])
  hits[pct1 >= min_pct & pct2 >= min_pct]
}

safe_fread <- function(path, ...) {
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s — returning NULL", path))
    return(NULL)
  }
  data.table::fread(path, ...)
}
