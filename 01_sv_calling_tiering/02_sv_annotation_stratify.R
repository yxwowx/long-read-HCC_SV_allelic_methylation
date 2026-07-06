#!/usr/bin/env Rscript
# 02_sv_annotation_stratify.R
# SV breakpoint annotation and stratification pipeline.
# Produces:
#   sv_tad_ctcf_annotation.csv.gz    -- per-breakpoint TAD/CTCF/CNV stratification
#   sv_bp.tier_nearest_gene.csv.gz   -- wide: one row per sv_id, nearest TSS
#   sv_bp.tier_with_gene.csv.gz      -- long: one row per sv_id × gene (±50 kb)
#
# Usage:
#   mamba run -n renv Rscript pipeline/02_sv_annotation_stratify.R \
#     [--outdir /node200data/kachungk/hcc_data/DMR_SVs] \
#     [--tad_bed <path>] [--ctcf_bed <path>] \
#     [--prom_rds <path>] [--gene_rds <path>] \
#     2>&1 | tee logs/02_sv_annotation_stratify.log

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(StructuralVariantAnnotation)
  library(VariantAnnotation)
  library(optparse)
  library(tidyr)
})
source(file.path(dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
)), "shared_utils.R"))

# ── Constants ─────────────────────────────────────────────────────────────────
HBV_JUNCTION_FLANK     <- 500L
COMPLEX_DETAILED_TYPES <- c(
  "Templated_ins", "Templated_ins_inv", "inv_tra", "complex_inv"
)

# ── Functions ─────────────────────────────────────────────────────────────────
simple_svtype <- function(gr) {
  partner_id  <- mcols(gr)$partner
  partner_idx <- match(partner_id, names(gr))

  valid_idx <- which(!is.na(partner_idx))
  dropped   <- names(gr)[is.na(partner_idx)]
  if (length(dropped) > 0)
    message("Dropped ", length(dropped), " unpaired bp: ",
            paste(head(dropped, 3), collapse = ", "))

  raw_svtype <- as.character(mcols(gr)$svtype)
  geom_type  <- raw_svtype

  det_type <- if ("DETAILED_TYPE" %in% names(mcols(gr)))
    as.character(mcols(gr)$DETAILED_TYPE)
  else
    rep(NA_character_, length(gr))

  bnd_valid <- valid_idx[raw_svtype[valid_idx] == "BND"]
  if (length(bnd_valid) > 0) {
    gr_bnd     <- gr[bnd_valid]
    gr_partner <- gr[partner_idx[bnd_valid]]
    same_chr <- as.character(seqnames(gr_bnd)) ==
                as.character(seqnames(gr_partner))
    str_v    <- as.character(strand(gr_bnd))
    str_p    <- as.character(strand(gr_partner))

    geom_type[bnd_valid][!same_chr]                       <- "TRA"
    geom_type[bnd_valid][same_chr & (str_v == str_p)]     <- "INV"
    geom_type[bnd_valid][same_chr & str_v == "+" & str_p == "-"] <- "DEL"
    geom_type[bnd_valid][same_chr & str_v == "-" & str_p == "+"] <- "DUP"

    geom_type[bnd_valid] <- dplyr::case_when(
      det_type[bnd_valid] %in% COMPLEX_DETAILED_TYPES           ~ "COM",
      geom_type[bnd_valid] == "INV" &
        det_type[bnd_valid] %in% c("foldback", "dup_inv_segment") ~ "DUP",
      geom_type[bnd_valid] == "INV" &
        det_type[bnd_valid] == "reciprocal_inv_del"              ~ "DEL",
      det_type[bnd_valid] == "Reciprocal_tra"                    ~ "TRA",
      TRUE ~ geom_type[bnd_valid]
    )
  }
  mcols(gr)$geom_type <- geom_type
  gr
}

get_purple_cn_at_sv <- function(sv_bp, seg_gr) {
  n_sv <- length(sv_bp)
  out  <- list(
    tumor_cn    = rep(NA_real_,      n_sv),
    minor_cn    = rep(NA_real_,      n_sv),
    is_reliable = rep(FALSE,         n_sv),
    germ_status = rep(NA_character_, n_sv)
  )
  if (length(seg_gr) == 0) return(out)

  seg_ok <- seg_gr[
    !is.na(seg_gr$depthWindowCount) &
      seg_gr$depthWindowCount >= MIN_DEPTH_WINDOW &
      !is.na(seg_gr$bafCount) &
      seg_gr$bafCount >= MIN_BAF_COUNT &
      !(seg_gr$germlineStatus %in%
          c("HET_DELETION", "HOM_DELETION", "AMPLIFICATION"))
  ]
  if (length(seg_ok) == 0) return(out)

  hits   <- findOverlaps(sv_bp, seg_ok, select = "first")
  hit_ok <- !is.na(hits)
  if (!any(hit_ok)) return(out)

  out$tumor_cn[hit_ok]    <- seg_ok$copyNumber[hits[hit_ok]]
  out$minor_cn[hit_ok]    <- seg_ok$minorAlleleCopyNumber[hits[hit_ok]]
  out$germ_status[hit_ok] <- seg_ok$germlineStatus[hits[hit_ok]]
  out$is_reliable[hit_ok] <- TRUE
  out
}

cn_to_cnv_class <- function(tumor_cn, minor_cn = NULL,
                             cn_range = CN_NORMAL_RANGE) {
  dplyr::case_when(
    is.na(tumor_cn)                                     ~ NA_character_,
    tumor_cn == 0                                       ~ "copy_losing",
    !is.na(minor_cn) & minor_cn < 0.2 &
      tumor_cn >= cn_range[1]                           ~ "copy_neutral",
    tumor_cn >= cn_range[1] & tumor_cn <= cn_range[2]   ~ "copy_neutral",
    tumor_cn > cn_range[2]                              ~ "copy_gaining",
    tumor_cn < cn_range[1]                              ~ "copy_losing",
    TRUE                                                ~ NA_character_
  )
}

classify_sv_cnv_final <- function(sv_gr, seg_gr = NULL) {
  n <- length(sv_gr)
  if (n == 0) return(sv_gr)
  if (!"geom_type" %in% names(mcols(sv_gr)))
    stop("sv_gr에 $geom_type 없음. simple_svtype() 먼저 실행하세요.")

  sv_mcols  <- mcols(sv_gr)
  svlen_col <- if ("svlen" %in% names(sv_mcols)) "svlen" else
               if ("svLen" %in% names(sv_mcols)) "svLen" else NULL
  svlen     <- if (!is.null(svlen_col)) as.numeric(sv_mcols[[svlen_col]]) else
               rep(NA_real_, n)
  if (length(svlen) != n) svlen <- rep(NA_real_, n)

  vaf_obs <- if ("VAF" %in% names(mcols(sv_gr)))
    as.numeric(mcols(sv_gr)$VAF) else rep(NA_real_, n)

  sv_bp  <- GRanges(seqnames(sv_gr), IRanges(start(sv_gr), start(sv_gr)))
  purple <- get_purple_cn_at_sv(sv_bp, seg_gr)

  if (length(purple$tumor_cn)    != n) purple$tumor_cn    <- rep(NA_real_, n)
  if (length(purple$minor_cn)    != n) purple$minor_cn    <- rep(NA_real_, n)
  if (length(purple$is_reliable) != n) purple$is_reliable <- rep(FALSE,   n)
  if (length(purple$germ_status) != n) purple$germ_status <- rep(NA_character_, n)

  purple$is_reliable <- purple$is_reliable & !is.na(svlen) & svlen >= MIN_SVLEN_FOR_CN
  purple_class <- ifelse(purple$is_reliable,
    cn_to_cnv_class(purple$tumor_cn, purple$minor_cn), NA_character_)

  vaf_concordant <- mapply(function(vo, tcn, thr = VAF_CONCORDANCE_THR) {
    if (is.na(vo) || is.na(tcn) || tcn <= 0) return(NA)
    abs(vo - 1 / tcn) <= thr
  }, vaf_obs, purple$tumor_cn, SIMPLIFY = TRUE)

  geom <- as.character(mcols(sv_gr)$geom_type)
  cnv_class <- dplyr::case_when(
    geom == "DEL"                                        ~ "copy_losing",
    geom == "DUP"                                        ~ "copy_gaining",
    geom == "INS"                                        ~ "insertion",
    geom %in% c("TRA", "INV") & purple$is_reliable &
      !is.na(purple_class) & purple_class != "copy_neutral" ~ purple_class,
    geom %in% c("TRA", "INV")                           ~ "copy_neutral",
    geom == "COM" & purple$is_reliable &
      !is.na(purple_class)                               ~ purple_class,
    geom == "COM"                                        ~ "COM",
    purple$is_reliable & !is.na(purple_class)            ~ purple_class,
    TRUE                                                 ~ "other"
  )

  cnv_source <- dplyr::case_when(
    geom %in% c("DEL", "DUP", "INS")                    ~ "geom_type",
    geom %in% c("TRA", "INV") &
      cnv_class != "copy_neutral"                        ~ "PURPLE_override",
    geom %in% c("TRA", "INV")                           ~ "geom_type",
    geom == "COM" & purple$is_reliable &
      !is.na(purple_class)                               ~ "PURPLE_CN",
    geom == "COM"                                        ~ "COM_no_cn",
    TRUE                                                 ~ "other"
  )

  mcols(sv_gr)$cnv_class        <- cnv_class
  mcols(sv_gr)$cnv_class_source <- cnv_source
  mcols(sv_gr)$tumor_cn         <- purple$tumor_cn
  mcols(sv_gr)$cn_reliable      <- purple$is_reliable
  mcols(sv_gr)$vaf_concordant   <- vaf_concordant
  mcols(sv_gr)$germline_status  <- purple$germ_status
  sv_gr
}

# ChromHMM 15-state active regulatory states (Roadmap Epigenomics)
CHROMHMM_ACTIVE_STATES <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "6_EnhG", "7_Enh")

annotate_sv_cre <- function(sv_gr, cre_liver_gr, enhancer_gr = NULL, buffer = 4000L) {
  n      <- length(sv_gr)
  sv_buf <- resize(sv_gr, width = buffer * 2L, fix = "center")

  ov_liver <- overlapsAny(sv_buf, cre_liver_gr)
  ov_enh   <- if (!is.null(enhancer_gr) && length(enhancer_gr) > 0)
                 overlapsAny(sv_buf, enhancer_gr)
               else rep(FALSE, n)

  cre_any  <- ov_liver | ov_enh

  cre_type <- dplyr::case_when(
    ov_liver ~ "chromhmm_liver_active",
    ov_enh   ~ "genehancer",
    TRUE     ~ NA_character_
  )

  # Sub-tier: meaningful only for TAD+CTCF disrupting SVs; NA elsewhere
  strat <- as.character(mcols(sv_gr)$stratification)
  tad_ctcf_cre_subtier <- dplyr::case_when(
    strat == "TAD+CTCF disrupting" &  cre_any ~ "TAD+CTCF+CRE",
    strat == "TAD+CTCF disrupting" & !cre_any ~ "TAD+CTCF-noCRE",
    TRUE                                       ~ NA_character_
  )

  mcols(sv_gr)$cre_any_overlap       <- cre_any
  mcols(sv_gr)$cre_chromhmm_liver    <- ov_liver
  mcols(sv_gr)$cre_type              <- cre_type
  mcols(sv_gr)$tad_ctcf_cre_subtier <- tad_ctcf_cre_subtier
  sv_gr
}


sv_stratification <- function(sv_gr, tad_body_gr, tad_bound_gr, ctcf_gr, hbv_genome_len) {
  n <- length(sv_gr)
  message("TAD/CTCF overlap annotation for ", unique(mcols(sv_gr)$sample)[1],
          " (", n, " breakpoints)...")

  partner_idx <- match(mcols(sv_gr)$partner, names(sv_gr))
  geom        <- as.character(mcols(sv_gr)$geom_type)
  sv_buffer   <- resize(sv_gr, width = 4000L, fix = "center")

  # ── Shared index sets ──────────────────────────────────────────────────────
  paired_idx     <- which(!is.na(partner_idx))
  intra_geom     <- c("DEL", "DUP", "INV", "COM")
  same_chr_mask  <- logical(n)
  if (length(paired_idx) > 0)
    same_chr_mask[paired_idx] <- as.character(seqnames(sv_gr))[paired_idx] ==
      as.character(seqnames(sv_gr))[partner_idx[paired_idx]]
  intra_idx <- which(!is.na(partner_idx) & geom %in% intra_geom & same_chr_mask)

  # ── Pre-compute which boundary regions contain CTCF peaks ─────────────────
  bound_has_ctcf <- overlapsAny(tad_bound_gr, ctcf_gr)

  # ── Condition A: BP1 and BP2 belong to different TAD bodies ───────────────
  # Priority: strongest evidence — SV must span ≥1 boundary
  tad_A <- logical(n)
  if (length(intra_idx) > 0) {
    b1 <- sv_gr[intra_idx]
    b2 <- sv_gr[partner_idx[intra_idx]]
    tad_b1 <- findOverlaps(b1, tad_body_gr, select = "first")
    tad_b2 <- findOverlaps(b2, tad_body_gr, select = "first")
    cond_A  <- !is.na(tad_b1) & !is.na(tad_b2) & tad_b1 != tad_b2
    tad_A[intra_idx[cond_A]]              <- TRUE
    tad_A[partner_idx[intra_idx[cond_A]]] <- TRUE
  }

  # ── Condition B: BP directly overlaps a TAD boundary gap region ───────────
  # boundary region = gap between TAD_i.end and TAD_{i+1}.start (or ±25kb fallback)
  tad_B_bp <- overlapsAny(sv_gr, tad_bound_gr)
  tad_B    <- tad_B_bp
  if (length(paired_idx) > 0)
    tad_B[paired_idx] <- tad_B_bp[paired_idx] | tad_B_bp[partner_idx[paired_idx]]

  # ── Condition C: SV span crosses a boundary region ────────────────────────
  tad_C <- logical(n)
  if (length(intra_idx) > 0) {
    b1s <- sv_gr[intra_idx]
    b2s <- sv_gr[partner_idx[intra_idx]]
    sv_spans <- GRanges(
      seqnames = seqnames(b1s),
      ranges   = IRanges(
        start = pmin(start(b1s), start(b2s)),
        end   = pmax(start(b1s), start(b2s))
      )
    )
    span_ov <- overlapsAny(sv_spans, tad_bound_gr)
    tad_C[intra_idx]              <- span_ov
    tad_C[partner_idx[intra_idx]] <- span_ov
  }
  # inter-chromosomal paired and unpaired SVs: buffer overlap
  other_idx <- setdiff(seq_len(n), intra_idx)
  if (length(other_idx) > 0)
    tad_C[other_idx] <- overlapsAny(sv_buffer[other_idx], tad_bound_gr)

  # ── Combined TAD disruption (priority A ≥ B > C) ─────────────────────────
  ov_tad <- tad_A | tad_B | tad_C
  tad_condition <- dplyr::case_when(
    tad_A                              ~ "A",
    tad_B & !tad_A                     ~ "B",
    tad_C & !tad_A & !tad_B           ~ "C",
    TRUE                               ~ NA_character_
  )

  # ── CTCF within disrupted boundary region ─────────────────────────────────
  # "TAD+CTCF disrupting" requires CTCF peak inside the affected boundary gap,
  # not just at the breakpoint; "CTCF-only" retains per-breakpoint check.
  ctcf_in_boundary <- logical(n)
  if (any(ov_tad)) {
    # For intra-chr SVs: query = span ± small buffer to capture boundary-B hits
    if (length(intra_idx) > 0) {
      b1s <- sv_gr[intra_idx]
      b2s <- sv_gr[partner_idx[intra_idx]]
      span_buf <- GRanges(
        seqnames = seqnames(b1s),
        ranges   = IRanges(
          start = pmax(pmin(start(b1s), start(b2s)) - 2000L, 1L),
          end   = pmax(start(b1s), start(b2s)) + 2000L
        )
      )
      hits <- as.data.table(findOverlaps(span_buf, tad_bound_gr))
      if (nrow(hits) > 0) {
        hits[, has_ctcf := bound_has_ctcf[subjectHits]]
        agg <- hits[, .(ctcf = any(has_ctcf)), by = queryHits]
        ctcf_in_boundary[intra_idx[agg$queryHits]] <- agg$ctcf
        ctcf_in_boundary[partner_idx[intra_idx[agg$queryHits]]] <- agg$ctcf
      }
    }
    # For inter-chr paired and other disrupted SVs: use buffer
    other_disrupted <- setdiff(which(ov_tad), intra_idx)
    if (length(other_disrupted) > 0) {
      hits2 <- as.data.table(findOverlaps(sv_buffer[other_disrupted], tad_bound_gr))
      if (nrow(hits2) > 0) {
        hits2[, has_ctcf := bound_has_ctcf[subjectHits]]
        agg2 <- hits2[, .(ctcf = any(has_ctcf)), by = queryHits]
        ctcf_in_boundary[other_disrupted[agg2$queryHits]] <- agg2$ctcf
      }
    }
  }

  # ── Per-breakpoint CTCF overlap (used only for the CTCF-only tier) ─────────
  ov_ctcf_bp <- overlapsAny(sv_buffer, ctcf_gr)
  if (length(paired_idx) > 0)
    ov_ctcf_bp[paired_idx] <- ov_ctcf_bp[paired_idx] |
                               ov_ctcf_bp[partner_idx[paired_idx]]

  # ── Distance to nearest feature ────────────────────────────────────────────
  dist_to_nearest <- function(query_gr, subject_gr) {
    if (length(subject_gr) == 0) return(rep(NA_real_, length(query_gr)))
    hits <- distanceToNearest(query_gr, subject_gr, ignore.strand = TRUE)
    out  <- rep(NA_real_, length(query_gr))
    out[queryHits(hits)] <- mcols(hits)$distance
    out
  }
  mcols(sv_gr)$dist_to_TAD  <- dist_to_nearest(sv_gr, tad_bound_gr)
  mcols(sv_gr)$dist_to_CTCF <- dist_to_nearest(sv_gr, ctcf_gr)

  # ── Tier classification ────────────────────────────────────────────────────
  mcols(sv_gr)$stratification <- dplyr::case_when(
    sv_gr$is_hbv                                            ~ "HBV-associated",
    ov_tad & ctcf_in_boundary                               ~ "TAD+CTCF disrupting",
    ov_tad & !ctcf_in_boundary                              ~ "TAD-only",
    !ov_tad & ov_ctcf_bp                                    ~ "CTCF-only",
    mcols(sv_gr)$cnv_class == "copy_neutral" &
      !ov_tad & !ov_ctcf_bp                                 ~ "Copy-neutral",
    TRUE                                                    ~ "Non-boundary"
  )

  mcols(sv_gr)$sv_tier <- match(
    sv_gr$stratification,
    c("HBV-associated", "TAD+CTCF disrupting", "TAD-only",
      "CTCF-only", "Copy-neutral", "Non-boundary")
  )
  mcols(sv_gr)$hbv_associated    <- mcols(sv_gr)$is_hbv
  mcols(sv_gr)$tad_ctcf_overlap  <- ov_tad | ov_ctcf_bp
  mcols(sv_gr)$tad_condition     <- tad_condition
  mcols(sv_gr)$ctcf_in_boundary  <- ctcf_in_boundary

  if (any(mcols(sv_gr)$is_hbv)) {
    hbv_chrom_idx <- which(as.character(seqnames(sv_gr)) == "HBV")
    mcols(sv_gr)$flag_hbv_junction <- FALSE
    if (length(hbv_chrom_idx) > 0) {
      pos <- start(sv_gr[hbv_chrom_idx])
      mcols(sv_gr)$flag_hbv_junction[hbv_chrom_idx] <-
        pos < HBV_JUNCTION_FLANK | pos > (hbv_genome_len - HBV_JUNCTION_FLANK)
      n_flagged <- sum(mcols(sv_gr)$flag_hbv_junction[hbv_chrom_idx])
      if (n_flagged > 0)
        message(sprintf("  %d HBV breakpoints flagged near linearisation junction — manual review recommended", n_flagged))
    }
  }
  sv_gr
}

# ── CLI options ───────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-t", "--tad_bed"),  type = "character",
    default = "/node200data/kachungk/reference/GRCh38/3Dgenomebrowser/HepG2-Control_Merged_MicroC_GSE278978_tad.bed",
    help = "HepG2 Micro-C TAD boundary BED (hg38) [default: %default]"),
  make_option(c("-c", "--ctcf_bed"), type = "character",
    default = "/node200data/kachungk/reference/GRCh38/ensembl/HepG2_ChIP_optpeaks_ENCFF543WTP.bed.gz",
    help = "ENCODE HepG2 CTCF peak BED (hg38) [default: %default]"),
  make_option(c("-o", "--outdir"),   type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/",
    help = "Output directory [default: %default]"),
  make_option(c("--hbv_genome_len"), type = "integer", default = 3215L,
    help = "HBV contig length in bp [default: %default]"),
  make_option(c("--prom_rds"),       type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/canonical_promoters_hg38.gencode_v49.rds",
    help = "GENCODE v49 canonical promoters RDS [default: %default]"),
  make_option(c("--gene_rds"),       type = "character",
    default = "/node200data/kachungk/hcc_data/DMR_SVs/canonical_genes_hg38.gencode_v49.rds",
    help = "GENCODE v49 canonical genes RDS [default: %default]"),
  make_option(c("-g", "--gene_gtf"),      type = "character",
    default = "/node200data/kachungk/reference/GRCh38/gencode.v49.basic.annotation.gtf.gz",
    help = "GTF for gene annotation used if --prom_rds or --gene_rds not found [default: %default]"),
  make_option(c("--out_suffix"),          type = "character", default = "",
    help = "Optional suffix appended before .csv.gz in all three output filenames [default: none]"),
  make_option(c("--chromhmm_liver"), type = "character",
    default = "/node200data/kachungk/reference/GRCh38/chromHMM/E066_15_coreMarks_hg38lift_mnemonics.bed.gz",
    help = "Roadmap E066 liver ChromHMM 15-state BED.gz (active states filtered) [default: %default]"),
  make_option(c("--chromhmm_hepg2"), type = "character",
    default = "/node200data/kachungk/reference/GRCh38/chromHMM/E118_15_coreMarks_hg38lift_mnemonics.bed.gz",
    help = "Roadmap E118 HepG2 ChromHMM 15-state BED.gz (combined with --chromhmm_liver) [default: %default]"),
  make_option(c("--enhancer_bed"),   type = "character",
    default = "/node200data/kachungk/reference/GRCh38/genomic_element/hg38_genehancer_enhancer.bed",
    help = "GeneHancer enhancer BED (liver-annotated) [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$outdir))
  dir.create(opt$outdir, recursive = TRUE)
dir.create(file.path(opt$outdir, "logs"), showWarnings = FALSE, recursive = TRUE)

DEC_LOG <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"

# ── 1. Load SV VCFs ───────────────────────────────────────────────────────────
setwd("/node200data/kachungk/hcc_data/hg38+HBV/somatic_sv")

info_cols <- c("MATE_ID", "MAPQ", "DETAILED_TYPE", "INSIDE_VNTR",
               "LOW_COV_IN", "HP", "PHASESETID")
sv_files  <- list.files(pattern = "*somaticSVs.vcf.gz$", full.names = TRUE)
message(sprintf("Loading %d SV VCF files ...", length(sv_files)))

sv_list <- lapply(sv_files, function(x) {
  vcf   <- readVcf(x)
  bp_df <- c(
    breakpointRanges(vcf, inferMissingBreakends = TRUE, info_columns = info_cols),
    breakendRanges(vcf, info_columns = info_cols)
  ) %>% as.data.frame()

  bnd_df <- bp_df %>%
    filter(svtype == "BND", str_detect(partner, "svrecord")) %>%
    mutate(partner = MATE_ID) %>%
    dplyr::select(-MATE_ID)
  non_bnd_df <- bp_df %>%
    filter(svtype != "BND") %>%
    dplyr::select(-MATE_ID)

  gr <- bind_rows(non_bnd_df, bnd_df) %>%
    distinct() %>%
    mutate(sample = str_remove(basename(x), "\\.severus_somaticSVs.vcf.gz$")) %>%
    GRanges()
  mcols(gr)$bp_id <- names(gr)

  fmt_df <- data.table(
    ID  = as.character(names(rowRanges(vcf))),
    VAF = as.numeric(geno(vcf)$VAF),
    DR  = as.integer(geno(vcf)$DR),
    DV  = as.integer(geno(vcf)$DV)
  )
  if ("hVAF" %in% names(geno(vcf))) {
    hv <- geno(vcf)$hVAF[, 1, , drop = FALSE]
    fmt_df[, `:=`(hVAF_H0 = as.numeric(hv[,,1]),
                  hVAF_H1 = as.numeric(hv[,,2]),
                  hVAF_H2 = as.numeric(hv[,,3]))]
  }
  idx <- match(as.character(mcols(gr)$sourceId), fmt_df$ID)
  mcols(gr)$VAF <- fmt_df$VAF[idx]
  mcols(gr)$DR  <- fmt_df$DR[idx]
  mcols(gr)$DV  <- fmt_df$DV[idx]
  if ("hVAF_H0" %in% colnames(fmt_df)) {
    mcols(gr)$hVAF_H0 <- fmt_df$hVAF_H0[idx]
    mcols(gr)$hVAF_H1 <- fmt_df$hVAF_H1[idx]
    mcols(gr)$hVAF_H2 <- fmt_df$hVAF_H2[idx]
  }
  tibble::as_tibble(gr)
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  mutate(is_hbv = str_detect(ALT, "HBV")) %>%
  dplyr::rename(patient_name = sample, sample = patient_code) %>%
  GRanges() %>%
  split(mcols(.)$sample) %>%
  lapply(function(gr) { names(gr) <- mcols(gr)$bp_id; gr }) %>%
  endoapply(simple_svtype)

PATIENT_IDS <- names(sv_list)
message(sprintf("Loaded SVs for %d patients: %s",
                length(PATIENT_IDS), paste(PATIENT_IDS, collapse = ", ")))

# ── 2. Load PURPLE CNV segments ───────────────────────────────────────────────
cnv_files <- list.files(
  "../../cnv_deepsomatic.out_hg38/purple",
  pattern = "*tumor.purple.segment.tsv$",
  full.names = TRUE
)

cnv_segments <- lapply(cnv_files, function(x) {
  df <- fread(x)
  if (!"minorAlleleCopyNumber" %in% colnames(df))
    df$minorAlleleCopyNumber <- NA_real_
  df %>%
    dplyr::select(chromosome, start, end, tumorCopyNumber,
                  minorAlleleCopyNumber, bafCount, observedBAF,
                  germlineStatus, depthWindowCount) %>%
    dplyr::rename(copyNumber = tumorCopyNumber) %>%
    mutate(sample = str_remove(basename(x), "\\_tumor.purple.segment.tsv$"))
}) %>%
  bind_rows() %>%
  anonym_sample() %>%
  dplyr::rename(patient_name = sample, sample = patient_code) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  split(mcols(.)$sample)

# ── 3. CNV classification ─────────────────────────────────────────────────────
sv_list_classified <- lapply(PATIENT_IDS, function(pt) {
  classified <- classify_sv_cnv_final(sv_list[[pt]], cnv_segments[[pt]])
  cat(sprintf("\n[%s] CNV class 분포:\n", pt))
  print(table(mcols(classified)$cnv_class,
              mcols(classified)$cnv_class_source,
              dnn = c("cnv_class", "source")))
  classified
}) %>% setNames(PATIENT_IDS)

rm(sv_list, cnv_segments); gc()

# ── 4. Load TAD / CTCF references ────────────────────────────────────────────
ctcf_gr <- fread(opt$ctcf_bed) %>%
  makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")

TAD_BOUNDARY_FALLBACK_BP <- 25000L  # ±25 kb when adjacent TADs have no gap

tad_result <- local({
  tad_bodies <- fread(opt$tad_bed) %>%
    dplyr::rename(seqnames = V1, start = V2, end = V3) %>%
    dplyr::arrange(seqnames, start)

  tad_body_gr <- makeGRangesFromDataFrame(tad_bodies)

  # Boundary region = gap between TAD_i.end and TAD_{i+1}.start.
  # When gap is zero (adjacent TADs), fall back to ±25 kb around the junction
  # to reflect Micro-C resolution limits.
  bound_df <- lapply(split(tad_bodies, tad_bodies$seqnames), function(cdf) {
    m <- nrow(cdf)
    if (m < 2) return(NULL)
    gap_start <- cdf$end[-m]
    gap_end   <- cdf$start[-1]
    has_gap   <- gap_end > gap_start
    data.frame(
      seqnames = cdf$seqnames[-m],
      start    = pmax(ifelse(has_gap, gap_start,
                             gap_start - TAD_BOUNDARY_FALLBACK_BP), 1L),
      end      = ifelse(has_gap, gap_end,
                        gap_start + TAD_BOUNDARY_FALLBACK_BP),
      has_gap  = has_gap
    )
  }) %>% bind_rows()

  list(
    body  = tad_body_gr,
    bound = makeGRangesFromDataFrame(bound_df, keep.extra.columns = TRUE)
  )
})
tad_body_gr  <- tad_result$body
tad_bound_gr <- tad_result$bound

n_gap    <- sum(tad_bound_gr$has_gap)
n_fb     <- sum(!tad_bound_gr$has_gap)
message(sprintf("TAD bodies: %d | boundary regions: %d (%d gap-based, %d ±%dkb fallback) across %d chr",
                length(tad_body_gr), length(tad_bound_gr), n_gap, n_fb,
                TAD_BOUNDARY_FALLBACK_BP / 1000L,
                length(unique(as.character(seqnames(tad_bound_gr))))))

# ── 4b. Load CRE references (liver ChromHMM + GeneHancer) ────────────────────
load_chromhmm_active <- function(path) {
  if (is.null(path) || !file.exists(path)) return(GRanges())
  dt <- data.table::fread(path, col.names = c("seqnames", "start", "end", "state"))
  dt <- dt[state %in% CHROMHMM_ACTIVE_STATES]
  if (nrow(dt) == 0L) return(GRanges())
  makeGRangesFromDataFrame(as.data.frame(dt))
}

cre_liver_gr <- GenomicRanges::reduce(
  c(load_chromhmm_active(opt$chromhmm_liver),
    load_chromhmm_active(opt$chromhmm_hepg2))
)

enhancer_gr <- tryCatch({
  if (!is.null(opt$enhancer_bed) && file.exists(opt$enhancer_bed)) {
    dt <- data.table::fread(opt$enhancer_bed, skip = 1L,
                            col.names = c("seqnames", "start", "end", "element", "type"))
    makeGRangesFromDataFrame(as.data.frame(dt[, 1:3]))
  } else NULL
}, error = function(e) { message("GeneHancer load failed: ", e$message); NULL })

message(sprintf("CRE references loaded: ChromHMM liver+HepG2 active = %d regions; GeneHancer = %d regions",
                length(cre_liver_gr),
                if (!is.null(enhancer_gr)) length(enhancer_gr) else 0L))

# ── 5. TAD / CTCF stratification + CRE annotation ────────────────────────────
sv_list_strat <- lapply(
  sv_list_classified,
  sv_stratification,
  tad_body_gr = tad_body_gr, tad_bound_gr = tad_bound_gr,
  ctcf_gr = ctcf_gr, hbv_genome_len = opt$hbv_genome_len
)

sv_list_strat <- lapply(
  sv_list_strat, annotate_sv_cre,
  cre_liver_gr = cre_liver_gr, enhancer_gr = enhancer_gr
)

sv_df <- lapply(sv_list_strat, as.data.frame) %>% dplyr::bind_rows()
cat("=== SV tier distribution ===\n")
print(sv_df %>% dplyr::count(stratification) %>% as.data.frame())
cat("=== CRE sub-tier (TAD+CTCF disrupting only) ===\n")
print(sv_df %>% dplyr::filter(stratification == "TAD+CTCF disrupting") %>%
        dplyr::count(tad_ctcf_cre_subtier) %>% as.data.frame())

# 6. Write sv_tad_ctcf_annotation.csv.gz =====================================
out_strat <- file.path(opt$outdir, paste0("sv_tad_ctcf_annotation", opt$out_suffix, ".csv.gz"))
fwrite(sv_df, out_strat)
message("Stratification saved to: ", out_strat)

#  7. Gene annotation (GENCODE v49) ===========================================
if (!file.exists(opt$prom_rds) || !file.exists(opt$gene_rds)) {
  message("RDS files for gene annotation not found.\n Create promoter/gene rds file with given GTF file")
  annot <- rtracklayer::import(opt$gene_gtf) #nolint
  keep_annot <- annot$type %in% c("gene", "transcript") &
    annot$gene_type %in% c("protein_coding", "lncRNA") &
    !is.na(annot$gene_name) &
    nzchar(annot$gene_name)
  annot <- annot[keep_annot]

  tx <- annot[annot$type == "transcript"] %>%
    as.data.frame() %>%
    mutate(
      tx.score = case_when(
        tag == "MANE_Select"        ~ 1L,
        tag == "MANE_Plus_Clinical" ~ 2L,
        tag == "appris_principal_1" ~ 3L,
        tag == "basic"              ~ 4L,
        TRUE                        ~ 5L
      )
    ) %>%
    group_by(gene_id) %>%
    arrange(tx.score) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    GRanges()

  canonical_prom_gr <- promoters(tx, upstream = 2000, downstream = 200)
  saveRDS(canonical_prom_gr, opt$prom_rds)

  # Additional cache for gene name
  mcols(tx) <- DataFrame(gene_id = tx$gene_id, gene_name = tx$gene_name, gene_type = tx$gene_type)
  saveRDS(tx, opt$gene_rds)

  prom_gr <- canonical_prom_gr
  gene_gr <- tx
} else {
  message("\nLoading GENCODE v49 annotations ...")
  prom_gr <- readRDS(opt$prom_rds)
  gene_gr <- readRDS(opt$gene_rds)
}

message(sprintf("  promoters: %d ranges", length(prom_gr)))
message(sprintf("  genes:     %d ranges", length(gene_gr)))

sv_dt <- data.table::as.data.table(sv_df)

sv_gr <- GRanges(
  seqnames = sv_dt$seqnames,
  ranges   = IRanges(start = sv_dt$start, width = 1L)
)
mcols(sv_gr) <- DataFrame(
  sv_id            = sv_dt$bp_id,
  svtype           = sv_dt$svtype,
  sample           = sv_dt$sample,
  sv_tier          = sv_dt$sv_tier,
  hbv_associated   = sv_dt$hbv_associated,
  tad_ctcf_overlap = sv_dt$tad_ctcf_overlap,
  stratification   = sv_dt$stratification
)

# Nearest TSS
message("Computing nearest TSS ...")
strand_v <- as.character(strand(prom_gr))
tss_pos  <- ifelse(strand_v == "+", start(prom_gr) + 2000L, end(prom_gr) - 2000L)
tss_gr   <- GRanges(
  seqnames  = seqnames(prom_gr),
  ranges    = IRanges(start = tss_pos, width = 1L),
  strand    = strand(prom_gr),
  gene_name = prom_gr$gene_name,
  gene_id   = prom_gr$gene_id,
  gene_type = prom_gr$gene_type
)

nn          <- distanceToNearest(sv_gr, tss_gr, ignore.strand = TRUE)
sv_pos_v    <- start(sv_gr)[queryHits(nn)]
tss_pos_v   <- start(tss_gr)[subjectHits(nn)]
gstrand_v   <- as.character(strand(tss_gr))[subjectHits(nn)]
signed_dist <- ifelse(gstrand_v == "+", sv_pos_v - tss_pos_v, tss_pos_v - sv_pos_v)

nearest_dt <- data.table(
  sv_id             = mcols(sv_gr)$sv_id[queryHits(nn)],
  nearest_gene      = tss_gr$gene_name[subjectHits(nn)],
  nearest_gene_id   = tss_gr$gene_id[subjectHits(nn)],
  nearest_gene_type = tss_gr$gene_type[subjectHits(nn)],
  dist_to_tss       = mcols(nn)$distance,
  signed_dist_tss   = signed_dist,
  tss_relative_pos  = fcase(
    mcols(nn)$distance == 0L,               "at_TSS",
    signed_dist < -200L,                    "upstream",
    signed_dist >= -200L & signed_dist < 0, "proximal_upstream",
    signed_dist >= 0L & signed_dist <= 200L,"proximal_downstream",
    default                                 = "downstream"
  )
)
message(sprintf("  Nearest-gene: %d / %d SVs matched", nrow(nearest_dt), length(sv_gr)))

# Promoter overlap
message("Finding promoter overlaps ...")
prom_hits <- findOverlaps(sv_gr, prom_gr, ignore.strand = TRUE)
prom_ov <- unique(data.table(
  sv_id     = mcols(sv_gr)$sv_id[queryHits(prom_hits)],
  gene_name = prom_gr$gene_name[subjectHits(prom_hits)]
))
prom_flag <- prom_ov[, .(promoter_overlap = TRUE), by = sv_id][, unique(.SD)]
message(sprintf("  Breakpoints overlapping a promoter: %d", uniqueN(prom_ov$sv_id)))

# Gene-body overlap
message("Finding gene-body overlaps ...")
gene_hits <- findOverlaps(sv_gr, gene_gr, ignore.strand = TRUE)
gene_flag <- unique(data.table(
  sv_id = mcols(sv_gr)$sv_id[queryHits(gene_hits)]
))[, gene_body_overlap := TRUE]
message(sprintf("  Breakpoints inside a gene body: %d", nrow(gene_flag)))

# All TSS within ±50 kb (long format)
message("Finding all TSS within ±50 kb ...")
sv_50k <- sv_gr
ranges(sv_50k) <- IRanges(
  start = pmax(start(sv_gr) - 50000L, 1L),
  end   = start(sv_gr) + 50000L
)
nearby_hits <- findOverlaps(sv_50k, tss_gr, ignore.strand = TRUE)
nearby_dt   <- unique(data.table(
  sv_id     = mcols(sv_gr)$sv_id[queryHits(nearby_hits)],
  gene_name = toupper(tss_gr$gene_name[subjectHits(nearby_hits)]),
  gene_id   = tss_gr$gene_id[subjectHits(nearby_hits)],
  gene_type = tss_gr$gene_type[subjectHits(nearby_hits)]
))
sv_pos_nb  <- start(sv_gr)[queryHits(nearby_hits)]
tss_pos_nb <- start(tss_gr)[subjectHits(nearby_hits)]
nearby_dt[, dist_to_sv := abs(sv_pos_nb - tss_pos_nb)]
message(sprintf("  SV–gene pairs ±50 kb: %d (unique genes: %d)",
                nrow(nearby_dt), uniqueN(nearby_dt$gene_name)))

# ── 8. Compile and write gene annotation outputs ─────────────────────────────
sv_meta <- sv_dt[, .(
  sv_id    = bp_id,
  seqnames, start,
  svtype, sample, patient_name,
  sv_tier, hbv_associated, tad_ctcf_overlap, stratification,
  dist_to_TAD, dist_to_CTCF,
  cre_any_overlap, cre_chromhmm_liver, cre_type, tad_ctcf_cre_subtier
)]

near_out <- merge(sv_meta,   nearest_dt, by = "sv_id", all.x = TRUE)
near_out <- merge(near_out,  prom_flag,  by = "sv_id", all.x = TRUE)
near_out <- merge(near_out,  gene_flag,  by = "sv_id", all.x = TRUE)
near_out[is.na(promoter_overlap),  promoter_overlap  := FALSE]
near_out[is.na(gene_body_overlap), gene_body_overlap := FALSE]

long_out <- merge(
  nearby_dt,
  sv_meta[, .(sv_id, svtype, sample, sv_tier)],
  by = "sv_id"
)
long_out <- merge(
  long_out,
  prom_ov[, .(sv_id, gene_name = toupper(gene_name), promoter_overlap = TRUE)],
  by = c("sv_id", "gene_name"), all.x = TRUE
)
long_out[is.na(promoter_overlap), promoter_overlap := FALSE]
setorder(long_out, sv_id, dist_to_sv)

NEAR_OUT <- file.path(opt$outdir, paste0("sv_bp.tier_nearest_gene", opt$out_suffix, ".csv.gz"))
LONG_OUT <- file.path(opt$outdir, paste0("sv_bp.tier_with_gene",    opt$out_suffix, ".csv.gz"))
fwrite(near_out, NEAR_OUT, compress = "gzip")
fwrite(long_out, LONG_OUT, compress = "gzip")

message(sprintf("\nNear-gene output : %d rows → %s", nrow(near_out), NEAR_OUT))
message(sprintf("Long-format output: %d rows → %s", nrow(long_out), LONG_OUT))

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n=== Annotation summary ===")
message(sprintf("Total breakpoints    : %d", nrow(near_out)))
message(sprintf("Promoter overlap     : %d (%.1f%%)",
                sum(near_out$promoter_overlap), 100 * mean(near_out$promoter_overlap)))
message(sprintf("Gene-body overlap    : %d (%.1f%%)",
                sum(near_out$gene_body_overlap), 100 * mean(near_out$gene_body_overlap)))
message(sprintf("Unique nearest genes : %d", uniqueN(near_out$nearest_gene)))
message(sprintf("SV-gene pairs ±50 kb : %d (unique genes: %d)",
                nrow(long_out), uniqueN(long_out$gene_name)))
message("\n=== By tier ===")
print(near_out[, .(
  n_sv        = .N,
  n_in_prom   = sum(promoter_overlap),
  n_in_gene   = sum(gene_body_overlap),
  n_uniq_gene = uniqueN(nearest_gene)
), by = sv_tier])

# ── Decision log ──────────────────────────────────────────────────────────────
if (!dir.exists(dirname(DEC_LOG)))
  dir.create(dirname(DEC_LOG), recursive = TRUE, showWarnings = FALSE)
cat(sprintf(
  "[%s] 02_sv_annotation_stratify.R: %d SV bp stratified; %d uniq nearest genes; %d SV-gene pairs ±50 kb\n",
  format(Sys.Date()), nrow(near_out), uniqueN(near_out$nearest_gene), nrow(long_out)
), file = DEC_LOG, append = TRUE)

message("\nDone.")
