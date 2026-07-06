#!/usr/bin/env Rscript
# P1-C: Common Fragile Sites (CFS) overlap with SV breakpoints and aDMR loci
#
# Purpose:
#   Test whether non-boundary SV breakpoints and Gold/Silver aDMR loci
#   are enriched at known Common Fragile Sites (CFS).
#   Validates the "fragility co-localization" model (C7/C14) with an
#   orthogonal, well-established genomic instability feature.
#
# CFS catalog: curated from Le Tallec 2013, Helmrich 2011, Durkin 2007,
#   and UCSC cytoBand hg38 coordinates.
#
# Run: mamba run -n renv Rscript post_processing/cfs_overlap.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(rtracklayer)
})

set.seed(42)

# ── Paths ─────────────────────────────────────────────────────────────────────
SV_FRAG  <- "/node200data/kachungk/hcc_data/DMR_SVs/result/sv_fragility_annotation.csv"
GOLD_CSV <- "/node200data/kachungk/hcc_data/DMR_SVs/04.final_candidate/gold_tier_final.csv"
SEGDUP   <- "/node200data/kachungk/reference/GRCh38/LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed"
FAI      <- "/node200data/kachungk/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"
CACHE    <- "/node200data/kachungk/hcc_data/DMR_SVs/external_validation_cache"
OUT_DIR  <- "/node200data/kachungk/hcc_data/DMR_SVs/result"
FIG_DIR  <- file.path(OUT_DIR, "figures")
LOG_FILE <- "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"
dir.create(FIG_DIR, showWarnings = FALSE)

# ── 1. CFS catalog (hg38) ─────────────────────────────────────────────────────
# Curated from:
#   Le Tallec et al. 2013 Nat Struct Mol Biol (APH-induced, top expressed in HepG2)
#   Helmrich et al. 2011 Cell (FHIT/WWOX)
#   Durkin & Glover 2007 Annu Rev Genet
#   Chromosomal coordinates from UCSC hg38 cytoBand track (cytoband midpoints ± flanking)

# Download hg38 cytoBand if needed
cytoband_cache <- file.path(CACHE, "hg38_cytoBand.txt")
if (!file.exists(cytoband_cache)) {
  message("Downloading hg38 cytoBand...")
  download.file(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    paste0(cytoband_cache, ".gz"),
    method = "curl", quiet = TRUE, extra = "-L"
  )
  system(paste("gunzip -f", shQuote(paste0(cytoband_cache, ".gz"))))
}

cytoband <- fread(cytoband_cache,
  col.names = c("chr", "start", "end", "name", "stain"),
  data.table = FALSE
) |> filter(grepl("^chr[0-9XY]+$", chr))

# Known CFS loci: cytoband → gene (literature-curated)
# Sources: Le Tallec 2013 Table S1, Helmrich 2011, Durkin 2007 review
cfs_catalog <- data.frame(
  cfs_name = c(
    "FRA3B",  "FRA16D", "FRA6E",  "FRA7H",  "FRA1H",
    "FRA2G",  "FRA4D",  "FRA6F",  "FRA7E",  "FRA7G",
    "FRA9C",  "FRA13A", "FRA14B", "FRA18D", "FRA20E",
    "FRA1B",  "FRA4F",  "FRA5C",  "FRA8C",  "FRA10D",
    "FRA11F", "FRA12A", "FRA15A", "FRA17A"
  ),
  cytoband = c(
    "3p14.2",  "16q23.2", "6q26",    "7q31.2",  "1p36.1",
    "2q31.1",  "4q22.1",  "6q21",    "7p22.1",  "7q11.23",
    "9p21.1",  "13q13.2", "14q24.1", "18q21.2", "20p12.2",
    "1p36.2",  "4q26",    "5q23.1",  "8q22.3",  "10q23.3",
    "11q14.2", "12q23.3", "15q22.2", "17p12"
  ),
  gene = c(
    "FHIT", "WWOX", "PARK2", "CTNNA2", "RFN213",
    "LYPD6B", "GRID2", "LAMA2", "SDK1", "AUTS2",
    "CNTNAP2_adj", "SMAD9_adj", "RAD51B", "LIPG_adj", "MACROD2",
    "RCC1_adj", "DPPA3_adj", "FAT2_adj", "KCNK9_adj", "PTEN",
    "FOLH1_adj", "TBX3_adj", "HERC2_adj", "MAP2K4_adj"
  ),
  stringsAsFactors = FALSE
)

# Map cytobands to hg38 coordinates
get_cytoband_coords <- function(band_name, cytoband_df) {
  chr_part <- sub("^([0-9XY]+)[pq].*", "chr\\1", band_name)  # "3p14.2" → "chr3"
  arm_part <- sub("^[0-9XY]+", "", band_name)                 # "3p14.2" → "p14.2"
  rows <- cytoband_df[cytoband_df$chr == chr_part &
                       cytoband_df$name == arm_part, ]
  if (nrow(rows) == 0) return(NULL)
  data.frame(chr = rows$chr[1], start = min(rows$start), end = max(rows$end))
}

cfs_coords <- lapply(seq_len(nrow(cfs_catalog)), function(i) {
  coords <- get_cytoband_coords(cfs_catalog$cytoband[i], cytoband)
  if (is.null(coords)) return(NULL)
  cbind(cfs_catalog[i, ], coords)
}) |> bind_rows()

cat(sprintf("CFS catalog: %d loci mapped to hg38 coordinates\n", nrow(cfs_coords)))
print(cfs_coords |> select(cfs_name, gene, chr, start, end))

cfs_gr <- GRanges(
  seqnames = cfs_coords$chr,
  ranges   = IRanges(cfs_coords$start, cfs_coords$end),
  cfs_name = cfs_coords$cfs_name,
  gene     = cfs_coords$gene
)

fwrite(cfs_coords, file.path(OUT_DIR, "cfs_hg38_catalog.csv"))
message(sprintf("CFS catalog saved: %d loci", nrow(cfs_coords)))

# ── 2. Load SV breakpoints ────────────────────────────────────────────────────
message("\nLoading SV breakpoints...")
sv <- fread(SV_FRAG, data.table = FALSE) |>
  filter(!is.na(seqnames), grepl("^chr[0-9XY]+$", seqnames)) |>
  distinct(bp_id, .keep_all = TRUE)

cat(sprintf("SV breakpoints: %d total\n", nrow(sv)))

sv_gr <- GRanges(seqnames = sv$seqnames,
                 ranges   = IRanges(sv$start, sv$start),
                 bp_id    = sv$bp_id,
                 tier     = sv$sv_tier,
                 tier_grp = sv$tier_group,
                 segdup   = sv$segdup_overlap)

# Non-boundary SVs (Tier 3–6)
sv_nonbound_gr <- sv_gr[sv_gr$tier_grp == "Non-boundary"]
sv_bound_gr    <- sv_gr[sv_gr$tier_grp != "Non-boundary"]
cat(sprintf("Non-boundary: %d | Boundary: %d\n",
            length(sv_nonbound_gr), length(sv_bound_gr)))

# ── 3. Load aDMR loci (Gold + Silver) ─────────────────────────────────────────
message("Loading Gold aDMR loci...")
gold <- fread(GOLD_CSV, data.table = FALSE)

find_col <- function(df, cands) {
  m <- intersect(cands, names(df)); if (length(m)) m[[1]] else stop("col not found")
}
col_chr <- find_col(gold, c("admr_chr", "chr", "seqnames"))
col_s   <- find_col(gold, c("admr_start", "start"))
col_e   <- find_col(gold, c("admr_end", "end"))
col_tier <- find_col(gold, c("tier_class", "tier", "gold_bonus"))

admr_gr <- GRanges(
  seqnames  = gold[[col_chr]],
  ranges    = IRanges(gold[[col_s]], gold[[col_e]]),
  tier_class = gold[[col_tier]]
)
admr_gr <- admr_gr[grepl("^chr[0-9XY]+$", as.character(seqnames(admr_gr)))]
gold_gr  <- admr_gr[admr_gr$tier_class == "Gold"]
cat(sprintf("aDMR: %d Gold loci\n", length(gold_gr)))

# ── 4. Genome background (random regions, chr-matched) ────────────────────────
message("Generating random background...")
chrom_sizes <- fread(FAI, col.names = c("chr","len","x","y","z"),
                     data.table = FALSE) |>
  filter(grepl("^chr[0-9XY]+$", chr), !grepl("_", chr)) |>
  select(chr, len)

rand_regions <- function(n, chr_dist = NULL, sizes = chrom_sizes, width = 1L) {
  if (is.null(chr_dist)) {
    chr_dist <- setNames(rep(1L, nrow(sizes)), sizes$chr)
  }
  chr_probs <- chr_dist / sum(chr_dist)
  chrs <- sample(names(chr_probs), n, replace = TRUE, prob = chr_probs)
  tab  <- table(chrs)
  grl  <- lapply(names(tab), function(ch) {
    len <- sizes$len[sizes$chr == ch]
    pos <- sample.int(len - width, tab[[ch]], replace = TRUE)
    GRanges(ch, IRanges(pos, pos + width - 1L))
  })
  do.call(c, Filter(Negate(is.null), grl))
}

# Match chr distribution of non-boundary SVs for background
sv_chr_tab <- table(as.character(seqnames(sv_nonbound_gr)))
n_bg_sv    <- length(sv_nonbound_gr) * 10L
bg_sv_gr   <- rand_regions(n_bg_sv, chr_dist = sv_chr_tab)

# Match chr distribution of Gold aDMRs for background
admr_chr_tab <- table(as.character(seqnames(gold_gr)))
n_bg_admr    <- length(gold_gr) * 100L
bg_admr_gr   <- rand_regions(n_bg_admr, chr_dist = admr_chr_tab,
                              width = 500L)  # ~aDMR width

# ── 5. CFS overlap analysis ───────────────────────────────────────────────────
message("\nComputing CFS overlap...")

fisher_cfs <- function(query_gr, bg_gr, cfs_ref_gr, label) {
  q_hit <- overlapsAny(query_gr, cfs_ref_gr)
  b_hit <- overlapsAny(bg_gr,    cfs_ref_gr)

  tbl <- matrix(c(sum(q_hit),  sum(!q_hit),
                  sum(b_hit),  sum(!b_hit)),
                nrow = 2,
                dimnames = list(c("CFS_overlap", "CFS_no_overlap"),
                                c("Query", "Background")))
  fish <- fisher.test(tbl, alternative = "greater")

  data.frame(
    group         = label,
    n_query       = length(query_gr),
    n_bg          = length(bg_gr),
    n_query_cfs   = sum(q_hit),
    n_bg_cfs      = sum(b_hit),
    pct_query_cfs = mean(q_hit) * 100,
    pct_bg_cfs    = mean(b_hit) * 100,
    OR            = fish$estimate,
    CI_lo         = fish$conf.int[1],
    CI_hi         = fish$conf.int[2],
    p             = fish$p.value,
    stringsAsFactors = FALSE
  )
}

results <- bind_rows(
  fisher_cfs(sv_nonbound_gr, bg_sv_gr,   cfs_gr, "Non-boundary SV"),
  fisher_cfs(sv_bound_gr,    bg_sv_gr,   cfs_gr, "Boundary SV"),
  fisher_cfs(gold_gr,        bg_admr_gr, cfs_gr, "Gold aDMR"),
  fisher_cfs(sv_gr,          bg_sv_gr,   cfs_gr, "All SV")
) |> mutate(
  p_adj = p.adjust(p, method = "BH"),
  sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                    p_adj < 0.05  ~ "*",   TRUE ~ "ns"),
  OR_lab = sprintf("%.2f [%.2f–%.2f] %s", OR, CI_lo, CI_hi, sig)
)

cat("\n=== CFS Overlap Results ===\n")
print(results |> select(group, n_query, pct_query_cfs, pct_bg_cfs, OR, p_adj, sig, OR_lab))

fwrite(results, file.path(OUT_DIR, "cfs_overlap.csv"))
message("Saved: cfs_overlap.csv")

# ── 6. Per-CFS locus overlap table ────────────────────────────────────────────
message("\nPer-CFS locus breakdown...")

per_locus <- lapply(seq_len(length(cfs_gr)), function(i) {
  locus_gr   <- cfs_gr[i]
  sv_hit     <- sum(overlapsAny(sv_nonbound_gr, locus_gr))
  gold_hit   <- sum(overlapsAny(gold_gr, locus_gr))
  segdup_pct <- mean(overlapsAny(locus_gr,
    import(SEGDUP, format = "BED") |>
      (\(x) { seqlevelsStyle(x) <- "UCSC"; x })()
  )) * 100

  data.frame(
    cfs_name    = cfs_gr$cfs_name[i],
    gene        = cfs_gr$gene[i],
    chr         = as.character(seqnames(locus_gr)),
    start       = start(locus_gr),
    end         = end(locus_gr),
    n_nb_sv     = sv_hit,
    n_gold_admr = gold_hit,
    has_sv      = sv_hit > 0,
    has_admr    = gold_hit > 0,
    both        = sv_hit > 0 & gold_hit > 0,
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

cat(sprintf("\nCFS loci with non-boundary SV: %d/%d\n",
            sum(per_locus$has_sv), nrow(per_locus)))
cat(sprintf("CFS loci with Gold aDMR: %d/%d\n",
            sum(per_locus$has_admr), nrow(per_locus)))
cat(sprintf("CFS loci with BOTH SV + aDMR: %d/%d\n",
            sum(per_locus$both), nrow(per_locus)))

print(per_locus |> arrange(desc(both), desc(n_nb_sv + n_gold_admr)))

fwrite(per_locus, file.path(OUT_DIR, "cfs_per_locus.csv"))
message("Saved: cfs_per_locus.csv")

# ── 7. Figure ─────────────────────────────────────────────────────────────────
message("\nGenerating figures...")

theme_hcc <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

# Panel A: forest plot of CFS enrichment ORs
p_forest <- results |>
  filter(group != "All SV") |>
  mutate(group = factor(group, levels = rev(c("Non-boundary SV", "Boundary SV", "Gold aDMR")))) |>
  ggplot(aes(x = OR, xmin = CI_lo, xmax = CI_hi, y = group, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), width = 0.25,
                orientation = "y") +
  geom_point(size = 4) +
  geom_text(aes(x = pmax(CI_hi, OR) * 1.1, label = OR_lab),
            hjust = 0, size = 3) +
  scale_color_manual(values = c("***" = "#E24B4A", "**" = "#E07B39",
                                "*"   = "#F0BA4E", "ns"  = "#888780")) +
  scale_x_log10(limits = c(0.3, 30)) +
  labs(
    title    = "Common Fragile Site Enrichment",
    subtitle = sprintf("CFS catalog: %d loci (Le Tallec 2013 + Helmrich 2011)", nrow(cfs_coords)),
    x        = "Odds Ratio vs random genome (log scale)",
    y        = NULL, color = "FDR"
  ) +
  theme_hcc + theme(legend.position = "right")

# Panel B: per-locus heatmap (SV / aDMR presence at each CFS)
heat_df <- per_locus |>
  select(cfs_name, gene, n_nb_sv, n_gold_admr) |>
  tidyr::pivot_longer(c(n_nb_sv, n_gold_admr),
                      names_to = "type", values_to = "count") |>
  mutate(
    type = recode(type,
                  n_nb_sv     = "Non-boundary SV",
                  n_gold_admr = "Gold aDMR"),
    label = paste0(cfs_name, "\n(", gene, ")"),
    label = factor(label, levels = rev(unique(label[order(per_locus$n_nb_sv + per_locus$n_gold_admr)])))
  )

p_heat <- ggplot(heat_df, aes(x = type, y = label, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(count > 0, count, "")),
            size = 3, color = "white") +
  scale_fill_gradient(low = "grey90", high = "#E24B4A", name = "Count") +
  labs(
    title = "SV & aDMR at CFS Loci",
    x = NULL, y = NULL
  ) +
  theme_hcc +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size = 7))

ggsave(file.path(FIG_DIR, "fig_s8a_cfs_forest.png"), p_forest,
       width = 10, height = 4, dpi = 150)
ggsave(file.path(FIG_DIR, "fig_s8b_cfs_heatmap.png"), p_heat,
       width = 8, height = 8, dpi = 150)
message("Saved: fig_s8a_cfs_forest.png + fig_s8b_cfs_heatmap.png")

# ── 8. Summary ────────────────────────────────────────────────────────────────
nb_row <- results |> filter(group == "Non-boundary SV")
gd_row <- results |> filter(group == "Gold aDMR")

cat(sprintf("\n=== P1-C Summary ===\n"))
cat(sprintf("Non-boundary SV at CFS: OR=%.2f [%.2f–%.2f] %s (p_adj=%.3g)\n",
            nb_row$OR, nb_row$CI_lo, nb_row$CI_hi, nb_row$sig, nb_row$p_adj))
cat(sprintf("Gold aDMR at CFS:       OR=%.2f [%.2f–%.2f] %s (p_adj=%.3g)\n",
            gd_row$OR, gd_row$CI_lo, gd_row$CI_hi, gd_row$sig, gd_row$p_adj))

cat(
  sprintf("[%s] P1-C cfs_overlap: Non-boundary SV OR=%.2f %s; Gold aDMR OR=%.2f %s; CFS loci with both=%d/%d\n",
          Sys.Date(),
          nb_row$OR, nb_row$sig,
          gd_row$OR, gd_row$sig,
          sum(per_locus$both), nrow(per_locus)),
  file = LOG_FILE, append = TRUE
)

message("\nDone.")
