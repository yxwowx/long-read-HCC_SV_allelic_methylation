# 09_hbv_integration

Manuscript: Methods > HBV integration; HBV integration per-read and
allele-anchored analyses. Results > HBV integration constitutes a second,
architecture-independent pathway.

## Files
- `12_hbv_analysis.R` <- `remodeled_constitutional_AMR/pipeline/12_hbv_analysis.R` (HBV integration x aDMR association tables) — **copied**
- `hbv_perread_meth.py` <- `remodeled_constitutional_AMR/post_processing/hbv_perread_meth.py` — **copied**
  (C32, Level A — per-read HBV-junction vs reference-read beta, 3 distance bins)
- `hbv_allele_anchored.R` <- `remodeled_constitutional_AMR/post_processing/hbv_allele_anchored.R` (C33, Level B — earlier pairing, superseded counts) — **copied** (kept for provenance; do not use for reported numbers, see `06_hbv_allele_specific.R` below)
- `06_hbv_allele_specific.R` <- `somatic_AMR/06_hbv_allele_specific.R` — **copied**
  (**current** — supersedes hbv_allele_anchored.R's pairing counts; produces the 510/25,725/26,235-pair pool and current LMM cited in text)
