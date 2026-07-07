# 05_fragility_enrichment

Manuscript: Results > Shared structural fragility associates SVs and aDMRs with SegDup/LAD/B-compartment features, independent of proximity to each other.

## Files
- `fragility_glm.R` — C13 (SV -> SegDup), C14 (somatic aDMR -> SegDup), C22 (+CpG), C23 (SV-far) weighted GLMs, window + merged-region granularity
- `fragility_multivariate.R` — P0-C: multivariate logistic regression, independent fragility effects on SV breakpoints
- `admr_independent_fragility.R` — P0-D: aDMR-own fragility enrichment (SV-independent), Gold vs Silver tier breakdown
- `admr_sv_exclusion.R` — A3: SegDup enrichment after excluding aDMRs near SVs (10kb/50kb thresholds)
- `agold1_proximity_decircularize.R` — A_GOLD1: recurrence-gradient analysis (Gold*/Gold*_R3/Gold*_R5); primary source of the manuscript's Fig 4B recurrence-scaling claim
- `aeq1_or_equivalence.R` — A_EQ1: formal Wald/TOST equivalence test, SegDup OR(SV) vs OR(aDMR) across 3 covariate levels
- `cpg_adjusted_glm.R` — A1: CpG-density adjusted C14 GLM
- `cfs_overlap.R` — P1-C: Common Fragile Site (CFS) overlap for non-boundary SVs and Gold aDMR
- `nahr_microhomology.R` — P0-B: NAHR breakpoint signature proxies (SV size, VNTR context, PC1) since Severus reports HOMLEN=0
- `repeatmask.R` — Figure S6 RepeatMasker stratification (negative control); re-exports a cached RDS if present, otherwise writes a placeholder (the from-scratch computation is not included in this repo)
- `replication_fragility_annotation.R` — LAD/SegDup/PC1/repeat-density annotation of SV breakpoints, boundary vs non-boundary
- `replication_timing.R` — P1-5: ENCODE Repli-seq + common fragile site enrichment for non-boundary SVs
- `replication_timing_admr.R` — MicroC PC1 vs ENCODE Repli-seq concordance at aDMR/SV loci, probe-level RT x is_admr interaction
- `segdup_permutation.R` / `segdup_permutation_control.R` — circular-shift permutation nulls for SegDup enrichment (C24/A2)
- `sensitivity.R` — C23 distance-cutoff sensitivity (1/10/50/100kb)
- `somatic_vs_constitutional.R` — three-way SegDup OR comparison: somatic vs recurrently remodeled constitutional vs normal aDMR
- `download_encode_repliseq.sh` — downloads ENCODE2 HepG2 Repli-seq phase bigWigs (hg19) used by `replication_timing_admr.R`