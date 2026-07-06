# 05_fragility_enrichment

Manuscript: Methods > Genomic fragility annotation; Results > Segmental
duplication exhibited distinct enrichment pattern in genomic fragile sites;
De novo somatic aDMRs show no independent SegDup enrichment; Constitutional
aDMR hotspots scale with recurrence at SegDup loci; RepeatMasker stratification.

This is the core claim-stack folder (C13/C14/C22/C23 GLMs, C24 permutation,
the constitutional recurrence gradient, and all sensitivity/robustness checks).

## Planned files
- `replication_fragility_annotation.R` <- `remodeled_constitutional_AMR/post_processing/replication_fragility_annotation.R`
  (run first — annotates SV breakpoints with SegDup/LAD/Micro-C PC1/CFS overlap)
- `02_fragility_glm.R` <- `somatic_AMR/02_fragility_glm.R` (current C13/C14/C22/C23 weighted binomial GLM)
- `02b_c23_sensitivity.R` <- `somatic_AMR/02b_c23_sensitivity.R` (C23 distance-cutoff sensitivity: 1/10/50/100kb)
- `03_segdup_permutation.R` <- `somatic_AMR/03_segdup_permutation.R` (circular-shift null, current)
- `05_somatic_vs_constitutional.R` <- `somatic_AMR/05_somatic_vs_constitutional.R`
  (3-way SegDup OR: somatic / recurrently-remodeled-constitutional / recurrently-remodeled-normal; stacked interaction model + bootstrap delta-OR)
- `agold1_proximity_decircularize.R` <- `remodeled_constitutional_AMR/post_processing/agold1_proximity_decircularize.R`
  (A_GOLD1 — constitutional aDMR SegDup OR by recurrence depth All/R3/R5 = 1.32/1.85/4.22; Fig 4B — key recurrence-gradient script)
- `fragility_multivariate.R` <- `remodeled_constitutional_AMR/post_processing/fragility_multivariate.R` (SV multivariate GLM: SegDup+LAD+PC1+RepeatMasker)
- `cpg_adjusted_glm.R` <- `remodeled_constitutional_AMR/post_processing/cpg_adjusted_glm.R` (A1 — verify vs 02_fragility_glm.R, may be superseded)
- `admr_sv_exclusion.R` <- `remodeled_constitutional_AMR/post_processing/admr_sv_exclusion.R` (A3 sensitivity)
- `segdup_permutation_control.R` <- `remodeled_constitutional_AMR/post_processing/segdup_permutation_control.R` (A2 — verify vs 03_segdup_permutation.R, may be superseded)
- `aeq1_or_equivalence.R` <- `remodeled_constitutional_AMR/post_processing/aeq1_or_equivalence.R` (A_EQ1, formal OR-equivalence test SV vs aDMR)
- `admr_independent_fragility.R` <- `remodeled_constitutional_AMR/post_processing/admr_independent_fragility.R`
- `nahr_microhomology.R` <- `remodeled_constitutional_AMR/post_processing/nahr_microhomology.R` (NAHR proxy: SV size, Hi-C PC1 compartment shift)
- `cfs_overlap.R` <- `remodeled_constitutional_AMR/post_processing/cfs_overlap.R` (Common Fragile Site null control)
- `figS6_repeatmask.R` <- `remodeled_constitutional_AMR/viz/v4/figS6_repeatmask.R`
  (analysis+figure combined: repeat-density quartile OR, Cochran-Q, LINE depletion — kept here rather than figures/ since it computes the stats, not just plots them)
- `replication_timing_admr.R` <- `remodeled_constitutional_AMR/post_processing/replication_timing_admr.R` (Repli-seq phase + Micro-C PC1 at aDMR/SV loci)
- `11_replication_timing.R` <- `remodeled_constitutional_AMR/pipeline/11_replication_timing.R`
- `download_encode_repliseq.sh` <- `remodeled_constitutional_AMR/post_processing/download_encode_repliseq.sh`
