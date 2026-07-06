# 10_external_validation_tcga

Manuscript: Methods > External validation in TCGA-LIHC.

**Open issue #3 — RESOLVED** (see `docs/OPEN_ISSUES.md`): the 14
Gold-promoter-gene survival analysis (log-rank, BH-FDR) is computed by
`fig8_tcga_validation.R` Module 4c only — confirmed the only script with
`survfit`/`survdiff` logic anywhere in the source tree, and its cached output
(`tcga_survival_summary.csv`) matches the manuscript's "14 genes, none
FDR<0.05" exactly (14 rows, min FDR=0.367). `viz/v1/fig5_tcga_validation.R`
and `viz/v2/fig5_tcga_validation.R` are downstream visualization-only
consumers of that same cached CSV, not independent computations — both
remain excluded (`EXCLUDED.md`).

## Files
- `fig8_tcga_validation.R` <- `remodeled_constitutional_AMR/viz/v1/fig8_tcga_validation.R` — **copied**. Kept whole (not just Module 4c) since the survival block depends on TCGA download/methylation-processing objects (Modules 1-2) built earlier in the same script; its other panels (A/B/D/E) duplicate `viz/v4` coverage and can be trimmed later if desired, but aren't the reason this file is here.

## Planned files
- `tcga_lihc.R` <- `remodeled_constitutional_AMR/external_validation/tcga_lihc.R`
- `tcga_lihc_shared_fragility.R` <- `remodeled_constitutional_AMR/external_validation/tcga_lihc_shared_fragility.R`
- `tcga_circularity_sensitivity.R` <- `remodeled_constitutional_AMR/external_validation/tcga_circularity_sensitivity.R`
- `tcga_scna_meth_comparison.R` <- `remodeled_constitutional_AMR/post_processing/tcga_scna_meth_comparison.R`
- `hm450_segdup_probe_density.R` <- `remodeled_constitutional_AMR/post_processing/hm450_segdup_probe_density.R` (C17, HM450 probe-depletion audit)
- `admr_normal_tissue_variance.R` <- `remodeled_constitutional_AMR/post_processing/admr_normal_tissue_variance.R` (C21, TCGA normal-liver inter-sample SD)
- `a_dup1_dup_power_tcga.R` <- `remodeled_constitutional_AMR/post_processing/a_dup1_dup_power_tcga.R` (A_DUP1, power calc + TCGA CN-gain distance-bin test)
