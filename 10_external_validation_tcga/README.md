# 10_external_validation_tcga

Manuscript: Results > External validation of SV-SegDup fragility and
allele-specific methylation associations using TCGA-LIHC.

## Files
- `tcga_validation.R` — Fig 8: TCGA-LIHC promoter methylation validation of Gold-tier genes (Tumor vs Normal, meth-expr correlation, survival, allelic concordance)
- `tcga_lihc.R` — P0-3: Gold-tier loci × TCGA-LIHC HM450, SV+ (CNA proxy) vs SV- sample comparison
- `tcga_lihc_shared_fragility.R` — V1: two-axis replication of shared SegDup fragility model (TCGA CNA breakpoints + 450K T-N DMR vs SegDup)
- `tcga_scna_meth_comparison.R` — Per Gold-tier locus, SCNA+ vs SCNA- TCGA tumor promoter β comparison
- `tcga_circularity_sensitivity.R` — Reviewer-response sensitivity analysis addressing circularity between array/NGS detection sensitivity and SegDup enrichment
- `a_dup1_dup_power_tcga.R` — A_DUP1: formal power calculation for DUP/INV distance-decay + TCGA CN-gain breakpoint proximity test
- `admr_normal_tissue_variance.R` — P1-D: normal-liver inter-sample β-variance at Gold/Silver aDMR loci (tests constitutive vs cancer-induced instability)
- `hm450_segdup_probe_density.R` — C17 probe-density audit quantifying HM450 probe depletion in SegDup regions

Some scripts (`hm450_segdup_probe_density.R`, `tcga_circularity_sensitivity.R`) additionally use `REFERENCE_DIR` for SegDup/LAD/PC1/FASTA-index reference files.
