# 04_admr_definition

Manuscript: Methods > CpG methylation extraction and DMR calling;
Recurrently remodeled constitutional aDMR loci and somatic aDMR definition.

Two distinct, non-mutually-exclusive populations built from different parent
sets — do not compare counts across them directly:
- `constitutional/`: tumor aDMR requiring bulk tumor-normal DMR overlap (>=30%
  reciprocal, >=100bp); further stratified by cross-patient recurrence depth.
- `somatic/`: tumor aDMR with <10% overlap to any matched-normal aDMR (constitutional
  signal explicitly excluded).

## Planned files (repo root of this folder)
- `01_dmr_recurrence_analysis.R` <- `remodeled_constitutional_AMR/pipeline/01_dmr_recurrence_analysis.R`
  (3-stage filter: per-patient confidence -> recurrence n>=3 -> `confident_dmr_per_patient.csv.gz` / `consensus_dmrs_per_patient.csv.gz`)
- `admr_hp_coverage_symmetry.py` <- `remodeled_constitutional_AMR/post_processing/admr_hp_coverage_symmetry.py`
  (HP1/HP2 coverage-bias check at SegDup aDMR loci, Fig S11)
- `segdup_coverage_bias.R` <- `remodeled_constitutional_AMR/post_processing/segdup_coverage_bias.R` (mosdepth SegDup coverage QC)
- `segdup_admr_cpg_variability.R` <- `remodeled_constitutional_AMR/post_processing/segdup_admr_cpg_variability.R`

## constitutional/
- (Gold/Silver recurrence-depth staircase lives in `05_fragility_enrichment/`, see that folder's README — `agold1_proximity_decircularize.R`)

## somatic/
- `00_somatic_admr_annotate.R` <- `somatic_AMR/00_somatic_admr_annotate.R`
  (ov_bulk / segdup / lad / b_compartment / log10_nCG / n_patients annotation -> `somatic_admr_annotated.csv.gz`, ~2.96M loci)
- `00b_tvn_validation.R` <- `somatic_AMR/00b_tvn_validation.R` (tumor vs normal HP |delta beta|, Fig 2A/2C support)
- `01_distribution.R` <- `somatic_AMR/01_distribution.R` (1Mb bin density, SV x CNV cross-table)
- `const_amr_threshold_sensitivity.R` <- `somatic_AMR/post_processing/const_amr_threshold_sensitivity.R`
  (5-50% constitutional-overlap threshold sensitivity, confirms null is not threshold-dependent, Fig 4C)
- `somatic_bidir_check.R` <- `somatic_AMR/somatic_bidir_check.R` (bidirectional HP shift breakdown; supports the 54-fold yield-advantage discussion)
- `a_bulk1_allele_dilution.R` <- `remodeled_constitutional_AMR/post_processing/a_bulk1_allele_dilution.R`
  (C31: allele-specific vs bulk-simulated detection power, Fig 1a panels I/J)
