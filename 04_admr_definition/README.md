# 04_admr_definition

Recurrently remodeled constitutional aDMR loci and somatic aDMR definition.

Two distinct, non-mutually-exclusive populations built from different parent
sets — do not compare counts across them directly:
- `constitutional/`: tumor aDMR requiring bulk tumor-normal DMR overlap (>=30%
  reciprocal, >=100bp); further stratified by cross-patient recurrence depth.
  (Directory currently empty in this checkout — populate before running the
  somatic/ scripts, which depend on its outputs.)
- `somatic/`: tumor aDMR with <10% overlap to any matched-normal aDMR (constitutional
  signal explicitly excluded).

## Files
- `dmr_recurrence_analysis.R` — builds confident/consensus DMR calls and cross-patient recurrence (feeds `constitutional/`)
- `segdup_admr_cpg_variability.R` — P1-B: CpG density and allelic methylation variability, SegDup vs non-SegDup Gold aDMR
- `segdup_coverage_bias.R` — mosdepth-based mapping-bias check at SegDup regions
- `admr_hp_coverage_symmetry.py` — HP1/HP2 coverage symmetry at aDMR loci (mappability-confound check)
- `somatic/00_somatic_admr_annotate.R` — annotates somatic aDMR with ov_bulk, SegDup/LAD/B-compartment, recurrence
- `somatic/00b_tvn_validation.R` — Tumor-vs-Normal |HP Δβ| table (Fig 2A)
- `somatic/01_distribution.R` — 1Mb bin density, SV type x CNV cross-table, aDMR width/nCG summary
- `somatic/a_bulk1_allele_dilution.R` — A_BULK1: allele-specific vs simulated-bulk detection power comparison
- `somatic/const_amr_threshold_sensitivity.R` — constitutional-overlap threshold sensitivity + normal aDMR stats
- `somatic/somatic_bidir_check.R` — HP1/HP2 direction breakdown and bulk-detectability tau sweep for somatic aDMR