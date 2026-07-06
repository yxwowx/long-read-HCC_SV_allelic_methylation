# Script mapping (old path -> new path)

All old paths are relative to `~/script/SV-DMR/`. Status: `ok` = confirmed
current/authoritative; `verify` = flagged in `OPEN_ISSUES.md`, confirm before
treating as final.

| New path | Old path | Status |
|---|---|---|
| `shared/shared_utils.R` | `shared_file/pipeline/shared_utils.R` | ok |
| `01_sv_calling_tiering/02_sv_annotation_stratify.R` | `remodeled_constitutional_AMR/pipeline/02_sv_annotation_stratify.R` | ok |
| `02_sv_ccf_svclone/prep_svclone_longread.py` | `remodeled_constitutional_AMR/post_processing/prep_svclone_longread.py` | ok |
| `02_sv_ccf_svclone/run_svclone_ccf.sh` | `remodeled_constitutional_AMR/post_processing/run_svclone_ccf.sh` | ok |
| `02_sv_ccf_svclone/svclone_config.ini` | `remodeled_constitutional_AMR/post_processing/svclone_config.ini` | ok |
| `02_sv_ccf_svclone/sv_clonality.R` | `remodeled_constitutional_AMR/post_processing/sv_clonality.R` | ok |
| `02_sv_ccf_svclone/sv_clonality_somatic.R` | `somatic_AMR/post_processing/sv_clonality_somatic.R` | ok (current) |
| `03_phasing_qc/phasing_quality.R` | `remodeled_constitutional_AMR/post_processing/phasing_quality.R` | ok |
| `03_phasing_qc/imprinting_hp_validation.py` | `remodeled_constitutional_AMR/post_processing/imprinting_hp_validation.py` | ok |
| `03_phasing_qc/data/imprinted_dmrs_hg38.bed` | `remodeled_constitutional_AMR/post_processing/data/imprinted_dmrs_hg38.bed` | ok |
| `04_admr_definition/01_dmr_recurrence_analysis.R` | `remodeled_constitutional_AMR/pipeline/01_dmr_recurrence_analysis.R` | ok |
| `04_admr_definition/admr_hp_coverage_symmetry.py` | `remodeled_constitutional_AMR/post_processing/admr_hp_coverage_symmetry.py` | ok |
| `04_admr_definition/segdup_coverage_bias.R` | `remodeled_constitutional_AMR/post_processing/segdup_coverage_bias.R` | ok |
| `04_admr_definition/segdup_admr_cpg_variability.R` | `remodeled_constitutional_AMR/post_processing/segdup_admr_cpg_variability.R` | ok |
| `04_admr_definition/somatic/00_somatic_admr_annotate.R` | `somatic_AMR/00_somatic_admr_annotate.R` | ok |
| `04_admr_definition/somatic/00b_tvn_validation.R` | `somatic_AMR/00b_tvn_validation.R` | ok |
| `04_admr_definition/somatic/01_distribution.R` | `somatic_AMR/01_distribution.R` | ok |
| `04_admr_definition/somatic/const_amr_threshold_sensitivity.R` | `somatic_AMR/post_processing/const_amr_threshold_sensitivity.R` | ok |
| `04_admr_definition/somatic/somatic_bidir_check.R` | `somatic_AMR/somatic_bidir_check.R` | ok |
| `04_admr_definition/somatic/a_bulk1_allele_dilution.R` | `remodeled_constitutional_AMR/post_processing/a_bulk1_allele_dilution.R` | ok |
| `05_fragility_enrichment/replication_fragility_annotation.R` | `remodeled_constitutional_AMR/post_processing/replication_fragility_annotation.R` | ok |
| `05_fragility_enrichment/02_fragility_glm.R` | `somatic_AMR/02_fragility_glm.R` | ok (current) |
| `05_fragility_enrichment/02b_c23_sensitivity.R` | `somatic_AMR/02b_c23_sensitivity.R` | ok |
| `05_fragility_enrichment/03_segdup_permutation.R` | `somatic_AMR/03_segdup_permutation.R` | ok (current) |
| `05_fragility_enrichment/05_somatic_vs_constitutional.R` | `somatic_AMR/05_somatic_vs_constitutional.R` | ok |
| `05_fragility_enrichment/agold1_proximity_decircularize.R` | `remodeled_constitutional_AMR/post_processing/agold1_proximity_decircularize.R` | ok |
| `05_fragility_enrichment/fragility_multivariate.R` | `remodeled_constitutional_AMR/post_processing/fragility_multivariate.R` | ok |
| `05_fragility_enrichment/cpg_adjusted_glm.R` | `remodeled_constitutional_AMR/post_processing/cpg_adjusted_glm.R` | verify (may be superseded by 02_fragility_glm.R) |
| `05_fragility_enrichment/admr_sv_exclusion.R` | `remodeled_constitutional_AMR/post_processing/admr_sv_exclusion.R` | ok |
| `05_fragility_enrichment/segdup_permutation_control.R` | `remodeled_constitutional_AMR/post_processing/segdup_permutation_control.R` | verify (may be superseded by 03_segdup_permutation.R) |
| `05_fragility_enrichment/aeq1_or_equivalence.R` | `remodeled_constitutional_AMR/post_processing/aeq1_or_equivalence.R` | ok |
| `05_fragility_enrichment/admr_independent_fragility.R` | `remodeled_constitutional_AMR/post_processing/admr_independent_fragility.R` | ok |
| `05_fragility_enrichment/nahr_microhomology.R` | `remodeled_constitutional_AMR/post_processing/nahr_microhomology.R` | ok |
| `05_fragility_enrichment/cfs_overlap.R` | `remodeled_constitutional_AMR/post_processing/cfs_overlap.R` | ok |
| `05_fragility_enrichment/figS6_repeatmask.R` | `remodeled_constitutional_AMR/viz/v4/figS6_repeatmask.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/03_sv_dmr_enrichment.R` | `remodeled_constitutional_AMR/pipeline/03_sv_dmr_enrichment.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/06_tad_ctcf_validation.R` | `remodeled_constitutional_AMR/pipeline/06_tad_ctcf_validation.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/06b_compute_insulation.py` | `remodeled_constitutional_AMR/pipeline/06b_compute_insulation.py` | ok |
| `06_spatial_cooccurrence_tad_ctcf/06c_insulation_sv_dmr.R` | `remodeled_constitutional_AMR/pipeline/06c_insulation_sv_dmr.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/13_insulation_lme.R` | `remodeled_constitutional_AMR/pipeline/13_insulation_lme.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/07_trans_negative_control.R` | `remodeled_constitutional_AMR/pipeline/07_trans_negative_control.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/trans_negative_control.R` | `remodeled_constitutional_AMR/post_processing/trans_negative_control.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/power_calibration.R` | `remodeled_constitutional_AMR/post_processing/power_calibration.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/compare_window_runs.R` | `remodeled_constitutional_AMR/post_processing/compare_window_runs.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/plot_window_enrichment.R` | `remodeled_constitutional_AMR/post_processing/plot_window_enrichment.R` | ok |
| `06_spatial_cooccurrence_tad_ctcf/run_window_enrichment.sh` | `remodeled_constitutional_AMR/post_processing/run_window_enrichment.sh` | ok |
| `06_spatial_cooccurrence_tad_ctcf/run_confident_dmr_enrichment.sh` | `remodeled_constitutional_AMR/post_processing/run_confident_dmr_enrichment.sh` | ok |
| `06_spatial_cooccurrence_tad_ctcf/run_compare_50kb_100kb.sh` | `remodeled_constitutional_AMR/post_processing/run_compare_50kb_100kb.sh` | ok |
| `07_phased_pairing_concordance/04_phaseblock_pairing.R` | `somatic_AMR/04_phaseblock_pairing.R` | ok (current) |
| `07_phased_pairing_concordance/04_haplotype_sv_admr_analysis.R` | `remodeled_constitutional_AMR/pipeline/04_haplotype_sv_admr_analysis.R` | ok |
| `07_phased_pairing_concordance/sv_admr_hp_concordance.R` | `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R` | **verify** (open issue #1) |
| `07_phased_pairing_concordance/concordance_distance_diagnostic.R` | `remodeled_constitutional_AMR/post_processing/concordance_distance_diagnostic.R` | verify (depends on #1) |
| `07_phased_pairing_concordance/sv_admr_hp_concordance_clonal.R` | `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance_clonal.R` | verify (depends on #1) |
| `07_phased_pairing_concordance/compute_within_block_hp_delta.R` | `remodeled_constitutional_AMR/post_processing/compute_within_block_hp_delta.R` | ok |
| `07_phased_pairing_concordance/concordance_power_icc.R` | `remodeled_constitutional_AMR/post_processing/concordance_power_icc.R` | ok |
| `08_locus_matched_lme/locus_matched_sv_lme.py` | `remodeled_constitutional_AMR/post_processing/locus_matched_sv_lme.py` | ok |
| `09_hbv_integration/12_hbv_analysis.R` | `remodeled_constitutional_AMR/pipeline/12_hbv_analysis.R` | ok |
| `09_hbv_integration/hbv_perread_meth.py` | `remodeled_constitutional_AMR/post_processing/hbv_perread_meth.py` | ok |
| `09_hbv_integration/hbv_allele_anchored.R` | `remodeled_constitutional_AMR/post_processing/hbv_allele_anchored.R` | superseded-but-kept (see 06_hbv_allele_specific.R) |
| `09_hbv_integration/06_hbv_allele_specific.R` | `somatic_AMR/06_hbv_allele_specific.R` | ok (current) |
| `10_external_validation_tcga/tcga_lihc.R` | `remodeled_constitutional_AMR/external_validation/tcga_lihc.R` | ok |
| `10_external_validation_tcga/tcga_lihc_shared_fragility.R` | `remodeled_constitutional_AMR/external_validation/tcga_lihc_shared_fragility.R` | ok |
| `10_external_validation_tcga/tcga_circularity_sensitivity.R` | `remodeled_constitutional_AMR/external_validation/tcga_circularity_sensitivity.R` | ok |
| `10_external_validation_tcga/tcga_scna_meth_comparison.R` | `remodeled_constitutional_AMR/post_processing/tcga_scna_meth_comparison.R` | ok |
| `10_external_validation_tcga/hm450_segdup_probe_density.R` | `remodeled_constitutional_AMR/post_processing/hm450_segdup_probe_density.R` | ok |
| `10_external_validation_tcga/admr_normal_tissue_variance.R` | `remodeled_constitutional_AMR/post_processing/admr_normal_tissue_variance.R` | ok |
| `10_external_validation_tcga/a_dup1_dup_power_tcga.R` | `remodeled_constitutional_AMR/post_processing/a_dup1_dup_power_tcga.R` | ok |
| `10_external_validation_tcga/fig8_tcga_validation.R` (name TBD) | `remodeled_constitutional_AMR/viz/v1/fig8_tcga_validation.R` OR `viz/v2/fig5_tcga_validation.R` | **verify** (open issue #3) |
| `05_fragility_enrichment/replication_timing_admr.R` | `remodeled_constitutional_AMR/post_processing/replication_timing_admr.R` | ok, supporting/Repli-seq compartment context |
| `05_fragility_enrichment/11_replication_timing.R` | `remodeled_constitutional_AMR/pipeline/11_replication_timing.R` | ok |
| `05_fragility_enrichment/download_encode_repliseq.sh` | `remodeled_constitutional_AMR/post_processing/download_encode_repliseq.sh` | ok |
| `figures/main/*`, `figures/supplementary/*` | `remodeled_constitutional_AMR/viz/v4/*` + `somatic_AMR/viz/v1/*` | ok, but re-map Fig/FigS numbers by hand |
| `figures/.../vignette_single_locus_zoom.R` | `remodeled_constitutional_AMR/post_processing/vignette_single_locus_zoom.R` | ok (Fig S7) |
