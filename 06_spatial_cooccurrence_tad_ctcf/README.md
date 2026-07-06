# 06_spatial_cooccurrence_tad_ctcf

Manuscript: Results > Broad SV-methylation co-occurrence: aDMRs enrich near SV
breakpoints at regional rather than focal scales; TAD/CTCF-disrupting SVs do
not show an allelic methylation magnitude gradient.

## Planned files
- `03_sv_dmr_enrichment.R` <- `remodeled_constitutional_AMR/pipeline/03_sv_dmr_enrichment.R`
  (Layer 2 window-permutation enrichment, 10kb-1Mb sweep, per-patient FDR)
- `06_tad_ctcf_validation.R` <- `remodeled_constitutional_AMR/pipeline/06_tad_ctcf_validation.R` (tier-magnitude Jonckheere-Terpstra tests)
- `06b_compute_insulation.py` <- `remodeled_constitutional_AMR/pipeline/06b_compute_insulation.py` (HepG2 Micro-C insulation score)
- `06c_insulation_sv_dmr.R` <- `remodeled_constitutional_AMR/pipeline/06c_insulation_sv_dmr.R`
- `13_insulation_lme.R` <- `remodeled_constitutional_AMR/pipeline/13_insulation_lme.R`
- `07_trans_negative_control.R` <- `remodeled_constitutional_AMR/pipeline/07_trans_negative_control.R`
- `trans_negative_control.R` <- `remodeled_constitutional_AMR/post_processing/trans_negative_control.R` (window + distance-band cis/mid/trans decay)
- `power_calibration.R` <- `remodeled_constitutional_AMR/post_processing/power_calibration.R` (per-patient power, P3/P11/P12 limiting cases)
- `compare_window_runs.R`, `plot_window_enrichment.R`, `run_window_enrichment.sh`, `run_confident_dmr_enrichment.sh`, `run_compare_50kb_100kb.sh`
  <- `remodeled_constitutional_AMR/post_processing/{same names}` (window-size sweep orchestration, 50kb vs 1Mb comparison)
