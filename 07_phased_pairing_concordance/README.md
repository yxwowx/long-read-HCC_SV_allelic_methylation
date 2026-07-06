# 07_phased_pairing_concordance

Manuscript: Methods > Phased SV-aDMR pair construction and haplotype
concordance; Results > Within SegDup co-localization, the SV-bearing allele
is not preferentially hypomethylated.

**Before copying these in, resolve `docs/OPEN_ISSUES.md` issue #1** — the
concordance pairing may need to be re-run against the current
`confident_dmr_per_patient.csv.gz` (regenerated 2026-06-10, after the last
`sv_admr_hp_concordance.R` output on 2026-06-02).

## Planned files
- `04_phaseblock_pairing.R` <- `somatic_AMR/04_phaseblock_pairing.R`
  (current canonical n=39,212 Gold/Silver/Bronze pairing; Fig 3B)
- `04_haplotype_sv_admr_analysis.R` <- `remodeled_constitutional_AMR/pipeline/04_haplotype_sv_admr_analysis.R`
  (`infer_sv_hp_map()` / `get_sv_hp_beta()` helpers, reused across the repo)
- `sv_admr_hp_concordance.R` <- `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R` (A_HP1 — verify/re-run, see open issue #1)
- `concordance_distance_diagnostic.R` <- `remodeled_constitutional_AMR/post_processing/concordance_distance_diagnostic.R`
  (A_HP3 — distance-bin reversal diagnostic; A vs B compartment stratification)
- `sv_admr_hp_concordance_clonal.R` <- `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance_clonal.R` (clonal CCF>=0.8 restriction)
- `compute_within_block_hp_delta.R` <- `remodeled_constitutional_AMR/post_processing/compute_within_block_hp_delta.R`
- `concordance_power_icc.R` <- `remodeled_constitutional_AMR/post_processing/concordance_power_icc.R` (ICC/DEFF/n_eff/MDES power calc)
