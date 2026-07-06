# Open issues — resolve before treating the repo as manuscript-final

## 1. Concordance pairing may be stale (confirmed by file timestamps)

- `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R` output
  `sv_admr_hp_concordance_pairs.csv.gz` was last generated **2026-06-02 18:25**.
- Its input `DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz` was
  regenerated **2026-06-10 09:25** — eight days *after* the concordance script
  last ran.
- `viz/v4/fig3_allele_specific_null.R`'s own header comment cites "51.5%,
  GLMER p=0.146" for this claim, but the manuscript text (v5.3) reports
  "51.9%, GLMER OR=1.088 [1.015-1.166], p=0.017."
- **Action:** re-run `sv_admr_hp_concordance.R` (and its dependents
  `concordance_distance_diagnostic.R`, `sv_admr_hp_concordance_clonal.R`)
  against the current `confident_dmr_per_patient.csv.gz` and confirm the
  51.9%/OR=1.088/p=0.017 numbers reproduce before including these in the
  final repo as-is.

## 2. Two coexisting Gold/Silver tiering systems

- `remodeled_constitutional_AMR/pipeline/05_sv_dmr_final_candidates.R` builds
  Gold/Silver/Bronze against the constitutional-inclusive `confident_dmr`
  set (old n=8,616 pairing — explicitly superseded per manuscript text).
- `somatic_AMR/04_phaseblock_pairing.R` builds a *different* Gold/Silver/Bronze
  against the somatic-only aDMR set (current n=39,212 — matches Fig 3B and
  the manuscript's Results text exactly).
- Several `post_processing/` scripts under `remodeled_constitutional_AMR`
  still read the *old* `04.final_candidate/gold_tier_final.csv` /
  `silver_tier.csv` path (e.g. `admr_normal_tissue_variance.R`,
  `nahr_microhomology.R`, `cfs_overlap.R`, `vignette_single_locus_zoom.R`,
  `aeq1_or_equivalence.R`, `agold1_proximity_decircularize.R`,
  `a_bulk1_allele_dilution.R`). This is very likely *correct* — these serve
  the "recurrently remodeled constitutional aDMR" claim stack, a legitimately
  separate population from the somatic aDMR Gold/Silver/Bronze in Fig 3 — but
  it must not be assumed; confirm the Gold/Silver population intended by each
  claim (constitutional vs. somatic) matches the file it actually reads.
- **Action:** before copying `post_processing/` scripts in, annotate each
  with which Gold/Silver system (constitutional-pipeline/05 vs.
  somatic-AMR/04) it consumes, to avoid silently mixing populations.

## 3. Survival analysis (14 Gold promoter genes) — script version unclear

- Manuscript text: "Survival analysis of 14 Gold-tier aDMR promoter genes in
  TCGA-LIHC found no gene reaching FDR < 0.05."
- The only scripts containing `survfit`/`survdiff` logic are
  `remodeled_constitutional_AMR/viz/v1/fig8_tcga_validation.R` and
  `viz/v2/fig5_tcga_validation.R` — both otherwise-superseded viz versions.
  Neither `viz/v4/fig5_validation_mechanism.R` nor any `external_validation/`
  script contains this analysis.
- **Action:** confirm which of the two (v1 or v2) actually produced the
  "14 genes, none FDR<0.05" number quoted in the manuscript, then pull only
  that one file into `10_external_validation_tcga/` (the survival-specific
  portion, not the whole viz script if it also contains superseded panels).

## Also worth a quick look

- `05_fragility_enrichment/cpg_adjusted_glm.R` (A1) and
  `segdup_permutation_control.R` (A2) look like earlier iterations of what
  `somatic_AMR/02_fragility_glm.R` and `03_segdup_permutation.R` now do.
  Compare outputs before deciding whether to keep both (e.g. for an
  appendix/history) or drop the older pair.
