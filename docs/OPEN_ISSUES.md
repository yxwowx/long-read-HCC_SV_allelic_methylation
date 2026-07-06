# Open issues — resolve before treating the repo as manuscript-final

## 1. RESOLVED — concordance numbers reproduce; the confusion was which script/file to use

**Root cause found:** `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R`
builds its *own* pairing from `sv_tad_ctcf_annotation.v2.csv.gz` +
`confident_dmr_per_patient.csv.gz` (the constitutional-inclusive DMR set) —
its output `sv_admr_hp_concordance_pairs.csv.gz` has exactly **8,616 rows**,
confirmed by direct row count. That is precisely the manuscript's explicitly
-named superseded pairing ("an earlier n=8,616 pairing ... that is no longer
used"). Three sibling scripts (`concordance_distance_diagnostic.R`,
`sv_admr_hp_concordance_clonal.R`, `concordance_power_icc.R`) all read this
same superseded file — none of them reproduce the manuscript.

The **correct current pairing** is `SV_aDMR/phaseblock_pairs.csv` (39,213
lines incl. header = 39,212 pairs, produced by `somatic_AMR/04_phaseblock_pairing.R`,
confirmed by direct row count). The correct current concordance script is
`somatic_AMR/viz/v1/figS6_haplotype_concordance.R`, which reads
`phaseblock_pairs.csv` and correctly recomputes `dmatch = sv_minus_wt < 0`
(bypassing the buggy `direction_match` column) — its *code* is right, even
though its own header comment and `figure_captions.md` cite stale cached
numbers from an earlier run (48.2%->54.9% "distance paradox" instead of the
current flat 51.2-52.5%).

**Verification performed:** re-ran the full statistical chain (overall
GLMER, distance-bin Clopper-Pearson + GLMER, HepG2 Micro-C PC1 compartment
join + GLMER, ICC/DEFF/n_eff/MDES) directly against `phaseblock_pairs.csv`.
Every number reproduces the manuscript exactly or within rounding noise of a
few pairs out of tens of thousands (see the table in
`07_phased_pairing_concordance/README.md`). No compartment-stratified or
power/ICC script existed anywhere in the source tree pointed at the correct
file — both were written fresh (`compartment_stratified_concordance.R`,
`concordance_power_icc.R` in `07_phased_pairing_concordance/`) as forks of
the logic in the superseded originals, re-pointed at the current file.

**Action for the source repo (not yet done, needs your OK):** the four
superseded scripts in `remodeled_constitutional_AMR/post_processing/` still
point at the wrong file and will keep producing misleading output if anyone
reruns them. Worth patching their `PAIRS_FILE`/`ADM_FILE` constants in place,
or at minimum adding a header warning, so the *source* repo doesn't silently
diverge from the manuscript too.

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
