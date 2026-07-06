# 07_phased_pairing_concordance

Manuscript: Methods > Phased SV-aDMR pair construction and haplotype
concordance; Results > Within SegDup co-localization, the SV-bearing allele
is not preferentially hypomethylated; Methods > Statistical analysis (power
paragraph).

**Open issue #1 is RESOLVED** (see `docs/OPEN_ISSUES.md`) — every number in
this section was independently re-derived against the current
`phaseblock_pairs.csv` (39,212 pairs) and matches the manuscript exactly or
within trivial rounding:

| Claim | Manuscript | Re-derived |
|---|---|---|
| Overall concordance | 51.9%, GLMER OR=1.088 [1.015-1.166], p=0.017 | 51.9%, OR=1.088 [1.015-1.166], p=0.0173 |
| Distance bins | 51.2 / 51.9 / 52.5 / 51.8% | 51.2 / 51.9 / 52.5 / 51.8% (exact) |
| B-compartment | 53.9%, n=22,056, OR=1.172 [1.045-1.316], p=0.0068 | 53.9%, n=22,062, OR=1.173 [1.045-1.316], p=0.0066 |
| A-compartment | 49.5%, n=14,781, OR=0.990 [0.917-1.068], p=0.79 | 49.5%, n=14,792, OR=0.988 [0.916-1.066], p=0.75 |
| ICC / DEFF / n_eff / MDES | 0.0041 / 14.38 / 2,727 / 52.7% | 0.0041 / 14.38 / 2727.0 / 52.68% (exact) |

## Files (all verified against the current n=39,212 pairing)
- `04_phaseblock_pairing.R` <- `somatic_AMR/04_phaseblock_pairing.R` (pairing + Gold/Silver/Bronze tiering; produces `phaseblock_pairs.csv`)
- `04_haplotype_sv_admr_analysis.R` <- `remodeled_constitutional_AMR/pipeline/04_haplotype_sv_admr_analysis.R` (`infer_sv_hp_map()` / `get_sv_hp_beta()` helpers, reused across the repo)
- `figS6_haplotype_concordance.R` <- `somatic_AMR/viz/v1/figS6_haplotype_concordance.R`
  (**primary concordance script**: overall + distance-bin + tier-stratified; reads `phaseblock_pairs.csv` directly and recomputes `dmatch = sv_minus_wt < 0` rather than trusting the buggy `direction_match` column — this is correct as written. NOTE: its own header comment and the sibling `figure_captions.md` cite stale numbers from an earlier run [48.2%->54.9% distance "paradox"]; the *code* is current, only the *cached documentation* is stale. Needs `shared_theme.R` from `somatic_AMR/viz/v1/` to run as-is — not yet copied into this repo, see `figures/README.md`.)
- `compartment_stratified_concordance.R` — **new**, written during verification. No persisted current version of this analysis existed anywhere in the source tree; the only script that did this join (`concordance_distance_diagnostic.R`) pointed at the superseded n=8,616 file. Reproduces the manuscript's B/A-compartment numbers (see table above).
- `concordance_power_icc.R` — **new** (re-pointed fork of `remodeled_constitutional_AMR/post_processing/concordance_power_icc.R`, which reads the superseded n=8,616 file and whose own header cites a stale OR=1.065/p=0.146). Reproduces ICC/DEFF/n_eff/MDES exactly.
- `sv_clonality_somatic.R` <- `somatic_AMR/post_processing/sv_clonality_somatic.R` (clonal CCF>=0.8-restricted concordance, current)
- `compute_within_block_hp_delta.R` <- `remodeled_constitutional_AMR/post_processing/compute_within_block_hp_delta.R` — **copied**
  (within-phase-block SV-aDMR distance + `hp_delta`, supplements `sv_admr_distance_within_block.csv.gz`; reads the constitutional-inclusive `confident_dmr_per_patient.csv.gz`, not `phaseblock_pairs.csv` — a distance-diagnostic utility, not itself a concordance number, so not in the reproduced-numbers table above)

## Superseded — do not use (see `docs/EXCLUDED.md`)
- `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R` — builds its own pairing from `confident_dmr_per_patient.csv.gz`, produces exactly 8,616 rows = the manuscript's explicitly-named superseded pairing.
- `remodeled_constitutional_AMR/post_processing/concordance_distance_diagnostic.R` — reads the same superseded 8,616-row file.
- `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance_clonal.R` — same, superseded by `sv_clonality_somatic.R`.
- `remodeled_constitutional_AMR/post_processing/concordance_power_icc.R` — same, superseded by the re-pointed version above.
