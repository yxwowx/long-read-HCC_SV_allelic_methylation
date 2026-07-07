# 07_phased_pairing_concordance

| Claim | Manuscript | Re-derived |
|---|---|---|
| Overall concordance | 51.9%, GLMER OR=1.088 [1.015-1.166], p=0.017 | 51.9%, OR=1.088 [1.015-1.166], p=0.0173 |
| Distance bins | 51.2 / 51.9 / 52.5 / 51.8% | 51.2 / 51.9 / 52.5 / 51.8% (exact) |
| B-compartment | 53.9%, n=22,056, OR=1.172 [1.045-1.316], p=0.0068 | 53.9%, n=22,062, OR=1.173 [1.045-1.316], p=0.0066 |
| A-compartment | 49.5%, n=14,781, OR=0.990 [0.917-1.068], p=0.79 | 49.5%, n=14,792, OR=0.988 [0.916-1.066], p=0.75 |
| ICC / DEFF / n_eff / MDES | 0.0041 / 14.38 / 2,727 / 52.7% | 0.0041 / 14.38 / 2727.0 / 52.68% (exact) |

## Files
- `phaseblock_pairing.R` — builds SV-aDMR phase-block pairs (Gold/Silver/Bronze tiers)
- `haplotype_concordance.R` — HYPO-concordance per patient / distance / boundary tier (Fig S6); reads
  `phaseblock_pairs.csv` (confirmed the correct current pairing per `docs/OPEN_ISSUES.md` issue #1)
- `compute_within_block_hp_delta.R` — within-block HP |Δβ| computation
- `compartment_stratified_concordance.R` — A/B compartment (PC1)-stratified concordance
- `concordance_power_icc.R` — ICC/DEFF/effective-n power calculation for the concordance test
- `haplotype_sv_admr_analysis.R` — SV-aDMR haplotype analysis (CRE/promoter-annotated)
- `sv_clonality_somatic.R` — somatic-context CCF x HYPO-concordance cross-analysis (SVclone CCF
  from `02_sv_ccf_svclone` joined with this folder's phase-block pairs); previously duplicated
  in `02_sv_ccf_svclone/`, kept here as the single canonical copy since it is primarily a
  concordance analysis