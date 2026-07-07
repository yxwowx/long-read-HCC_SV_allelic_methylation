# 09_hbv_integration

Manuscript: Results > HBV integration co-occurs with allele-specific methylation
disruption near the integration site.

## Files
- `hbv_analysis.R` — HBV integration characterization: BED clustering, somatic filter, SV-HBV positional match, phase-block HP |Δβ|, genome integration map, SV-tier enrichment, distance-decay
- `hbv_allele_anchored.R` — Level B (C33): allele-anchored deviation from normal methylation for HBV-proximal phase blocks, + A-III covariate-adjusted LMM and matched-resample sensitivity
- `hbv_allele_specific.R` — Level A (existing per-read output check) + Level B re-anchored to somatic aDMR with normal_beta from the genome-wide DSS-smoothed haplotype model
- `hbv_perread_meth.py` — Level A (C32): per-read 5mC methylation at HBV junctions, HBV-carrying vs reference reads across distance bins