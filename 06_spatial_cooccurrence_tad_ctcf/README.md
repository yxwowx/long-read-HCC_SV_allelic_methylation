# 06_spatial_cooccurrence_tad_ctcf

Manuscript: Results > Broad SV-methylation co-occurrence: aDMRs enrich near SV
breakpoints at regional rather than focal scales; TAD/CTCF-disrupting SVs do
not show an allelic methylation magnitude gradient.

## Files
- `sv_dmr_enrichment.R` — Layer 2 window-permutation enrichment, 10kb-1Mb sweep, per-patient FDR
- `tad_ctcf_validation.R` — Analyses A/B/C mechanistic follow-up (CTCF-anchor precision, TAD-level enrichment, distance-to-boundary trend)
- `compute_insulation.py` — HepG2 Micro-C insulation score computation
- `insulation_sv_dmr.R` — Analysis D: reference insulation score at SV breakpoints and correlation with nearby DMR count
- `insulation_lme.R` — P2-7: continuous insulation-score LME vs the categorical tier model (robustness check)
- `trans_negative_control.R` — **canonical**: window-scale enrichment decay (50kb vs 1Mb paired Wilcoxon, one-sample at 1Mb) + aDMR-level proximity-zone check. Reproduces the manuscript's reported trans-negative-control numbers (p=0.098, p=1.5e-4).
- `07_trans_negative_control.R` — **superseded, kept for provenance only**. An earlier KS/JT-based distance-class test on a different input file; does not match the manuscript-reported numbers. Do not use for reported figures — see `trans_negative_control.R` instead.
- `power_calibration.R` — per-patient power calibration for non-significant patients (P3/P11/P12 limiting cases)
- `compare_window_runs.R` — combines/compares two fixed enrichment runs (e.g. primary_50kb vs primary_100kb)
- `plot_window_enrichment.R` — generic plot of enrichment_ratio ~ window_kb across all runs found in a directory
- `run_enrichment.sh` — parameterized runner for `sv_dmr_enrichment.R` (modes: `window`, `compare`, `confident`); consolidates what were previously three separate scripts (`run_window_enrichment.sh`, `run_compare_50kb_100kb.sh`, `run_confident_dmr_enrichment.sh`)
