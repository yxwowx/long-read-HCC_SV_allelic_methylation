# 03_phasing_qc

Manuscript: Methods > Haplotype phasing and haplotagging; Phasing quality assessment.
Block N50 / phased fraction (whatshap stats), reference-frame switch rate
(whatshap compare, hg38 vs hg38+HBV), and biological validation at 16 canonical
germline-imprinted control regions (ICRs).

## Files
- `phasing_quality.R` — A-I: whatshap stats (block N50/phased fraction) + whatshap compare (hg38 vs hg38+HBV switch rate)
- `imprinting_hp_validation.py` — A-II: HP1/HP2 |Δβ| at imprinted DMRs vs random background (biological positive control)
- `data/imprinted_dmrs_hg38.bed` — 16 canonical germline-imprinted control regions (hg38)
