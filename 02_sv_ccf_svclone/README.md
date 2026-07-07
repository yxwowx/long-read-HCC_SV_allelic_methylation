# 02_sv_ccf_svclone

SVclone workflow (prepare -> filter -> cluster -> postassign) with PURPLE-supplied
purity/ploidy/CN segmentation.

## Files
- `prep_svclone_longread.py` — Severus SV VCF + PURPLE CNV/purity -> SVclone input files (long-read adaptation)
- `run_svclone_ccf.sh` — runs prep -> filter -> cluster -> postassign for all patients, then collects CCF results
- `svclone_config.ini` — SVclone parameters (long-read adapted: `skip_norm_adjust`, `filter_outliers` disabled, etc.)
- `sv_clonality.R` — SV clonality analysis using SVclone CCF estimates (SegDup vs non-SegDup CCF/timing)

The somatic-context CCF x HYPO-concordance cross-analysis (combining this
folder's CCF output with phase-block pairs) lives in
`07_phased_pairing_concordance/sv_clonality_somatic.R` rather than here, since
it is primarily a concordance analysis; it was previously duplicated in both
folders and has been consolidated to a single copy.
