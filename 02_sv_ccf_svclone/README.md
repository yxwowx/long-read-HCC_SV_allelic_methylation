# 02_sv_ccf_svclone

Manuscript: Methods > SV cancer cell fraction.
SVclone workflow (prepare -> filter -> cluster -> postassign) with PURPLE-supplied
purity/ploidy/CN segmentation. Downstream subclonality claim: SegDup-overlapping
SVs are late/subclonal (median CCF 0.261 vs 0.348 non-SegDup).

## Files
- `prep_svclone_longread.py` <- `remodeled_constitutional_AMR/post_processing/prep_svclone_longread.py` — **copied**
- `run_svclone_ccf.sh` <- `remodeled_constitutional_AMR/post_processing/run_svclone_ccf.sh` — **copied**
- `svclone_config.ini` <- `remodeled_constitutional_AMR/post_processing/svclone_config.ini` — **copied**
- `sv_clonality.R` <- `remodeled_constitutional_AMR/post_processing/sv_clonality.R` (in-house CCF baseline) — **copied**
- `sv_clonality_somatic.R` <- `somatic_AMR/post_processing/sv_clonality_somatic.R` (current: somatic-context CCF + clonal-restricted concordance, CCF>=0.8) — **copied** (also present in `07_phased_pairing_concordance/`, which reuses it for the clonal-restricted concordance check)
