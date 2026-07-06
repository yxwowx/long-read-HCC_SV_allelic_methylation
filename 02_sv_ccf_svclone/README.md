# 02_sv_ccf_svclone

Manuscript: Methods > SV cancer cell fraction.
SVclone workflow (prepare -> filter -> cluster -> postassign) with PURPLE-supplied
purity/ploidy/CN segmentation. Downstream subclonality claim: SegDup-overlapping
SVs are late/subclonal (median CCF 0.261 vs 0.348 non-SegDup).

## Planned files
- `prep_svclone_longread.py` <- `remodeled_constitutional_AMR/post_processing/prep_svclone_longread.py`
- `run_svclone_ccf.sh` <- `remodeled_constitutional_AMR/post_processing/run_svclone_ccf.sh`
- `svclone_config.ini` <- `remodeled_constitutional_AMR/post_processing/svclone_config.ini`
- `sv_clonality.R` <- `remodeled_constitutional_AMR/post_processing/sv_clonality.R` (in-house CCF baseline)
- `sv_clonality_somatic.R` <- `somatic_AMR/post_processing/sv_clonality_somatic.R` (current: somatic-context CCF + clonal-restricted concordance, CCF>=0.8)
