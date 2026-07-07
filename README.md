# Allele-Specific Epigenetic Associations with Somatic Structural Variants in HBV+ HCC

Analysis code for: *Allele-Specific Epigenetic Associations with Somatic Structural
Variants in Hepatocellular Carcinoma Using Long-Read Sequencing*

## Status

Directory structure mirrors the Methods/Results flow of the manuscript.
Analysis scripts for sections `01`-`10` (and `shared/`) are copied in; each
folder's `README.md` documents provenance and notes superseded/verify-flagged
files. Two pieces remain outstanding:
- `figures/` — figure-generation scripts not yet populated; still need to be
  copied in and manually re-mapped to final Fig/FigS numbers.
- `04_admr_definition/constitutional/` — no script maps directly into this
  subfolder; the relevant analysis (Gold/Silver recurrence-depth staircase)
  currently lives in `05_fragility_enrichment/` instead (see that folder's
  README).

## Layout

| Folder | Manuscript section |
|---|---|
| `shared/` | Common utilities sourced by all pipeline scripts |
| `01_sv_calling_tiering/` | Somatic Mutation Calling — SV topological tiers, CRE classification |
| `02_sv_ccf_svclone/` | SV cancer cell fraction (SVclone) |
| `03_phasing_qc/` | Haplotype phasing and haplotagging — phasing quality assessment |
| `04_admr_definition/` | Recurrently remodeled constitutional aDMR loci and somatic aDMR definition |
| `05_fragility_enrichment/` | Genomic fragility annotation (SegDup/LAD/B-compartment GLMs, permutation) |
| `06_spatial_cooccurrence_tad_ctcf/` | Broad SV-aDMR co-occurrence; TAD/CTCF topological tiers |
| `07_phased_pairing_concordance/` | Phased SV-aDMR pair construction and haplotype concordance |
| `08_locus_matched_lme/` | Genome-wide locus-matched LME |
| `09_hbv_integration/` | HBV integration per-read and allele-anchored analyses |
| `10_external_validation_tcga/` | External validation in TCGA-LIHC |
| `figures/` | Main (Fig 1-5) and supplementary (Fig S1-S8) figure scripts |

## Known open issues

Three provenance items need to be resolved against the live pipeline before
these scripts can be treated as fully reproducing the manuscript's reported
numbers (concordance pairing staleness, dual Gold/Silver tiering systems,
survival-analysis script version).
