# Allele-Specific Epigenetic Associations with Somatic Structural Variants in HBV+ HCC

Analysis code for: *Allele-Specific Epigenetic Associations with Somatic Structural
Variants in Hepatocellular Carcinoma Using Long-Read Sequencing*

## Status

Directory structure mirrors the Methods/Results flow of the manuscript.
Analysis scripts for sections `01`-`10` (and `shared/`) are copied in; each
folder's `README.md` documents provenance and notes superseded/verify-flagged
files.

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

## AI Usage Disclosure

During the development of this repository, Anthropic's Claude was utilized as an AI coding assistant to improve overall code quality and maintainability. 

* **Refactoring & File Organization:** Claude provided assistance in refactoring Python and R scripts, optimizing the directory structure, and enhancing code readability.
* **Automated Commits:** The appearance of the AI bot in the GitHub contributor list is strictly due to the automated generation and modification of project configuration files by the AI-integrated development environment. 
