# Deliberately excluded from this repo

| Item | Reason |
|---|---|
| `SV-DMR/denovo_SD/` (entire directory) | Separate project (assembly-diff-based de novo segmental duplication detection, T2T-CHM13). Not referenced anywhere in manuscript v5.3. |
| `remodeled_constitutional_AMR/viz/v1/`, `viz/v2/` | Superseded by `viz/v4/`, except the survival-analysis exception tracked in `OPEN_ISSUES.md` #3. |
| `remodeled_constitutional_AMR/pipeline/10_ase_admr_integration.R` | Manuscript Limitations explicitly states: "this study lacks transcriptomic (RNA-seq) data" — this RNA-seq/ASE linkage analysis is not part of the final manuscript. |
| `remodeled_constitutional_AMR/post_processing/ase_admr_link.R` | Same reason as above. |
| `remodeled_constitutional_AMR/post_processing/sv_cre_directionality.R` / `.py` | No explicit citation found in manuscript text; re-check if you believe this feeds a specific claim before dropping it permanently. |
| `remodeled_constitutional_AMR/post_processing/segdup_class_cooccurrence.R` | No explicit citation found in manuscript text (appears to be an earlier/parallel exploratory analysis). |
| `remodeled_constitutional_AMR/post_processing/admr_ccf_estimation.R` + `extract_admr_reads_meth.py` | aDMR-own CCF from per-read methylation — not explicitly cited in the manuscript text reviewed (the CCF claim in text is about *SV* CCF, covered by `sv_clonality.R` / `sv_clonality_somatic.R`). |
| `.claude/` directories, `logs/*.log` (except optionally `claude_decisions.log` for provenance) | Agent tooling config / run logs, not analysis code. |
| `Rplots.pdf`, `last.dump.rda`, `somatic_AMR/exploratory/hbv_recurrent_admr_memo.md` | Scratch/debug artifacts and an exploratory memo, not part of the citable pipeline. |
| `remodeled_constitutional_AMR/pipeline/diag_admr_overlap_matrix.R` | Ad-hoc diagnostic scratch script. |
| `remodeled_constitutional_AMR/pipeline/05_sv_dmr_final_candidates.R` | Superseded — produced the old n=8,616 SV-aDMR pairing that the manuscript explicitly says is "no longer used," replaced by `somatic_AMR/04_phaseblock_pairing.R` (n=39,212). Keep only if you want pre-reorg provenance/history. |

If you find a manuscript claim not covered by anything in `SCRIPT_MAPPING.md` and not listed here, flag it — it likely means a script exists but wasn't identified, or the analysis was done ad-hoc and never persisted as a script.
