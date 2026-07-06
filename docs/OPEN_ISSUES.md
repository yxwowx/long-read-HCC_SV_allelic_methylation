# Open issues — resolve before treating the repo as manuscript-final

## 1. RESOLVED — concordance numbers reproduce; the confusion was which script/file to use

**Root cause found:** `remodeled_constitutional_AMR/post_processing/sv_admr_hp_concordance.R`
builds its *own* pairing from `sv_tad_ctcf_annotation.v2.csv.gz` +
`confident_dmr_per_patient.csv.gz` (the constitutional-inclusive DMR set) —
its output `sv_admr_hp_concordance_pairs.csv.gz` has exactly **8,616 rows**,
confirmed by direct row count. That is precisely the manuscript's explicitly
-named superseded pairing ("an earlier n=8,616 pairing ... that is no longer
used"). Three sibling scripts (`concordance_distance_diagnostic.R`,
`sv_admr_hp_concordance_clonal.R`, `concordance_power_icc.R`) all read this
same superseded file — none of them reproduce the manuscript.

The **correct current pairing** is `SV_aDMR/phaseblock_pairs.csv` (39,213
lines incl. header = 39,212 pairs, produced by `somatic_AMR/04_phaseblock_pairing.R`,
confirmed by direct row count). The correct current concordance script is
`somatic_AMR/viz/v1/figS6_haplotype_concordance.R`, which reads
`phaseblock_pairs.csv` and correctly recomputes `dmatch = sv_minus_wt < 0`
(bypassing the buggy `direction_match` column) — its *code* is right, even
though its own header comment and `figure_captions.md` cite stale cached
numbers from an earlier run (48.2%->54.9% "distance paradox" instead of the
current flat 51.2-52.5%).

**Verification performed:** re-ran the full statistical chain (overall
GLMER, distance-bin Clopper-Pearson + GLMER, HepG2 Micro-C PC1 compartment
join + GLMER, ICC/DEFF/n_eff/MDES) directly against `phaseblock_pairs.csv`.
Every number reproduces the manuscript exactly or within rounding noise of a
few pairs out of tens of thousands (see the table in
`07_phased_pairing_concordance/README.md`). No compartment-stratified or
power/ICC script existed anywhere in the source tree pointed at the correct
file — both were written fresh (`compartment_stratified_concordance.R`,
`concordance_power_icc.R` in `07_phased_pairing_concordance/`) as forks of
the logic in the superseded originals, re-pointed at the current file.

**Action for the source repo — DONE (2026-07-06):** patched all five affected
scripts in `~/script/SV-DMR/remodeled_constitutional_AMR/post_processing/`
(`sv_admr_hp_concordance.R`, `concordance_distance_diagnostic.R`,
`concordance_power_icc.R`, `sv_admr_hp_concordance_clonal.R`, and
`segdup_stratified_hp_concordance.R` — the last one shares the same root
cause but wasn't in the original 3-issue list) with a prominent
"SUPERSEDED" header block pointing to the current replacement script/file,
plus a `message()` printed at runtime so a naive rerun can't silently
produce misleading numbers without a warning. Did not silently rewrite
their internals to the new pairing methodology in place — `04_phaseblock_pairing.R`'s
pairing/tiering logic differs enough (different aDMR input set, `infer_sv_hp_map`/
`get_sv_hp_beta` helpers, continuous `sv_minus_wt` instead of binary
`sv_hp==admr_hypo_hp`, Gold/Silver/Bronze tiers) that an in-place rewrite would
just be duplicating — with drift risk — the already-verified scripts in this
repo's `07_phased_pairing_concordance/`. `segdup_stratified_hp_concordance.R`
has no verified current replacement anywhere yet — it's flagged, not fixed
at the analysis level.

## 2. RESOLVED — two coexisting Gold/Silver tiering systems are legitimate, not a bug

**Confirmed:** the two tiering systems are correctly kept separate and no
script silently mixes them.

- `remodeled_constitutional_AMR/pipeline/05_sv_dmr_final_candidates.R` builds
  Gold/Silver/Bronze from `confident_dmr_per_patient.csv.gz` (constitutional-
  inclusive aDMR set) -> `04.final_candidate/gold_tier_final.csv` /
  `silver_tier.csv`. This is a *different, still-current* population from the
  old n=8,616 concordance pairing (issue #1) — pipeline/05's own pairing
  count is unrelated to that number and this script is not superseded.
- `somatic_AMR/04_phaseblock_pairing.R` builds a separate Gold/Silver/Bronze
  from the somatic-only aDMR set (current n=39,212 — Fig 3B, Gold=80/
  Silver=7,096/Bronze=32,036).
- Every `post_processing/` script that reads `04.final_candidate/gold_tier_final.csv`
  /`silver_tier.csv` (`admr_normal_tissue_variance.R`, `cfs_overlap.R`,
  `vignette_single_locus_zoom.R`, `aeq1_or_equivalence.R`,
  `agold1_proximity_decircularize.R`) does so correctly for the
  constitutional-aDMR claim stack — verified by tracing each number back to
  manuscript_v5.3.md:
  - `admr_independent_fragility.R`: Gold+Silver combined SegDup OR = 2.25
    (p=4.7e-6) == manuscript's "SV-proximity-selected Gold+Silver tier aDMR
    subset" number exactly (Results: Constitutional aDMR hotspots).
  - `agold1_proximity_decircularize.R`'s `a_gold1_or_table.csv` (re-read
    directly): All aDMR OR=1.3208, Gold\*_R3 OR=1.8479, Gold\*_R5 OR=4.2213
    [3.8187-4.6664] — these are, digit-for-digit, the manuscript's
    "1.32 / 1.85 / 4.22 [3.82-4.67], p<10⁻¹⁷⁴" recurrence-gradient numbers
    (Fig 4B / repo-internal Fig2 panel D, `C28`). Confirmed by grepping
    `fig2_segdup_coexistence.R` panel D, which reads this exact file.
  - `admr_normal_tissue_variance.R`'s `admr_normal_variance.csv`, Gold-aDMR
    row: fc=1.6108, wilcox_p=1.40e-8, perm_p=0.01 — exact match to
    manuscript's "1.61-fold higher... Wilcoxon p=1.4e-8; permutation p=0.010".
  - `cfs_overlap.R`'s `cfs_overlap.csv` + `cfs_hg38_catalog.csv` (21 CFS loci,
    confirmed by row count): Non-boundary-SV OR=0.948->0.95 ns, Gold-aDMR
    row shows 0/88 loci overlapping any CFS — matches manuscript's "OR=0.95,
    ns; 0/21 CFS loci" (the sentence loosely bundles the SV-side OR and the
    aDMR-side zero-count together, but both numbers check out).
- `nahr_microhomology.R` doesn't read either Gold/Silver system at all (pure
  SV-side SegDup/LAD/PC1 analysis) — it was a false positive in the original
  flagged list, harmless either way.

**Non-obvious finding worth flagging to the user:** `agold1_proximity_decircularize.R`'s
own header/purpose ("A_GOLD1: does Gold OR=3.80 reflect proximity selection
bias?") targets a Gold-tier-only number that does **not** appear anywhere in
manuscript_v5.3.md (grepped for "3.80"/"OR=3.8" — no hits). That original
sensitivity-check framing appears to have been dropped from the manuscript.
But the same script's *other* output rows (the recurrence-only "Gold\*"
subsets, built to answer that sensitivity question) are — via `fig2_segdup_coexistence.R`
panel D — the actual, current primary source of the manuscript's headline
"Constitutional aDMR hotspots scale with recurrence" result (Fig 4B, R>=5
OR=4.22, p<10⁻¹⁷⁴, cited in Abstract/Results/Discussion/Conclusions). In
other words: this script is far more load-bearing for the final manuscript
than its own header comment suggests. `05_fragility_enrichment/README.md`
already flags this ("Fig 4B — key recurrence-gradient script"); worth
double-checking the source script's header comment too before deposition,
since a future reader relying only on the header would come away thinking
this is a minor robustness check rather than the source of a headline
result. Not fixed in the source repo (unlike issue #1, this isn't wrong
output, just an undersold header comment) — flagging for your call on
whether to touch it.

## 3. RESOLVED — survival analysis (14 Gold promoter genes) traced to `viz/v1/fig8_tcga_validation.R`

**Confirmed:** `remodeled_constitutional_AMR/viz/v1/fig8_tcga_validation.R`
Module 4c (lines ~367-439) is the *sole* script anywhere in the source tree
that actually computes the KM/log-rank survival statistics. Checked every
candidate:

- `viz/v1/fig8_tcga_validation.R` — computes `survfit`/`survdiff` per gene,
  `p.adjust(logrank_p, method="BH")`, and writes
  `result/tcga_survival_summary.csv` (`fwrite(surv_summary, ...)`, line 438).
- `viz/v1/fig5_tcga_validation.R` — only reads
  `SURV_SUM <- .../tcga_survival_summary.csv` via `fread(SURV_SUM)` for
  plotting; no `survfit`/`survdiff` call anywhere in the file.
- `viz/v2/fig5_tcga_validation.R` — same: reads the cached
  `tcga_survival_summary.csv` in `make_panelC()`, no recomputation (its only
  `p.adjust` call is for the unrelated SCNA-meth panel D, not survival).
- `viz/v4/fig5_validation_mechanism.R` and everything under
  `external_validation/` — grepped, no `survfit`/`survdiff` at all.

**Data check:** `result/tcga_survival_summary.csv` (mtime 2026-05-18, i.e.
produced by an actual run of `fig8_tcga_validation.R`) has exactly 14 data
rows, all with `logrank_fdr > 0.05` (min FDR = 0.367, SMIM24) — an exact
match to "Survival analysis of 14 Gold-tier aDMR promoter genes in
TCGA-LIHC found no gene reaching FDR < 0.05." The 14 genes come from 15
unique promoter-source genes in `genes_gold_silver_admr.csv` (one dropped
for <20 matched samples, per the `nrow(gb) < 20` filter in Module 4c).

**Action taken:** copied `viz/v1/fig8_tcga_validation.R` into
`10_external_validation_tcga/fig8_tcga_validation.R` in this repo as the
source-of-record for the survival number (see `SCRIPT_MAPPING.md`). Kept
the whole script rather than extracting just Module 4c, since 4c depends on
`tumor_meth`/`clin_surv` objects built earlier in the same script (Modules
1-2: TCGA download/cache, methylation processing) — extracting only the
survival block would silently orphan those dependencies. The script's other
panels (A/B/D/E: T-vs-N boxplots, meth-expr correlation, concordance tiles)
duplicate what `viz/v4` already covers and are not the reason this file is
here; if you want to trim it down to just the survival module for the public
repo, that's a follow-up, not a correctness issue. `viz/v1/fig5_tcga_validation.R`
and `viz/v2/fig5_tcga_validation.R` remain excluded (`EXCLUDED.md`) — both are
pure downstream visualization consumers of `fig8_tcga_validation.R`'s output,
not independent analyses.

## Also worth a quick look

- `05_fragility_enrichment/cpg_adjusted_glm.R` (A1) and
  `segdup_permutation_control.R` (A2) look like earlier iterations of what
  `somatic_AMR/02_fragility_glm.R` and `03_segdup_permutation.R` now do.
  Compare outputs before deciding whether to keep both (e.g. for an
  appendix/history) or drop the older pair.
