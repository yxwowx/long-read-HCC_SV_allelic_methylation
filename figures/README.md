# figures

Figure-generation scripts, current version only. Manuscript figure numbering
(Fig 1-5 main, Fig S1-S8 supplementary) does not map 1:1 onto the source
scripts' internal version numbering — expect to re-map filenames to final
Fig/FigS numbers when copying these in.

## main/ (Fig 1-5)
- `remodeled_constitutional_AMR/viz/v4/*.R` + `README.md` + `shared_theme.R` (Fig 1a/1b, Fig 2-5)
- `somatic_AMR/viz/v1/*.R` + `figure_captions.md` + `shared_theme.R` (this is the only/current viz set for somatic_AMR — there is no v2+)

## supplementary/ (Fig S1-S8)
Same two source directories as above cover the supplementary panels
(figS1-figS10, figS_hbv, figS_phasing, figS_gviz_vignettes in v4; figS1-figS9
+ figure_captions.md in somatic_AMR/viz/v1). Sort into main/ vs supplementary/
by cross-checking against the final Figure Legends section of the manuscript,
not by the scripts' own figS-numbering.

Also note: `vignette_single_locus_zoom.R` (post_processing, six highest-recurrence
Gold aDMR hotspots, Fig S7) is analysis+figure combined — decide whether it lives
here or in `05_fragility_enrichment/` alongside its sibling analysis scripts.
