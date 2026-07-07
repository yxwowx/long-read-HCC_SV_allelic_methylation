library(data.table)
REPO_ROOT <- local({
  d <- dirname(normalizePath(sub("--file=", "",
    grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
  while (!dir.exists(file.path(d, "shared")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(REPO_ROOT, "shared", "shared_utils.R"))

d <- fread(file.path(Sys.getenv("HCC_DATA_DIR"), "SV_aDMR/somatic_admr_annotated.csv.gz"), nThread=4)
d <- d[abs(diff.Methy) >= 0.2 & nCG >= 5]
tot <- nrow(d)
cat(sprintf('Total somatic aDMR: %d\n', tot))

d[, bulk_mean  := (HP1.Methy + HP2.Methy) / 2]
d[, bulk_delta := abs(bulk_mean - 0.5)]

d[, direction := ifelse(HP1.Methy < 0.5 & HP2.Methy < 0.5, 'both_hypo',
                 ifelse(HP1.Methy > 0.5 & HP2.Methy > 0.5, 'both_hyper', 'bidirectional'))]

cat('\n=== Direction breakdown ===\n')
dir_tab <- d[, .N, by=direction][order(-N)]
dir_tab[, pct := round(N/tot*100, 2)]
print(dir_tab)

cat('\n=== tau sweep (normal.Methy approx = 0.5) ===\n')
for (tau in c(0.10, 0.15, 0.20, 0.25, 0.30)) {
  bulk_n <- sum(d$bulk_delta >= tau)
  excl   <- tot - bulk_n
  cat(sprintf('tau=%.2f  Bulk-detectable=%d (%.1f%%)  AS-only=%d (%.1f%%)\n',
    tau, bulk_n, bulk_n/tot*100, excl, excl/tot*100))
}

cat('\n=== ov_bulk x direction ===\n')
print(d[, .N, by=.(direction, ov_bulk)][order(direction, ov_bulk)])
