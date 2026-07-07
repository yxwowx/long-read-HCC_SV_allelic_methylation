#!/usr/bin/env python3
"""
admr_hp_coverage_symmetry.py

Tests whether allele-specific aDMRs, particularly those at SegDup loci, are
confounded by asymmetric HP1/HP2 haplotype read coverage (a mappability bias
that could generate spurious methylation differences).

For each aDMR × patient: fetch per-CpG coverage from HP1 and HP2 tabix-indexed
pb-CpG-tools beds, compute coverage symmetry metrics (ratio, log2-ratio), then
compare SegDup-overlapping vs non-SegDup aDMRs by Wilcoxon test.
If ratio ≈ 1 and no SegDup-vs-non-SegDup difference, the mappability-confound
hypothesis is not supported.

Run: mamba run -n hifiasm python post_processing/admr_hp_coverage_symmetry.py
"""

import os
import sys
import csv
import gzip
import random
import numpy as np
import pandas as pd
from multiprocessing import Pool
from pathlib import Path
from scipy import stats

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.paths import HCC_DATA_DIR, REFERENCE_DIR, PATIENT_CODE_MAP  # noqa: E402

# Paths ==========
HAP_DIR     = str(HCC_DATA_DIR / "DMR_minimap2.out_hg38/cpg_sites")
MAPPING_CSV = str(PATIENT_CODE_MAP)
GOLD_CSV    = str(HCC_DATA_DIR / "DMR_SVs/04.final_candidate/gold_tier_final.csv")
SILVER_CSV  = str(HCC_DATA_DIR / "DMR_SVs/04.final_candidate/silver_tier.csv")
SEGDUP_BED  = str(REFERENCE_DIR / "LOLACore_180423/hg38/ucsc_features/regions/genomicSuperDups.bed")
FAI         = str(REFERENCE_DIR / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai")
OUT_DIR     = str(HCC_DATA_DIR / "DMR_SVs/result")
FIG_DIR     = str(HCC_DATA_DIR / "DMR_SVs/figs/v2")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Parameters ==========
MIN_COV      = 5    # min per-CpG coverage per haplotype (matches A-II)
MIN_CPG      = 3    # min CpGs passing filter per locus
N_RANDOM     = 500  # matched random autosomal intervals per patient
N_WORKERS    = 3    # parallelism across patients (consistent with A-II)

CHR_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422,"chr11": 135086622,"chr12": 133275309,
    "chr13": 114364328,"chr14": 107043718,"chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468,
}


def fetch_coverage(tabix_fh, chrom, start, end):
    """Return mean per-CpG coverage for CpGs with cov>=MIN_COV, else None."""
    covs = []
    try:
        for row in tabix_fh.fetch(chrom, start, end):
            cols = row.split("\t")
            c = int(cols[5])
            if c >= MIN_COV:
                covs.append(c)
    except (ValueError, KeyError):
        pass
    if len(covs) < MIN_CPG:
        return None, len(covs)
    return float(np.mean(covs)), len(covs)


def load_segdup(bed_path):
    """Load SegDup intervals as a dict of chrom → sorted list of (start, end)."""
    seg = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            seg.setdefault(chrom, []).append((s, e))
    return {c: sorted(v) for c, v in seg.items()}


def overlaps_segdup(chrom, start, end, segdup):
    """Binary search whether interval overlaps any SegDup entry."""
    import bisect
    intervals = segdup.get(chrom, [])
    if not intervals:
        return False
    starts = [iv[0] for iv in intervals]
    idx = bisect.bisect_right(starts, end) - 1
    if idx < 0:
        return False
    for i in range(max(0, idx), min(len(intervals), idx + 5)):
        if intervals[i][0] <= end and intervals[i][1] >= start:
            return True
    return False


def process_patient(args):
    """Worker: compute HP1/HP2 coverage symmetry at aDMR loci for one patient."""
    pt_code, pt_name, admr_loci, segdup, sample_type = args

    hap1_path = os.path.join(HAP_DIR, f"{pt_name}_{sample_type}.meth.hap1.bed.gz")
    hap2_path = os.path.join(HAP_DIR, f"{pt_name}_{sample_type}.meth.hap2.bed.gz")

    if not (os.path.exists(hap1_path) and os.path.exists(hap2_path)):
        print(f"  [SKIP] {pt_code}: beds not found ({sample_type})", flush=True)
        return []

    try:
        tbx1 = pysam.TabixFile(hap1_path)
        tbx2 = pysam.TabixFile(hap2_path)
    except Exception as e:
        print(f"  [SKIP] {pt_code}: tabix open failed — {e}", flush=True)
        return []

    rows = []
    for locus_type, chrom, start, end, locus_id, is_segdup in admr_loci:
        cov1, n1 = fetch_coverage(tbx1, chrom, start, end)
        cov2, n2 = fetch_coverage(tbx2, chrom, start, end)
        if cov1 is None or cov2 is None:
            continue
        ratio = min(cov1, cov2) / max(cov1, cov2) if max(cov1, cov2) > 0 else np.nan
        log2r = np.log2(cov1 / cov2) if cov2 > 0 else np.nan
        rows.append({
            "patient_id": pt_code,
            "locus_type": locus_type,
            "chrom": chrom, "start": start, "end": end,
            "locus_id": locus_id,
            "is_segdup": is_segdup,
            "mean_cov_hp1": cov1, "n_cpg_hp1": n1,
            "mean_cov_hp2": cov2, "n_cpg_hp2": n2,
            "cov_ratio": ratio,     # min/max: 1=perfect symmetry, 0=complete asymmetry
            "log2_hp1_hp2": log2r,  # 0=symmetric; >0 HP1 more covered
        })

    tbx1.close(); tbx2.close()
    print(f"  {pt_code}: {len(rows)} loci processed", flush=True)
    return rows


def main():
    # Load patient mapping ==========
    pmap = pd.read_csv(MAPPING_CSV)
    code2name = dict(zip(pmap["patient_code"], pmap["Samples_ID"]))
    print(f"Patients: {len(code2name)}")

    # Load SegDup ==========
    print("Loading SegDup …")
    segdup = load_segdup(SEGDUP_BED)

    # Load Gold+Silver aDMR loci ==========
    print("Loading aDMR loci …")
    gold   = pd.read_csv(GOLD_CSV)
    silver = pd.read_csv(SILVER_CSV)
    gold["locus_type"]   = "Gold"
    silver["locus_type"] = "Silver"
    admr = pd.concat([gold, silver], ignore_index=True)
    admr = admr[admr["admr_chr"].str.match(r"^chr[0-9XY]+$")]
    admr = admr.drop_duplicates(subset=["admr_chr", "admr_start", "admr_end"])

    # Annotate SegDup
    admr["is_segdup"] = admr.apply(
        lambda r: overlaps_segdup(r["admr_chr"], int(r["admr_start"]),
                                   int(r["admr_end"]), segdup), axis=1
    )
    print(f"aDMR loci: {len(admr)} total; "
          f"SegDup-overlapping: {admr['is_segdup'].sum()}")

    # Build per-patient locus list: (locus_type, chrom, start, end, id, is_segdup)
    # Use union across all patients (per-locus entry once; coverage tested per patient)
    all_loci = [
        (r["locus_type"],
         r["admr_chr"], int(r["admr_start"]), int(r["admr_end"]),
         f"{r['admr_chr']}:{r['admr_start']}-{r['admr_end']}",
         bool(r["is_segdup"]))
        for _, r in admr.iterrows()
    ]
    print(f"Total aDMR intervals to query: {len(all_loci)}")

    # Run per patient ==========
    tasks = [
        (code, code2name[code], all_loci, segdup, "tumor")
        for code in sorted(code2name.keys())
    ]
    print(f"Running {len(tasks)} patients (Pool={N_WORKERS}) …")
    with Pool(N_WORKERS) as pool:
        results = pool.map(process_patient, tasks)

    all_rows = [r for sublist in results for r in sublist]
    print(f"Total records: {len(all_rows)}")

    if not all_rows:
        print("ERROR: no data collected — check bed file paths", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(all_rows)

    # Statistical comparisons ==========
    print("\n=== HP Coverage Symmetry ===")

    # 1. Overall summary
    for grp, gdf in df.groupby("is_segdup"):
        label = "SegDup" if grp else "non-SegDup"
        print(f"  {label}: n={len(gdf)}, "
              f"median cov_ratio={gdf['cov_ratio'].median():.3f}, "
              f"median |log2(HP1/HP2)|={gdf['log2_hp1_hp2'].abs().median():.3f}")

    # 2. Wilcoxon: cov_ratio SegDup vs non-SegDup
    seg_ratio    = df[df["is_segdup"]  == True]["cov_ratio"].dropna()
    nonseg_ratio = df[df["is_segdup"]  == False]["cov_ratio"].dropna()
    if len(seg_ratio) > 5 and len(nonseg_ratio) > 5:
        w, p = stats.mannwhitneyu(seg_ratio, nonseg_ratio, alternative="two-sided")
        print(f"\n  Wilcoxon (MWU) cov_ratio SegDup vs non-SegDup:")
        print(f"    SegDup median={seg_ratio.median():.4f}  "
              f"non-SegDup median={nonseg_ratio.median():.4f}  p={p:.4g}")
        wilcox_ratio_p = p
    else:
        wilcox_ratio_p = np.nan

    # 3. One-sample test: log2(HP1/HP2) ≈ 0 at SegDup aDMRs
    seg_log2 = df[df["is_segdup"] == True]["log2_hp1_hp2"].dropna()
    if len(seg_log2) > 5:
        t, p_one = stats.ttest_1samp(seg_log2, 0)
        print(f"\n  One-sample t-test log2(HP1/HP2)=0 at SegDup aDMRs:")
        print(f"    mean={seg_log2.mean():.4f}  t={t:.3f}  p={p_one:.4g}")
        ttest_p = p_one
    else:
        ttest_p = np.nan

    # Summary CSV ==========
    summary_rows = []
    for is_sd in [True, False]:
        sub = df[df["is_segdup"] == is_sd]
        if sub.empty:
            continue
        summary_rows.append({
            "group": "SegDup" if is_sd else "non-SegDup",
            "n_locus_patient": len(sub),
            "n_loci_unique": sub["locus_id"].nunique(),
            "median_cov_ratio": sub["cov_ratio"].median(),
            "iqr_cov_ratio_lo": sub["cov_ratio"].quantile(0.25),
            "iqr_cov_ratio_hi": sub["cov_ratio"].quantile(0.75),
            "median_abs_log2": sub["log2_hp1_hp2"].abs().median(),
            "pct_ratio_gt0.8": (sub["cov_ratio"] > 0.8).mean() * 100,
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df["wilcox_segdup_vs_nonsegdup_p"] = wilcox_ratio_p
    summary_df["ttest_log2_zero_segdup_p"]     = ttest_p
    print("\nSummary:")
    print(summary_df.to_string(index=False))

    # Figures ==========
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        # Panel A: cov_ratio distribution SegDup vs non-SegDup
        ax = axes[0]
        colors = {True: "#d73027", False: "#4393c3"}
        labels = {True: f"SegDup (n={len(seg_ratio)})",
                  False: f"non-SegDup (n={len(nonseg_ratio)})"}
        for is_sd, c in colors.items():
            sub = df[df["is_segdup"] == is_sd]["cov_ratio"].dropna()
            ax.hist(sub, bins=50, alpha=0.6, color=c, label=labels[is_sd],
                    density=True)
        ax.axvline(1.0, color="k", linestyle="--", lw=1)
        ax.set_xlabel("Coverage ratio (min/max haplotype)")
        ax.set_ylabel("Density")
        ax.set_title("HP1/HP2 coverage symmetry\nat aDMR loci")
        ax.legend(fontsize=8)
        if not np.isnan(wilcox_ratio_p):
            ax.text(0.05, 0.95, f"MWU p={wilcox_ratio_p:.3g}",
                    transform=ax.transAxes, va="top", fontsize=9)

        # Panel B: log2(HP1/HP2) at SegDup aDMRs
        ax2 = axes[1]
        seg_log2_clean = seg_log2.clip(-3, 3)
        ax2.hist(seg_log2_clean, bins=60, color="#d73027", alpha=0.7)
        ax2.axvline(0, color="k", linestyle="--", lw=1)
        ax2.set_xlabel("log₂(HP1 cov / HP2 cov)")
        ax2.set_ylabel("Count")
        ax2.set_title(f"Coverage asymmetry direction\nSegDup aDMRs (n={len(seg_log2)})")
        if not np.isnan(ttest_p):
            ax2.text(0.05, 0.95,
                     f"mean={seg_log2.mean():.3f}\nt-test p={ttest_p:.3g}",
                     transform=ax2.transAxes, va="top", fontsize=9)

        plt.tight_layout()
        fig_path = os.path.join(FIG_DIR, "figS_hp_coverage_symmetry.png")
        plt.savefig(fig_path, dpi=150)
        plt.close()
        print(f"\nFigure saved: {fig_path}")
    except Exception as e:
        print(f"Figure generation failed: {e}")

    # Save CSVs ==========
    per_locus_path = os.path.join(OUT_DIR, "admr_hp_coverage_symmetry.csv")
    df.to_csv(per_locus_path, index=False)

    summary_path = os.path.join(OUT_DIR, "admr_hp_coverage_symmetry_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    print(f"\nSaved: {per_locus_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()
