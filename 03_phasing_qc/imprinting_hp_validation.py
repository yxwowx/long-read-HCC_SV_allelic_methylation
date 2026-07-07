#!/usr/bin/env python3
"""
imprinting_hp_validation.py — A-II: Biological positive control for phasing fidelity.

At canonical germline imprinted DMRs one parental allele is constitutively
methylated and the other unmethylated.  Correctly phased samples therefore show
large |HP1 − HP2| (≫ background random regions).  12/12 concordance confirms
that the HP1/HP2 labels separate parental alleles across the cohort.

Run: mamba run -n hifiasm python post_processing/imprinting_hp_validation.py
"""

import os
import sys
import random
import numpy as np
import pandas as pd
from multiprocessing import Pool
from pathlib import Path
from scipy import stats

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.paths import HCC_DATA_DIR, PATIENT_CODE_MAP  # noqa: E402

# Paths ==========
HAP_DIR     = str(HCC_DATA_DIR / "DMR_minimap2.out_hg38/cpg_sites")
MAPPING_CSV = str(PATIENT_CODE_MAP)
ICR_BED     = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "data/imprinted_dmrs_hg38.bed")
OUT_DIR     = str(HCC_DATA_DIR / "DMR_SVs/result")
FIG_DIR     = str(HCC_DATA_DIR / "DMR_SVs/figs/v2")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
MIN_COV      = 5      # min per-CpG coverage for both haplotypes
MIN_CPG      = 3      # min CpGs passing filter per locus per haplotype
N_BGND       = 1000   # random background regions per patient
DELTA_THRESH = 0.25   # |Δβ| threshold to count as "large"
N_WORKERS    = 3

# hg38 autosomal chromosome sizes (for background sampling)
CHR_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
    "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
    "chr22": 50818468,
}


def fetch_beta(tabix_fh, chrom, start, end):
    """Return (mean_beta, n_cpg) for CpGs with cov>=MIN_COV; None if <MIN_CPG pass."""
    vals = []
    try:
        for row in tabix_fh.fetch(chrom, start, end):
            cols = row.split("\t")
            if int(cols[5]) >= MIN_COV:
                vals.append(float(cols[3]) / 100.0)
    except (ValueError, KeyError):
        pass
    if len(vals) < MIN_CPG:
        return None, len(vals)
    return float(np.mean(vals)), len(vals)


def load_icr_bed(bed_path):
    loci = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            loci.append((cols[0], int(cols[1]), int(cols[2]),
                         cols[3] if len(cols) > 3 else f"{cols[0]}:{cols[1]}-{cols[2]}",
                         cols[4] if len(cols) > 4 else "unknown"))
    return loci


def random_autosomal_regions(icr_loci, n=N_BGND, seed=42):
    """Draw n random regions matched on width distribution from ICR widths."""
    rng = random.Random(seed)
    widths = [end - start for chrom, start, end, *_ in icr_loci]
    chroms = list(CHR_SIZES.keys())
    regions = []
    for _ in range(n):
        w      = rng.choice(widths)
        chrom  = rng.choice(chroms)
        maxs   = CHR_SIZES[chrom] - w
        start  = rng.randint(0, max(0, maxs))
        regions.append((chrom, start, start + w))
    return regions


def process_patient(args):
    name, patient_code, icr_loci, bg_regions = args

    hap1_path = os.path.join(HAP_DIR, f"{name}_HCC_normal.meth.hap1.bed.gz")
    hap2_path = os.path.join(HAP_DIR, f"{name}_HCC_normal.meth.hap2.bed.gz")

    if not os.path.exists(hap1_path) or not os.path.exists(hap2_path):
        print(f"  SKIP {patient_code}: per-CpG BED files not found", flush=True)
        return None, None

    tbx1 = pysam.TabixFile(hap1_path)
    tbx2 = pysam.TabixFile(hap2_path)

    # ICR loci
    icr_rows = []
    for chrom, start, end, name_locus, parent in icr_loci:
        b1, n1 = fetch_beta(tbx1, chrom, start, end)
        b2, n2 = fetch_beta(tbx2, chrom, start, end)
        if b1 is None or b2 is None:
            continue
        icr_rows.append({
            "patient_code": patient_code,
            "locus":        name_locus,
            "parent":       parent,
            "chrom":        chrom,
            "start":        start,
            "end":          end,
            "beta_hap1":    round(b1, 4),
            "beta_hap2":    round(b2, 4),
            "abs_delta":    round(abs(b1 - b2), 4),
            "n_cpg_hap1":   n1,
            "n_cpg_hap2":   n2,
        })

    # Background
    bg_deltas = []
    for chrom, start, end in bg_regions:
        b1, _ = fetch_beta(tbx1, chrom, start, end)
        b2, _ = fetch_beta(tbx2, chrom, start, end)
        if b1 is not None and b2 is not None:
            bg_deltas.append(abs(b1 - b2))

    tbx1.close()
    tbx2.close()

    if not icr_rows:
        print(f"  {patient_code}: 0 ICR loci passed coverage filter", flush=True)
        return None, None

    icr_df = pd.DataFrame(icr_rows)
    icr_arr = icr_df["abs_delta"].values

    mwu_p = None
    if len(bg_deltas) >= 10:
        mwu = stats.mannwhitneyu(icr_arr, np.array(bg_deltas), alternative="greater")
        mwu_p = float(mwu.pvalue)

    n_large = int((icr_arr >= DELTA_THRESH).sum())
    summary = {
        "patient_code":     patient_code,
        "n_icr_loci":       len(icr_df),
        "median_icr_delta": float(np.median(icr_arr)),
        "n_large_delta":    n_large,
        "pct_large_delta":  round(100 * n_large / len(icr_arr), 1),
        "n_bg_regions":     len(bg_deltas),
        "median_bg_delta":  float(np.median(bg_deltas)) if bg_deltas else None,
        "mwu_p":            mwu_p,
    }
    bg_str = f"{summary['median_bg_delta']:.3f}" if summary['median_bg_delta'] is not None else "NA"
    p_str  = f"{mwu_p:.4f}" if mwu_p is not None else "NA"
    print(f"  {patient_code}: {len(icr_df)} ICRs, "
          f"median |Δ|={summary['median_icr_delta']:.3f} vs bg={bg_str}, p={p_str}",
          flush=True)
    return icr_df, summary


if __name__ == "__main__":
    pm = pd.read_csv(MAPPING_CSV)
    pm["name"] = pm["Samples_ID"].str.replace("_HCC", "", regex=False)

    icr_loci   = load_icr_bed(ICR_BED)
    bg_regions = random_autosomal_regions(icr_loci, n=N_BGND)
    print(f"Loaded {len(icr_loci)} ICR loci; {N_BGND} background regions", flush=True)

    patient_args = [
        (row["name"], row["patient_code"], icr_loci, bg_regions)
        for _, row in pm.iterrows()
    ]
    print(f"Processing {len(patient_args)} patients ({N_WORKERS} workers)...", flush=True)

    with Pool(N_WORKERS) as pool:
        results = pool.map(process_patient, patient_args)

    icr_dfs   = [r[0] for r in results if r[0] is not None]
    summaries = [r[1] for r in results if r[1] is not None]

    if not icr_dfs:
        print("ERROR: no results produced", file=sys.stderr)
        sys.exit(1)

    all_icr = pd.concat(icr_dfs, ignore_index=True)
    all_icr.to_csv(os.path.join(OUT_DIR, "imprinting_hp_validation.csv"), index=False)

    summary_df = pd.DataFrame(summaries).sort_values("patient_code")
    summary_df.to_csv(os.path.join(OUT_DIR, "imprinting_hp_summary.csv"), index=False)

    print("\nCohort summary:", flush=True)
    print(summary_df[["patient_code", "n_icr_loci",
                       "median_icr_delta", "median_bg_delta", "mwu_p"]].to_string(index=False))

    n_concordant = int(
        (summary_df["median_icr_delta"] > 2 * summary_df["median_bg_delta"]).sum()
    )
    print(f"\n{n_concordant}/{len(summary_df)} patients: ICR median |Δ| > 2× background")

    # Figure ==========
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n_pts = len(summaries)
        ncols = 4
        nrows = (n_pts + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.5 * nrows))
        axes = axes.flatten()

        sorted_summaries = sorted(summaries, key=lambda x: x["patient_code"])
        for idx, s in enumerate(sorted_summaries):
            ax   = axes[idx]
            pt   = s["patient_code"]
            data = sorted(
                all_icr[all_icr["patient_code"] == pt]["abs_delta"].values,
                reverse=True
            )
            ax.scatter(range(len(data)), data, s=30, color="#d73027",
                       alpha=0.85, zorder=3, label="ICR |Δβ|")
            bg_med = s["median_bg_delta"] or 0
            ax.axhline(bg_med, color="#2166ac", linestyle="--",
                       linewidth=1.2, label="bg median")
            ax.axhline(DELTA_THRESH, color="#888888", linestyle=":",
                       linewidth=1.0, label=f"|Δ|={DELTA_THRESH}")
            p_str = f"p={s['mwu_p']:.3f}" if s["mwu_p"] is not None else "p=NA"
            ax.set_title(f"{pt}\n{p_str}", fontsize=9)
            ax.set_ylim(0, 1.05)
            ax.set_ylabel("|HP1−HP2|", fontsize=7)
            ax.set_xlabel("ICR rank", fontsize=7)

        for idx in range(len(sorted_summaries), len(axes)):
            axes[idx].set_visible(False)

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", fontsize=9, frameon=True)
        fig.suptitle(
            "A-II: Imprinting-locus haplotype validation (phasing biological positive control)",
            fontsize=11, fontweight="bold"
        )
        plt.tight_layout()
        fig_path = os.path.join(FIG_DIR, "figS_imprinting_validation.png")
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        print(f"\nSaved figure: {fig_path}")
    except Exception as e:
        print(f"Figure generation skipped: {e}")

    print("Done.")
