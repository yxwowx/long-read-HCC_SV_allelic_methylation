#!/usr/bin/env python3
"""
hbv_perread_meth.py — Level A: per-read HBV-junction methylation (C32)

For each somatic HBV integration locus, partitions reads at the genomic
flanking sequence into:
  HBV-carrying : read_id appears in chimeric/INS/clipped breakpoint BED files
  Reference    : all other reads overlapping the same locus

Compares per-read 5mC beta between the two groups at three non-overlapping
distance bins from the junction: 0-1 kb, 1-2 kb, 2-5 kb.

Prediction (cis): HBV-carrying reads are hypomethylated near the integration
site (beta_HBV < beta_ref); effect decays with distance.

Parallelism: multiprocessing.Pool(3) — one worker per patient, each opens
its haplotagged BAM once and processes all loci for that patient.

Usage: mamba run -n hifiasm python post_processing/hbv_perread_meth.py
"""

import gzip
import glob
import os
import sys
from collections import defaultdict
from multiprocessing import Pool
from datetime import date

import pysam
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ── Paths ──────────────────────────────────────────────────────────────────────
HBV_LOCI_CSV  = "/node200data/kachungk/hcc_data/DMR_SVs/12.HBV_analysis/hbv_v1_somatic_hbv_loci.csv"
BED_DIR_CHIM  = "/node200data/kachungk/hcc_data/hg38+HBV/HBV/breakpoint"
BED_DIR_INS   = "/node200data/kachungk/hcc_data/hg38+HBV/HBV/ins/breakpoint"
MAPPING_CSV   = os.path.expanduser("~/patient_code_mapping.csv")
BAM_TMPL      = ("/node200data/kachungk/hcc_data/hg38+HBV/minimap2.out/"
                 "haplotagged_bam/{name}_HCC_tumor.haplotagged.bam")
OUT_DIR       = "/node200data/kachungk/hcc_data/DMR_SVs/result/hbv_perread"
OUT_TSV       = os.path.join(OUT_DIR, "hbv_perread_meth.tsv.gz")
OUT_LOCUS_CSV = "/node200data/kachungk/hcc_data/DMR_SVs/result/hbv_perread_locus_stats.csv"
OUT_POOL_CSV  = "/node200data/kachungk/hcc_data/DMR_SVs/result/hbv_perread_pooled.csv"
FIG_DIR       = "/node200data/kachungk/hcc_data/DMR_SVs/figs/v2"
LOG_FILE      = "/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/logs/claude_decisions.log"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# ── Constants ──────────────────────────────────────────────────────────────────
N_WORKERS      = 3
MIN_CPG        = 3       # minimum CpG calls per read to include
PROB_SCALE     = 255.0
MIN_PER_GROUP  = 3       # minimum reads per group (HBV/ref) for per-locus MWU
BED_COLS       = ["chrom", "start", "end", "hbv_loc", "sampleid", "readname"]

# Non-overlapping distance bins from junction position
# (lo, hi, label): for bin_lo==0, use symmetric window [pos-hi, pos+hi];
#                  otherwise, two flanking sub-windows [pos-hi, pos-lo] and [pos+lo, pos+hi]
DIST_BINS = [
    (0,    1000, "0-1kb"),
    (1000, 2000, "1-2kb"),
    (2000, 5000, "2-5kb"),
]


# ── Per-read 5mC extractor for one genomic window ─────────────────────────────
def _fetch_window(bam, chrom, win_start, win_end, hbv_readnames):
    """
    Fetch reads overlapping [win_start, win_end] and compute per-read mean
    5mC beta using only CpGs within that window.
    Returns dict: read_id -> {read_id, hp, n_cpg, beta_read, is_hbv_read}.
    Reads with <MIN_CPG CpGs in window are skipped.
    """
    win_start = max(0, int(win_start))
    win_end   = int(win_end)
    if win_start >= win_end:
        return {}

    results = {}
    try:
        fetch_iter = bam.fetch(chrom, win_start, win_end)
    except (ValueError, KeyError):
        return {}

    for read in fetch_iter:
        if read.is_unmapped or read.is_supplementary or read.is_secondary:
            continue
        hp = read.get_tag("HP") if read.has_tag("HP") else None
        try:
            mods = read.modified_bases
        except Exception:
            continue
        cpg_key = ("C", 0, "m")
        if cpg_key not in mods:
            continue

        q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)
               if r is not None}
        ref_probs = []
        for q_pos, prob in mods[cpg_key]:
            if q_pos in q2r:
                r_pos = q2r[q_pos]
                if win_start <= r_pos < win_end:
                    ref_probs.append(prob)

        if len(ref_probs) < MIN_CPG:
            continue

        beta = sum(p / PROB_SCALE for p in ref_probs) / len(ref_probs)
        rid  = read.query_name
        entry = {
            "read_id":     rid,
            "hp":          hp,
            "n_cpg":       len(ref_probs),
            "beta_read":   round(beta, 4),
            "is_hbv_read": (rid in hbv_readnames),
        }
        # If same read spans two sub-windows, keep entry with more CpG evidence
        if rid not in results or entry["n_cpg"] > results[rid]["n_cpg"]:
            results[rid] = entry
    return results


def _fetch_bin(bam, chrom, pos, bin_lo, bin_hi, hbv_readnames):
    """
    Fetch and compute beta for one non-overlapping distance bin.
    bin_lo==0: symmetric window [pos-bin_hi, pos+bin_hi].
    Otherwise: two flanking sub-windows merged (dedup by read_id).
    """
    if bin_lo == 0:
        d = _fetch_window(bam, chrom, pos - bin_hi, pos + bin_hi, hbv_readnames)
    else:
        left  = _fetch_window(bam, chrom, pos - bin_hi, pos - bin_lo, hbv_readnames)
        right = _fetch_window(bam, chrom, pos + bin_lo, pos + bin_hi, hbv_readnames)
        d = left.copy()
        for rid, entry in right.items():
            if rid not in d or entry["n_cpg"] > d[rid]["n_cpg"]:
                d[rid] = entry
    return list(d.values())


# ── Per-patient worker (module-level for multiprocessing) ─────────────────────
def process_patient(args):
    """
    Worker: opens one BAM, processes all loci for one patient, returns row list.
    """
    patient_code, name, loci_list, hbv_readnames = args
    bam_path = BAM_TMPL.format(name=name)

    if not os.path.exists(bam_path):
        print(f"  [WARN] BAM not found: {bam_path}", flush=True)
        return []

    print(f"\n[{patient_code}] opening {os.path.basename(bam_path)} "
          f"({len(loci_list)} loci, {len(hbv_readnames)} HBV readnames)", flush=True)

    bam  = pysam.AlignmentFile(bam_path, "rb")
    rows = []

    for chrom, pos, locus_id, hbv_loc in loci_list:
        for bin_lo, bin_hi, bin_label in DIST_BINS:
            reads = _fetch_bin(bam, chrom, pos, bin_lo, bin_hi, hbv_readnames)
            n_hbv = sum(1 for r in reads if r["is_hbv_read"])
            n_ref = sum(1 for r in reads if not r["is_hbv_read"])
            print(f"  [{patient_code}] {locus_id} {bin_label}: "
                  f"n_hbv={n_hbv} n_ref={n_ref}", flush=True)
            for r in reads:
                rows.append({
                    "patient":     patient_code,
                    "locus_id":    locus_id,
                    "hbv_loc":     hbv_loc,
                    "chrom":       chrom,
                    "pos":         pos,
                    "dist_bin":    bin_label,
                    **r,
                })

    bam.close()
    print(f"[{patient_code}] done — {len(rows)} read-rows", flush=True)
    return rows


# ── Stats helpers ──────────────────────────────────────────────────────────────
def _mwu(a, b):
    if len(a) < MIN_PER_GROUP or len(b) < MIN_PER_GROUP:
        return float("nan"), float("nan")
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return float(stat), float(p)


def _wsr(deltas):
    if len(deltas) < 5:
        return float("nan"), float("nan")
    try:
        stat, p = wilcoxon(deltas, alternative="two-sided")
        return float(stat), float(p)
    except Exception:
        return float("nan"), float("nan")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    # Patient code → sample name
    mapping = pd.read_csv(MAPPING_CSV)
    mapping["name"] = mapping["Samples_ID"].str.replace("_HCC", "", regex=False)
    code_to_name = dict(zip(mapping["patient_code"], mapping["name"]))
    name_to_code = {v + "_HCC": k for k, v in code_to_name.items()}

    # Somatic HBV loci — already T_−N_ filtered
    hbv_loci = pd.read_csv(HBV_LOCI_CSV)
    # Keep only tumor-side records (is_tumor may be bool or string)
    if "is_tumor" in hbv_loci.columns:
        hbv_loci = hbv_loci[hbv_loci["is_tumor"].astype(str).str.upper().isin(["TRUE", "1"])]
    hbv_loci = hbv_loci[hbv_loci["pcode"].isin(code_to_name)].copy()
    print(f"Somatic HBV loci: {len(hbv_loci)} rows, "
          f"{hbv_loci['pcode'].nunique()} patients")

    # HBV read names from BED files → per patient frozenset
    hbv_reads_by_patient = defaultdict(set)
    bed_patterns = [
        os.path.join(BED_DIR_CHIM, "T_*_chimeric_breakpoint_hm_hbv.bed"),
        os.path.join(BED_DIR_INS,  "T_*_INS_breakpoint_hm_hbv.bed"),
        os.path.join(BED_DIR_INS,  "T_*_clipped_breakpoint_hm_hbv.bed"),
    ]
    for pattern in bed_patterns:
        for bed_path in glob.glob(pattern):
            try:
                bed = pd.read_csv(bed_path, sep="\t", header=None, names=BED_COLS,
                                  dtype=str)
            except Exception as exc:
                print(f"  [WARN] Could not read {bed_path}: {exc}", flush=True)
                continue
            # sampleid like "T_JJT_HCC" → strip T_/N_ prefix → "JJT_HCC" → lookup pcode
            bed["sid_norm"] = bed["sampleid"].str.replace(r"^[TN]_", "", regex=True)
            for _, row in bed.iterrows():
                pcode = name_to_code.get(row["sid_norm"])
                if pcode:
                    hbv_reads_by_patient[pcode].add(row["readname"])

    print("HBV readname counts per patient:")
    for pc in sorted(hbv_reads_by_patient):
        print(f"  {pc}: {len(hbv_reads_by_patient[pc])} reads")

    # Build per-patient worker args
    patient_args = []
    for pcode, grp in hbv_loci.groupby("pcode"):
        name = code_to_name.get(pcode)
        if name is None:
            continue
        loci_list = [
            (str(row["chrom"]), int(row["pos"]),
             f"{row['chrom']}:{int(row['pos'])}", str(row["hbv_loc"]))
            for _, row in grp.iterrows()
        ]
        readnames = frozenset(hbv_reads_by_patient.get(pcode, set()))
        patient_args.append((pcode, name, loci_list, readnames))

    print(f"\nDispatching {len(patient_args)} patients → Pool(N_WORKERS={N_WORKERS})")

    with Pool(N_WORKERS) as pool:
        results = pool.map(process_patient, patient_args)

    all_rows = [r for patient_rows in results for r in patient_rows]

    if not all_rows:
        print("ERROR: no reads extracted — check BAM paths and loci.", flush=True)
        sys.exit(1)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    print(f"\nSaved per-read table: {OUT_TSV}  ({len(df)} rows)")

    # ── Per-locus MWU stats ────────────────────────────────────────────────────
    locus_rows = []
    for (patient, locus_id, dist_bin), grp in df.groupby(
            ["patient", "locus_id", "dist_bin"]):
        hbv_b = grp.loc[ grp["is_hbv_read"], "beta_read"].to_numpy()
        ref_b = grp.loc[~grp["is_hbv_read"], "beta_read"].to_numpy()
        mwu_stat, mwu_p = _mwu(hbv_b, ref_b)
        med_hbv = float(pd.Series(hbv_b).median()) if len(hbv_b) else float("nan")
        med_ref = float(pd.Series(ref_b).median()) if len(ref_b) else float("nan")
        locus_rows.append({
            "patient":      patient,
            "locus_id":     locus_id,
            "dist_bin":     dist_bin,
            "n_hbv":        len(hbv_b),
            "n_ref":        len(ref_b),
            "median_hbv":   med_hbv,
            "median_ref":   med_ref,
            "median_delta": med_hbv - med_ref,
            "mwu_stat":     mwu_stat,
            "mwu_p":        mwu_p,
        })

    locus_df = pd.DataFrame(locus_rows)
    locus_df.to_csv(OUT_LOCUS_CSV, index=False)
    print(f"Saved locus stats: {OUT_LOCUS_CSV}  ({len(locus_df)} rows)")

    # ── Pooled stats per distance bin ─────────────────────────────────────────
    pooled_rows = []
    for dist_bin, grp in locus_df.groupby("dist_bin"):
        valid   = grp.dropna(subset=["median_delta"])
        deltas  = valid["median_delta"].to_numpy()
        n_loci  = len(deltas)
        n_hypo  = int((deltas < 0).sum())
        wsr_stat, wsr_p = _wsr(deltas)

        # Per-patient sign test
        pat_signs = valid.groupby("patient")["median_delta"].apply(
            lambda x: (x < 0).mean()
        )
        n_pat_hypo = int((pat_signs > 0.5).sum())

        pooled_rows.append({
            "dist_bin":         dist_bin,
            "n_loci":           n_loci,
            "median_delta":     round(float(pd.Series(deltas).median()), 4) if n_loci else float("nan"),
            "n_loci_hypo":      n_hypo,
            "pct_loci_hypo":    round(100 * n_hypo / n_loci, 1) if n_loci else float("nan"),
            "wsr_stat":         wsr_stat,
            "wsr_p":            wsr_p,
            "n_patients":       int(valid["patient"].nunique()),
            "n_pat_majority_hypo": n_pat_hypo,
        })

    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df.to_csv(OUT_POOL_CSV, index=False)
    print(f"\n=== Pooled results ===")
    print(pooled_df.to_string(index=False))

    # ── Figure ─────────────────────────────────────────────────────────────────
    if HAS_MPL and len(df) > 0:
        _make_figure(df, locus_df)

    # ── Log ────────────────────────────────────────────────────────────────────
    n_loci_0_1kb = locus_df[locus_df["dist_bin"] == "0-1kb"].shape[0]
    with open(LOG_FILE, "a") as f:
        f.write(f"[{date.today()}] hbv_perread_meth.py (C32 Level A, Pool={N_WORKERS}): "
                f"{len(df)} per-read rows, {n_loci_0_1kb} loci at 0-1kb; "
                f"see hbv_perread_pooled.csv\n")
    print("\nDone.")


def _make_figure(df, locus_df):
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: violin per dist_bin (HBV vs ref, pooled)
    ax = axes[0]
    bin_order = ["0-1kb", "1-2kb", "2-5kb"]
    positions = []
    tick_pos, tick_lab = [], []
    colors = {"HBV-carrying": "#d73027", "Reference": "#4292c6"}
    for i, dbin in enumerate(bin_order):
        sub = df[df["dist_bin"] == dbin]
        hbv_b = sub.loc[ sub["is_hbv_read"], "beta_read"].to_numpy()
        ref_b = sub.loc[~sub["is_hbv_read"], "beta_read"].to_numpy()
        for j, (data, label, col) in enumerate([
                (hbv_b, "HBV-carrying", "#d73027"),
                (ref_b, "Reference", "#4292c6")]):
            pos = i * 2.5 + j * 0.9
            vp = ax.violinplot([data] if len(data) else [[float("nan")]],
                               positions=[pos], widths=0.7,
                               showmedians=True, showextrema=False)
            for pc in vp["bodies"]:
                pc.set_facecolor(col)
                pc.set_alpha(0.65)
            vp["cmedians"].set_color("black")
            vp["cmedians"].set_linewidth(1.5)
            positions.append((pos, label))
        tick_pos.append(i * 2.5 + 0.45)
        tick_lab.append(dbin)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, fontsize=11)
    ax.set_ylabel("Per-read 5mC β", fontsize=11)
    ax.set_title("HBV-carrying vs reference reads\nby distance bin", fontsize=11, fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=c, label=l, alpha=0.7)
                        for l, c in colors.items()], fontsize=9)
    ax.set_ylim(-0.05, 1.05)

    # Panel B: per-locus median delta at 0-1kb, colored by patient
    ax = axes[1]
    sub0 = locus_df[locus_df["dist_bin"] == "0-1kb"].dropna(subset=["median_delta"])
    if len(sub0) > 0:
        patients = sorted(sub0["patient"].unique())
        pal = plt.cm.tab10.colors
        for k, pt in enumerate(patients):
            pts = sub0[sub0["patient"] == pt]
            ax.scatter(range(len(pts)), pts["median_delta"].to_numpy(),
                       label=pt, color=pal[k % 10], s=60, alpha=0.8, zorder=3)
        ax.axhline(0, color="grey", linewidth=1, linestyle="--")
        ax.set_xlabel("Locus index (0-1 kb)", fontsize=11)
        ax.set_ylabel("Median Δβ (HBV − ref)", fontsize=11)
        ax.set_title("Per-locus methylation difference\n(HBV-carrying minus reference, 0-1 kb)",
                     fontsize=11, fontweight="bold")
        ax.legend(fontsize=8, title="Patient")

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, "fig_hbv_perread.png")
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f"Saved figure: {fig_path}")


if __name__ == "__main__":
    main()
