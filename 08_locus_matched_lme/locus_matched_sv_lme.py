#!/usr/bin/env python3
"""
Analysis #1 (C8-extended): Genome-wide locus-matched cross-patient LME
=======================================================================
For every aDMR locus present in ≥MIN_PATIENTS patients, test whether
patients with a somatic SV within ±WINDOW_KB show significantly higher
|HP Δβ| using a linear mixed-effects model (locus as random intercept).

Hypothesis (cis-induction): β_sv_present > 0, p < 0.05
Null (shared-fragility):    β_sv_present ≈ 0

Data:
  confident_dmr_per_patient.csv.gz  → all phased aDMR loci + HP betas
  sv_tad_ctcf_annotation.v2.csv.gz  → all SV breakpoints per patient

Run:
  mamba run -n renv python post_processing/locus_matched_sv_lme.py
  python post_processing/locus_matched_sv_lme.py  (needs pandas, statsmodels, scipy, matplotlib)
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime

warnings.filterwarnings("ignore", category=FutureWarning)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from shared.paths import HCC_DATA_DIR  # noqa: E402

# Paths ==========
ROOT     = HCC_DATA_DIR / "DMR_SVs"
DMR_FILE = ROOT / "01.DMR_recurrence/confident_dmr_per_patient.csv.gz"
SV_FILE  = ROOT / "02.sv_dmr_enrichment/sv_tad_ctcf_annotation.v2.csv.gz"
OUTDIR   = ROOT / "result"
OUTDIR.mkdir(exist_ok=True)

# Parameters ==========
BIN_KB       = 50          # genomic bin size for locus definition (matches enrichment window)
MIN_PATIENTS = 2           # minimum patients with aDMR at locus for LME
WINDOWS_KB   = [10, 50, 100]   # SV–aDMR proximity windows for sensitivity analysis
SV_TYPES     = ["DEL", "DUP", "INS", "BND"]

def log(msg):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


# 1. Load aDMR data ==========
log("Loading aDMR data (confident_dmr_per_patient.csv.gz)...")
dmr = pd.read_csv(
    DMR_FILE,
    usecols=["seqnames", "start", "end", "HP1.Methy", "HP2.Methy",
             "patient_code", "n_patients"],
    dtype={"seqnames": str, "patient_code": str},
)
dmr.rename(columns={"seqnames": "chr", "n_patients": "locus_recurrence"}, inplace=True)
dmr["hp_abs_diff"] = (dmr["HP1.Methy"] - dmr["HP2.Methy"]).abs()
dmr = dmr.dropna(subset=["hp_abs_diff", "chr", "start", "end", "patient_code"])
dmr["center"] = (dmr["start"] + dmr["end"]) // 2
log(f"  {len(dmr):,} phased aDMR records; {dmr['patient_code'].nunique()} patients")

# Assign 50kb bin locus_id
dmr["locus_id"] = (
    dmr["chr"] + ":"
    + ((dmr["center"] // (BIN_KB * 1000)) * BIN_KB).astype(str) + "kb"
)

# Aggregate to one row per (patient, locus): mean hp_abs_diff
# (a patient can have multiple DMR segments mapping to the same 50kb bin)
dmr_locus = (
    dmr.groupby(["locus_id", "patient_code"], observed=True)
    .agg(
        hp_abs_diff=("hp_abs_diff", "mean"),
        locus_recurrence=("locus_recurrence", "max"),
        chr=("chr", "first"),
        center=("center", "median"),
    )
    .reset_index()
)
log(f"  {len(dmr_locus):,} (locus × patient) records; "
    f"{dmr_locus['locus_id'].nunique():,} unique loci")


# 2. Load SV breakpoints ==========
log("Loading SV breakpoints (sv_tad_ctcf_annotation.v2.csv.gz)...")
sv = pd.read_csv(
    SV_FILE,
    usecols=["seqnames", "start", "end", "geom_type", "sample",
             "svtype", "cnv_class"],
    dtype={"seqnames": str, "sample": str},
)
# 'sample' column == patient_code (P1-P12)
sv.rename(columns={
    "seqnames": "sv_chr",
    "start":    "sv_start",
    "end":      "sv_end",
    "sample":   "patient_code",
    "geom_type": "sv_type",
}, inplace=True)
sv["sv_mid"] = (sv["sv_start"] + sv["sv_end"]) // 2
sv = sv.dropna(subset=["patient_code", "sv_chr", "sv_start"])
sv["sv_chr"] = sv["sv_chr"].astype(str)
log(f"  {len(sv):,} SV breakpoints; {sv['patient_code'].nunique()} patients")


# 3. Flag SV presence per (locus, patient, window) ==========
def flag_sv_presence(dmr_locus_df, sv_df, window_kb):
    """
    For each (locus, patient) row, flag True if a same-patient SV breakpoint
    falls within ±window_kb of the locus center.
    Uses sorted searchsorted for O(n log n) performance.
    """
    window_bp = window_kb * 1000
    result = np.zeros(len(dmr_locus_df), dtype=bool)
    idx_arr = dmr_locus_df.index.to_numpy()

    sv_by_pat = {pat: grp for pat, grp in sv_df.groupby("patient_code")}

    for pat, d_pat in dmr_locus_df.groupby("patient_code"):
        if pat not in sv_by_pat:
            continue
        s_pat = sv_by_pat[pat]

        for chrom, d_chr in d_pat.groupby("chr"):
            s_chr = np.sort(
                s_pat.loc[s_pat["sv_chr"] == chrom, "sv_mid"].dropna().values
            )
            if len(s_chr) == 0:
                continue

            centers = d_chr["center"].values
            local_idx = np.where(np.isin(idx_arr, d_chr.index))[0]

            # searchsorted: for each center, closest SV on left and right
            ir = np.searchsorted(s_chr, centers, side="right")
            il = ir - 1

            right_dist = np.where(
                ir < len(s_chr), s_chr[np.clip(ir, 0, len(s_chr) - 1)] - centers, np.inf
            )
            left_dist = np.where(
                il >= 0, centers - s_chr[np.clip(il, 0, len(s_chr) - 1)], np.inf
            )
            within = (np.minimum(right_dist, left_dist) <= window_bp)
            result[local_idx] = within

    return pd.Series(result, index=dmr_locus_df.index, name=f"sv_{window_kb}kb")


log(f"Flagging SV presence at windows: {WINDOWS_KB} kb...")
for wkb in WINDOWS_KB:
    col = f"sv_{wkb}kb"
    dmr_locus[col] = flag_sv_presence(dmr_locus, sv, wkb).values
    n_pos = dmr_locus[col].sum()
    log(f"  ±{wkb:>3}kb: {n_pos:,} / {len(dmr_locus):,} records with SV "
        f"({100*n_pos/len(dmr_locus):.1f}%)")


# 4. LME - main analysis (50 kb window) ==========
log("Fitting linear mixed-effects models...")

try:
    import statsmodels.formula.api as smf
    HAS_STATSMODELS = True
except ImportError:
    log("  WARNING: statsmodels not found — skipping LME, outputting summary stats only")
    HAS_STATSMODELS = False


def run_lme(df, sv_col, label, min_pts=MIN_PATIENTS):
    """
    Fit hp_abs_diff ~ sv_present + (1 | locus_id).
    Returns dict of results.
    """
    # Filter to loci with ≥ min_pts patients
    locus_counts = df.groupby("locus_id")["patient_code"].nunique()
    keep = locus_counts[locus_counts >= min_pts].index
    sub = df[df["locus_id"].isin(keep)].copy()
    sub[sv_col] = sub[sv_col].astype(int)  # 0/1 for formula

    if len(sub) < 20:
        return {"label": label, "n_loci": len(keep), "n_obs": len(sub),
                "beta": np.nan, "se": np.nan, "z": np.nan, "p": np.nan,
                "ci_lo": np.nan, "ci_hi": np.nan, "converged": False}

    # Summary stats (always computed)
    sv_pos = sub[sub[sv_col] == 1]["hp_abs_diff"]
    sv_neg = sub[sub[sv_col] == 0]["hp_abs_diff"]
    mwu_stat, mwu_p = stats.mannwhitneyu(sv_pos, sv_neg, alternative="greater")
    log(f"  [{label}] n_loci={len(keep)}, n_obs={len(sub)}, "
        f"median SV+={sv_pos.median():.4f}, SV−={sv_neg.median():.4f}, "
        f"MWU p={mwu_p:.4f}")

    if not HAS_STATSMODELS:
        return {"label": label, "n_loci": len(keep), "n_obs": len(sub),
                "median_sv_pos": sv_pos.median(), "median_sv_neg": sv_neg.median(),
                "mwu_p": mwu_p, "converged": False}

    formula = f"hp_abs_diff ~ {sv_col}"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mdl = smf.mixedlm(formula, data=sub, groups=sub["locus_id"])
            fit = mdl.fit(reml=True, method="lbfgs")
        coef_key = sv_col
        beta  = fit.fe_params[coef_key]
        se    = fit.bse_fe[coef_key]
        z_val = fit.tvalues[coef_key]
        pval  = fit.pvalues[coef_key]
        ci    = fit.conf_int().loc[coef_key]
        log(f"    LME β={beta:.5f} (SE={se:.5f}), z={z_val:.3f}, p={pval:.4f}, "
            f"95%CI=[{ci[0]:.5f},{ci[1]:.5f}]")
        return {"label": label, "n_loci": len(keep), "n_obs": len(sub),
                "median_sv_pos": sv_pos.median(), "median_sv_neg": sv_neg.median(),
                "mwu_p": mwu_p,
                "beta": beta, "se": se, "z": z_val, "p": pval,
                "ci_lo": ci[0], "ci_hi": ci[1], "converged": True}
    except Exception as e:
        log(f"    LME failed ({e}); reporting MWU only")
        return {"label": label, "n_loci": len(keep), "n_obs": len(sub),
                "median_sv_pos": sv_pos.median(), "median_sv_neg": sv_neg.median(),
                "mwu_p": mwu_p, "beta": np.nan, "converged": False}


# 4a. Main model (50kb, all aDMRs) ==========
results = []

log("  --- Main model (±50 kb, all phased aDMRs, min_patients≥2) ---")
results.append(run_lme(dmr_locus, "sv_50kb", "All aDMRs | ±50kb"))

# 4b. Sensitivity: window size ==========
log("  --- Sensitivity: window size ---")
for wkb in WINDOWS_KB:
    col = f"sv_{wkb}kb"
    results.append(run_lme(dmr_locus, col, f"All aDMRs | ±{wkb}kb"))

# 4c. Sensitivity: recurrent loci only (n_patients >= 3) ==========
log("  --- Sensitivity: recurrent loci (locus_recurrence ≥ 3) ---")
recurrent = dmr_locus[dmr_locus["locus_recurrence"] >= 3].copy()
results.append(run_lme(recurrent, "sv_50kb", "Recurrent loci (n≥3) | ±50kb",
                       min_pts=2))

# 4d. Sensitivity: by SV type ==========
log("  --- Sensitivity: by SV type ---")
# For each SV type, flag sv_present only for that type
sv_by_type = {t: sv[sv["sv_type"] == t] for t in SV_TYPES}
for sv_t, sv_sub in sv_by_type.items():
    if len(sv_sub) == 0:
        continue
    col_t = f"sv_{sv_t}_50kb"
    dmr_locus[col_t] = flag_sv_presence(dmr_locus, sv_sub, 50).values
    n = dmr_locus[col_t].sum()
    if n < 10:
        log(f"  [{sv_t}] n={n} — too few, skipping")
        continue
    results.append(run_lme(dmr_locus, col_t, f"{sv_t} SVs | ±50kb"))


# 5. Per-locus summary table ==========
log("Building per-locus summary...")
per_locus = (
    dmr_locus.groupby("locus_id")
    .apply(lambda g: pd.Series({
        "n_patients_admr": g["patient_code"].nunique(),
        "median_hp_abs_diff_all": g["hp_abs_diff"].median(),
        "median_hp_abs_diff_sv_pos": g.loc[g["sv_50kb"], "hp_abs_diff"].median()
            if g["sv_50kb"].any() else np.nan,
        "median_hp_abs_diff_sv_neg": g.loc[~g["sv_50kb"], "hp_abs_diff"].median()
            if (~g["sv_50kb"]).any() else np.nan,
        "n_sv_pos": g["sv_50kb"].sum(),
        "n_sv_neg": (~g["sv_50kb"]).sum(),
        "locus_recurrence": g["locus_recurrence"].max(),
        "chr": g["chr"].iloc[0],
        "center": int(g["center"].median()),
    }), include_groups=False)
    .reset_index()
)
per_locus["delta_median"] = (
    per_locus["median_hp_abs_diff_sv_pos"] - per_locus["median_hp_abs_diff_sv_neg"]
)

per_locus_out = OUTDIR / "locus_matched_sv_lme_per_locus.csv"
per_locus.to_csv(per_locus_out, index=False)
log(f"  Per-locus table → {per_locus_out}")


# 6. Results table ==========
results_df = pd.DataFrame(results)
res_out = OUTDIR / "locus_matched_sv_lme_results.csv"
results_df.to_csv(res_out, index=False)
log(f"  LME results → {res_out}")
print("\n" + results_df.to_string(index=False) + "\n")


# 7. Figures ==========
log("Generating figures...")
FIG_DPI = 180

# Fig A: Violin — SV+ vs SV- |HP Δβ| (main 50kb window, all loci) ==========
plot_df = dmr_locus[
    dmr_locus["locus_id"].isin(
        dmr_locus.groupby("locus_id")["patient_code"]
        .nunique()[lambda x: x >= MIN_PATIENTS].index
    )
].copy()
plot_df["group"] = np.where(plot_df["sv_50kb"], "SV+ (≤50 kb)", "SV− (>50 kb)")

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Violin
ax = axes[0]
groups = ["SV+ (≤50 kb)", "SV− (>50 kb)"]
colors = ["#E24B4A", "#5B9BD5"]
data_parts = [plot_df[plot_df["group"] == g]["hp_abs_diff"].dropna().values for g in groups]
parts = ax.violinplot(data_parts, positions=[1, 2], showmedians=True, widths=0.6)
for pc, col in zip(parts["bodies"], colors):
    pc.set_facecolor(col); pc.set_alpha(0.7)
for med in [parts["cmedians"]]:
    med.set_color("black"); med.set_linewidth(2)
ax.set_xticks([1, 2]); ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel("|HP Δβ|", fontsize=11)
ax.set_title("A.  |HP Δβ| by SV proximity (±50 kb)", fontweight="bold", fontsize=11)
# Annotate MWU p
main_res = next((r for r in results if r.get("label") == "All aDMRs | ±50kb"), None)
if main_res and "mwu_p" in main_res:
    plab = (f"MWU p={main_res['mwu_p']:.3f}"
            if main_res["mwu_p"] >= 0.001 else f"MWU p={main_res['mwu_p']:.2e}")
    ax.text(0.98, 0.97, plab, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color="black")
    if HAS_STATSMODELS and not np.isnan(main_res.get("beta", np.nan)):
        blab = (f"LME β={main_res['beta']:.4f}\n"
                f"p={main_res['p']:.3f}" if main_res["p"] >= 0.001
                else f"LME β={main_res['beta']:.4f}\np={main_res['p']:.2e}")
        ax.text(0.98, 0.84, blab, transform=ax.transAxes, ha="right", va="top",
                fontsize=9, color="#555")
for n_pts, x_pos in zip(
    [len(d) for d in data_parts], [1, 2]
):
    ax.text(x_pos, ax.get_ylim()[0] + 0.01, f"n={n_pts}",
            ha="center", va="bottom", fontsize=8, color="grey")
ax.grid(axis="y", alpha=0.3, linestyle="--")

# Fig B: Forest plot — sensitivity analyses ==========
ax2 = axes[1]
if HAS_STATSMODELS:
    forest_df = results_df[results_df["converged"] == True].copy()
    if len(forest_df) > 0:
        y_pos = range(len(forest_df) - 1, -1, -1)
        cols = ["#E24B4A" if b > 0 else "#5B9BD5" for b in forest_df["beta"]]
        ax2.axvline(0, color="black", linewidth=0.8, linestyle="--")
        for i, (_, row) in enumerate(forest_df.iterrows()):
            y = len(forest_df) - 1 - i
            ax2.plot([row["ci_lo"], row["ci_hi"]], [y, y],
                     color=cols[i], linewidth=1.8, solid_capstyle="round")
            ax2.plot(row["beta"], y, "o", color=cols[i], markersize=6, zorder=5)
            pstr = f"p={row['p']:.3f}" if row["p"] >= 0.001 else f"p={row['p']:.2e}"
            ax2.text(
                ax2.get_xlim()[1] if ax2.get_xlim()[1] > 0 else 0.02,
                y, f"  {pstr}", va="center", fontsize=7.5
            )
        ax2.set_yticks(list(y_pos))
        ax2.set_yticklabels(forest_df["label"].tolist()[::-1], fontsize=8.5)
        ax2.set_xlabel("LME β (SV-present effect on |HP Δβ|)", fontsize=10)
        ax2.set_title("B.  Forest plot — sensitivity analyses", fontweight="bold", fontsize=11)
        ax2.grid(axis="x", alpha=0.3, linestyle="--")
    else:
        ax2.text(0.5, 0.5, "No converged LME models", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=10, color="grey")
        ax2.set_title("B.  Forest plot (no converged models)", fontsize=11)
else:
    # MWU-based forest using median differences
    mwu_df = results_df[~results_df["mwu_p"].isna()].copy()
    y_vals = range(len(mwu_df))
    ax2.barh(list(y_vals),
             mwu_df["median_sv_pos"].values - mwu_df["median_sv_neg"].values,
             color=["#E24B4A" if v > 0 else "#5B9BD5"
                    for v in mwu_df["median_sv_pos"] - mwu_df["median_sv_neg"]],
             alpha=0.75, height=0.6)
    ax2.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax2.set_yticks(list(y_vals))
    ax2.set_yticklabels(mwu_df["label"].tolist(), fontsize=8.5)
    ax2.set_xlabel("Δ median |HP Δβ| (SV+ − SV−)", fontsize=10)
    ax2.set_title("B.  Sensitivity: Δ median |HP Δβ|", fontweight="bold", fontsize=11)
    ax2.grid(axis="x", alpha=0.3, linestyle="--")

fig.suptitle(
    "Genome-wide locus-matched analysis: does SV proximity predict aDMR magnitude?\n"
    "Model: |HP Δβ| ~ sv_present + (1 | locus_id);  loci with ≥2 patients",
    fontsize=10, y=1.01
)
plt.tight_layout()
fig_out = OUTDIR / "locus_matched_sv_lme_plot.png"
fig.savefig(fig_out, dpi=FIG_DPI, bbox_inches="tight")
plt.close()
log(f"  Figure → {fig_out}")

# Fig C: per-locus delta scatter ==========
pl_plot = per_locus.dropna(subset=["delta_median"]).copy()
if len(pl_plot) > 0:
    fig2, ax3 = plt.subplots(figsize=(6, 5))
    ax3.axhline(0, color="grey", linewidth=0.8, linestyle="--")
    sc = ax3.scatter(
        pl_plot["locus_recurrence"],
        pl_plot["delta_median"],
        c=pl_plot["delta_median"],
        cmap="RdBu_r", vmin=-0.15, vmax=0.15,
        s=30 + 10 * pl_plot["n_sv_pos"], alpha=0.6, edgecolors="none"
    )
    plt.colorbar(sc, ax=ax3, label="Δ median |HP Δβ| (SV+ − SV−)")
    ax3.set_xlabel("Locus recurrence (n patients with aDMR)", fontsize=10)
    ax3.set_ylabel("Δ median |HP Δβ| (SV+ − SV−)", fontsize=10)
    ax3.set_title(
        "C.  Per-locus effect of SV proximity (±50 kb)\n"
        "Positive = SV-bearing patients show larger allelic imbalance",
        fontsize=10
    )
    ax3.grid(alpha=0.3, linestyle="--")
    fig2.tight_layout()
    fig2_out = OUTDIR / "locus_matched_sv_lme_per_locus.png"
    fig2.savefig(fig2_out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()
    log(f"  Per-locus scatter → {fig2_out}")


log("Done.")
print(f"\nOutputs in {OUTDIR}:")
for p in [res_out, per_locus_out, fig_out]:
    print(f"  {p}")
