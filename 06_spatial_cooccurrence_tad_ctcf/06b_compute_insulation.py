#!/usr/bin/env python
"""
Compute genome-wide insulation score from HepG2 Micro-C mcool using cooltools.

Output (TSV):
  chrom  start  end  is_bad_bin  log2_insulation_score_*  boundary_strength_*  ...

Multi-window analysis at resolutions=10kb (closest available: 8000 bp).
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import cooler
import cooltools

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mcool", required=True)
    ap.add_argument("--resolution", type=int, default=8000)
    ap.add_argument("--windows", type=int, nargs="+", default=[80000, 240000, 480000])
    ap.add_argument("--out", required=True, help="Output TSV(.gz) path")
    ap.add_argument("--chromsizes", default=None, help="Optional chromsizes TSV")
    ap.add_argument("--nproc", type=int, default=4)
    args = ap.parse_args()

    cool_uri = f"{args.mcool}::resolutions/{args.resolution}"
    clr = cooler.Cooler(cool_uri)
    print(f"[mcool] {cool_uri}")
    print(f"[mcool] bins: {len(clr.bins())}  resolution: {clr.binsize}  chroms: {clr.chromnames[:5]}...")

    # Restrict to standard autosomes + chrX for speed
    keep_chroms = [c for c in clr.chromnames if c in
                   [f"chr{i}" for i in range(1, 23)] + ["chrX"]]
    view_df = pd.DataFrame({
        "chrom": keep_chroms,
        "start": [0] * len(keep_chroms),
        "end":   [clr.chromsizes[c] for c in keep_chroms],
        "name":  keep_chroms
    })

    print(f"[insulation] windows = {args.windows}  on {len(view_df)} chroms")
    ins = cooltools.insulation(
        clr,
        window_bp=args.windows,
        view_df=view_df,
        clr_weight_name="weight",
        ignore_diags=2,
        min_dist_bad_bin=2,
        nproc=args.nproc,
        verbose=True
    )

    print(f"[output] shape={ins.shape}  cols={list(ins.columns)[:8]}...")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    ins.to_csv(args.out, sep="\t", index=False, compression="gzip" if args.out.endswith(".gz") else None)
    print(f"[done] saved → {args.out}")

if __name__ == "__main__":
    main()
