#!/usr/bin/env python3
"""
Prepare SVclone input files from Severus SV VCF + PURPLE CNV/purity outputs.

Long-read (PacBio HiFi) adaptation:
  (1) DR/DV from Severus VCF used directly as norm/support — no BAM recounting.
  (2) Normal read count adjustment for DNA gains is omitted (skip_norm_adjust=True
      in svclone_config.ini); each long read represents a single DNA molecule, so
      coverage at duplicated regions is not overestimated.

Outputs per patient in <out_dir>/<sample>/:
  {sample}_svinfo.txt     — SVclone count-step equivalent (tab-separated)
  purity_ploidy.txt       — SVclone purity/ploidy file
  {sample}_cnv.tsv        — CNV in Battenberg-like format (chr, startpos, endpos, gtype)
  {sample}_snv.vcf.gz     — Symlink to ClairS PASS SNV VCF for co-clustering

Usage:
  mamba run -n svclone python prep_svclone_longread.py \\
      --sample JJT \\
      --sv_vcf  $HCC_DATA_DIR/severus_minimap2.out_hg38/JJT_HCC.severus.somatic.vcf.gz \\
      --purple_purity $HCC_DATA_DIR/cnv_deepsomatic.out_hg38/purple/JJT_HCC_tumor.purple.purity.tsv \\
      --purple_cnv    $HCC_DATA_DIR/cnv_deepsomatic.out_hg38/purple/JJT_HCC_tumor.purple.cnv.somatic.tsv \\
      --snv_vcf  $HCC_DATA_DIR/clairS_minimap2.out_hg38/JJT_HCC.somatic.merged.PASS.vcf.gz \\
      --out_dir  $HCC_DATA_DIR/DMR_SVs/result/svclone
"""

import argparse
import csv
import gzip
import math
import os
import re
import sys


# Severus VCF → SVclone direction mapping ==========
# Severus STRANDS encodes orientation of both breakend reads.
# SVclone uses +/- to indicate which side of the breakpoint the read extends.
STRANDS_TO_DIRS = {
    "+-": ("+", "-"),   # deletion-like
    "-+": ("-", "+"),   # duplication/tandem-like
    "++": ("+", "+"),   # inversion type 1
    "--": ("-", "-"),   # inversion type 2
}

# Severus SVTYPE → SVclone classification
# Must match svDetectFuncs.getResultType() strings: DEL, DUP, INV, INTDUP, TRX, INS, INTRX
def sv_classification(svtype, chr1, chr2, strands):
    if svtype == "DEL":
        return "DEL"
    if svtype == "DUP":
        return "DUP"
    if svtype == "INV":
        return "INV"
    if svtype == "INS":
        return "INS"
    if svtype == "BND":
        if chr1 != chr2:
            return "INTRX"
        # Same-chromosome BND: infer from strands
        if strands in ("+-", "-+"):
            return "DEL" if strands == "+-" else "DUP"
        return "TRX"
    return "Unknown"


# VCF helpers ==========
def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info_str):
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def parse_format(fmt_keys, fmt_vals):
    keys = fmt_keys.split(":")
    vals = fmt_vals.split(":")
    return dict(zip(keys, vals))


def parse_bnd_partner_pos(alt):
    """Extract chr2, pos2 from BND ALT notation like N[chr1:7745317[ or ]chr2:100[N."""
    m = re.search(r"[\[\]]([^[\]:]+):(\d+)[\[\]]", alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


# Severus VCF → svinfo.txt ==========
def parse_severus_vcf(vcf_path):
    """
    Returns list of SV dicts with keys:
      chr1, pos1, dir1, chr2, pos2, dir2, classification,
      norm (=DR), support (=DV), original_id
    BND pairs are deduplicated: only the first record of each MATE pair is kept.
    """
    svs = []
    seen_mate = set()

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, sv_id, ref, alt, qual, filt = cols[:7]
            info_str = cols[7]
            fmt_keys = cols[8]
            fmt_vals = cols[9]

            if filt not in ("PASS", "."):
                continue

            info = parse_info(info_str)
            fmt  = parse_format(fmt_keys, fmt_vals)

            svtype  = info.get("SVTYPE", "")
            strands = info.get("STRANDS", "")
            pos1    = int(pos)
            chr1    = chrom

            # Read counts from Severus FORMAT
            try:
                dr = int(fmt.get("DR", 0))
                dv = int(fmt.get("DV", 0))
            except (ValueError, TypeError):
                dr, dv = 0, 0

            if dr + dv == 0:
                continue  # skip zero-depth records

            # Handle BND pairs: keep only first of each pair
            mate_id = info.get("MATE_ID", "")
            if svtype == "BND":
                if sv_id in seen_mate:
                    continue
                if mate_id:
                    seen_mate.add(mate_id)
                chr2, pos2 = parse_bnd_partner_pos(alt)
                if chr2 is None:
                    continue
            else:
                chr2 = chr1
                end  = info.get("END")
                pos2 = int(end) if end else pos1 + abs(int(info.get("SVLEN", 1)))

            dir1, dir2 = STRANDS_TO_DIRS.get(strands, ("+", "-"))
            cls = sv_classification(svtype, chr1, chr2, strands)

            svs.append({
                "original_id": sv_id,
                "chr1": chr1, "pos1": pos1, "dir1": dir1,
                "chr2": chr2, "pos2": pos2, "dir2": dir2,
                "classification": cls,
                "norm": dr,
                "support": dv,
            })

    return svs


# Write SVclone svinfo.txt ==========
SVINFO_HEADER = [
    "ID", "chr1", "pos1", "dir1", "chr2", "pos2", "dir2", "classification",
    "split_norm1", "norm_olap_bp1", "span_norm1", "win_norm1",
    "split1", "sc_bases1", "total_reads1",
    "split_norm2", "norm_olap_bp2", "span_norm2", "win_norm2",
    "split2", "sc_bases2", "total_reads2",
    "anomalous", "spanning", "norm1", "norm2", "support",
    "vaf1", "vaf2",
    "original_ID", "original_pos1", "original_pos2",
]

def write_svinfo(svs, out_path):
    """
    Long-read mapping:
      norm1 = norm2 = DR  (reference reads at breakpoint)
      support = spanning = DV  (all variant reads treated as spanning)
      split1 = split2 = 0  (no short-read split distinction)
      span_norm1 = span_norm2 = DR
    """
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SVINFO_HEADER)
        for idx, sv in enumerate(svs):
            norm    = sv["norm"]
            support = sv["support"]
            total   = norm + support
            vaf     = round(support / total, 4) if total > 0 else 0.0

            w.writerow([
                idx,                       # ID
                sv["chr1"], sv["pos1"], sv["dir1"],
                sv["chr2"], sv["pos2"], sv["dir2"],
                sv["classification"],
                0,    # split_norm1
                0,    # norm_olap_bp1
                norm, # span_norm1 — normal reads spanning bp1
                0,    # win_norm1
                0,    # split1
                0,    # sc_bases1
                total,# total_reads1
                0,    # split_norm2
                0,    # norm_olap_bp2
                norm, # span_norm2
                0,    # win_norm2
                0,    # split2
                0,    # sc_bases2
                total,# total_reads2
                0,    # anomalous
                support, # spanning = DV (all variant reads)
                norm,    # norm1
                norm,    # norm2
                support, # support
                vaf, vaf,
                sv["original_id"], sv["pos1"], sv["pos2"],
            ])
    print(f"  Wrote {len(svs)} SVs to {out_path}")


# PURPLE purity/ploidy → SVclone purity_ploidy.txt ==========
def write_purity_ploidy(purple_purity_tsv, sample, out_dir):
    purity = ploidy = None
    with open(purple_purity_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            purity = float(row["purity"])
            ploidy = float(row["ploidy"])
            break
    if purity is None:
        raise ValueError(f"Could not read purity from {purple_purity_tsv}")

    out_path = os.path.join(out_dir, "purity_ploidy.txt")
    with open(out_path, "w") as fh:
        fh.write("purity\tploidy\n")
        fh.write(f"{purity:.4f}\t{ploidy:.4f}\n")
    print(f"  Purity={purity:.3f}, Ploidy={ploidy:.3f} → {out_path}")
    return purity, ploidy


# PURPLE CNV → SVclone Battenberg-like CNV TSV ==========
def write_cnv(purple_cnv_tsv, out_path):
    """
    PURPLE cnv.somatic.tsv → SVclone Battenberg format.

    Columns written: chr, startpos, endpos, nMaj1_A, nMin1_A, frac1_A
    SVclone's load_cnvs() recognises the nMaj1_A header as Battenberg input.
    Copy numbers are rounded to nearest integer (minimum 0).
    Clonal fraction = 1.0 (PURPLE segments are clonal by definition).
    """
    with open(purple_cnv_tsv) as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        # Include subclonal columns (nMaj2_A etc.) as NaN — required by SVclone's
        # Battenberg loader even when all segments are clonal (frac1_A = 1.0).
        writer.writerow([
            "chr", "startpos", "endpos",
            "nMaj1_A", "nMin1_A", "frac1_A",
            "nMaj2_A", "nMin2_A", "frac2_A",
        ])
        n = 0
        for row in reader:
            chrom = row["chromosome"]
            start = int(row["start"])
            end   = int(row["end"])
            major = max(0, round(float(row["majorAlleleCopyNumber"])))
            minor = max(0, round(float(row["minorAlleleCopyNumber"])))
            # frac1_A = 1.0 → clonal; subclonal columns = NaN (not accessed)
            writer.writerow([chrom, start, end, major, minor, 1.0, "NA", "NA", "NA"])
            n += 1
    print(f"  Wrote {n} CNV segments to {out_path}")


# Main ==========
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample",        required=True, help="Sample name (e.g. JJT)")
    ap.add_argument("--sv_vcf",        required=True, help="Severus somatic SV VCF (.vcf.gz)")
    ap.add_argument("--purple_purity", required=True, help="PURPLE purity TSV (*purple.purity.tsv)")
    ap.add_argument("--purple_cnv",    required=True, help="PURPLE somatic CNV TSV (*purple.cnv.somatic.tsv)")
    ap.add_argument("--snv_vcf",       default="",    help="ClairS somatic SNV VCF (.vcf.gz) for co-clustering")
    ap.add_argument("--out_dir",       required=True, help="Output root directory")
    args = ap.parse_args()

    sample_dir = os.path.join(args.out_dir, args.sample)
    os.makedirs(sample_dir, exist_ok=True)

    print(f"[{args.sample}] Parsing Severus VCF ...")
    svs = parse_severus_vcf(args.sv_vcf)
    print(f"  {len(svs)} somatic SVs loaded")

    svinfo_path = os.path.join(sample_dir, f"{args.sample}_svinfo.txt")
    write_svinfo(svs, svinfo_path)

    print(f"[{args.sample}] Writing purity/ploidy ...")
    write_purity_ploidy(args.purple_purity, args.sample, sample_dir)

    print(f"[{args.sample}] Converting PURPLE CNV ...")
    cnv_path = os.path.join(sample_dir, f"{args.sample}_cnv.tsv")
    write_cnv(args.purple_cnv, cnv_path)

    if args.snv_vcf and os.path.exists(args.snv_vcf):
        snv_link = os.path.join(sample_dir, f"{args.sample}_snvs.vcf.gz")
        if os.path.islink(snv_link):
            os.remove(snv_link)
        os.symlink(os.path.abspath(args.snv_vcf), snv_link)
        # Also symlink the index if it exists
        for ext in (".tbi", ".csi"):
            idx = args.snv_vcf + ext
            if os.path.exists(idx):
                link = snv_link + ext
                if os.path.islink(link):
                    os.remove(link)
                os.symlink(os.path.abspath(idx), link)
        print(f"  SNV VCF symlinked → {snv_link}")

    print(f"[{args.sample}] Done. Output in {sample_dir}")


if __name__ == "__main__":
    main()
