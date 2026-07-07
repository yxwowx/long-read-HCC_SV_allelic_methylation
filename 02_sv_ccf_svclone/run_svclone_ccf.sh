#!/usr/bin/env bash
# run_svclone_ccf.sh — SVclone CCF estimation for all HCC patients (long-read mode)
#
# Long-read adaptations (PacBio HiFi):
#   (1) DR/DV from Severus VCF used directly; no BAM recounting (annotate/count skipped)
#   (2) Normal read count adjustment for DNA gains omitted (skip_norm_adjust=True in config)
#   (3) Purity/ploidy and CNV from PURPLE; SNVs from ClairS for co-clustering
#
# Steps per patient:
#   1. prep  — Severus VCF + PURPLE → SVclone input files
#   2. filter — SVclone filter (CNV annotation + depth/size filter)
#   3. cluster — CCF estimation with co-clustering of SNVs
#   4. postassign — assign unclustered SVs to nearest cluster
#
# Output: $HCC_DATA_DIR/DMR_SVs/result/svclone/<sample>/
# Final CCF table: $HCC_DATA_DIR/DMR_SVs/result/svclone_ccf_all.csv
#
# Usage: bash run_svclone_ccf.sh [SAMPLE1 SAMPLE2 ...]
#        (no args = run all 12 patients)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
set -a
source "$REPO_ROOT/.env"
set +a

# Paths
BASE=${HCC_DATA_DIR}
SV_DIR=${BASE}/severus_minimap2.out_hg38
PURPLE_DIR=${BASE}/cnv_deepsomatic.out_hg38/purple
CLAIRS_DIR=${BASE}/clairS_minimap2.out_hg38
OUT_DIR=${BASE}/DMR_SVs/result/svclone
CFG="${SCRIPT_DIR}/svclone_config.ini"

mkdir -p "${OUT_DIR}"

# Patient list 
ALL_SAMPLES=(JJT KIS KSJ LHS LSS MSB NSH PJS PSY WSY YJS YMS)
if [[ $# -gt 0 ]]; then
    SAMPLES=("$@")
else
    SAMPLES=("${ALL_SAMPLES[@]}")
fi

# Per-patient function 
run_patient() {
    local SAMPLE=$1
    local SDIR=${OUT_DIR}/${SAMPLE}

    # Input files 
    local SV_VCF=${SV_DIR}/${SAMPLE}_HCC.severus.somatic.vcf.gz
    local PURITY_TSV=${PURPLE_DIR}/${SAMPLE}_HCC_tumor.purple.purity.tsv
    local CNV_TSV=${PURPLE_DIR}/${SAMPLE}_HCC_tumor.purple.cnv.somatic.tsv

    # ClairS: prefer hg38 merged PASS VCF; fallback to deepSomatic
    local SNV_VCF=""
    if [[ -f ${CLAIRS_DIR}/${SAMPLE}_HCC.somatic.merged.PASS.vcf.gz ]]; then
        SNV_VCF=${CLAIRS_DIR}/${SAMPLE}_HCC.somatic.merged.PASS.vcf.gz
    else
        echo "  [WARN] No ClairS SNV VCF for ${SAMPLE}; running without SNV co-clustering"
    fi

    for f in "${SV_VCF}" "${PURITY_TSV}" "${CNV_TSV}"; do
        if [[ ! -f "$f" ]]; then
            echo "  [ERROR] Missing: $f — skipping ${SAMPLE}"
            return 1
        fi
    done

    # Step 1: prepare input files 
    echo "  [1/4] Preparing input files ..."
    mamba run -n svclone python "${SCRIPT_DIR}/prep_svclone_longread.py" \
        --sample    "${SAMPLE}" \
        --sv_vcf    "${SV_VCF}" \
        --purple_purity "${PURITY_TSV}" \
        --purple_cnv    "${CNV_TSV}" \
        ${SNV_VCF:+--snv_vcf "${SNV_VCF}"} \
        --out_dir   "${OUT_DIR}" 

    local SVINFO=${SDIR}/${SAMPLE}_svinfo.txt
    local PP=${SDIR}/purity_ploidy.txt
    local CNV=${SDIR}/${SAMPLE}_cnv.tsv

    # Read params: set long-read appropriate values (rlen=15000, insert=0 → not used)
    cat > "${SDIR}/read_params.txt" << 'EOF'
read_len	insert_mean	insert_std
15000	0	0
EOF

    # Step 2: filter 
    echo "  [2/4] Running svclone filter ..."
    local FILTER_ARGS=(
        -s "${SAMPLE}"
        -i "${SVINFO}"
        -c "${CNV}"
        -p "${PP}"
        -o "${SDIR}"
        -cfg "${CFG}"
        --params "${SDIR}/read_params.txt"
    )
    if [[ -f "${SDIR}/${SAMPLE}_snvs.vcf.gz" ]]; then
        FILTER_ARGS+=(--snvs "${SDIR}/${SAMPLE}_snvs.vcf.gz" --snv_format mutect)
    fi
    mamba run -n svclone svclone filter "${FILTER_ARGS[@]}"

    local FILT_SVS=${SDIR}/${SAMPLE}_filtered_svs.tsv
    if [[ ! -f "${FILT_SVS}" ]]; then
        echo "  [ERROR] Filter step produced no output for ${SAMPLE}" 
        return 1
    fi

    # Step 3: cluster 
    echo "  [3/4] Running svclone cluster ..."
    local CLUS_ARGS=(
        -s "${SAMPLE}"
        -o "${SDIR}"
        -i "${FILT_SVS}"
        -p "${PP}"
        -cfg "${CFG}"
        --params "${SDIR}/read_params.txt"
    )
    local FILT_SNVS=${SDIR}/${SAMPLE}_filtered_snvs.tsv
    if [[ -f "${FILT_SNVS}" ]]; then
        CLUS_ARGS+=(--snvs "${FILT_SNVS}")
    fi
    mamba run -n svclone svclone cluster "${CLUS_ARGS[@]}" 

    # Step 4: postassign 
    echo "  [4/4] Running svclone postassign ..." 
    local RDATA_SVS=${SDIR}/ccube_out/${SAMPLE}_ccube_sv_results.RData
    local RDATA_SNVS=${SDIR}/ccube_out/snvs/${SAMPLE}_ccube_snv_results.RData

    local PA_ARGS=(-s "${SAMPLE}" -o "${SDIR}" --svs "${RDATA_SVS}")
    if [[ -f "${RDATA_SNVS}" ]]; then
        PA_ARGS+=(-j --snvs "${RDATA_SNVS}")
    fi
    mamba run -n svclone svclone postassign "${PA_ARGS[@]}"
        echo "  [WARN] postassign step skipped or failed (non-fatal)" 

    echo "  [OK] ${SAMPLE} complete."
}

# Run all patients 
echo "Running SVclone CCF estimation for: ${SAMPLES[*]}"
echo "Output directory: ${OUT_DIR}"
echo ""

FAILED=()
for SAMPLE in "${SAMPLES[@]}"; do
    run_patient "${SAMPLE}" || FAILED+=("${SAMPLE}")
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "[WARN] Failed samples: ${FAILED[*]}"
fi

#  Collect results across patients 
echo ""
echo "Collecting CCF results ..."
mamba run -n renv Rscript - << 'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR <- file.path(Sys.getenv("HCC_DATA_DIR"), "DMR_SVs/result/svclone")
SAMPLES <- c("JJT","KIS","KSJ","LHS","LSS","MSB","NSH","PJS","PSY","WSY","YJS","YMS")

# Read per-SV cluster certainty (postassign output: chr1/pos1/../most_likely_assignment/average_proportion1/2)
# and join with subclonal structure to get CCF (proportion) per cluster.
collect_ccf <- function(samp) {
  cert_f <- file.path(OUT_DIR, samp, paste0(samp, "_cluster_certainty.txt"))
  sub_f  <- file.path(OUT_DIR, samp, paste0(samp, "_subclonal_structure.txt"))
  if (!file.exists(cert_f) || !file.exists(sub_f)) return(NULL)

  cert <- fread(cert_f)
  sub  <- fread(sub_f)
  # sub columns: cluster, n_ssms, proportion (CCF)
  setnames(sub, "cluster", "most_likely_assignment")
  cert <- merge(cert, sub[, .(most_likely_assignment, ccf = proportion)],
                by = "most_likely_assignment", all.x = TRUE)
  cert[, sample := samp]
  cert
}

all_ccf <- rbindlist(lapply(SAMPLES, collect_ccf), fill = TRUE)
if (nrow(all_ccf) > 0) {
  fwrite(all_ccf, file.path(OUT_DIR, "svclone_ccf_raw_all.csv"))
  cat(sprintf("Collected %d SV CCF estimates across %d patients\n",
              nrow(all_ccf), length(unique(all_ccf$sample))))
} else {
  cat("No CCF results found.\n")
}
RSCRIPT

echo ""
echo "Done. Run mutationtimer_sv_timing.R next for timing categories."
echo "[$(date +%F)] run_svclone_ccf.sh complete: ${#SAMPLES[@]} patients, failed=${#FAILED[@]}"
