#!/usr/bin/env bash
# Download ENCODE2 HepG2 Repli-seq phase bigWigs (hg19) to reference dir.
# Experiment: ENCSR000CXG (Repli-seq, HepG2, hg19, 6 phases)
# Run: bash post_processing/download_encode_repliseq.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/.env"

OUT_DIR="${HCC_DATA_DIR}/DMR_SVs/external_validation_cache/repliseq"
mkdir -p "${OUT_DIR}"

declare -A FILES=(
  ["ENCFF001GPK_S1_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPK/@@download/ENCFF001GPK.bigWig"
  ["ENCFF001GPO_S2_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPO/@@download/ENCFF001GPO.bigWig"
  ["ENCFF001GPU_S3_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPU/@@download/ENCFF001GPU.bigWig"
  ["ENCFF001GPX_S4_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPX/@@download/ENCFF001GPX.bigWig"
  ["ENCFF001GPC_G1b_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPC/@@download/ENCFF001GPC.bigWig"
  ["ENCFF001GPF_G2_hg19.bigWig"]="https://www.encodeproject.org/files/ENCFF001GPF/@@download/ENCFF001GPF.bigWig"
)

for fname in "${!FILES[@]}"; do
  dest="${OUT_DIR}/${fname}"
  if [[ -f "${dest}" ]]; then
    echo "Already exists: ${fname}"
    continue
  fi
  echo "Downloading ${fname}..."
  wget -q --show-progress -O "${dest}" "${FILES[$fname]}"
  echo "  → saved: ${dest}"
done

echo "Done. Files in ${OUT_DIR}:"
ls -lh "${OUT_DIR}"
