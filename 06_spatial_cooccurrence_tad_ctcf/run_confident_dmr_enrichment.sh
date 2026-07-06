#!/bin/bash
# Two-phase enrichment run using confident DMRs
#
# Phase 1 (default): nperm=100, --no_plot — window selection scan
# Phase 2:           nperm=1000, with plots — run after reviewing Phase 1 output
#
# Usage:
#   Phase 1:  bash run_confident_dmr_enrichment.sh           [runs nperm=100]
#   Phase 2:  bash run_confident_dmr_enrichment.sh 1000 50000,100000
#             (arg1 = n_perm, arg2 = comma-separated primary windows to run)

set -euo pipefail

N_PERM="${1:-100}"
# Phase 2: pass specific primary windows; Phase 1: scan all candidates
PHASE2_PRIMARIES="${2:-}"

SCRIPT="/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/pipeline/03_sv_dmr_enrichment.R"
DMR_FILE="/node200data/kachungk/hcc_data/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz"
SV_FILE="/node200data/kachungk/hcc_data/DMR_SVs/sv_tad_ctcf_annotation.csv.gz"
WINDOWS="10000,25000,50000,100000,250000,500000,1000000"
BASE_OUTDIR="/node200data/kachungk/hcc_data/DMR_SVs/02.sv_dmr_enrichment/window_enrich/confident_dmr"

NO_PLOT_FLAG=""
if [[ "$N_PERM" -le 100 ]]; then
  NO_PLOT_FLAG="--no_plot"
  echo "[INFO] Phase 1 mode: nperm=${N_PERM}, plots disabled"
else
  echo "[INFO] Phase 2 mode: nperm=${N_PERM}, plots enabled"
fi

# Determine primary windows to run
if [[ -n "$PHASE2_PRIMARIES" ]]; then
  IFS=',' read -ra PRIMARY_WINDOWS <<< "$PHASE2_PRIMARIES"
else
  PRIMARY_WINDOWS=(25000 50000 100000 250000 500000)
fi

GROUP_MODES=("cnv_class" "tier")
MIN_SV=5

pids=()
labels=()

for GROUP in "${GROUP_MODES[@]}"; do
  OUTDIR="${BASE_OUTDIR}/${GROUP}_nperm${N_PERM}"
  mkdir -p "$OUTDIR"
  MASTER_LOG="${OUTDIR}/master_nperm${N_PERM}.log"
  echo "[$(date '+%F %T')] START — group_by=${GROUP}, n_perm=${N_PERM}" > "$MASTER_LOG"

  for PW in "${PRIMARY_WINDOWS[@]}"; do
    LBL="primary_$(( PW / 1000 ))kb"
    RUN_ID="${LBL}"
    LOG="${OUTDIR}/${RUN_ID}.log"

    echo "[$(date '+%F %T')] ▶ ${GROUP}/${RUN_ID}" | tee -a "$MASTER_LOG"

    mamba run -n renv Rscript "$SCRIPT" \
      --primary_window "$PW" \
      --windows        "$WINDOWS" \
      --n_perm         "$N_PERM" \
      --outdir         "$OUTDIR" \
      --run_id         "$RUN_ID" \
      --dmr_file       "$DMR_FILE" \
      --sv_strat_file  "$SV_FILE" \
      --group_by       "$GROUP" \
      --min_sv_per_group "$MIN_SV" \
      $NO_PLOT_FLAG \
      > "$LOG" 2>&1 &

    pids+=($!)
    labels+=("${GROUP}/${LBL}")
  done
done

echo "[INFO] ${#pids[@]} processes launched — waiting..."

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "[$(date '+%F %T')] ✓ ${labels[$i]}"
  else
    echo "[$(date '+%F %T')] ✗ ${labels[$i]} FAILED (exit $?)"
    failed=1
  fi
done

if [[ $failed -eq 0 ]]; then
  echo "[$(date '+%F %T')] All jobs completed successfully"
  if [[ "$N_PERM" -le 100 ]]; then
    echo ""
    echo "Next step — review enrichment_ratio and p_perm in:"
    echo "  ${BASE_OUTDIR}/cnv_class_nperm${N_PERM}/"
    echo "  ${BASE_OUTDIR}/tier_nperm${N_PERM}/"
    echo ""
    echo "Then run Phase 2 with selected windows, e.g.:"
    echo "  bash run_confident_dmr_enrichment.sh 1000 50000,100000"
  fi
else
  echo "[$(date '+%F %T')] Some jobs failed — check logs in ${BASE_OUTDIR}"
  exit 1
fi
