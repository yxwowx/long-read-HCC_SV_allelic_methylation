#!/bin/bash
# Run sv_dmr_enrichment.R for multiple PRIMARY_WINDOW values in parallel
# Usage: bash run_window_enrichment.sh [n_perm]

N_PERM="${1:-100}"
OUTDIR="/node200data/kachungk/hcc_data/DMR_SVs/window_enrichment"
WINDOWS="10000,50000,100000,500000,1000000"
SCRIPT="/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/pipeline/03_sv_dmr_enrichment.R"

PRIMARY_WINDOWS=(25000 50000 100000 500000)
LABELS=(25kb 50kb 100kb 500kb)

mkdir -p "$OUTDIR"
MASTER_LOG="$OUTDIR/master_nperm${N_PERM}.log"
echo "[$(date '+%F %T')] START — n_perm=${N_PERM}" > "$MASTER_LOG"
echo "출력 디렉토리: $OUTDIR" >> "$MASTER_LOG"
echo "===================================================" >> "$MASTER_LOG"

pids=()
for i in "${!PRIMARY_WINDOWS[@]}"; do
  PW="${PRIMARY_WINDOWS[$i]}"
  LBL="${LABELS[$i]}"
  RUN_ID="primary_${LBL}"
  LOG="$OUTDIR/${RUN_ID}_nperm${N_PERM}.log"

  echo "[$(date '+%F %T')] ▶ ${RUN_ID}  (primary_window=${PW}, n_perm=${N_PERM})" | tee -a "$MASTER_LOG"

  mamba run -n renv Rscript "$SCRIPT" \
    --primary_window "$PW" \
    --windows "$WINDOWS" \
    --n_perm "$N_PERM" \
    --outdir "$OUTDIR" \
    --run_id "$RUN_ID" \
    > "$LOG" 2>&1 &
  pids+=($!)
done

echo "[$(date '+%F %T')] 4개 프로세스 시작 완료 (PIDs: ${pids[*]})" | tee -a "$MASTER_LOG"

# Wait for all to finish
failed=0
for i in "${!pids[@]}"; do
  pid="${pids[$i]}"
  LBL="${LABELS[$i]}"
  if wait "$pid"; then
    echo "[$(date '+%F %T')] ✓ primary_${LBL} 완료" | tee -a "$MASTER_LOG"
  else
    echo "[$(date '+%F %T')] ✗ primary_${LBL} 실패 (exit $?)" | tee -a "$MASTER_LOG"
    failed=1
  fi
done

if [[ $failed -eq 0 ]]; then
  echo "[$(date '+%F %T')] 전체 완료" | tee -a "$MASTER_LOG"
else
  echo "[$(date '+%F %T')] 일부 실패 — 로그 확인 필요" | tee -a "$MASTER_LOG"
  exit 1
fi
