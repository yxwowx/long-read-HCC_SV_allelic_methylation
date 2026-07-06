#!/bin/bash
# Run sv_dmr_enrichment.R with primary_50kb and primary_100kb in parallel,
# then compare results.

set -euo pipefail

OUTDIR="/node200data/kachungk/hcc_data/DMR_SVs/window_enrichment/nperm_1000"
WINDOWS="10000,50000,100000,500000,1000000"
N_PERM=1000
SCRIPT="/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/pipeline/03_sv_dmr_enrichment.R"
COMPARE_SCRIPT="/home/kachungk/script/SV-DMR/remodeled_constitutional_AMR/post_processing/compare_window_runs.R"

mkdir -p "$OUTDIR"
MASTER_LOG="$OUTDIR/master_compare.log"
echo "[$(date '+%F %T')] START — n_perm=${N_PERM}" | tee "$MASTER_LOG"

# ── 실행 1: primary_50kb ──────────────────────────────────────────────────────
LOG50="$OUTDIR/primary_50kb.log"
echo "[$(date '+%F %T')] ▶ primary_50kb 시작" | tee -a "$MASTER_LOG"
mamba run -n renv Rscript "$SCRIPT" \
  --primary_window 50000 \
  --windows        "$WINDOWS" \
  --n_perm         "$N_PERM" \
  --outdir         "$OUTDIR" \
  --run_id         "primary_50kb" \
  > "$LOG50" 2>&1 &
PID50=$!

# ── 실행 2: primary_100kb ─────────────────────────────────────────────────────
LOG100="$OUTDIR/primary_100kb.log"
echo "[$(date '+%F %T')] ▶ primary_100kb 시작" | tee -a "$MASTER_LOG"
mamba run -n renv Rscript "$SCRIPT" \
  --primary_window 100000 \
  --windows        "$WINDOWS" \
  --n_perm         "$N_PERM" \
  --outdir         "$OUTDIR" \
  --run_id         "primary_100kb" \
  > "$LOG100" 2>&1 &
PID100=$!

echo "[$(date '+%F %T')] 두 프로세스 시작 (PIDs: $PID50 $PID100)" | tee -a "$MASTER_LOG"

# ── 완료 대기 ─────────────────────────────────────────────────────────────────
FAILED=0
if wait "$PID50"; then
  echo "[$(date '+%F %T')] ✓ primary_50kb 완료" | tee -a "$MASTER_LOG"
else
  echo "[$(date '+%F %T')] ✗ primary_50kb 실패 — 로그: $LOG50" | tee -a "$MASTER_LOG"
  FAILED=1
fi

if wait "$PID100"; then
  echo "[$(date '+%F %T')] ✓ primary_100kb 완료" | tee -a "$MASTER_LOG"
else
  echo "[$(date '+%F %T')] ✗ primary_100kb 실패 — 로그: $LOG100" | tee -a "$MASTER_LOG"
  FAILED=1
fi

if [[ $FAILED -ne 0 ]]; then
  echo "[$(date '+%F %T')] 일부 실패 — compare 스크립트 실행 중단" | tee -a "$MASTER_LOG"
  exit 1
fi

# ── 비교 스크립트 실행 ─────────────────────────────────────────────────────────
echo "[$(date '+%F %T')] ▶ compare_window_runs.R 실행 중..." | tee -a "$MASTER_LOG"
COMPARE_LOG="$OUTDIR/compare.log"
mamba run -n renv Rscript "$COMPARE_SCRIPT" \
  2>&1 | tee "$COMPARE_LOG" | tee -a "$MASTER_LOG"

echo "[$(date '+%F %T')] 전체 완료 — 결과: $OUTDIR" | tee -a "$MASTER_LOG"
