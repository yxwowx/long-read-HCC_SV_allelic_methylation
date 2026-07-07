#!/bin/bash
# Run sv_dmr_enrichment.R across one or more PRIMARY_WINDOW values in parallel.
# Consolidates three formerly separate scripts into three modes of one script:
#   window     - quick multi-window scan with default inputs (was run_window_enrichment.sh)
#   compare    - fixed 50kb vs 100kb comparison, then compare_window_runs.R (was run_compare_50kb_100kb.sh)
#   confident  - two-phase, dual-grouping (cnv_class/tier) run using confident DMRs
#                (was run_confident_dmr_enrichment.sh)
#
# Usage:
#   bash run_enrichment.sh window    [n_perm]
#   bash run_enrichment.sh compare   [n_perm]
#   bash run_enrichment.sh confident [n_perm] [primary_windows_csv]
#
# Any mode's defaults can be overridden with flags:
#   --outdir DIR --windows LIST --primary_windows LIST
#   --dmr_file PATH --sv_file PATH --group_by MODE[,MODE...] --min_sv_per_group N
#   --no_plot

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/.env"

SCRIPT="$SCRIPT_DIR/sv_dmr_enrichment.R"
COMPARE_SCRIPT="$SCRIPT_DIR/compare_window_runs.R"

MODE="${1:-}"
if [[ -z "$MODE" || "$MODE" != "window" && "$MODE" != "compare" && "$MODE" != "confident" ]]; then
  echo "Usage: bash run_enrichment.sh {window|compare|confident} [n_perm] [primary_windows_csv]" >&2
  exit 1
fi
shift

N_PERM="${1:-100}"
[[ $# -gt 0 ]] && shift
PHASE2_PRIMARIES="${1:-}"
[[ $# -gt 0 ]] && shift

# ── Mode defaults ─────────────────────────────────────────────────────────────
DMR_FILE=""
SV_FILE=""
GROUP_BY=""
MIN_SV_PER_GROUP=5
NO_PLOT_FLAG=""

case "$MODE" in
  window)
    OUTDIR="$HCC_DATA_DIR/DMR_SVs/window_enrichment"
    WINDOWS="10000,50000,100000,500000,1000000"
    PRIMARY_WINDOWS=(25000 50000 100000 500000)
    LABELS=(25kb 50kb 100kb 500kb)
    ;;
  compare)
    N_PERM=1000
    OUTDIR="$HCC_DATA_DIR/DMR_SVs/window_enrichment/nperm_1000"
    WINDOWS="10000,50000,100000,500000,1000000"
    PRIMARY_WINDOWS=(50000 100000)
    LABELS=(50kb 100kb)
    ;;
  confident)
    DMR_FILE="$HCC_DATA_DIR/DMR_SVs/01.DMR_recurrence/confident_dmr_per_patient.csv.gz"
    SV_FILE="$HCC_DATA_DIR/DMR_SVs/sv_tad_ctcf_annotation.csv.gz"
    WINDOWS="10000,25000,50000,100000,250000,500000,1000000"
    BASE_OUTDIR="$HCC_DATA_DIR/DMR_SVs/02.sv_dmr_enrichment/window_enrich/confident_dmr"
    GROUP_BY="cnv_class,tier"
    if [[ -n "$PHASE2_PRIMARIES" ]]; then
      IFS=',' read -ra PRIMARY_WINDOWS <<< "$PHASE2_PRIMARIES"
    else
      PRIMARY_WINDOWS=(25000 50000 100000 250000 500000)
    fi
    LABELS=()
    for pw in "${PRIMARY_WINDOWS[@]}"; do LABELS+=("$(( pw / 1000 ))kb"); done
    if [[ "$N_PERM" -le 100 ]]; then
      NO_PLOT_FLAG="--no_plot"
      echo "[INFO] Phase 1 mode: nperm=${N_PERM}, plots disabled"
    else
      echo "[INFO] Phase 2 mode: nperm=${N_PERM}, plots enabled"
    fi
    ;;
esac

# ── Flag overrides (apply to any mode) ────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --windows) WINDOWS="$2"; shift 2 ;;
    --primary_windows)
      IFS=',' read -ra PRIMARY_WINDOWS <<< "$2"
      LABELS=()
      for pw in "${PRIMARY_WINDOWS[@]}"; do LABELS+=("$(( pw / 1000 ))kb"); done
      shift 2 ;;
    --dmr_file) DMR_FILE="$2"; shift 2 ;;
    --sv_file) SV_FILE="$2"; shift 2 ;;
    --group_by) GROUP_BY="$2"; shift 2 ;;
    --min_sv_per_group) MIN_SV_PER_GROUP="$2"; shift 2 ;;
    --no_plot) NO_PLOT_FLAG="--no_plot"; shift ;;
    *) echo "Unknown flag: $1" >&2; exit 1 ;;
  esac
done

IFS=',' read -ra GROUP_MODES <<< "${GROUP_BY:-__none__}"

mkdir -p "${BASE_OUTDIR:-$OUTDIR}"

run_one_group() {
  local group="$1"
  local outdir="$2"
  mkdir -p "$outdir"
  local master_log="$outdir/master_nperm${N_PERM}.log"
  echo "[$(date '+%F %T')] START mode=$MODE group=${group:-none} n_perm=${N_PERM}" | tee "$master_log"

  local pids=() labels=()
  for i in "${!PRIMARY_WINDOWS[@]}"; do
    local pw="${PRIMARY_WINDOWS[$i]}"
    local lbl="${LABELS[$i]}"
    local run_id="primary_${lbl}"
    local log="$outdir/${run_id}.log"

    echo "[$(date '+%F %T')] Launching ${group:+$group/}${run_id} (primary_window=${pw})" | tee -a "$master_log"

    local args=(--primary_window "$pw" --windows "$WINDOWS" --n_perm "$N_PERM"
                --outdir "$outdir" --run_id "$run_id")
    [[ -n "$DMR_FILE" ]] && args+=(--dmr_file "$DMR_FILE")
    [[ -n "$SV_FILE" ]] && args+=(--sv_strat_file "$SV_FILE")
    [[ -n "$group" ]] && args+=(--group_by "$group" --min_sv_per_group "$MIN_SV_PER_GROUP")
    [[ -n "$NO_PLOT_FLAG" ]] && args+=("$NO_PLOT_FLAG")

    mamba run -n renv Rscript "$SCRIPT" "${args[@]}" > "$log" 2>&1 &
    pids+=($!)
    labels+=("${group:+$group/}${run_id}")
  done

  local failed=0
  for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
      echo "[$(date '+%F %T')] OK   ${labels[$i]}" | tee -a "$master_log"
    else
      echo "[$(date '+%F %T')] FAIL ${labels[$i]} (log: $outdir/${labels[$i]##*/}.log)" | tee -a "$master_log"
      failed=1
    fi
  done
  return $failed
}

OVERALL_FAILED=0
for group in "${GROUP_MODES[@]}"; do
  if [[ "$group" == "__none__" ]]; then
    run_one_group "" "$OUTDIR" || OVERALL_FAILED=1
  else
    run_one_group "$group" "$BASE_OUTDIR/${group}_nperm${N_PERM}" || OVERALL_FAILED=1
  fi
done

if [[ "$OVERALL_FAILED" -ne 0 ]]; then
  echo "[$(date '+%F %T')] Some jobs failed — see logs above" >&2
  exit 1
fi

echo "[$(date '+%F %T')] All jobs completed successfully"

if [[ "$MODE" == "compare" ]]; then
  echo "[$(date '+%F %T')] Running compare_window_runs.R..."
  compare_log="$OUTDIR/compare.log"
  mamba run -n renv Rscript "$COMPARE_SCRIPT" 2>&1 | tee "$compare_log"
  echo "[$(date '+%F %T')] Done — results in $OUTDIR"
fi

if [[ "$MODE" == "confident" && "$N_PERM" -le 100 ]]; then
  echo ""
  echo "Next step — review enrichment_ratio and p_perm in:"
  for group in "${GROUP_MODES[@]}"; do
    echo "  $BASE_OUTDIR/${group}_nperm${N_PERM}/"
  done
  echo ""
  echo "Then run Phase 2 with selected windows, e.g.:"
  echo "  bash run_enrichment.sh confident 1000 50000,100000"
fi
