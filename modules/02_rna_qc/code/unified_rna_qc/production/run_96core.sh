#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/nfsdata/apadetect/analysis_projects/reprogramming_unified_rna_qc_v3.2_release
SNAKEMAKE=/root/miniconda3/envs/single-cell-qc/bin/snakemake
PYTHON=/root/miniconda3/envs/single-cell-qc/bin/python
LOG="$ROOT/production_qc_v3.2.log"
PID_FILE="$ROOT/production_qc_v3.2.pid"
MODE=${1:-status}

preflight() {
  cd "$ROOT"
  "$SNAKEMAKE" --snakefile production/Snakefile \
    --config run_mode=resume \
    --cores 1 preflight --printshellcmds
}

dryrun() {
  cd "$ROOT"
  "$SNAKEMAKE" --snakefile production/Snakefile \
    --config run_mode=resume \
    --cores 96 --jobs 96 \
    --rerun-incomplete --latency-wait 120 \
    --printshellcmds --dry-run \
    --resources mem_mb=512000
}

case "$MODE" in
  preflight)
    preflight
    ;;
  dryrun)
    preflight
    dryrun
    ;;
  start)
    if [[ -s "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
      echo "QC production is already running: PID $(cat "$PID_FILE")" >&2
      exit 1
    fi
    preflight
    dryrun
    cd "$ROOT"
    nohup "$SNAKEMAKE" --snakefile production/Snakefile \
      --config run_mode=resume \
      --cores 96 --jobs 96 \
      --rerun-incomplete --latency-wait 120 \
      --printshellcmds \
      --resources mem_mb=512000 \
      > "$LOG" 2>&1 &
    echo "$!" > "$PID_FILE"
    echo "Started QC production: PID $(cat "$PID_FILE")"
    echo "Log: $LOG"
    ;;
  status)
    if [[ -s "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
      echo "RUNNING PID $(cat "$PID_FILE")"
    else
      echo "NOT_RUNNING"
    fi
    "$PYTHON" "$ROOT/production/check_acceptance.py" --root "$ROOT"
    if [[ -f "$LOG" ]]; then
      tail -n 30 "$LOG"
    fi
    ;;
  acceptance)
    "$PYTHON" "$ROOT/production/check_acceptance.py" --root "$ROOT" --strict
    ;;
  *)
    echo "Usage: $0 {preflight|dryrun|start|status|acceptance}" >&2
    exit 2
    ;;
esac
