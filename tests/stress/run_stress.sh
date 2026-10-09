#!/usr/bin/env bash
# Stress test: compare dubsteppy with R DUBStepR on 16 generated datasets/parameter sets.
# Usage (from the repo root): tests/stress/run_stress.sh [workdir]   (default: a new dir in $TMPDIR)
# Runs inside container/dubsteppy.sif; generated inputs/outputs go to workdir (keep it off GPFS).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SIF=${SIF:-$HERE/../../container/dubsteppy.sif}
WORK=${1:-$(mktemp -d "${TMPDIR:-/tmp}/dubsteppy_stress.XXXX")}
mkdir -p "$WORK"
echo "workdir: $WORK"
apptainer exec --cleanenv "$SIF" python -I "$HERE/make_jobs.py" "$WORK"
apptainer exec --cleanenv "$SIF" Rscript "$HERE/run_r.R" "$WORK" 2>&1 | tr '\r' '\n' | grep -E "done in|ERROR"
apptainer exec --cleanenv "$SIF" python -I "$HERE/compare.py" "$WORK" 2>&1 | tee "$WORK/summary.txt"
apptainer exec --cleanenv "$SIF" Rscript "$HERE/irlba_first_step_check.R" "$WORK" 2>&1 | grep "step=" | tee -a "$WORK/summary.txt"
