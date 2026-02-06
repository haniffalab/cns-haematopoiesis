#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scvi_integrate_any.sh /path/to/input.h5ad [out_dir] [logs_dir]
#
# Requires:
#   - You must provide the input .h5ad path (first arg).
#   - Optionally provide an output dir (default: same as input dir).
#   - Optionally provide a logs dir (default: ./logs).

# ── Args ─────────────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 /path/to/input.h5ad [out_dir] [logs_dir]" >&2
  exit 1
fi
H5="$1"
OUT_DIR="${2:-$(dirname "$H5")}"
LOG_DIR="${3:-logs}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ── LSF resources (tweak as needed) ──────────────────────────────────────
q=gpu-normal
mem=50000
ncpu=4

# ── Python script path ───────────────────────────────────────────────────
PY="/nfs/team298/sm54/BoneAtlasProject/src/integration/13_standard_reintegation_script.py"

# ── Submit single job ────────────────────────────────────────────────────
job_name="scvi_integrate_any_$(basename "$H5" .h5ad)"
echo "Submitting: $job_name  ($H5 → $OUT_DIR)"

bsub -q "$q" \
     -n "$ncpu" \
     -M "$mem" \
     -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
     -gpu "mode=shared:j_exclusive=yes:gmem=20:num=1" \
     -G cellulargenetics-priority \
     -o "${LOG_DIR}/${job_name}.%J.out" \
     "module load cellgen/conda; module load cellgen/scvi; \
      python -u \"$PY\" \
        --input \"$H5\" \
        --out_dir \"$OUT_DIR\""
