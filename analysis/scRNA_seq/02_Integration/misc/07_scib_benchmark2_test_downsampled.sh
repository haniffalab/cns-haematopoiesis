#!/bin/bash

# Usage:
#   ./07_scib_benchmark2_test_downsampled.sh /path/to/adata_file.h5ad /path/to/logs
# Example Usage 
#bash /nfs/team298/sm54/BoneAtlasProject/src/integration/07_scib_benchmark2_test_downsampled.sh /lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_test_int.h5ad /nfs/team298/sm54/BoneAtlasProject/src/integration 
#Job <295000>
# or set LOG_DIR in the environment:
#   export LOG_DIR=/path/to/logs
#   ./submit_scib_one.sh /path/to/adata_file.h5ad

# ─── Parse args ───────────────────────────────────────────────────────────
INTEGRATED_FILE="${1:?Must pass path to one .h5ad}"   # e.g. /…/adata_nlatent15_nlayers4_…hvg7500.h5ad
LOG_DIR="${2:-${LOG_DIR:-logs}}"
mkdir -p "$LOG_DIR"

# ─── Job params ────────────────────────────────────────────────────────────
q=gpu-normal
mem=200000
ncpu=4

# ─── Paths ─────────────────────────────────────────────────────────────────
SCIB_SCRIPT=/nfs/team298/sm54/BoneAtlasProject/src/integration/07_scib_benchmark_2.py
SCIB_CONTAINER=/nfs/cellgeni/singularity/images/scib-1.1.5.sif
PCA=/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_test_unintegrated_pca.parquet
OUTDIR=/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_output
mkdir -p "$OUTDIR"

# ─── Derive job name ─────────────────────────────────────────────────────────
base=$(basename "${INTEGRATED_FILE%.h5ad}")
job_name="scib_${base}_downsampled_new_method"
echo "Submitting single‐dataset job: $job_name"

# ─── Submit ──────────────────────────────────────────────────────────────────
bsub -q "$q" \
     -n "$ncpu" \
     -J "$job_name" \
     -M "$mem" \
     -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
     -gpu "mode=shared:j_exclusive=yes:gmem=10:num=1" \
     -G cellulargenetics-priority \
     -o "${LOG_DIR}/ooo.${job_name}.%J.out" \
   "module load cellgen/singularity; \
     singularity exec --cleanenv --home /tmp/$(whoami) --bind /lustre,/nfs /nfs/cellgeni/singularity/images/scib-1.1.5.sif \
    python -u $SCIB_SCRIPT \
      --input    $INTEGRATED_FILE \
      --pca      $PCA \
      --output   $OUTDIR/${job_name}_scib.csv \
      --batch-key Technology \
      --label-key majority_voting_Level2 \
      --neighbors 15 \
      --jobs     $ncpu \
  "
