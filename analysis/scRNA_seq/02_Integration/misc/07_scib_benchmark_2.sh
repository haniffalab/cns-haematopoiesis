#!/bin/bash

# Usage:
#   ./submit_scib.sh /path/to/logs
# or set LOG_DIR in the environment:
#   export LOG_DIR=/path/to/logs
#   ./submit_scib.sh

# ─── Log directory ────────────────────────────────────────────────────────
LOG_DIR="${1:-${LOG_DIR:-logs}}"
mkdir -p "$LOG_DIR"

# ─── Job params ──────────────────────────────────────────────────────────
q=gpu-normal
mem=300000
ncpu=8

# ─── Paths ────────────────────────────────────────────────────────────────
SCIB_SCRIPT=/nfs/team298/sm54/BoneAtlasProject/src/processing/07_scib_benchmark_2.py
SCIB_CONTAINER=/nfs/cellgeni/singularity/images/scib-1.1.5.sif
INTEGRATED_DIR=/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/processed_adata
PCA=/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata_pca_embedding.parquet
OUTDIR=/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_output
mkdir -p "$OUTDIR"

# ─── Loop over all integrated Adatas ─────────────────────────────────────
for integrated in "${INTEGRATED_DIR}"/*.h5ad; do
  base=$(basename "${integrated%.h5ad}")
  job_name="scib_${base}"
  echo "Submitting job: $job_name"

  bsub -q "$q" \
       -n "$ncpu" \
       -M "$mem" \
       -J scib_faster \
       -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
       -gpu "mode=shared:j_exclusive=yes:gmem=20:num=1" \
       -G cellulargenetics-priority \
       -o "${LOG_DIR}/ooo.${job_name}.%J.out" \
       -e "${LOG_DIR}/eee.${job_name}.%J.err" \
      "module load cellgen/singularity; \
     singularity exec --cleanenv --home /tmp/$(whoami) --bind /lustre,/nfs /nfs/cellgeni/singularity/images/scib-1.1.5.sif \
      python $SCIB_SCRIPT \
        --input    $integrated \
        --pca      $PCA \
        --output   $OUTDIR/${job_name}_scib.csv \
        --batch-key Technology \
        --label-key majority_voting_Level1 \
        --neighbors 15 \
        --jobs     $ncpu \
    "
done
