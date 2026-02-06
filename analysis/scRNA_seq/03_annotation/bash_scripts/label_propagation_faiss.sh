#!/bin/bash
# -*- coding: utf-8 -*-

# ────────────────────────────── USAGE ───────────────────────────────
# Run with defaults:
#   ./submit_label_prop_gpu.sh
#
# Or override specific arguments:
#   ./label_propagation_faiss.sh /my/logs /my/data/adata.h5ad published_anno X_scVI
# ───────────────────────────────────────────────────────────────────

# ─── Arguments with defaults ────────────────────────────────────────
LOG_DIR="${1:-logs}"   # default: logs/
H5AD_PATH="${2:-/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_best_embedding_with_hvg_label_propagation.h5ad}"
LABEL_COL="${3:-published_anno}"
EMBED_KEY="${4:-X_scVI}"

mkdir -p "$LOG_DIR"

# ─── Fixed job parameters (edit if needed) ──────────────────────────
QUEUE="gpu-normal"
MEM=300000        # MB (200 GB)   ← comment matches value now
NCPU=8
GROUP="cellulargenetics-priority"
ENV_NAME="faiss_new"

# ─── GPU params ─────────────────────────────────────────────────────
GPUMODE="mode=shared:j_exclusive=yes:gmem=40:num=1"

# ─── Job name ──────────────────────────────────────────────────────
JOB_NAME="label_prop_gpu"

# ─── Submit job ────────────────────────────────────────────────────
bsub -q "$QUEUE" \
     -n "$NCPU" \
     -M "$MEM" \
     -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" \
     -gpu "$GPUMODE" \
     -G "$GROUP" \
     -o "${LOG_DIR}/${JOB_NAME}.%J.out" \
     "module load cellgen/conda; conda activate $ENV_NAME; \
      python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/label_propagation_faiss.py \
        --h5ad $H5AD_PATH \
        --embed-key $EMBED_KEY \
        --label-col $LABEL_COL \
        --k 25 \
        --min-score 0.6 \
        --batch-size 200000 \
        --index ivf \
        --nlist 16384 \
        --nprobe 32 \
        --gpu-id 0 \
        --out /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_propagated_labels_gpu_faiss_28Nov_0_6_Conf.csv"
