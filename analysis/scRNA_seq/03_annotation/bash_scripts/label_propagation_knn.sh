#!/bin/bash
# -*- coding: utf-8 -*-

# ────────────────────────────── USAGE ───────────────────────────────
# Run with defaults:
#   ./submit_label_prop_cpu.sh
#
# Or override specific arguments:
#   ./label_propagation_knn.sh /my/logs /my/data/adata.h5ad published_anno X_scVI
# ───────────────────────────────────────────────────────────────────

# ─── Arguments with defaults ────────────────────────────────────────
LOG_DIR="${1:-logs}"   # default: logs/
H5AD_PATH="${2:-/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_best_embedding_with_hvg_label_propagation.h5ad}"
LABEL_COL="${3:-published_anno}"
EMBED_KEY="${4:-X_scVI}"

mkdir -p "$LOG_DIR"

# ─── Fixed job parameters (edit if needed) ──────────────────────────
QUEUE="long"
MEM=300000      # MB (200 GB)
NCPU=8
GROUP="cellulargenetics-priority"
ENV_NAME="boneatlas"

# ─── Job name ──────────────────────────────────────────────────────
JOB_NAME="label_prop_cpu"

# ─── Submit job ────────────────────────────────────────────────────
bsub -q "$QUEUE" \
     -n "$NCPU" \
     -M "$MEM" \
     -R "select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" \
     -G "$GROUP" \
     -o "${LOG_DIR}/${JOB_NAME}.%J.out" \
     -e "${LOG_DIR}/${JOB_NAME}.%J.err" \
     "module load cellgen/conda; conda activate $ENV_NAME; \
      python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/label_propagation_knn.py \
        --h5ad $H5AD_PATH \
        --embed-key $EMBED_KEY \
        --label-col $LABEL_COL \
        --k 25 \
        --min-score 0.4 \
        --batch-size 200000 \
        --workers -1 \
        --out /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/boneatlas_propagated_labels_knn_25Dec_0_4_Conf.csv"
