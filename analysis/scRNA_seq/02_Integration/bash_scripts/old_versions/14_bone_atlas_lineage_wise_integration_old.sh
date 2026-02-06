#!/bin/bash
# -*- coding: utf-8 -*-
# Usage: ./submit_scvi_nonimmune.sh /path/to/logs

# ─── Parse log directory ───────────────────────────────────────────────
LOG_DIR="${1:-${LOG_DIR:-logs}}"
mkdir -p "$LOG_DIR"

# ─── Job parameters ─────────────────────────────────────────────────────
q=gpu-normal
mem=150000
ncpu=4

# ─── SCVI hyperparameters ───────────────────────────────────────────────
n_latent_list=(30)
n_layers_list=(2)
n_hidden_list=(128)
n_hvg_list=(3000)

# ─── Directory with non-immune h5ads ───────────────────────────────────
DATA_DIR="/nfs/team298/sm54/BoneAtlasProject/data/bone_atlas_anndatas_non_immune"

# ─── Loop over all non-immune .h5ad files ───────────────────────────────
for h5ad_file in "$DATA_DIR"/bone_*_filtered_with_metadata_scvi.h5ad; do
    # Extract cell type name from filename
    cell_type=$(basename "$h5ad_file" | sed 's/bone_//' | sed 's/_filtered_with_metadata_scvi.h5ad//')

    for n_latent in "${n_latent_list[@]}"; do
        for n_layers in "${n_layers_list[@]}"; do
            for n_hidden in "${n_hidden_list[@]}"; do
                for n_hvg in "${n_hvg_list[@]}"; do

                    job_name="scvi_${cell_type}_latent${n_latent}_layers${n_layers}_hidden${n_hidden}_hvg${n_hvg}"
                    echo "Submitting job: $job_name"

                    bsub -q "$q" \
                         -n "$ncpu" \
                         -M "$mem" \
                         -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
                         -gpu "mode=shared:j_exclusive=yes:gmem=20:num=1" \
                         -G cellulargenetics-priority \
                         -o "${LOG_DIR}/scvi_${cell_type}.%J.out" \
                         "module load cellgen/conda; module load cellgen/scvi; \
                          python -u /nfs/team298/sm54/BoneAtlasProject/src/integration/14_bone_atlas_lineage_wise_integration.py \
                            --n_latent $n_latent \
                            --n_layers $n_layers \
                            --n_hidden $n_hidden \
                            --n_hvg $n_hvg \
                            --data_dir $DATA_DIR"
                done
            done
        done
    done
done
