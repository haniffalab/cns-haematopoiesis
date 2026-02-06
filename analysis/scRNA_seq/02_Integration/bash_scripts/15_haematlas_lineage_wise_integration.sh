#!/bin/bash
# -*- coding: utf-8 -*-
# Usage: ./submit_scvi_immune.sh /path/to/logs

LOG_DIR="${1:-${LOG_DIR:-logs}}"
mkdir -p "$LOG_DIR"

q=gpu-normal
mem=50000
ncpu=4

n_latent_list=(10 20 30)
n_layers_list=(2)
n_hidden_list=(128)
n_hvg_list=(3000 5000 7500)

DATA_DIR="/nfs/team298/sm54/BoneAtlasProject/data/bone_atlas_anndatas_immune_subsets"

for h5ad_file in "$DATA_DIR"/bone_*_filtered_with_metadata_scvi.h5ad; do
    cell_type=$(basename "$h5ad_file" | sed 's/bone_//' | sed 's/_filtered_with_metadata_scvi.h5ad//')

    for n_latent in "${n_latent_list[@]}"; do
      for n_layers in "${n_layers_list[@]}"; do
        for n_hidden in "${n_hidden_list[@]}"; do
          for n_hvg in "${n_hvg_list[@]}"; do

            job_name="haem_atlas_scvi_${cell_type}_latent${n_latent}_layers${n_layers}_hidden${n_hidden}_hvg${n_hvg}"
            echo "Submitting job: $job_name"

            bsub -q "$q" \
                 -n "$ncpu" \
                 -M "$mem" \
                 -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
                 -gpu "mode=shared:j_exclusive=yes:gmem=20:num=1" \
                 -G cellulargenetics-priority \
                 -o "${LOG_DIR}/scvi_${cell_type}.%J.out" \
                 "module load cellgen/conda; module load cellgen/scvi; \
                  python -u /nfs/team298/sm54/BoneAtlasProject/src/integration/15_haem_atlas_lineage_wise_integration.py \
                    --h5ad \"$h5ad_file\" \
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
