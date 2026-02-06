#!/bin/bash
# -*- coding: utf-8 -*-

# Usage:
#   ./submit_scvi.sh /path/to/logs
# or set LOG_DIR in the environment before calling:
#   export LOG_DIR=/path/to/logs
#   ./submit_scvi.sh

# ─── Parse out the log directory (first arg or env var, default "logs") ───
LOG_DIR="${1:-${LOG_DIR:-logs}}"
mkdir -p "$LOG_DIR"

# ─── Job parameters ───────────────────────────────────────────────────────
q=gpu-normal
mem=100000
ncpu=8

# ─── Parameter lists ──────────────────────────────────────────────────────
n_latent_list=(10 20 30)
n_layers_list=(3)
n_hidden_list=(128)
n_hvg_list=(3000 5000)

# ─── Loop over combinations ────────────────────────────────────────────────
for n_latent in "${n_latent_list[@]}"; do
  for n_layers in "${n_layers_list[@]}"; do
    for n_hidden in "${n_hidden_list[@]}"; do
      for n_hvg in "${n_hvg_list[@]}"; do

        job_name="scvi_best_integration${n_latent}_${n_layers}_${n_hidden}_${n_hvg}"
        echo "Submitting job: $job_name"

        bsub -q "$q" \
             -n "$ncpu" \
             -M "$mem" \
             -R "select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
             -gpu "mode=shared:j_exclusive=yes:gmem=20:num=1" \
             -G cellulargenetics-priority \
             -o "${LOG_DIR}/scvi_integration.${job_name}.%J_permissive.out" \
             "module load cellgen/conda; module load cellgen/scvi; \
              python -u /nfs/team298/sm54/BoneAtlasProject/src/integration/10_3_batch_correction_scvi_haem_best_integration_stroma_filtered_permissive_doublet_removed.py  \
                --n_latent $n_latent \
                --n_layers $n_layers \
                --n_hidden $n_hidden \
                --n_hvg $n_hvg"

      done
    done
  done
done
