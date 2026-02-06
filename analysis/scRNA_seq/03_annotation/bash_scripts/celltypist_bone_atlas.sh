#!/bin/bash

# ────────────────────────────────────────────────
# LSF settings
q=normal
mem=350000        # memory per job in MB
ncpu=1           # number of cores
# ────────────────────────────────────────────────

# ─────── Configure these three paths ────────────
INPUT_H5AD="/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/best_integration/processed_adata/Best_Integration_bone_atlas_nlatent30_nlayers3_nhidden128_nhvg3000.h5ad"
MODEL_PKL="/nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/Human_WholeEmbryo_Public.pkl"
OUTPUT_DIR="/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/best_integration/celltypist_results"
# ────────────────────────────────────────────────

bsub -q "$q" \
     -n"${ncpu}" \
     -M"${mem}" \
     -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
     -G team298 \
     -o "%J.celltypist_bone_atlas_integrated.out" \
"module load cellgen/conda; \
 conda activate boneatlas; \
 mkdir -p ${OUTPUT_DIR}; \
 python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/celltypist_bone_atlas.py \
    --input_h5ad  "$INPUT_H5AD" \
    --model_pkl   "$MODEL_PKL" \
    --output_dir  "$OUTPUT_DIR"
"
