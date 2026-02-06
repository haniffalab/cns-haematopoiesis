#!/bin/bash
set -euo pipefail

UNINTEGRATED_ADATA="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_filtered_with_metadata_scvi.h5ad"
INTEGRATED_DIR="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/processed_adata"
OUTPUT_DIR="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_output"
SCRIPT="06_evaluate_batch_correction.py"
BATCH='Technology'
LABEL='majority_voting_Level2'
CPU=8
MEM=650000
GROUP="team298"
QUEUE="long"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

for integrated_adata in "$INTEGRATED_DIR"/*.h5ad; do
    filename=$(basename "$integrated_adata" .h5ad)
    output_csv="${OUTPUT_DIR}/${filename}_evaluation_metrics.csv"

    echo "Submitting test job for: $integrated_adata"
    bsub -n $CPU \
         -M $MEM \
         -J scib_test \
         -R "select[mem>${MEM}] rusage[mem=${MEM}]" \
         -G $GROUP \
         -q $QUEUE \
         -o "logs/scib.${filename}_test.%J.out" \
         "module load cellgen/singularity; \
          singularity exec --cleanenv --home /tmp/$(whoami) --bind /lustre,/nfs \
              /nfs/cellgeni/singularity/images/scib-1.1.5.sif \
              python $SCRIPT --adata_uninteg $UNINTEGRATED_ADATA \
                             --adata_int $integrated_adata \
                             --output_csv $output_csv \
                             --batch_key $BATCH \
                             --label_key $LABEL"

    # break out after the first file so only one test job is submitted
    break
done
