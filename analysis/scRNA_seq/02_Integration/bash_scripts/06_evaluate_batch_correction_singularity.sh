#!/bin/bash

# Script to submit multiple batch correction evaluation jobs
set -euo pipefail

# Define constants
UNINTEGRATED_ADATA="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_filtered_with_metadata_scvi.h5ad"

INTEGRATED_DIR="/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_scvi_runs/processed_adata"
OUTPUT_DIR="/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_scvi_runs/scib_output"
SCRIPT="06_evaluate_batch_correction.py"
BATCH='Technology'
LABEL='majority_voting_Level2'

# Job Submission Settings
CPU=8                           # Number of CPUs
MEM=650000                       # Memory in MB
GROUP="cellulargenetics-priority"                # Group name
QUEUE="long"              # Queue to submit to


# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Iterate over all integrated .h5ad files in the directory
index=0
for integrated_adata in "$INTEGRATED_DIR"/*.h5ad; do
    # Extract filename without path and extension for naming output files
    filename=$(basename "$integrated_adata" .h5ad)
    output_csv="${OUTPUT_DIR}/bone_atlas_${filename}_evaluation_metrics.csv"

    # Submit the job
    echo "Submitting job for: $integrated_adata"
    bsub -n $CPU \
         -M $MEM \
         -R "select[mem>${MEM}] rusage[mem=${MEM}]" \
         -G $GROUP \
         -q $QUEUE \
         -o "logs/bone_atlas_scib_${filename}.%J.out" \
         "module load cellgen/singularity; \
     singularity exec --cleanenv --home /tmp/$(whoami) --bind /lustre,/nfs /nfs/cellgeni/singularity/images/scib-1.1.5.sif \
          python -u /nfs/team298/sm54/BoneAtlasProject/src/integration/$SCRIPT --adata_uninteg $UNINTEGRATED_ADATA \
                         --adata_int $integrated_adata \
                         --output_csv $output_csv \
                         --batch_key $BATCH \
                         --label_key $LABEL"

    # Increment index for file tracking
    #index=$((index + 1))
done
