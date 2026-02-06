#!/bin/bash
set -euo pipefail

# paths to your test datasets
UNINTEGRATED_ADATA="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_test_unintegrated.h5ad"
INTEGRATED_ADATA="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_test_int.h5ad"

OUTPUT_DIR="/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_scvi_runs/scib_output"
SCRIPT="06_evaluate_batch_correction.py"
BATCH="Technology"
LABEL="majority_voting_Level2"
CPU=4
MEM=400000
GROUP="team298"
QUEUE="long"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# derive a job‐friendly name
filename=$(basename "$INTEGRATED_ADATA" .h5ad)
output_csv="${OUTPUT_DIR}/${filename}_evaluation_metrics.csv"

echo "Submitting single test job for: $INTEGRATED_ADATA"
bsub -n $CPU \
     -M $MEM \
     -R "select[mem>${MEM}] rusage[mem=${MEM}]" \
     -G $GROUP \
     -q $QUEUE \
     -o "logs/ooo.downsampled_test_${MEM}_GB.%J.out" \
  "module load cellgen/singularity; \
   singularity exec --cleanenv --home /tmp/$(whoami) --bind /lustre,/nfs \
     /nfs/cellgeni/singularity/images/scib-1.1.5.sif \
     python -u $SCRIPT \
       --adata_uninteg $UNINTEGRATED_ADATA \
       --adata_int     $INTEGRATED_ADATA \
       --output_csv   $output_csv \
       --batch_key    $BATCH \
       --label_key    $LABEL"
