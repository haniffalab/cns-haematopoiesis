#!/bin/bash
# submit_generate_velocity.sh
#
# EXAMPLE USAGE:
#   Edit the variables below, then run:
#       bash submit_generate_velocity.sh
#
# This submits one job that runs:
#   python generate_velocity_h5ads.py --input_dir ... --output_dir ... --merged_name ...

#########################
# USER-EDITABLE PARAMETERS
#########################

INPUT_DIR="/nfs/cellgeni/tickets/tic-3721/actions/data/STARsolo"
OUTPUT_DIR="/lustre/scratch124/cellgen/haniffa/users/sm54/data/velocity"
MERGED_NAME="velocity_merged.h5ad"

QUEUE="normal"
MEM=650000
NCPU=8
JOB="generate_velocity"


# path to your script
SCRIPT="/nfs/team298/sm54/BoneAtlasProject/src/processing/generate_velocity_h5ads.py"

#########################
# SUBMISSION
#########################

bsub -J "$JOB" \
     -q "$QUEUE" \
     -n "$NCPU" \
     -M "$MEM" \
     -R "select[mem>$MEM] rusage[mem=$MEM]" \
     -G cellulargenetics-priority \
     -o "%J.${JOB}.out" \
  "module load cellgen/conda;
   conda activate boneatlas;
   python -u $SCRIPT \
       --input_dir $INPUT_DIR \
       --output_dir $OUTPUT_DIR \
       --merged_name $MERGED_NAME"
