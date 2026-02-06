#!/bin/bash
# 03_submit_doublets_bone_atlas.sh
# Single LSF job to run Scrublet on all samples

q=long
mem=450000
ncpu=2
group=team298
JOB=doublet_bulk_filtered_linnarson

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G $group \
     -o %J.${JOB}.out \
  "module load cellgen/conda; conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python 03_scrublet_doublet_detection.py \
     --merged_h5ad /lustre/scratch124/cellgen/haniffa/users/sm54/data/Linnarson/Linnarson_filtered_with_metadata.h5ad \
     --adata_status filtered_linnarson \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Linnarson/ \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/processing/logs"
