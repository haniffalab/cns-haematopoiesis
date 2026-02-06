#!/bin/bash
# submit_generate_h5ad_bone_atlas.sh
# One job to generate all per-sample h5ads & merged h5ad

q=normal
mem=300000
ncpu=8
JOB=generate_h5ad_linnarson

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python -u 01_generate_h5ad_linnarson.py \
     --metadata /nfs/team298/sm54/BoneAtlasProject/metadata/sample_metadata/Linnarson_Metadata_Final.csv \
     --dataset Linnarson \
     --raw_dir /lustre/scratch127/cellgen/cellgeni/tickets/tic-3721/work/Cellbender/cellbender-results/ \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Linnarson \
     --merged_name Linnarson_unfiltered.h5ad \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/raw/logs"
