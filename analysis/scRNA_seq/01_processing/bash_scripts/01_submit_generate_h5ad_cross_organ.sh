#!/bin/bash
# submit_generate_h5ad_bone_atlas.sh
# One job to generate all per-sample h5ads & merged h5ad

q=normal
mem=300000
ncpu=8
JOB=generate_h5ad_crossorgan

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -e %J.${JOB}.err -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python 01_generate_h5ad.py \
     --metadata /nfs/team298/sm54/BoneAtlasProject/metadata/sample_metadata/COMBINED_BONE_CROSS_ORGAN_METADATA.csv \
     --dataset Cross_Organ \
     --raw_dir /nfs/cellgeni/tickets/tic-3721/Cellbender/cellbender-results/ \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Cross_Organ/ \
     --merged_name cross_organ_unfiltered.h5ad \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/raw/logs"
