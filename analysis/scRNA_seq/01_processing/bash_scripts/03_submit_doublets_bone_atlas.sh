#!/bin/bash
# 03_submit_doublets_bone_atlas.sh
# Single LSF job to run Scrublet on all samples

q=long
mem=250000
ncpu=2
group=team298
JOB=doublet_bulk_bone_atlas

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G $group \
     -e %J.${JOB}.err -o %J.${JOB}.out \
  "module load cellgen/conda; conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python 03_scrublet_doublet_detection.py \
     --merged_h5ad /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata.h5ad \
     --adata_status filtered_bone_atlas \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/ \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/processing/logs"
