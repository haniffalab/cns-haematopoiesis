#!/bin/bash
# submit_celltypist_bulk.sh
# Single LSF job to run CellTypist on all samples

q=long
mem=300000
ncpu=8
group=cellulargenetics-priority
JOB=celltypist_bone_atlas

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G $group \
     -e %J.${JOB}.err -o %J.${JOB}.out \
  "module load cellgen/conda; conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing/;
   python 02_celltypist_annotation.py \
     --merged_h5ad /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata.h5ad \
     --model_path /nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/Human_WholeEmbryo_Public.pkl \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/ \
     --adata_status filtered_bone_atlas \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/src/processing/logs"


# Done for unfiltered already 
# bsub -J $JOB \
#      -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
#      -G $group \
#      -e %J.${JOB}.err -o %J.${JOB}.out \
#   "module load cellgen/conda; conda activate boneatlas;
#    cd /nfs/team298/sm54/BoneAtlasProject/src/processing/;
#    python 02_celltypist_annotation.py \
#      --merged_h5ad /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_unfiltered.h5ad \
#      --model_path /nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/Human_WholeEmbryo_Public.pkl \
#      --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/ \
#      --adata_status unfiltered \
#      --log_level INFO \
#      --log_dir /nfs/team298/sm54/BoneAtlasProject/src/processing/logs"
