#!/bin/bash
# submit_celltypist_bulk.sh
# Single LSF job to run CellTypist on all samples

q=long
mem=300000
ncpu=8
group=cellulargenetics-priority
JOB=celltypist_bulk_cross_organ_filtered

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G $group \
     -o %J.${JOB}.out \
  "module load cellgen/conda; conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing/;
   python -u 02_celltypist_annotation.py \
     --merged_h5ad /lustre/scratch124/cellgen/haniffa/users/sm54/data/Cross_Organ/cross_organ_filtered_with_metadata.h5ad \
     --model_path /nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/Human_WholeEmbryo_Public.pkl \
     --output_dir /lustre/scratch124/cellgen/haniffa/users/sm54/data/Cross_Organ/ \
     --adata_status filtered \
     --log_level INFO \
     --log_dir /nfs/team298/sm54/BoneAtlasProject/src/processing/logs"
