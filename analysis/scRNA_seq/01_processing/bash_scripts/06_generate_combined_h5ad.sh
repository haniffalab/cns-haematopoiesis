#!/bin/bash
# 06_generate_combined_h5ad.sh
# One job to generate all per-sample h5ads & merged h5ad

# Set default values for variables
q=normal
mem=450000
ncpu=2
JOB=combine_all_objects

# Paths for the input files and output file
adata1="/lustre/scratch124/cellgen/haniffa/users/sm54/data/combined_adata_for_scvi.h5ad"
adata2="/lustre/scratch124/cellgen/haniffa/users/sm54/data/Linnarson/Linnarson_filtered_with_metadata.h5ad"
output_file="/lustre/scratch124/cellgen/haniffa/users/sm54/data/combined_adata_bone_cross_organ_brain_for_scvi.h5ad"

# Submit the job
bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python 06_generate_combined_h5ad.py --adata1 $adata1 --adata2 $adata2 --output $output_file"
