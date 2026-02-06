#!/bin/bash
# submit_generate_h5ad_bone_atlas.sh
# One job to generate all per-sample h5ads & merged h5ad

q=normal
mem=650000
ncpu=8
JOB=cross_organ_qc_filtering

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python -u 05_basic_qc_filtering_add_metadata.py \
    --adata-file /lustre/scratch124/cellgen/haniffa/users/sm54/data/Cross_Organ/cross_organ_unfiltered.h5ad \
     --metadata-file /nfs/team298/sm54/BoneAtlasProject/metadata/sample_metadata/COMBINED_BONE_CROSS_ORGAN_METADATA.csv \
     --output-file /lustre/scratch124/cellgen/haniffa/users/sm54/data/Cross_Organ/cross_organ_filtered_with_metadata.h5ad "
