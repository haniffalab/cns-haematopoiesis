#!/bin/bash
# submit_generate_h5ad_bone_atlas.sh
# One job to generate all per-sample h5ads & merged h5ad
# runs non-integrated basic pca analysis 
q=normal
mem=300000
ncpu=2
JOB=run_pca_bone_atlas

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python run_pca.py"
