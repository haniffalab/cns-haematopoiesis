#!/bin/bash
# submit_generate_h5ad_bone_atlas.sh
# One job to generate all per-sample h5ads & merged h5ad

q=normal
mem=300000
ncpu=2
JOB=qc_bone_atlas

bsub -J $JOB \
     -q $q -n $ncpu -M $mem -R "select[mem>$mem] rusage[mem=$mem]" \
     -G cellulargenetics-priority \
     -e %J.${JOB}.err -o %J.${JOB}.out \
  "module load cellgen/conda;
   conda activate boneatlas;
   cd /nfs/team298/sm54/BoneAtlasProject/src/processing;
   python 04_run_cell_cycle.py"
