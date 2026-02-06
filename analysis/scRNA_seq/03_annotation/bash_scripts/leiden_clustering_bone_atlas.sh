#!/bin/bash

q=normal
mem=350000
ncpu=4


bsub -q $q -n${ncpu}  -M${mem} -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" \
	-G team298 \
	 -o %J.leiden_clustering_bone_atlas.out \
"module load cellgen/conda; conda activate boneatlas; python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/leiden_clustering_bone_atlas.py"