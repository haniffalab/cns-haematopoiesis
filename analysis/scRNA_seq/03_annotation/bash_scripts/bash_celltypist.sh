q=normal
mem=400000
ncpu=2


bsub -q $q -n${ncpu}  -M${mem} -R"select[mem>${mem}] rusage[mem=${mem}]" \
	-G team298 \
	 -o %J.celltypist_model_training.out \
        "module load cellgen/conda; conda activate boneatlas; python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/cell_typist_model_training_myeloid_all_organs.py"
