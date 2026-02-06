q=long
mem=200000
ncpu=4


bsub -q $q -n${ncpu}  -M${mem} -R"select[mem>${mem}] rusage[mem=${mem}]" \
	-G cellulargenetics-priority \
	 -o %J.Bone_Atlaslabel_propagation.out \
        "module load cellgen/conda; conda activate CellRank; python /nfs/team298/sm54/BoneAtlasProject/src/annotation/Published_Anno_Label_Propagation.py"
