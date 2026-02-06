INPUT_DIR="/nfs/team298/sm54/BoneAtlasProject/data/celltypist_training_data/first_trimester/"
OUTPUT_DIR="/nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/ft"

for INPUT_FILE in ${INPUT_DIR}/*.h5ad
do
    BASENAME=$(basename "$INPUT_FILE" .h5ad)
    echo "Submitting CellTypist training for ${BASENAME}"

    bsub -q normal -n 2 -M 200000 -R"select[mem>200000] rusage[mem=200000]" \
        -G cellulargenetics-priority \
        -o logs/${BASENAME}.celltypist_training.%J.out \
        "module load cellgen/conda; conda activate boneatlas; \
        python -u /nfs/team298/sm54/BoneAtlasProject/src/annotation/cell_typist_model_training.py \
        ${INPUT_FILE} ${OUTPUT_DIR}"
done

