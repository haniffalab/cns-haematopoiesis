
All the commands and files for creating and maintaining your conda env.

conda env create -f environment.yml

# create a kernel 
module load cellgen/conda
conda activate /software/cellgen/team298/sm54/envs/boneatlas

conda install ipykernel
python -m ipykernel install \
  --user \
  --name boneatlas \
  --display-name "boneatlas"
  

saved at location:  /software/cellgen/team298/sm54/envs/boneatlas