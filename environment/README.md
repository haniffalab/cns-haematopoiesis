Use this directory to store environment information

I have used the following steps 


conda create -n cns_haem python=3 scanpy scvi


conda activate cns_haem


python -m ipykernel install --user --name myenv --display-name "cns_haem"



conda env export -n myenv > environment_cns_haem.yml