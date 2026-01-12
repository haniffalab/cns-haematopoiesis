print("Importing libraries")

# import packages 
import os
import tempfile

#import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
import time

print("Reading in datasets")


#  import datasets 

#adata_ref=sc.read_h5ad("/nfs/team298/sm54/BoneAtlasProject/data/macs_liver.h5ad")

adata_ref= sc.read_h5ad('/nfs/team298/sm54/BoneAtlasProject/data/MYELOID_panfetal.h5ad')
# adata_ref.X= adata_ref.layers["raw_counts"].copy()

# adata_ref.X.max()

print("normalising")

sc.pp.normalize_total(adata_ref, target_sum = 1e4)
sc.pp.log1p(adata_ref)


print("check that per cell count is 10000")

np.expm1(adata_ref.X).sum(axis=1)

print("Training with SGD for feature selection...")




# Use `celltypist.train` to quickly train a rough CellTypist model.
# You can also set `mini_batch = True` to enable mini-batch training.
t_start = time.time()
model_fs = celltypist.train(adata_ref, 'cell.labels', n_jobs = 10, max_iter = 6, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")



#  Model training with downsampled genes 

gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -600, axis = 1)[:, -600:]
gene_index = np.unique(gene_index)


print(f"Number of genes selected: {len(gene_index)}")


print("Training with reduced features...")

# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(adata_ref[:, gene_index], 'cell.labels', check_expression = False, n_jobs = 10, max_iter = 100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")

model.write('/nfs/team298/sm54/BoneAtlasProject/metadata/celltypist_models/macs_.pkl')
