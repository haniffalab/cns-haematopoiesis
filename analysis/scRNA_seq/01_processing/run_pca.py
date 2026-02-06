import scanpy as sc 
import pandas as pd
pd.set_option('display.max_rows', 70)
pd.set_option('display.max_columns', 70)  # show all columns
pd.set_option('display.width', None)        # don’t wrap columns
adata= sc.read_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata_scvi.h5ad')
adata.layers['raw_counts']= adata.X.copy()
# Preprocessing
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key="10X_Chemistry")  # Use n_hvg argument

adata_hvg = adata[:, adata.var['highly_variable']].copy()
adata_hvg.X= adata_hvg.layers['raw_counts'].copy()

sc.pp.normalize_per_cell(adata_hvg, counts_per_cell_after=1e4)
sc.pp.log1p(adata_hvg)

sc.tl.pca(adata_hvg)
sc.pp.neighbors(adata_hvg)
sc.tl.umap(adata_hvg)


adata_hvg.write_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata_pca.h5ad')
