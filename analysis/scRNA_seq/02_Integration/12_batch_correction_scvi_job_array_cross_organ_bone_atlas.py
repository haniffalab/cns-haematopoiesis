import scanpy as sc
import scvi
import os
import matplotlib.pyplot as plt
import argparse
import pandas as pd   # ← added
import numpy as np

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Train SCVI model with different parameters.")
parser.add_argument('--n_latent', type=int, required=True, help="Number of latent dimensions.")
parser.add_argument('--n_layers', type=int, required=True, help="Number of layers in the SCVI model.")
parser.add_argument('--n_hidden', type=int, required=True, help="Number of hidden units per layer.")
parser.add_argument('--n_hvg', type=int, required=True, help="Number of highly variable genes.")
args = parser.parse_args()

# Load AnnData
adata = sc.read_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/combined_adata_for_scvi.h5ad')
adata.layers['raw_counts']= adata.X.copy()

if np.issubdtype(adata.X.dtype, np.integer):
    print("adata.X is raw counts (integer matrix)")
else:
    print("adata.X is already float (likely normalized)")

adata.X = adata.X.astype(np.float64)

# Preprocessing
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key="Technology")
adata = adata[:, adata.var['highly_variable']].copy()

# Output directories
model_dir = "/lustre/scratch124/cellgen/haniffa/users/sm54/data/bone_atlas_cross_organ_scvi_runs/trained_scvi_models"
umap_dir  = "/lustre/scratch124/cellgen/haniffa/users/sm54/data/bone_atlas_cross_organ_scvi_runs/umap_plots"
data_dir  = "/lustre/scratch124/cellgen/haniffa/users/sm54/data/bone_atlas_cross_organ_scvi_runs/processed_adata"
embed_dir = "/lustre/scratch124/cellgen/haniffa/users/sm54/data/bone_atlas_cross_organ_scvi_runs/embed_parquet"

os.makedirs(model_dir, exist_ok=True)
os.makedirs(umap_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)
os.makedirs(embed_dir, exist_ok=True)

# Setup SCVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer='raw_counts',
    batch_key='Technology',
    categorical_covariate_keys=['Donor_ID','10X_Chemistry']
)

# Initialize and train SCVI model
model = scvi.model.SCVI(
    adata,
    n_hidden=args.n_hidden,
    n_layers=args.n_layers,
    n_latent=args.n_latent,
    gene_likelihood='nb',
    dispersion='gene-batch',
    use_observed_lib_size=False,
)
model.train()

# Get latent representation & attach
latent = model.get_latent_representation()
adata.obsm['X_scVI'] = latent

# ← one extra line to save the scVI embedding immediately
pd.DataFrame(latent,
             index=adata.obs['Cell_ID'].astype(str),
             columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]).to_parquet(
               os.path.join(embed_dir, 
                            f"scvi_embedding_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.parquet"
             ),
             engine="pyarrow",
             compression="snappy"
)

# Perform UMAP and downstream analysis
sc.pp.neighbors(adata, use_rep='X_scVI', key_added='X_scVI')
sc.tl.umap(adata, neighbors_key='X_scVI')

# Generate and save UMAP plots
for color, fname in [
    ('Technology', f"umap_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png"),
    ('majority_voting_Level1', f"umap_anno_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png")
]:
    plt.figure(figsize=(20, 20))
    sc.pl.umap(adata, color=color, show=False,
               title=f"UMAP ({color}) n_latent={args.n_latent}, n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}")
    plt.tight_layout()
    plt.savefig(os.path.join(umap_dir, fname), dpi=300)
    plt.close()

# Save the trained model and processed AnnData
model_filename = f"scvi_model_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}"
model.save(os.path.join(model_dir, model_filename), overwrite=True)

adata_filename = f"adata_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.h5ad"
adata.write_h5ad(os.path.join(data_dir, adata_filename))

print(f"Model, embeddings, and outputs saved for n_latent={args.n_latent}, "
      f"n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}")

