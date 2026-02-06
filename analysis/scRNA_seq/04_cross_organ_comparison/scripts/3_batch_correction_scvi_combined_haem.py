import scanpy as sc
import scvi
import os
import matplotlib.pyplot as plt
import argparse
import pandas as pd   # ← added

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Train SCVI model with different parameters.")
parser.add_argument('--n_latent', type=int, required=True, help="Number of latent dimensions.")
parser.add_argument('--n_layers', type=int, required=True, help="Number of layers in the SCVI model.")
parser.add_argument('--n_hidden', type=int, required=True, help="Number of hidden units per layer.")
parser.add_argument('--n_hvg', type=int, required=True, help="Number of highly variable genes.")
args = parser.parse_args()

# Load AnnData
adata = sc.read_h5ad("/nfs/team298/sm54/BoneAtlasProject/data/combined_haem_only_adata_bone_cross_organ_brain_fr_scvi.h5ad")
adata.X = adata.X.astype(float)
adata.layers['raw_counts']= adata.X.copy()

# Preprocessing
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key="Technology")
adata_hvg = adata[:, adata.var['highly_variable']].copy()

# Output directories
model_dir = "/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/trained_scvi_model"
umap_dir  = "/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/umap_plots"
data_dir  = "/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/processed_adata"
embed_dir = "/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/embedding"

os.makedirs(model_dir, exist_ok=True)
os.makedirs(umap_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)
os.makedirs(embed_dir, exist_ok=True)

# Setup SCVI
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer='raw_counts',
    batch_key='Technology',
    categorical_covariate_keys=['Donor_V3','10X_Chemistry_V2']
)

# Initialize and train SCVI model
model = scvi.model.SCVI(
    adata_hvg,
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
adata_hvg.obsm['X_scVI'] = latent

# # ← one extra line to save the scVI embedding immediately
# pd.DataFrame(latent,
#              index=adata_hvg.obs['Cell_ID'].astype(str),
#              columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]).to_parquet(
#                os.path.join(embed_dir, 
#                             f"scvi_embedding_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.parquet"
#              ),
#              engine="pyarrow",
#              compression="snappy"
# )
# Save scVI embedding using obs_names, not Cell_ID column
pd.DataFrame(
    latent,
    index=adata_hvg.obs_names,   # <---- use index, guaranteed to exist
    columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]
).to_parquet(
    os.path.join(embed_dir,
        f"scvi_embedding_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.parquet"
    ),
    engine="pyarrow",
    compression="snappy"
)
# Perform UMAP and downstream analysis
sc.pp.neighbors(adata_hvg, use_rep='X_scVI', key_added='X_scVI')
sc.tl.umap(adata_hvg, neighbors_key='X_scVI')

# Generate and save UMAP plots
for color, fname in [
    ('Technology', f"umap_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png"),
    ('published_propagated_anno_level2', f"umap_anno_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png")
]:
    plt.figure(figsize=(20, 20))
    sc.pl.umap(adata_hvg, color=color, show=False,
               title=f"UMAP ({color}) n_latent={args.n_latent}, n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}")
    plt.tight_layout()
    plt.savefig(os.path.join(umap_dir, fname), dpi=300)
    plt.close()

# Save the trained model and processed AnnData
model_filename = f"scvi_model_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}"
model.save(os.path.join(model_dir, model_filename), overwrite=True)

adata_filename = f"Best_Integration_bone_atlas_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.h5ad"
adata_hvg.write_h5ad(os.path.join(data_dir, adata_filename))

print(f"Model, embeddings, and outputs saved for n_latent={args.n_latent}, "
      f"n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}")

