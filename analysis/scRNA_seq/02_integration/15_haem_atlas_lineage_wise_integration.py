#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scvi

# --------------------------
# Parse command-line arguments
# --------------------------
parser = argparse.ArgumentParser(description="Train SCVI for a single immune dataset (HVGs for training; save results to FULL AnnData).")
parser.add_argument('--h5ad', type=str, required=True, help="Path to one input .h5ad file.")
parser.add_argument('--n_latent', type=int, required=True, help="Number of latent dimensions.")
parser.add_argument('--n_layers', type=int, required=True, help="Number of layers in the SCVI model.")
parser.add_argument('--n_hidden', type=int, required=True, help="Number of hidden units per layer.")
parser.add_argument('--n_hvg', type=int, required=True, help="Number of highly variable genes.")
parser.add_argument('--data_dir', type=str, default="/nfs/team298/sm54/BoneAtlasProject/data/bone_atlas_anndatas_immune_subsets",
                    help="Base directory for outputs (subfolders per cell type).")
args = parser.parse_args()

# --------------------------
# Resolve paths, names, folders
# --------------------------
adata_path = args.h5ad
if not os.path.isfile(adata_path):
    raise FileNotFoundError(f"Input h5ad not found: {adata_path}")

cell_type = os.path.basename(adata_path).replace("bone_", "").replace("_filtered_with_metadata_scvi.h5ad", "")
print(f"\nProcessing cell type: {cell_type}")

out_base = os.path.join(args.data_dir, cell_type)
model_dir = os.path.join(out_base, "trained_scvi_model")
umap_dir  = os.path.join(out_base, "umap_plots")
data_dir  = os.path.join(out_base, "processed_adata")
embed_dir = os.path.join(out_base, "embedding")
for d in [model_dir, umap_dir, data_dir, embed_dir]:
    os.makedirs(d, exist_ok=True)

# --------------------------
# Load FULL AnnData
# --------------------------
adata_full = sc.read_h5ad(adata_path)

# Ensure raw counts layer exists (do NOT overwrite X)
if 'raw_counts' not in adata_full.layers:
    adata_full.layers['raw_counts'] = adata_full.X.copy()

# --------------------------
# HVG selection on a temporary normalized copy (keep FULL object untouched)
# --------------------------
tmp = adata_full.copy()
sc.pp.normalize_per_cell(tmp, counts_per_cell_after=1e4)
sc.pp.log1p(tmp)

batch_key = "Technology" if "Technology" in tmp.obs.columns else None
if args.n_hvg > tmp.n_vars:
    raise ValueError(f"Requested n_hvg={args.n_hvg} > number of genes={tmp.n_vars}")

sc.pp.highly_variable_genes(tmp, n_top_genes=args.n_hvg, batch_key=batch_key)
if "highly_variable" not in tmp.var or tmp.var["highly_variable"].sum() == 0:
    raise RuntimeError("No HVGs found. Check inputs/batch_key.")

hvgs = tmp.var_names[tmp.var["highly_variable"]].tolist()

# --------------------------
# Subset FULL object to HVGs for training only
# --------------------------
adata_hvg = adata_full[:, hvgs].copy()

# --------------------------
# Setup SCVI on HVG view
# --------------------------
cat_covars_all = ['Donor_clean', 'phase', '10X_Chemistry']
cat_covars = [c for c in cat_covars_all if c in adata_hvg.obs.columns]
for c in cat_covars:
    if not pd.api.types.is_categorical_dtype(adata_hvg.obs[c]):
        adata_hvg.obs[c] = adata_hvg.obs[c].astype('category')

setup_kwargs = dict(adata=adata_hvg, layer='raw_counts')
if batch_key is not None:
    setup_kwargs['batch_key'] = batch_key
if cat_covars:
    setup_kwargs['categorical_covariate_keys'] = cat_covars

scvi.model.SCVI.setup_anndata(**setup_kwargs)

# --------------------------
# Train SCVI
# --------------------------
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

# --------------------------
# Get latent; store in FULL; save embedding parquet
# --------------------------
latent = model.get_latent_representation()
if latent.shape[0] != adata_full.n_obs:
    raise RuntimeError("Latent rows != n_obs; obs order mismatch.")

adata_full.obsm["X_scVI"] = latent

pd.DataFrame(
    latent,
    index=adata_full.obs_names.astype(str),
    columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]
).to_parquet(
    os.path.join(
        embed_dir,
        f"scvi_embedding_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.parquet"
    ),
    engine="pyarrow",
    compression="snappy"
)

# --------------------------
# Graph / UMAP / Leiden on FULL using X_scVI
# --------------------------
sc.pp.neighbors(adata_full, use_rep='X_scVI', key_added='X_scVI')
sc.tl.umap(adata_full, neighbors_key='X_scVI')
sc.tl.leiden(adata_full, resolution=3, key_added='leiden_res_3', neighbors_key='X_scVI')

# --------------------------
# DEG on a normalized/log1p COPY; store in full.uns
# --------------------------
ad_de = adata_full.copy()
sc.pp.normalize_per_cell(ad_de, counts_per_cell_after=1e4)
sc.pp.log1p(ad_de)
sc.tl.rank_genes_groups(ad_de, groupby='leiden_res_3', method='wilcoxon', key_added='DE_leiden_res_3')
adata_full.uns['DE_leiden_res_3'] = ad_de.uns['DE_leiden_res_3']

# --------------------------
# UMAP plots from FULL object
# --------------------------
for color, fname in [
    ('Technology', f"umap_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png"),
    ('Haem_Manual_Annotation_Level1_V5', f"umap_anno_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png")
]:
    if color in adata_full.obs.columns:
        plt.figure(figsize=(20, 20))
        sc.pl.umap(
            adata_full,
            color=color,
            show=False,
            title=f"UMAP ({color}) {cell_type} n_latent={args.n_latent}, n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}",
            neighbors_key='X_scVI'
        )
        plt.tight_layout()
        plt.savefig(os.path.join(umap_dir, fname), dpi=300)
        plt.close()

# --------------------------
# Save model and FULL AnnData
# --------------------------
model_filename = f"scvi_model_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}"
model.save(os.path.join(model_dir, model_filename), overwrite=True)

adata_filename = f"Best_Integration_{cell_type}_FULL_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.h5ad"
adata_full.write_h5ad(os.path.join(data_dir, adata_filename))

print(f"[DONE] Completed SCVI integration for {cell_type}. Saved: {adata_filename}")
