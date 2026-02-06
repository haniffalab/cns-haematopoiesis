#!/usr/bin/env python3
"""
Train an SCVI model on a single .h5ad file using HVGs for training only.
The latent embedding is written back to the full object, with neighbors/UMAP/Leiden
and DE results computed and saved.

Usage:
    python train_scvi_single.py \
        --h5ad input.h5ad \
        --out_dir output_directory \
        --n_latent 20 --n_layers 2 --n_hidden 128 --n_hvg 5000
"""

import os
import argparse
import pandas as pd
import scanpy as sc
import scvi
import matplotlib.pyplot as plt

# --------------------------
# Parse arguments
# --------------------------
parser = argparse.ArgumentParser(description="Train SCVI on one dataset (HVGs for training; results saved to full AnnData).")
parser.add_argument('--h5ad', type=str, required=True, help="Path to input .h5ad file.")
parser.add_argument('--out_dir', type=str, required=True, help="Directory to save all outputs.")
parser.add_argument('--n_latent', type=int, required=True)
parser.add_argument('--n_layers', type=int, required=True)
parser.add_argument('--n_hidden', type=int, required=True)
parser.add_argument('--n_hvg', type=int, required=True)
args = parser.parse_args()

# --------------------------
# Fixed keys
# --------------------------
BATCH_KEY = "Technology"
CATEGORICAL_COVARS = ['Donor_clean', '10X_Chemistry']
CONTINUOUS_COVARS = ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt"]

# --------------------------
# Load AnnData
# --------------------------
adata = sc.read_h5ad(args.h5ad)
os.makedirs(args.out_dir, exist_ok=True)
print(f"\n[INFO] Loaded {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Ensure raw counts layer exists
if 'raw_counts' not in adata.layers:
    adata.layers['raw_counts'] = adata.X.copy()

# --------------------------
# HVG selection
# --------------------------
tmp = adata.copy()
sc.pp.normalize_per_cell(tmp, counts_per_cell_after=1e4)
sc.pp.log1p(tmp)

if args.n_hvg > tmp.n_vars:
    raise ValueError(f"Requested {args.n_hvg} HVGs > number of genes ({tmp.n_vars}).")

sc.pp.highly_variable_genes(
    tmp,
    n_top_genes=args.n_hvg,
    batch_key=BATCH_KEY if BATCH_KEY in tmp.obs else None
)
hvgs = tmp.var_names[tmp.var['highly_variable']].tolist()
print(f"[INFO] Selected {len(hvgs):,} highly variable genes.")

adata_hvg = adata[:, hvgs].copy()

# --------------------------
# Setup SCVI
# --------------------------
cat_covars = [c for c in CATEGORICAL_COVARS if c in adata_hvg.obs.columns]
con_covars = [c for c in CONTINUOUS_COVARS if c in adata_hvg.obs.columns]
for c in cat_covars:
    if not pd.api.types.is_categorical_dtype(adata_hvg.obs[c]):
        adata_hvg.obs[c] = adata_hvg.obs[c].astype('category')

batch_key = BATCH_KEY if BATCH_KEY in adata_hvg.obs.columns else None

print("\n[INFO] SCVI setup summary:")
print(f"  Batch key: {batch_key}")
print(f"  Categorical covariates: {cat_covars if cat_covars else 'None'}")
print(f"  Continuous covariates: {con_covars if con_covars else 'None'}")

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer='raw_counts',
    batch_key=batch_key,
    categorical_covariate_keys=cat_covars if cat_covars else None,
    continuous_covariate_keys=con_covars if con_covars else None,
)

# --------------------------
# Train model
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
# Save latent embedding
# --------------------------
latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

latent_df = pd.DataFrame(
    latent,
    index=adata.obs_names,
    columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]
)
latent_path = os.path.join(args.out_dir, "scvi_embedding.parquet")
latent_df.to_parquet(latent_path, compression='snappy')
print(f"[SAVED] Latent embedding → {latent_path}")

# --------------------------
# Graph, UMAP, Leiden
# --------------------------
neighbors_key = "X_scVI"

# Remove any old neighbor graphs with same key
to_remove = [f"{neighbors_key}_connectivities", f"{neighbors_key}_distances"]
for key in to_remove:
    if key in adata.obsp:
        del adata.obsp[key]
        print(f"[INFO] Removed existing obsp['{key}']")

if neighbors_key in adata.uns:
    del adata.uns[neighbors_key]
    print(f"[INFO] Removed existing uns['{neighbors_key}']")

# Compute neighbors, UMAP, Leiden
sc.pp.neighbors(adata, use_rep='X_scVI', key_added=neighbors_key)
sc.tl.umap(adata, neighbors_key=neighbors_key)
sc.tl.leiden(adata, resolution=2, key_added='leiden_res_2', neighbors_key=neighbors_key)
sc.tl.leiden(adata, resolution=1, key_added='leiden_res_1', neighbors_key=neighbors_key)

sc.tl.leiden(adata, resolution=3, key_added='leiden_res_3', neighbors_key=neighbors_key)
# --------------------------
# Differential expression
# --------------------------
ad_de = adata.copy()
sc.pp.normalize_per_cell(ad_de, counts_per_cell_after=1e4)
sc.pp.log1p(ad_de)
sc.tl.rank_genes_groups(ad_de, groupby='leiden_res_1', method='wilcoxon')
adata.uns['leiden_res_1'] = ad_de.uns['rank_genes_groups']

# --------------------------
# Save UMAP plots (legend always shown in right margin)
# --------------------------
def safe_umap_plot(adata, color, out_path):
    """Save UMAP with legend always visible on right margin."""
    if color not in adata.obs:
        print(f"[WARN] '{color}' not found in obs; skipping.")
        return
    plt.figure(figsize=(10, 10))
    sc.pl.umap(
        adata,
        color=color,
        show=False,
        legend_loc='right margin',
        title=f"UMAP – {color}",
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SAVED] {out_path}")

for color in ["Technology", "Haem_Manual_Annotation_Level1_V7"]:
    out_path = os.path.join(args.out_dir, f"umap_{color}.png")
    safe_umap_plot(adata, color, out_path)

# --------------------------
# Save model and AnnData
# --------------------------
model.save(os.path.join(args.out_dir, "trained_scvi_model"), overwrite=True)
adata.write_h5ad(os.path.join(args.out_dir, "FULL_scvi.h5ad"))

print(f"\n[DONE] SCVI integration complete. Results saved to: {args.out_dir}")
