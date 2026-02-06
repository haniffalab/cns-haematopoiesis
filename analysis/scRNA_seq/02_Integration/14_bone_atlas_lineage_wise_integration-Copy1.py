import scanpy as sc
import scvi
import os
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import glob

# --------------------------
# Parse command-line arguments
# --------------------------
parser = argparse.ArgumentParser(description="Train SCVI model for multiple non-immune datasets.")
parser.add_argument('--n_latent', type=int, required=True, help="Number of latent dimensions.")
parser.add_argument('--n_layers', type=int, required=True, help="Number of layers in the SCVI model.")
parser.add_argument('--n_hidden', type=int, required=True, help="Number of hidden units per layer.")
parser.add_argument('--n_hvg', type=int, required=True, help="Number of highly variable genes.")
parser.add_argument('--data_dir', type=str, default="/nfs/team298/sm54/BoneAtlasProject/data/bone_atlas_anndatas_non_immune",
                    help="Directory containing non-immune h5ads")

args = parser.parse_args()

# --------------------------
# Find all non-immune h5ads
# --------------------------
adata_files = glob.glob(os.path.join(args.data_dir, "*.h5ad"))
if not adata_files:
    raise FileNotFoundError(f"No .h5ad files found in {args.data_dir}")

# --------------------------
# Loop over datasets
# --------------------------
for adata_path in adata_files:
    # Extract cell type name from filename
    cell_type = os.path.basename(adata_path).replace("bone_", "").replace("_filtered_with_metadata_scvi.h5ad", "")
    print(f"\nProcessing cell type: {cell_type}")

    # Create output subfolders for this cell type
    out_base = os.path.join(args.data_dir, cell_type)
    model_dir = os.path.join(out_base, "trained_scvi_model")
    umap_dir  = os.path.join(out_base, "umap_plots")
    data_dir  = os.path.join(out_base, "processed_adata")
    embed_dir = os.path.join(out_base, "embedding")
    for d in [model_dir, umap_dir, data_dir, embed_dir]:
        os.makedirs(d, exist_ok=True)

    # --------------------------
    # Load AnnData
    # --------------------------
    adata = sc.read_h5ad(adata_path)
    if 'raw_counts' not in adata.layers:
        adata.layers['raw_counts'] = adata.X.copy()

    # --------------------------
    # Preprocessing
    # --------------------------
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key="Technology")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # --------------------------
    # Setup SCVI
    # --------------------------
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer='raw_counts',
        batch_key='Technology',
        categorical_covariate_keys=['Donor_clean','phase','10X_Chemistry']
    )

    # --------------------------
    # Initialize and train SCVI model
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
    # Get latent representation & save as parquet
    # --------------------------
    latent = model.get_latent_representation()
    adata_hvg.obsm['X_scVI'] = latent

    pd.DataFrame(
        latent,
        index=adata_hvg.obs_names.astype(str),
        columns=[f"SCVI{i+1}" for i in range(latent.shape[1])]
    ).to_parquet(
        os.path.join(embed_dir,
                     f"scvi_embedding_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.parquet"),
        engine="pyarrow",
        compression="snappy"
    )

    # --------------------------
    # UMAP & neighbors
    # --------------------------
    sc.pp.neighbors(adata_hvg, use_rep='X_scVI', key_added='X_scVI')
    sc.tl.umap(adata_hvg, neighbors_key='X_scVI')

    # --------------------------
    # Generate UMAP plots
    # --------------------------
    for color, fname in [
        ('Technology', f"umap_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png"),
        ('major_celltype_manual', f"umap_anno_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.png")
    ]:
        plt.figure(figsize=(20, 20))
        sc.pl.umap(
            adata_hvg,
            color=color,
            show=False,
            title=f"UMAP ({color}) {cell_type} n_latent={args.n_latent}, n_layers={args.n_layers}, n_hidden={args.n_hidden}, n_hvg={args.n_hvg}"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(umap_dir, fname), dpi=300)
        plt.close()

    # --------------------------
    # Save model and processed AnnData
    # --------------------------
    model_filename = f"scvi_model_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}"
    model.save(os.path.join(model_dir, model_filename), overwrite=True)

    adata_filename = f"Best_Integration_{cell_type}_nlatent{args.n_latent}_nlayers{args.n_layers}_nhidden{args.n_hidden}_nhvg{args.n_hvg}.h5ad"
    adata_hvg.write_h5ad(os.path.join(data_dir, adata_filename))

    print(f"Completed SCVI integration for {cell_type}.\n")
