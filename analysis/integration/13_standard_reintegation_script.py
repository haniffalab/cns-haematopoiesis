#!/usr/bin/env python3
import os, sys, argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

# ---- Defaults (override with args if you want later) --------------------
N_LATENT = 25
N_LAYERS = 3
N_HIDDEN = 128
N_HVG    = 3000
RESOLUTION = 3.0

BATCH_KEY = "Technology"
CATEGORICAL_COVARS = ['Donor_clean','phase','10X_Chemistry']
CONTINUOUS_COVARS  = ['log1p_total_counts','pct_counts_mt']

def clear_previous(adata):
    # Remove common embeddings/graphs to avoid stale state
    for k in list(adata.obsm.keys()):
        lk = k.lower()
        if lk.startswith("x_scvi") or lk.startswith("x_umap") or lk.startswith("x_pca"):
            adata.obsm.pop(k, None)
    for k in ["distances","connectivities","X_scVI_distances","X_scVI_connectivities"]:
        if k in adata.obsp:
            adata.obsp.pop(k, None)
    for k in list(adata.uns.keys()):
        if k in {"neighbors","umap"} or "neighbors" in k or "umap" in k:
            adata.uns.pop(k, None)
    leiden_cols = [c for c in adata.obs.columns if str(c).startswith("leiden")]
    for c in leiden_cols:
        del adata.obs[c]

def main():
    ap = argparse.ArgumentParser(description="Integrate any AnnData with scVI, saving a full-gene .h5ad (no model/parquet).")
    ap.add_argument("--input", required=True, help="Path to input .h5ad")
    ap.add_argument("--out_dir", default=None, help="Output directory (default: same as input)")
    args = ap.parse_args()

    in_h5ad = args.input
    if not os.path.isfile(in_h5ad):
        raise FileNotFoundError(f"Input not found: {in_h5ad}")

    base = os.path.splitext(os.path.basename(in_h5ad))[0]
    out_dir = args.out_dir or os.path.dirname(in_h5ad)
    os.makedirs(out_dir, exist_ok=True)

    # ---- Load FULL AnnData ------------------------------------------------
    adata_full = sc.read_h5ad(in_h5ad)

    # Ensure raw counts layer exists
    if 'raw_counts' not in adata_full.layers:
        adata_full.layers['raw_counts'] = adata_full.X.copy()

    # ---- HVG selection on a temp normalized copy (keep FULL untouched) ----
    tmp = adata_full.copy()
    sc.pp.normalize_per_cell(tmp, counts_per_cell_after=1e4)
    sc.pp.log1p(tmp)

    if N_HVG > tmp.n_vars:
        raise ValueError(f"Requested n_hvg={N_HVG} > number of genes={tmp.n_vars}")

    batch_key = BATCH_KEY if BATCH_KEY in tmp.obs.columns else None
    sc.pp.highly_variable_genes(tmp, n_top_genes=N_HVG, batch_key=batch_key)
    if "highly_variable" not in tmp.var or tmp.var["highly_variable"].sum() == 0:
        raise RuntimeError("No HVGs found; check inputs/batch_key.")
    hvgs = tmp.var_names[tmp.var["highly_variable"]].tolist()

    # ---- Train scVI on HVGs only -----------------------------------------
    adata_hvg = adata_full[:, hvgs].copy()

    cat_covars = [c for c in CATEGORICAL_COVARS if c in adata_hvg.obs.columns]
    con_covars = [c for c in CONTINUOUS_COVARS  if c in adata_hvg.obs.columns]

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer='raw_counts',
        batch_key=batch_key,
        categorical_covariate_keys=cat_covars if cat_covars else None,
        continuous_covariate_keys=con_covars if con_covars else None,
    )

    model = scvi.model.SCVI(
        adata_hvg,
        n_hidden=N_HIDDEN,
        n_layers=N_LAYERS,
        n_latent=N_LATENT,
        gene_likelihood='nb',
        dispersion='gene-batch',
        use_observed_lib_size=False,
    )
    model.train()

    # ---- Latent → FULL object --------------------------------------------
    latent = model.get_latent_representation()
    if not np.array_equal(adata_hvg.obs_names.values, adata_full.obs_names.values):
        # align rows by obs_names
        df_lat = pd.DataFrame(latent, index=adata_hvg.obs_names.astype(str))
        df_lat = df_lat.reindex(adata_full.obs_names.astype(str))
        latent = df_lat.values
    clear_previous(adata_full)
    adata_full.obsm["X_scVI"] = latent

    # ---- Fresh neighbors/UMAP/Leiden on FULL using X_scVI -----------------
    sc.pp.neighbors(adata_full, use_rep='X_scVI', key_added='X_scVI')
    sc.tl.umap(adata_full, neighbors_key='X_scVI')
    sc.tl.leiden(adata_full, resolution=RESOLUTION, key_added='leiden_res_3', neighbors_key='X_scVI')

    # ---- DEGs on a normalized/log1p COPY; attach to full.uns --------------
    ad_de = adata_full.copy()
    sc.pp.normalize_per_cell(ad_de, counts_per_cell_after=1e4)
    sc.pp.log1p(ad_de)
    sc.tl.rank_genes_groups(ad_de, groupby='leiden_res_3', method='wilcoxon', key_added='DE_leiden_res_3')
    adata_full.uns['DE_leiden_res_3'] = ad_de.uns['DE_leiden_res_3']

    # ---- Save ONLY the full-gene integrated .h5ad -------------------------
    out_h5ad = os.path.join(out_dir, f"{base}.scvi_integrated_fullgenes.h5ad")
    adata_full.write_h5ad(out_h5ad)

    print(f"[OK] Saved: {out_h5ad}")
    print(f"Cells: {adata_full.n_obs} | Genes: {adata_full.n_vars}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

