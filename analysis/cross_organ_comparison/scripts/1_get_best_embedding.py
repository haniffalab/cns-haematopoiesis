#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import sys
from pathlib import Path

# ---- Setup paths ----
sys.path.append("/nfs/team298/sm54/BoneAtlasProject/src/annotation")
from helper_functions import save_obs, restore_obs

repo_root = Path("/nfs/team298/sm54/BoneAtlasProject")
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from metadata.marker_genes.marker_dict import (
    marker_dict,
    CATEGORY_MARKERS_fbm,
    celltypist_lymphoid,
    celltypist_progenitors,
    s_genes,
    g2m_genes,
    T_Cells,
    T_cell_type_genes,
    tcell_diff_markers,
)

# ---- Load data ----
# adata = sc.read_h5ad(
#     "/lustre/scratch124/cellgen/haniffa/users/sm54/data/combined_adata_bone_cross_organ_brain_for_scvi.h5ad"
# )
adata = sc.read_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/combined_adata_bone_cross_organ_brain_for_scvi_with_gene_names.h5ad')

adata_int = sc.read_h5ad(
    "/nfs/team298/sm54/BoneAtlasProject/src/cross_organ_comparison/scib_results/adata_embedding.h5ad"
)

# ---- Find common cells ----
adata.obs["Cell_ID"] = adata.obs["Cell_ID"].astype(str)
adata_int.obs_names = adata_int.obs_names.astype(str)

common_cells = pd.Index(adata.obs["Cell_ID"]).intersection(adata_int.obs_names)
print(f"✅ Found {len(common_cells)} overlapping cells.")

# ---- Create index map (so embeddings go to right cells) ----
cell_to_index = pd.Series(adata.obs.index, index=adata.obs["Cell_ID"]).to_dict()

# ---- Transfer embeddings directly ----
keys_to_transfer = [
    "scvi_params_2_hvg_3000",
    "umap_PCA_params_unintegrated_hvg_3000",
    "umap_scvi_params_2_hvg_3000",
]

for key in keys_to_transfer:
    if key in adata_int.obsm:
        emb = adata_int.obsm[key]
        print(f"Processing {key} ({emb.shape})...")

        # Match order between adata_int.obs_names and adata.obs["Cell_ID"]
        common_mask = adata_int.obs_names.isin(common_cells)
        emb_common = emb[common_mask, :]

        adata.obsm[key] = np.zeros((adata.n_obs, emb.shape[1]))
        adata.obsm[key][:] = np.nan  # fill with NaN to keep dimensions consistent

        idx = adata.obs.index[adata.obs["Cell_ID"].isin(common_cells)]
        adata.obsm[key][adata.obs["Cell_ID"].isin(common_cells)] = emb_common
        print(f"✅ Transferred {key} for {len(common_cells)} cells.")
    else:
        print(f"⚠️ Skipped {key}: not found in adata_int")

# ---- Save raw-count version ----
adata.write_h5ad(
    "/nfs/team298/sm54/BoneAtlasProject/data/combined_adata_bone_cross_organ_brain_with_integrated_embedding_raw_counts_and_gene_names.h5ad"
)
print("💾 Saved raw counts version.")

# ---- Normalise + log-transform ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write_h5ad(
    "/nfs/team298/sm54/BoneAtlasProject/data/combined_adata_bone_cross_organ_brain_with_integrated_embedding_log_norm_counts.h5ad"
)
print("✅ Finished successfully — both raw and log-normalised versions saved.")
