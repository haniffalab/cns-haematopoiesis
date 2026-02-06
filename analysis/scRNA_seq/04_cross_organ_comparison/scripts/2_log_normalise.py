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
adata = sc.read_h5ad("/nfs/team298/sm54/BoneAtlasProject/data/combined_adata_bone_cross_organ_brain_with_integrated_embedding_raw_counts.h5ad")

# ---- Normalise + log-transform ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write_h5ad(
    "/nfs/team298/sm54/BoneAtlasProject/data/combined_adata_bone_cross_organ_brain_with_integrated_embedding_log_norm_counts.h5ad"
)
print("✅ Finished successfully — both raw and log-normalised versions saved.")
