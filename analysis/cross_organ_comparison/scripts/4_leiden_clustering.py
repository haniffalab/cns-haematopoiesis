#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import sys

# ---- Load data ----
adata= sc.read_h5ad('/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/processed_adata/Best_Integration_bone_atlas_nlatent20_nlayers2_nhidden128_nhvg3000.h5ad')

sc.tl.leiden(adata, resolution= 3, key_added='leiden_res_3',neighbors_key='X_scVI')
sc.tl.rank_genes_groups(adata, groupby='leiden_res_3', method='wilcoxon', corr_method='benjamini-hochberg', groups='all', reference='rest', n_genes=500, use_raw=False, log_transformed=True)


adata.write_h5ad(
'/nfs/team298/sm54/BoneAtlasProject/data/cross_organ_haem/best_integration/processed_adata/Best_Integration_bone_atlas_nlatent20_nlayers2_nhidden128_nhvg3000_with_leiden_clustering.h5ad'
)
print("✅ Finished successfully — leiden clustering done.")
