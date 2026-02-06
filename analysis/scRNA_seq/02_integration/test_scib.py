#!/usr/bin/env python

import sys
import scanpy as sc
from scib.metrics.kbet import kBET

# 1) point this at one of your integrated files
fn = "/lustre/scratch126/cellgen/haniffa/sm54/BoneAtlasProject/data/bone_atlas_filtered_with_metadata_scvi.h5ad"
adata = sc.read_h5ad(fn)

# 2) take a small random subset so it’s fast
#    (here 500 cells; adjust as you like)
if adata.n_obs > 5000:
    adata = adata[adata.obs_names.to_list()[:5000]].copy()

# 3) run kBET
try:
    score = kBET(
        adata,
        batch_key="Technology",
        label_key="majority_voting_Level2",
        # if your alternative kBET has a different signature, call it here instead
    )
    print("✅ kBET ran successfully, score =", score)
except Exception as e:
    print("❌ kBET failed with:", repr(e))
    sys.exit(1)
