# #!/usr/bin/env python3
# # coding: utf-8

# import argparse
# import scanpy as sc
# import pandas as pd
# import helper_functions
# import numpy as np
# from scipy.stats import median_abs_deviation

# # Command-line argument parsing
# parser = argparse.ArgumentParser(description="Generate and merge h5ad files")
# parser.add_argument("--adata1", type=str, required=True, help="Path to the first AnnData file")
# parser.add_argument("--adata2", type=str, required=True, help="Path to the second AnnData file")
# parser.add_argument("--output", type=str, required=True, help="Path to the output combined h5ad file")

# args = parser.parse_args()

# # Load the datasets
# adata1 = sc.read_h5ad(args.adata1)
# adata2 = sc.read_h5ad(args.adata2)

# # Align by gene names before combining
# common_genes = adata1.var_names.intersection(adata2.var_names)
# adata1 = adata1[:, common_genes].copy()
# adata2 = adata2[:, common_genes].copy()

# adata_combined = sc.concat([adata1, adata2], join='inner', index_unique=None)

# # Ensure var_names are gene symbols
# adata_combined.var_names_make_unique()


# print(f"Combined h5ad file saved at: {args.output}")


#!/usr/bin/env python3
# coding: utf-8

import argparse
import scanpy as sc
import numpy as np

# --------------------------
# Command-line argument parsing
# --------------------------
parser = argparse.ArgumentParser(description="Check gene alignment and merge if identical")
parser.add_argument("--adata1", type=str, required=True, help="Path to the first AnnData file")
parser.add_argument("--adata2", type=str, required=True, help="Path to the second AnnData file")
parser.add_argument("--output", type=str, required=True, help="Path to the output combined h5ad file")

args = parser.parse_args()

# --------------------------
# Load metadata only (no .X) to save memory
# --------------------------
print("🔹 Reading AnnData var indices only (light mode)...")
adata1 = sc.read(args.adata1, backed="r")  # read only metadata
adata2 = sc.read(args.adata2, backed="r")

genes1 = adata1.var_names
genes2 = adata2.var_names

# --------------------------
# Check identity of gene sets and order
# --------------------------
if np.array_equal(genes1, genes2):
    print("✅ Gene names and order are identical between both AnnData files.")
    identical = True
else:
    identical = False
    n_common = len(genes1.intersection(genes2))
    print(f"⚠️ Genes differ! Common genes: {n_common} / {len(genes1)}")
    only_in_1 = len(set(genes1) - set(genes2))
    only_in_2 = len(set(genes2) - set(genes1))
    print(f"  - Only in adata1: {only_in_1}")
    print(f"  - Only in adata2: {only_in_2}")

# --------------------------
# Merge only if identical
# --------------------------
if identical:
    print("🔹 Reloading full objects for concatenation...")
    adata1_full = sc.read_h5ad(args.adata1)
    adata2_full = sc.read_h5ad(args.adata2)

    print("🔹 Merging...")
    adata_combined = sc.concat([adata1_full, adata2_full], join="inner", index_unique=None)
    adata_combined.var_names_make_unique()

    for col in adata_combined.obs.columns:   # Only convert problematic columns
        if adata_combined.obs[col].dtype == "object":# Ensure everything is a string (even numbers or lists)
            adata_combined.obs[col] = adata_combined.obs[col].astype(str)

    print("💾 Saving combined file...")
    adata_combined.write(args.output)
    print(f"✅ Combined h5ad file saved at: {args.output}")
else:
    print("❌ Aborting merge because gene sets or order differ.")

