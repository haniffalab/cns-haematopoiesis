#!/usr/bin/env python3
# coding: utf-8

import argparse
import scanpy as sc
import pandas as pd
import helper_functions
import numpy as np
from scipy.stats import median_abs_deviation

# Command-line argument parsing
parser = argparse.ArgumentParser(description="Generate and merge h5ad files")
parser.add_argument("--adata1", type=str, required=True, help="Path to the first AnnData file")
parser.add_argument("--adata2", type=str, required=True, help="Path to the second AnnData file")
parser.add_argument("--output", type=str, required=True, help="Path to the output combined h5ad file")

args = parser.parse_args()

# Load the datasets
adata1 = sc.read_h5ad(args.adata1)
adata2 = sc.read_h5ad(args.adata2)

# Combine datasets
adata_combined = sc.concat([adata1, adata2], join='inner', index_unique=None)

# Save the combined dataset
adata_combined.write(args.output)

print(f"Combined h5ad file saved at: {args.output}")
