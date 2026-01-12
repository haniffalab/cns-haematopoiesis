#!/usr/bin/env python
import os
import argparse
import numpy as np

import scanpy as sc
import celltypist
from celltypist import models


def main():
    parser = argparse.ArgumentParser(
        description="Annotate single-cell data with a custom CellTypist model"
    )
    parser.add_argument(
        "--input_h5ad", required=True,
        help="Path to input AnnData (.h5ad) file"
    )
    parser.add_argument(
        "--model_pkl", required=True,
        help="Path to the CellTypist model (.pkl) file"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Directory to write prediction outputs"
    )
    parser.add_argument(
        "--no_majority_voting", action="store_false", dest="majority_voting",
        help="Disable majority voting (default: enabled)"
    )
    args = parser.parse_args()

    # 1) Load your custom CellTypist model
    print(f"Loading model from {args.model_pkl}...")
    model = models.Model.load(args.model_pkl)

    # 2) Read the integrated AnnData
    print(f"Reading data from {args.input_h5ad}...")
    adata = sc.read_h5ad(args.input_h5ad)

    # 3) Restore raw counts and normalize
    if 'raw_counts' in adata.layers:
        adata.X = adata.layers['raw_counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 4) CHECK: after normalization, counts should sum to ~10000 per cell
    sums = np.asarray(adata.X.expm1().sum(axis=1)).flatten()
    if not np.allclose(sums, 1e4, atol=1e-2):
        raise ValueError(
            "Post-normalization counts do *not* sum to 1e4 for all cells:\n"
            f" min={sums.min():.2f}, max={sums.max():.2f}"
        )
    print("✔ Normalization check passed: all cells sum to 1e4.")

    # 5) Hook in the scVI embedding for CellTypist graph lookup
    #    so that celltypist.label will use your 'X_scVI' neighbors
    if 'X_scVI' in adata.obsm:
        adata.uns['neighbors'] = {}  # reset default
        adata.uns['neighbors']['params'] = {}
        adata.obsp['connectivities'] = adata.obsp['X_scVI_connectivities'].copy()
        adata.obsp['distances']     = adata.obsp['X_scVI_distances'].copy()
        print("✔ scVI embedding graph registered for CellTypist.")
    else:
        print("⚠️  No X_scVI embedding found in .obsm – CellTypist will compute its own neighbors.")

    # 6) Annotate with CellTypist
    print("Running CellTypist annotation...")
    result = celltypist.annotate(
        adata,
        model=model,
        majority_voting=args.majority_voting
    )

    # 7) Prepare output paths
    os.makedirs(args.output_dir, exist_ok=True)
    sample_name = os.path.splitext(os.path.basename(args.input_h5ad))[0]
    pred_csv     = os.path.join(args.output_dir, f"{sample_name}_predicted_labels.csv")
    prob_csv     = os.path.join(args.output_dir, f"{sample_name}_probability_matrix.csv")
    decision_csv = os.path.join(
        args.output_dir,
        f"{sample_name}_decision_matrix_{os.path.basename(args.model_pkl).replace('.pkl','')}.csv"
    )

    # 8) Save outputs
    print(f"Saving predicted labels to {pred_csv}")
    result.predicted_labels.to_csv(pred_csv)
    print(f"Saving probability matrix to {prob_csv}")
    result.probability_matrix.to_csv(prob_csv)
    print(f"Saving decision matrix to {decision_csv}")
    result.decision_matrix.to_csv(decision_csv)

    print("✅ All done.")


if __name__ == "__main__":
    main()

