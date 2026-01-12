import argparse
import os
import scanpy as sc
import pandas as pd
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics.cluster import normalized_mutual_info_score as NMI
import scib

def evaluate_batch_correction(adata_uninteg, adata_int, batch_key, label_key, output_csv):
    """
    Evaluate batch correction performance on a single AnnData object.
    """
    # 1) Transfer batch_key and label_key
    adata_int.obs[batch_key] = adata_uninteg.obs[batch_key].copy()
    adata_int.obs[label_key] = adata_uninteg.obs[label_key].copy()
    print(f"Transferred {batch_key} and {label_key} to integrated AnnData.", flush=True)

    # 2) Alias SCVI‐computed graph into SCIB‐expected keys
    if 'X_scVI' in adata_int.uns:
        # Move the dict
        adata_int.uns['neighbors'] = adata_int.uns.pop('X_scVI')
        # Rename distance/connectivity matrices
        if 'X_scVI_distances' in adata_int.obsp:
            adata_int.obsp['distances'] = adata_int.obsp.pop('X_scVI_distances')
        if 'X_scVI_connectivities' in adata_int.obsp:
            adata_int.obsp['connectivities'] = adata_int.obsp.pop('X_scVI_connectivities')
        print("Aliased SCVI graph to neighbors/distances/connectivities.", flush=True)
    else:
        print("No uns['X_scVI'] found; using existing neighbors keys.", flush=True)
    # Check for NaN in labels
    if adata_int.obs[label_key].isna().any():
        print(f"Warning: {label_key} contains NaN values. Dropping these cells...")
        adata_int = adata_int[~adata_int.obs[label_key].isna()].copy()

    # 3) Run SCIB metrics
    metrics = scib.metrics.metrics(
        adata_uninteg,
        adata_int,
        batch_key=batch_key,
        label_key=label_key,
        embed='X_scVI',
        cluster_key='cluster',
        ari_=True,
        nmi_=True,
        silhouette_=True,
        graph_conn_=True,
        kBET_=True,
         hvg_score_=True,
        isolated_labels_=False,
        trajectory_=False,
        lisi_graph_=False, 
        ilisi_=False,
        clisi_=False,
        subsample=0.4,
        n_cores=4,
        type_='embed',
        verbose=True,
       
    )

    # 4) Save
    metrics.to_csv(output_csv, index=True)
    print(f"Metrics saved to: {output_csv}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate batch correction for a single-cell integration."
    )
    parser.add_argument("--adata_uninteg", type=str, required=True,
                        help="Path to unintegrated .h5ad file.")
    parser.add_argument("--adata_int", type=str, required=True,
                        help="Path to integrated .h5ad file.")
    parser.add_argument("--output_csv", type=str, required=True,
                        help="Where to save evaluation CSV.")
    parser.add_argument("--batch_key", type=str, required=True,
                        help="Batch annotation key in AnnData.obs.")
    parser.add_argument("--label_key", type=str, required=True,
                        help="Label (cell type) key in AnnData.obs.")

    args = parser.parse_args()
    print("Reading unintegrated AnnData...", flush=True)
    adata_uninteg = sc.read_h5ad(args.adata_uninteg)
    print("Unintegrated data loaded.", flush=True)

    print("Reading integrated AnnData...", flush=True)
    adata_int = sc.read_h5ad(args.adata_int)
    print("Integrated data loaded.", flush=True)

    print("Starting evaluation...", flush=True)
    evaluate_batch_correction(
        adata_uninteg=adata_uninteg,
        adata_int=adata_int,
        batch_key=args.batch_key,
        label_key=args.label_key,
        output_csv=args.output_csv
    )

