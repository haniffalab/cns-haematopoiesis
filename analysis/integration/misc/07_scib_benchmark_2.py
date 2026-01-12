#!/usr/bin/env python3
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import anndata as ad
import faiss
from scib_metrics.nearest_neighbors import NeighborsResults
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection


def faiss_brute_force_nn(X: np.ndarray, k: int) -> NeighborsResults:
    """
    GPU‐accelerated brute‐force nearest neighbors using FAISS.
    Returns distances (L2) and indices as NeighborsResults.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.IndexFlatL2(X.shape[1])
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    # distances returned are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run full scIB large‐scale benchmark on a scVI‐integrated AnnData"
    )
    p.add_argument("-i", "--input",    required=True,
                   help="scVI‐integrated AnnData (.h5ad) with X_scVI in .obsm")
    p.add_argument("-p", "--pca",      required=True,
                   help="Pre‐integration PCA in Parquet (index=Cell_ID)")
    p.add_argument("-o", "--output",   required=True,
                   help="CSV path for final scIB results table")
    p.add_argument("--batch-key", default="sample",
                   help="obs column for batch labels")
    p.add_argument("--label-key", default="cell_type",
                   help="obs column for cell‐type labels")
    p.add_argument("--neighbors", type=int, default=15,
                   help="k for nearest‐neighbor metrics")
    p.add_argument("--jobs", type=int, default=1,
                   help="Number of CPU cores for parallel metrics")
    return p.parse_args()


def main():
    args = parse_args()

    # 1) load AnnData
    print(f"[1/6] Loading integrated AnnData from {args.input}", file=sys.stderr)
    adata = ad.read_h5ad(args.input)

    # 2) read and align PCA
    print(f"[2/6] Reading PCA from {args.pca}", file=sys.stderr)
    pca_df = pd.read_parquet(args.pca, engine="pyarrow").astype(float)

    print("[3/6] Aligning cells to PCA…", file=sys.stderr)
    ids = adata.obs["Cell_ID"].astype(str)
    pos = pca_df.index.get_indexer(ids)

    mask = pos >= 0
    n_missing = (~mask).sum()
    if n_missing:
        print(f"  → Dropping {n_missing} cells (no PCA)", file=sys.stderr)
        adata._inplace_subset_obs(mask)
        pos = pos[mask]
        ids = ids[mask]

    aligned = pca_df.values[pos, :]
    if np.isnan(aligned).any():
        bad = ids[np.isnan(aligned).any(axis=1)].unique().tolist()[:5]
        raise ValueError(f"PCA has NaNs for some cells; examples: {bad}")

    adata.obsm["X_unintegrated_pca"] = aligned

    # 3) prepare embeddings for scIB
    # copy into the keys that Benchmarker expects
    adata.obsm["Unintegrated"] = adata.obsm["X_unintegrated_pca"]
    if "X_scVI" not in adata.obsm:
        raise KeyError("adata.obsm['X_scVI'] not found")
    adata.obsm["scVI"] = adata.obsm["X_scVI"]

    # 4) instantiate Benchmarker
    print("[4/6] Setting up Benchmarker…", file=sys.stderr)
    biocons = BioConservation(isolated_labels=False)
    bm = Benchmarker(
        adata,
        batch_key=args.batch_key,
        label_key=args.label_key,
        embedding_obsm_keys=["Unintegrated", "scVI"],
        pre_integrated_embedding_obsm_key="X_unintegrated_pca",
        bio_conservation_metrics=biocons,
        batch_correction_metrics=BatchCorrection(),
        n_jobs=args.jobs,
    )

    # 5) run metrics
    print(f"[5/6] Computing neighbors (k={args.neighbors}) & benchmarking…", file=sys.stderr)
    bm.prepare(neighbor_computer=lambda X: faiss_brute_force_nn(X, args.neighbors))
    start = time.time()
    bm.benchmark()
    end = time.time()
    print(f"  → Benchmark completed in {int((end-start)//60)}m {int((end-start)%60)}s", file=sys.stderr)

    # 6) collect and save results
    print("[6/6] Saving results…", file=sys.stderr)
    df = bm.get_results(min_max_scale=False)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output)
    print(f"Saved full scIB table to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
