#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Label propagation via exact k-nearest neighbors on CPU (SciPy cKDTree).

Overview
--------
This script propagates annotations from "reference" cells to "query" cells
using majority vote among k nearest neighbors in a given embedding (e.g., scVI).
- "Reference" cells are those with a non-null value in a chosen obs column.
- "Query" cells are those with a null (NaN) in that column.

This version implements several performance fixes for very large datasets:
1) Works directly on contiguous float32 numpy arrays (no pandas in the hot path).
2) Streams query batches (no giant dense/sparse matrices).
3) Avoids get_dummies/one-hot matrices: uses per-row np.bincount on integer-encoded classes.
4) Uses cKDTree.query(..., workers=...) for multithreaded search (SciPy ≥1.6).

Inputs
------
- H5AD file with:
  - embedding in .obsm[EMB_KEY] (e.g., "X_scVI")
  - annotation column in .obs[LABEL_COL] (e.g., "published_anno")
- Query cells are those where .obs[LABEL_COL] is NA; reference cells are not NA.

Outputs
-------
- CSV with a single column 'propagated_anno' for query cells only.
  The row index matches the AnnData obs_names of query cells.

Complexity (roughly)
--------------------
- Build tree: O(N_ref log N_ref)
- Query:      O(N_qry log N_ref + N_qry * k)
- Memory:     O(N_ref * dim) + batch-sized temporaries

Example
-------
python label_propagation_knn.py \
    --h5ad /path/to/data.h5ad \
    --embed-key X_scVI \
    --label-col published_anno \
    --k 35 \
    --min-score 0.30 \
    --batch-size 200000 \
    --workers -1 \
    --out /path/to/propagated_labels.csv
"""

import argparse
import time
from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree


def make_obs_columns_unique(adata) -> None:
    """
    Ensure adata.obs has unique column names. Duplicate names in .obs can cause
    pandas to return a DataFrame instead of a Series and break downstream code.
    """
    if not adata.obs.columns.is_unique:
        seen = {}
        new_cols = []
        for c in adata.obs.columns.astype(str):
            if c in seen:
                seen[c] += 1
                new_cols.append(f"{c}__dup{seen[c]}")
            else:
                seen[c] = 0
                new_cols.append(c)
        adata.obs.columns = new_cols


def build_reference_and_query(
    adata, label_col: str, embed_key: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Split AnnData into reference and query subsets, and prepare all arrays.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with .obs[label_col] and .obsm[embed_key].
    label_col : str
        Column denoting known labels for reference cells. NaN => query cell.
    embed_key : str
        Key in .obsm containing the embedding (e.g., "X_scVI").

    Returns
    -------
    X_ref : float32 array, shape (N_ref, d)
    X_qry : float32 array, shape (N_qry, d)
    y_ref : int32 array, shape (N_ref,)
        Integer-encoded labels for reference cells.
    classes : np.ndarray of str
        The category names corresponding to y_ref codes.
    qry_names : np.ndarray of str
        obs_names for query cells (used as index in the output).
    """
    if embed_key not in adata.obsm:
        raise KeyError(f"Embedding '{embed_key}' not found in .obsm")

    if label_col not in adata.obs:
        raise KeyError(f"Label column '{label_col}' not found in .obs")

    # Ensure consistent .obs names and types
    make_obs_columns_unique(adata)

    # Prepare embedding as contiguous float32 for faster math
    X = np.asarray(adata.obsm[embed_key], dtype=np.float32, order="C")

    # Boolean masks for reference/query
    is_ref = adata.obs[label_col].notna().values
    is_qry = ~is_ref

    if is_ref.sum() == 0 or is_qry.sum() == 0:
        raise ValueError(
            "Need both reference and query cells. "
            f"Found ref={is_ref.sum()}, qry={is_qry.sum()}."
        )

    X_ref = X[is_ref]
    X_qry = X[is_qry]

    # Encode reference labels to integers (categories)
    ref_cats = adata.obs.loc[is_ref, label_col].astype("category")
    classes = ref_cats.cat.categories.to_numpy()
    y_ref = ref_cats.cat.codes.to_numpy().astype(np.int32)

    qry_names = adata.obs_names[is_qry].to_numpy()

    return X_ref, X_qry, y_ref, classes, qry_names


def propagate_knn_cpu(
    X_ref: np.ndarray,
    X_qry: np.ndarray,
    y_ref: np.ndarray,
    classes: np.ndarray,
    k: int = 25,
    min_score: float = 0.80,
    batch_size: int = 200_000,
    workers: int = -1,
) -> np.ndarray:
    """
    Propagate labels from reference to query via exact kNN majority vote.

    Implementation details for speed:
    - Builds a cKDTree on X_ref (float32 contiguous).
    - Queries are processed in batches to control RAM.
    - For each query, we map neighbor indices -> class ids -> bincount per row.
    - Probability per class = counts / k; low-confidence threshold applied.

    Parameters
    ----------
    X_ref, X_qry : float32 arrays
        Reference and query embeddings.
    y_ref : int32 array
        Integer-encoded labels for reference cells.
    classes : array-like of str
        Label names corresponding to the integer codes in y_ref.
    k : int
        Number of nearest neighbors.
    min_score : float
        Minimum winning class fraction (0..1) to accept the label; else "low_confidence".
    batch_size : int
        Query batch size (tune to RAM).
    workers : int
        Threads for cKDTree.query. -1 means "all cores" (SciPy ≥1.6).

    Returns
    -------
    out_labels : object array of shape (N_qry,)
        Predicted labels for all query cells (strings + "low_confidence" for low scores).
    """
    t0 = time.time()
    tree = cKDTree(X_ref)

    n_q = X_qry.shape[0]
    n_classes = len(classes)
    out = np.empty(n_q, dtype=object)

    print(
        f"[kNN-CPU] refs={X_ref.shape[0]:,}  qry={n_q:,}  dim={X_ref.shape[1]}  k={k}  "
        f"batch={batch_size:,}  workers={workers}"
    )

    # Batch over queries
    for s in range(0, n_q, batch_size):
        e = min(s + batch_size, n_q)

        # Query neighbors. cKDTree supports 'workers' in SciPy>=1.6.
        # If workers is not supported in your SciPy, this will raise TypeError;
        # we fallback to single-threaded query in that case.
        try:
            _, knn_idx = tree.query(X_qry[s:e], k=k, workers=workers)
        except TypeError:
            _, knn_idx = tree.query(X_qry[s:e], k=k)

        if k == 1:
            knn_idx = knn_idx[:, None]  # (B,) -> (B,1)

        # Map neighbor indices -> class ids
        knn_cls = y_ref[knn_idx]  # shape (B, k), int32

        # Per-row bincount to get class counts
        counts = np.zeros((knn_cls.shape[0], n_classes), dtype=np.int32)
        for i in range(knn_cls.shape[0]):
            counts[i] = np.bincount(knn_cls[i], minlength=n_classes)

        # Turn counts into probabilities and pick winners
        probs = counts / counts.sum(axis=1, keepdims=True).clip(min=1)
        best_id = probs.argmax(axis=1)
        best_sc = probs.max(axis=1)

        labels = classes[best_id].astype(object)
        labels[best_sc <= min_score] = "low_confidence"
        out[s:e] = labels

    print(f"[kNN-CPU] done in {time.time() - t0:.1f}s")
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Label propagation via exact kNN (CPU) with streaming and bincount."
    )
    ap.add_argument("--h5ad", required=True, help="Path to input .h5ad file.")
    ap.add_argument("--embed-key", required=True, help="Key in .obsm with embedding (e.g., X_scVI).")
    ap.add_argument("--label-col", required=True, help="obs column with reference labels (NaN marks query).")
    ap.add_argument("--k", type=int, default=25, help="Number of neighbors (default: 25).")
    ap.add_argument("--min-score", type=float, default=0.80, help="Min winning fraction to accept label (default: 0.80).")
    ap.add_argument("--batch-size", type=int, default=200_000, help="Query batch size (default: 200000).")
    ap.add_argument("--workers", type=int, default=-1, help="Threads for cKDTree.query (-1=all cores).")
    ap.add_argument("--out", required=True, help="Output CSV path for propagated labels (query cells only).")
    args = ap.parse_args()

    print(f"[load] {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    # Ensure embedding is float32 contiguous
    adata.obsm[args.embed_key] = np.asarray(adata.obsm[args.embed_key], dtype=np.float32, order="C")

    X_ref, X_qry, y_ref, classes, qry_names = build_reference_and_query(
        adata, label_col=args.label_col, embed_key=args.embed_key
    )

    out = propagate_knn_cpu(
        X_ref=X_ref,
        X_qry=X_qry,
        y_ref=y_ref,
        classes=classes,
        k=args.k,
        min_score=args.min_score,
        batch_size=args.batch_size,
        workers=args.workers,
    )

    df = pd.DataFrame({"propagated_anno": out}, index=qry_names)
    df.to_csv(args.out)
    print(f"[save] wrote: {args.out}  (rows={len(df):,})")


if __name__ == "__main__":
    main()
