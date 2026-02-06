#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Label propagation with FAISS (GPU if available, CPU fallback).

Overview
--------
This script uses FAISS for fast k-nearest neighbor search on large datasets.
It supports:
- Exact search (Flat L2) on GPU/CPU.
- Approximate search (IVF-Flat) with tunable nlist/nprobe for speed vs. recall.

Like the CPU script, it streams queries and uses per-row bincount voting.

Inputs
------
- H5AD with .obsm[EMBED_KEY] and .obs[LABEL_COL].

Outputs
-------
- CSV with 'propagated_anno' for query cells only.

FAISS indices
-------------
- flat : Exact L2 (fastest on GPU if fits; uses more memory).
- ivf  : Inverted file with Flat residuals (trainable; excellent for scale).
         Tune nlist (coarse clusters) and nprobe (clusters scanned at query).
         Higher nprobe => higher recall, slower.

Example
-------
python label_propagation_faiss.py \
  --h5ad /path/to/data.h5ad \
  --embed-key X_scVI \
  --label-col published_anno \
  --k 35 \
  --index ivf --nlist 16384 --nprobe 32 \
  --batch-size 200000 \
  --gpu-id 0 \
  --out /path/to/propagated_labels_faiss.csv
"""

import argparse
import time
from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc


def make_obs_columns_unique(adata) -> None:
    """Ensure unique column names in .obs to avoid pandas ambiguity."""
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
    Prepare reference/query arrays and integer-encoded labels (as in the CPU script).
    """
    if embed_key not in adata.obsm:
        raise KeyError(f"Embedding '{embed_key}' not found in .obsm")
    if label_col not in adata.obs:
        raise KeyError(f"Label column '{label_col}' not found in .obs")

    make_obs_columns_unique(adata)
    X = np.asarray(adata.obsm[embed_key], dtype=np.float32, order="C")

    is_ref = adata.obs[label_col].notna().values
    is_qry = ~is_ref
    if is_ref.sum() == 0 or is_qry.sum() == 0:
        raise ValueError(
            f"Need both reference and query cells. Found ref={is_ref.sum()}, qry={is_qry.sum()}."
        )

    X_ref = X[is_ref]
    X_qry = X[is_qry]

    ref_cats = adata.obs.loc[is_ref, label_col].astype("category")
    classes = ref_cats.cat.categories.to_numpy()
    y_ref = ref_cats.cat.codes.to_numpy().astype(np.int32)

    qry_names = adata.obs_names[is_qry].to_numpy()
    return X_ref, X_qry, y_ref, classes, qry_names


def faiss_search(
    X_ref: np.ndarray,
    X_qry: np.ndarray,
    k: int,
    index_type: str = "ivf",
    nlist: int = 16384,
    nprobe: int = 32,
    gpu_id: int = 0,
    force_cpu: bool = False,
) -> np.ndarray:
    """
    Build FAISS index and return neighbor indices for all queries.

    Parameters
    ----------
    index_type : {'ivf','flat'}
        'flat' => exact L2; 'ivf' => IVF-Flat (approximate).
    nlist : int
        # of coarse centroids for IVF (ignored for 'flat').
    nprobe : int
        # of centroids probed at search time (ignored for 'flat').
    gpu_id : int
        GPU id (0-based) to use when FAISS GPU is available.
    force_cpu : bool
        If True, use CPU even if FAISS GPU is present.

    Returns
    -------
    I : ndarray, shape (N_qry, k)
        Indices into X_ref of the k nearest neighbors per query.
    """
    import faiss

    d = X_ref.shape[1]
    use_gpu = (not force_cpu) and hasattr(faiss, "StandardGpuResources")

    if index_type == "flat":
        # Exact L2 index
        index = faiss.IndexFlatL2(d)
        if use_gpu:
            res = faiss.StandardGpuResources()
            gindex = faiss.index_cpu_to_gpu(res, gpu_id, index)
            gindex.add(X_ref)
            # Batch queries to be memory-safe
            bs = 200_000
            parts = []
            for s in range(0, X_qry.shape[0], bs):
                e = min(s + bs, X_qry.shape[0])
                _, I = gindex.search(X_qry[s:e], k)
                parts.append(I)
            return np.vstack(parts)
        else:
            index.add(X_ref)
            bs = 200_000
            parts = []
            for s in range(0, X_qry.shape[0], bs):
                e = min(s + bs, X_qry.shape[0])
                _, I = index.search(X_qry[s:e], k)
                parts.append(I)
            return np.vstack(parts)

    # IVF-Flat (trainable)
    quant = faiss.IndexFlatL2(d)
    index = faiss.IndexIVFFlat(quant, d, nlist, faiss.METRIC_L2)

    # FAISS requires training: use a random sample of references
    train_n = min(200_000, X_ref.shape[0])
    train_idx = np.random.choice(X_ref.shape[0], train_n, replace=False)
    index.train(X_ref[train_idx])

    if use_gpu:
        res = faiss.StandardGpuResources()
        gindex = faiss.index_cpu_to_gpu(res, gpu_id, index)
        gindex.nprobe = nprobe
        gindex.add(X_ref)
        bs = 200_000
        parts = []
        for s in range(0, X_qry.shape[0], bs):
            e = min(s + bs, X_qry.shape[0])
            _, I = gindex.search(X_qry[s:e], k)
            parts.append(I)
        return np.vstack(parts)
    else:
        # CPU IVF
        index.nprobe = nprobe
        index.add(X_ref)
        bs = 200_000
        parts = []
        for s in range(0, X_qry.shape[0], bs):
            e = min(s + bs, X_qry.shape[0])
            _, I = index.search(X_qry[s:e], k)
            parts.append(I)
        return np.vstack(parts)


def propagate_with_neighbors(
    I: np.ndarray,
    y_ref: np.ndarray,
    classes: np.ndarray,
    min_score: float,
    batch_size: int,
) -> np.ndarray:
    """
    Given neighbor indices and integer-encoded ref labels, produce final labels.
    """
    n_q = I.shape[0]
    n_classes = len(classes)
    out = np.empty(n_q, dtype=object)

    for s in range(0, n_q, batch_size):
        e = min(s + batch_size, n_q)
        knn_cls = y_ref[I[s:e]]  # (B, k)

        counts = np.zeros((knn_cls.shape[0], n_classes), dtype=np.int32)
        for i in range(knn_cls.shape[0]):
            counts[i] = np.bincount(knn_cls[i], minlength=n_classes)

        probs = counts / counts.sum(axis=1, keepdims=True).clip(min=1)
        best = probs.argmax(axis=1)
        score = probs.max(axis=1)

        lab = classes[best].astype(object)
        lab[score <= min_score] = "low_confidence"
        out[s:e] = lab

    return out


def main():
    ap = argparse.ArgumentParser(
        description="Label propagation with FAISS (GPU preferred, CPU fallback)."
    )
    ap.add_argument("--h5ad", required=True, help="Path to input .h5ad file.")
    ap.add_argument("--embed-key", required=True, help="Key in .obsm with embedding (e.g., X_scVI).")
    ap.add_argument("--label-col", required=True, help="obs column with reference labels (NaN marks query).")
    ap.add_argument("--k", type=int, default=35, help="Number of neighbors (default: 35).")
    ap.add_argument("--min-score", type=float, default=0.30, help="Min winning fraction to accept label (default: 0.30).")
    ap.add_argument("--batch-size", type=int, default=200_000, help="Batch size for voting (default: 200000).")
    ap.add_argument("--index", choices=["ivf", "flat"], default="ivf", help="FAISS index type (default: ivf).")
    ap.add_argument("--nlist", type=int, default=16384, help="IVF: number of coarse centroids (default: 16384).")
    ap.add_argument("--nprobe", type=int, default=32, help="IVF: centroids probed at search (default: 32).")
    ap.add_argument("--gpu-id", type=int, default=0, help="GPU id (default: 0).")
    ap.add_argument("--force-cpu", action="store_true", help="Force CPU even if FAISS GPU is present.")
    ap.add_argument("--out", required=True, help="Output CSV path for propagated labels (query cells only).")
    args = ap.parse_args()

    print(f"[load] {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    # Ensure embedding is float32 contiguous
    adata.obsm[args.embed_key] = np.asarray(adata.obsm[args.embed_key], dtype=np.float32, order="C")

    X_ref, X_qry, y_ref, classes, qry_names = build_reference_and_query(
        adata, label_col=args.label_col, embed_key=args.embed_key
    )

    print(
        f"[faiss] refs={X_ref.shape[0]:,}  qry={X_qry.shape[0]:,}  "
        f"dim={X_ref.shape[1]}  k={args.k}  index={args.index}"
    )

    t0 = time.time()
    I = faiss_search(
        X_ref=X_ref,
        X_qry=X_qry,
        k=args.k,
        index_type=args.index,
        nlist=args.nlist,
        nprobe=args.nprobe,
        gpu_id=args.gpu_id,
        force_cpu=args.force_cpu,
    )
    print(f"[faiss] neighbor search done in {time.time()-t0:.1f}s")

    out = propagate_with_neighbors(
        I=I,
        y_ref=y_ref,
        classes=classes,
        min_score=args.min_score,
        batch_size=args.batch_size,
    )

    df = pd.DataFrame({"propagated_anno": out}, index=qry_names)
    df.to_csv(args.out)
    print(f"[save] wrote: {args.out}  (rows={len(df):,})")


if __name__ == "__main__":
    main()
