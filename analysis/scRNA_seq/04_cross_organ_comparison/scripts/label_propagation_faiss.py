#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Label propagation with FAISS (GPU if available, CPU fallback) — robust & documented.

Goal
----
Propagate labels from "reference" cells (non-null in obs[label_col]) to "query" cells (null in obs[label_col])
by nearest-neighbour voting in a chosen embedding (e.g., scVI in `.obsm['X_scVI']`).

Why this version is robust
--------------------------
1) **Parameter safety for IVF**:
   - If `nlist` (number of IVF coarse centroids) is too big for the number of reference points, we
     automatically *downscale* it to satisfy FAISS's rough training rule: `train_points ≈ 40 * nlist`.
   - We also cap `nprobe` to never exceed `nlist`.
   - If the reference set is small (heuristically < 50k for typical 20–64D embeddings), we auto-switch to
     **exact** `IndexFlatL2` because IVF overhead is not helpful.

2) **Deterministic training**:
   - IVF requires k-means training on a subset of reference vectors. We sample deterministically
     (`np.random.seed(0)`) so runs are reproducible.

3) **Performance & memory**:
   - We keep the embedding as contiguous `float32`.
   - We stream queries in batches so we never allocate massive temporary arrays.
   - Voting is done with fast `numpy.bincount` per row (no huge one-hot matrices).

4) **Clear logging & summaries**:
   - Sanity info (ref/query counts, dim, index type).
   - IVF parameter adjustments are printed.
   - Final low-confidence fraction is reported.

Inputs
------
- H5AD with:
  - `.obsm[--embed-key]` (e.g. "X_scVI"): the embedding to run kNN on.
  - `.obs[--label-col]`: the known labels for reference cells. Null marks query cells.

Outputs
-------
- CSV written to `--out` with a single column `propagated_anno` for *query* cells only.
  The index is the cell `obs_names` for those query cells.

CLI (unchanged from your original)
----------------------------------
--h5ad PATH            : input AnnData file
--embed-key STR        : embedding key in .obsm (e.g., X_scVI)
--label-col STR        : obs column holding reference labels (NaN marks queries)
--k INT                : number of neighbors (default 35)
--min-score FLOAT      : min winning fraction to accept label; else "low_confidence"
--batch-size INT       : how many queries to process per batch (default 200k)
--index {ivf,flat}     : FAISS index type (default ivf)
--nlist INT            : IVF coarse centroids (default 16384)
--nprobe INT           : IVF clusters probed per query (default 32)
--gpu-id INT           : which GPU to use (default 0)
--force-cpu            : force CPU FAISS even if GPU is available
--out PATH             : output CSV path

Typical settings
----------------
- IVF (default) for big reference sets: nlist=8192–16384, nprobe=32–64.
- Exact FLAT if you want deterministic neighbors and it fits in memory.

"""

import argparse
import time
from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc


# ------------------------------- Utilities ------------------------------- #

def make_obs_columns_unique(adata) -> None:
    """
    Ensure .obs column names are unique.

    Why: Duplicate column names in adata.obs can make pandas return ambiguous
    DataFrames instead of Series and break downstream selection by label_col.
    """
    if not adata.obs.columns.is_unique:
        seen = {}
        new_cols = []
        for c in adata.obs.columns.astype(str):
            # if we've seen this name, append a suffix to make it unique
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
    Split the AnnData into reference/query sets and produce the arrays needed downstream.

    Parameters
    ----------
    adata : AnnData
        Input object with embedding in .obsm[embed_key] and labels in .obs[label_col].
    label_col : str
        Name of the obs column containing known labels (reference). NaN marks queries.
    embed_key : str
        Name of the embedding in .obsm to use for kNN (e.g., "X_scVI").

    Returns
    -------
    X_ref : (N_ref, d) float32
        Embedding rows for reference cells.
    X_qry : (N_qry, d) float32
        Embedding rows for query cells.
    y_ref : (N_ref,) int32
        Integer-encoded reference labels.
    classes : (C,) object/str
        The string labels corresponding to the integer codes in y_ref (category order).
    qry_names : (N_qry,) str
        AnnData obs_names for query cells (to use as CSV index).
    """
    # Basic presence checks
    if embed_key not in adata.obsm:
        raise KeyError(f"Embedding '{embed_key}' not found in .obsm")
    if label_col not in adata.obs:
        raise KeyError(f"Label column '{label_col}' not found in .obs")

    # Avoid .obs duplicate-column issues
    make_obs_columns_unique(adata)

    # Ensure the embedding is contiguous float32 (fast math, predictable memory)
    X = np.asarray(adata.obsm[embed_key], dtype=np.float32, order="C")
    if X.ndim != 2 or X.size == 0:
        raise ValueError(f"Embedding '{embed_key}' is empty or malformed: shape={X.shape}")

    # Define reference (non-null label) vs query (null label)
    is_ref = adata.obs[label_col].notna().values
    is_qry = ~is_ref
    n_ref, n_qry = int(is_ref.sum()), int(is_qry.sum())
    if n_ref == 0 or n_qry == 0:
        raise ValueError(f"Need both reference and query cells. Found ref={n_ref}, qry={n_qry}.")

    # Select rows for reference/query
    X_ref = X[is_ref]
    X_qry = X[is_qry]

    # Convert reference labels to categorical codes (int) for fast bincount voting
    ref_cats = adata.obs.loc[is_ref, label_col].astype("category")
    classes = ref_cats.cat.categories.to_numpy()        # ordered unique label names
    y_ref = ref_cats.cat.codes.to_numpy().astype(np.int32)  # integer codes aligned to `classes`

    # Keep query names to index the output CSV
    qry_names = adata.obs_names[is_qry].to_numpy()

    return X_ref, X_qry, y_ref, classes, qry_names


def _choose_index_type(requested: str, n_ref: int, dim: int) -> str:
    """
    Decide whether to use 'ivf' or 'flat'.

    Heuristic:
      - If the user requested 'flat' → use 'flat'.
      - If the reference set is small (heuristically < 50k in typical scVI dims),
        IVF overhead isn't worth it → use 'flat'.
      - Otherwise, use 'ivf'.

    Notes:
      - This is a pragmatic heuristic; you can adjust the 50k threshold if needed.
    """
    if requested == "flat":
        return "flat"
    if n_ref < 50_000:
        return "flat"
    return "ivf"


def _adjust_ivf_params(n_ref: int, nlist: int, nprobe: int):
    """
    Adjust IVF parameters to satisfy FAISS's rough training rule:
      number_of_training_points ≈ 40 * nlist

    What we do:
      - Compute desired training points (40 * nlist).
      - If we don't have enough references to meet that, *downscale nlist* so
        the rule can be satisfied.
      - Training size is max(200k, 40*nlist) but never larger than #refs.
      - Ensure nprobe ≤ nlist at the end.

    Returns
    -------
    nlist_adj : int
        Possibly reduced nlist.
    nprobe_adj : int
        nprobe capped at nlist_adj.
    train_n : int
        Number of reference points to sample for IVF training.
    """
    # Desired training size by rule-of-thumb
    desired_train = 40 * nlist

    # Start from a practical floor (200k) but never exceed available refs
    train_n = min(n_ref, max(200_000, desired_train))

    # If we cannot meet desired_train with available refs, reduce nlist
    if n_ref < desired_train:
        # new_nlist so that 40 * new_nlist <= n_ref
        new_nlist = max(1024, n_ref // 40)  # keep ≥1024 to avoid tiny nlist
        nlist = int(new_nlist)
        desired_train = 40 * nlist
        train_n = min(n_ref, max(200_000, desired_train))

    # Cap nprobe to not exceed the actual number of lists
    nprobe = min(nprobe, nlist) if nlist > 0 else nprobe

    return int(nlist), int(nprobe), int(train_n)


# ------------------------- FAISS neighbor search ------------------------- #

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
    Build a FAISS index (GPU if available) and return neighbor indices for all queries.

    Parameters
    ----------
    X_ref : (N_ref, d) float32
        Reference embeddings.
    X_qry : (N_qry, d) float32
        Query embeddings.
    k : int
        Number of neighbors to retrieve.
    index_type : {"ivf","flat"}
        FAISS index family to use.
    nlist, nprobe : int
        IVF hyperparameters (ignored for "flat").
    gpu_id : int
        Which GPU to use when FAISS-GPU is available.
    force_cpu : bool
        If True, force CPU path even if GPU is present.

    Returns
    -------
    I : (N_qry, k) int64
        For each query, indices of its k nearest neighbors in X_ref.
    """
    # Import FAISS and check if GPU is available
    try:
        import faiss
    except Exception as e:
        raise RuntimeError(
            "FAISS is not installed or failed to import. "
            "Install 'faiss-gpu' (preferred) or 'faiss-cpu' in your environment."
        ) from e

    d = X_ref.shape[1]  # embedding dimensionality
    use_gpu = (not force_cpu) and hasattr(faiss, "StandardGpuResources")

    # Decide effective index type (may override user request if refs are small)
    index_type_effective = _choose_index_type(index_type, n_ref=X_ref.shape[0], dim=d)

    # ------------------------- Exact FLAT path ------------------------- #
    if index_type_effective == "flat":
        index = faiss.IndexFlatL2(d)  # exact L2 distance

        if use_gpu:
            # Allocate FAISS GPU resources and move the CPU index to GPU
            res = faiss.StandardGpuResources()
            gindex = faiss.index_cpu_to_gpu(res, gpu_id, index)
            # Add all reference vectors to the GPU index
            gindex.add(X_ref)

            # Query in batches to limit memory
            bs = 200_000
            parts = []
            for s in range(0, X_qry.shape[0], bs):
                e = min(s + bs, X_qry.shape[0])
                _, I = gindex.search(X_qry[s:e], k)  # distances not needed
                parts.append(I)
            return np.vstack(parts)

        else:
            # CPU fallback: still exact, but slower without GPU
            print("[faiss] GPU not available; using CPU flat index.", flush=True)
            index.add(X_ref)
            bs = 200_000
            parts = []
            for s in range(0, X_qry.shape[0], bs):
                e = min(s + bs, X_qry.shape[0])
                _, I = index.search(X_qry[s:e], k)
                parts.append(I)
            return np.vstack(parts)

    # ------------------------- IVF-Flat path -------------------------- #
    # Adjust IVF params to be safe and effective given the available refs
    nlist_adj, nprobe_adj, train_n = _adjust_ivf_params(
        n_ref=X_ref.shape[0], nlist=nlist, nprobe=nprobe
    )
    if (nlist_adj != nlist) or (nprobe_adj != nprobe):
        print(
            f"[faiss][ivf] adjusted params: nlist={nlist}→{nlist_adj}, "
            f"nprobe={nprobe}→{nprobe_adj}, train_n={train_n}",
            flush=True,
        )
    else:
        print(
            f"[faiss][ivf] using nlist={nlist_adj}, nprobe={nprobe_adj}, train_n={train_n}",
            flush=True,
        )

    # IVF-Flat needs a coarse quantizer (flat L2) and the IVF index
    try:
        import faiss
    except Exception:
        # already handled above, but keeps linters happy
        raise
    quant = faiss.IndexFlatL2(d)
    index = faiss.IndexIVFFlat(quant, d, nlist_adj, faiss.METRIC_L2)

    # Deterministic sampling for k-means training on reference vectors
    np.random.seed(0)
    train_idx = np.random.choice(X_ref.shape[0], train_n, replace=False)
    index.train(X_ref[train_idx])

    if use_gpu:
        # Move trained IVF index to GPU
        res = faiss.StandardGpuResources()
        gindex = faiss.index_cpu_to_gpu(res, gpu_id, index)
        gindex.nprobe = nprobe_adj  # how many coarse lists to scan per query
        gindex.add(X_ref)           # add all reference vectors to the IVF index

        # Query in batches
        bs = 200_000
        parts = []
        for s in range(0, X_qry.shape[0], bs):
            e = min(s + bs, X_qry.shape[0])
            _, I = gindex.search(X_qry[s:e], k)
            parts.append(I)
        return np.vstack(parts)

    else:
        # CPU IVF fallback
        print("[faiss] GPU not available; using CPU IVF index.", flush=True)
        index.nprobe = nprobe_adj
        index.add(X_ref)
        bs = 200_000
        parts = []
        for s in range(0, X_qry.shape[0], bs):
            e = min(s + bs, X_qry.shape[0])
            _, I = index.search(X_qry[s:e], k)
            parts.append(I)
        return np.vstack(parts)


# --------------------------- Voting / Propagation -------------------------- #

def propagate_with_neighbors(
    I: np.ndarray,
    y_ref: np.ndarray,
    classes: np.ndarray,
    min_score: float,
    batch_size: int,
) -> np.ndarray:
    """
    Convert neighbor indices to labels via per-row bincount voting.

    Mechanism
    ---------
    - For each query, we have k neighbor indices of reference cells: I[q, :].
    - We map those indices to integer labels `y_ref`, then count how many times each
      class occurs (np.bincount).
    - The winning probability is max(counts)/k; if it is ≤ min_score, we emit "low_confidence".

    Returns
    -------
    out : (N_qry,) object
        Predicted label strings (or "low_confidence") for every query.
    """
    n_q = I.shape[0]
    n_classes = len(classes)
    out = np.empty(n_q, dtype=object)  # preallocate array for final labels

    # Process queries in batches to keep memory bounded
    for s in range(0, n_q, batch_size):
        e = min(s + batch_size, n_q)

        # y_ref maps reference index -> int code (0..n_classes-1)
        # I[s:e] are neighbor indices into the reference set
        knn_cls = y_ref[I[s:e]]  # shape (B, k)

        # For each query row, count class occurrences with a fixed-length bincount
        counts = np.zeros((knn_cls.shape[0], n_classes), dtype=np.int32)
        for i in range(knn_cls.shape[0]):
            counts[i] = np.bincount(knn_cls[i], minlength=n_classes)

        # Convert counts to probabilities and choose the argmax class
        probs = counts / counts.sum(axis=1, keepdims=True).clip(min=1)  # safe divide
        best = probs.argmax(axis=1)
        score = probs.max(axis=1)

        # Map integer class ids back to string labels
        lab = classes[best].astype(object)

        # Enforce minimum confidence threshold
        lab[score <= min_score] = "low_confidence"

        # Store into output slice
        out[s:e] = lab

    return out


# ---------------------------------- Main ---------------------------------- #

def main():
    # Parse command-line arguments (kept identical to your original)
    ap = argparse.ArgumentParser(
        description="Label propagation with FAISS (GPU preferred, CPU fallback)."
    )
    ap.add_argument("--h5ad", required=True, help="Path to input .h5ad file.")
    ap.add_argument("--embed-key", required=True, help="Key in .obsm with embedding (e.g., X_scVI).")
    ap.add_argument("--label-col", required=True, help="obs column with reference labels (NaN marks query).")
    ap.add_argument("--k", type=int, default=25, help="Number of neighbors (default: 35).")
    ap.add_argument("--min-score", type=float, default=0.80, help="Min winning fraction to accept label (default: 0.80).")
    ap.add_argument("--batch-size", type=int, default=200_000, help="Batch size for voting (default: 200000).")
    ap.add_argument("--index", choices=["ivf", "flat"], default="ivf", help="FAISS index type (default: ivf).")
    ap.add_argument("--nlist", type=int, default=16384, help="IVF: number of coarse centroids (default: 16384).")
    ap.add_argument("--nprobe", type=int, default=32, help="IVF: centroids probed at search (default: 32).")
    ap.add_argument("--gpu-id", type=int, default=0, help="GPU id (default: 0).")
    ap.add_argument("--force-cpu", action="store_true", help="Force CPU even if FAISS GPU is present.")
    ap.add_argument("--out", required=True, help="Output CSV path for propagated labels (query cells only).")
    args = ap.parse_args()

    # Load the data
    print(f"[load] {args.h5ad}", flush=True)
    adata = sc.read_h5ad(args.h5ad)
    # Convert to string first to avoid category mismatch
    adata.obs[args.label_col] = adata.obs[args.label_col].astype(str)
    adata.obs[args.label_col] = adata.obs[args.label_col].replace(['nan', 'NaN', 'None', ''], np.nan)

    # Ensure the embedding is in float32 C-order for performance
    adata.obsm[args.embed_key] = np.asarray(adata.obsm[args.embed_key], dtype=np.float32, order="C")

    # Build arrays for reference/query and encoded labels
    X_ref, X_qry, y_ref, classes, qry_names = build_reference_and_query(
        adata, label_col=args.label_col, embed_key=args.embed_key
    )

    # Log a clear sanity line up-front
    print(
        f"[sanity] ref={X_ref.shape[0]:,}  qry={X_qry.shape[0]:,}  dim={X_ref.shape[1]}  "
        f"embed_key={args.embed_key}  label_col={args.label_col}",
        flush=True,
    )

    # Decide effective index type (may override IVF→FLAT if refs are small)
    effective_index = _choose_index_type(args.index, n_ref=X_ref.shape[0], dim=X_ref.shape[1])
    print(
        f"[faiss] refs={X_ref.shape[0]:,}  qry={X_qry.shape[0]:,}  dim={X_ref.shape[1]}  k={args.k}  index={effective_index}",
        flush=True,
    )

    # Run FAISS neighbor search (GPU if available)
    t0 = time.time()
    I = faiss_search(
        X_ref=X_ref,
        X_qry=X_qry,
        k=args.k,
        index_type=effective_index,
        nlist=args.nlist,
        nprobe=args.nprobe,
        gpu_id=args.gpu_id,
        force_cpu=args.force_cpu,
    )
    print(f"[faiss] neighbor search done in {time.time()-t0:.1f}s", flush=True)

    # Vote / propagate labels
    out = propagate_with_neighbors(
        I=I,
        y_ref=y_ref,
        classes=classes,
        min_score=args.min_score,
        batch_size=args.batch_size,
    )

    # Save as a single-column CSV indexed by query cell names
    df = pd.DataFrame({"propagated_anno": out}, index=qry_names)
    df.to_csv(args.out)

    # Final summary
    lc = (df["propagated_anno"] == "low_confidence").mean()
    print(f"[save] wrote: {args.out}  (rows={len(df):,})", flush=True)
    print(f"[summary] low_confidence = {lc:.3%}", flush=True)


if __name__ == "__main__":
    main()
