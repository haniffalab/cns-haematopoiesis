#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report


def main(input_h5ad, organ_label, output_dir, n_genes):

    # ===========================================================
    # 1. LOAD DATA
    # ===========================================================
    print("\nLoading data...")
    adata = sc.read_h5ad(input_h5ad)
    adata.obs[organ_label] = adata.obs[organ_label].astype(str).str.strip()

    # ===========================================================
    # 2. SUBSAMPLE ORGANS >20k cells → keep 9k cells
    # ===========================================================
    print("Subsampling organs >20k cells down to 9000 cells...")
    keep = []
    for organ in adata.obs[organ_label].unique():
        idx = np.where(adata.obs[organ_label] == organ)[0]
        if len(idx) > 20000:
            idx = np.random.choice(idx, 9000, replace=False)
        keep.extend(idx)
    adata = adata[keep].copy()

    print(f"Final dataset shape: {adata.shape}")

    # ===========================================================
    # 3. NORMALIZATION (use raw_counts → normalized/log1p)
    # ===========================================================
    print("\nNormalizing data...")
    # start from raw counts layer
    adata.X = adata.layers["raw_counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ===========================================================
    # 3b. RESTRICT TO HIGHLY VARIABLE GENES ONLY
    # ===========================================================
    if "highly_variable" in adata.var.columns:
        print("Using existing HVG annotation in adata.var['highly_variable']...")
    else:
        print("Computing highly variable genes (Seurat v3 flavor)...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=5000,   # adjust if you want a different HVG size
        batch_key='10X_Chemistry_V2'
        )

    hvg_mask = adata.var["highly_variable"].values.astype(bool)
    print(f"Number of HVGs: {hvg_mask.sum()}")

    # subset to HVGs only
    adata = adata[:, hvg_mask].copy()
    print(f"Shape after HVG restriction: {adata.shape}")

    # ===========================================================
    # 4. BUILD X (expression) AND y (label) — HVGs ONLY
    # ===========================================================
    print("\nBuilding design matrix X (HVGs) and labels y...")
    X = adata.X.copy()
    X = X.A if sparse.issparse(X) else np.array(X)
    y = adata.obs[organ_label].values
    genes = np.array(adata.var_names)

    # ===========================================================
    # 5. TRAIN/TEST SPLIT for classifier evaluation (Method A)
    # ===========================================================
    print("\nSplitting train/test for evaluation...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )

    # ===========================================================
    # 6. METHOD A — Classifier Evaluation + ANOVA
    # ===========================================================
    print("\n===================================================")
    print(" METHOD A — Classification + Univariate ANOVA (HVGs only)")
    print("===================================================")

    clf_eval = LogisticRegression(
        penalty="l1",
        solver="saga",
        multi_class="multinomial",
        C=0.1,
        tol=1e-4,
        max_iter=2000,
        n_jobs=-1,
        random_state=42
    )

    print("\nFitting logistic regression on TRAINING SET...")
    clf_eval.fit(X_train, y_train)

    print("\n--- Classification Report (TEST SET) ---")
    print(classification_report(y_test, clf_eval.predict(X_test)))

    # --------------------------
    # Univariate ANOVA OVO markers
    # --------------------------
    mv_dir = os.path.join(output_dir, "univariate_ovo")
    os.makedirs(mv_dir, exist_ok=True)

    print("\nSelecting top genes using ANOVA (per organ, within HVGs)...")
    for organ in np.unique(y_train):
        y_bin = (y_train == organ).astype(int)

        sel = SelectKBest(score_func=f_classif, k=n_genes)
        sel.fit(X_train, y_bin)

        top_genes = genes[sel.get_support()]

        out_dir = os.path.join(mv_dir, organ.replace(" ", "_"))
        os.makedirs(out_dir, exist_ok=True)

        pd.Series(top_genes, name="gene").to_csv(
            os.path.join(out_dir, f"top_{n_genes}_genes_{organ}.csv"),
            index=False
        )

        print(f"[Method A] {organ}: saved {len(top_genes)} genes")

    # ===========================================================
    # 7. METHOD B — Multinomial L1 Logistic Regression (Feature Discovery)
    # ===========================================================
    print("\n===================================================")
    print(" METHOD B — Multinomial L1 Logistic Regression (HVGs only)")
    print("===================================================")

    clf_full = LogisticRegression(
        penalty="l1",
        solver="saga",
        multi_class="multinomial",
        C=0.1,
        tol=1e-4,
        max_iter=2000,
        n_jobs=-1,
        random_state=42
    )

    print("\nFitting logistic regression on FULL DATASET (HVGs only)...")
    clf_full.fit(X, y)

    lg_dir = os.path.join(output_dir, "logistic_l1")
    os.makedirs(lg_dir, exist_ok=True)

    print("\nSelecting top genes by coefficient magnitude...")
    for i, organ in enumerate(clf_full.classes_):
        coefs = clf_full.coef_[i]
        idx = np.argsort(np.abs(coefs))[::-1][:n_genes]
        top_genes = genes[idx]

        out_dir = os.path.join(lg_dir, organ.replace(" ", "_"))
        os.makedirs(out_dir, exist_ok=True)

        pd.Series(top_genes, name="gene").to_csv(
            os.path.join(out_dir, f"top_{n_genes}_genes_{organ}.csv"),
            index=False
        )

        print(f"[Method B] {organ}: saved {len(top_genes)} genes")

    print("\nDONE ✓✓✓")


# ===========================================================
# CLI
# ===========================================================
if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Organ-specific gene discovery using ANOVA + L1 Logistic Regression (HVGs only)"
    )
    p.add_argument("--input_h5ad", type=str, required=True)
    p.add_argument("--organ_label", type=str, default="Organ")
    p.add_argument("--output_dir", type=str, required=True)
    p.add_argument("--n_genes", type=int, default=100)
    args = p.parse_args()

    main(
        input_h5ad=args.input_h5ad,
        organ_label=args.organ_label,
        output_dir=args.output_dir,
        n_genes=args.n_genes
    )
