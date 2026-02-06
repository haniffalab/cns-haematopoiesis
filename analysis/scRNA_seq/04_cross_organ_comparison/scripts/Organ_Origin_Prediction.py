#!/usr/bin/env python3
import os
import argparse

import scanpy as sc
from scipy import sparse
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report

def main(input_h5ad, organ_label, output_dir, n_genes):
    # Load (assumes HVG filter already applied)
    adata = sc.read_h5ad(input_h5ad)
    adata.obs['Organ'] = adata.obs['Organ'].astype(str).str.strip()
    keep = []
    for organ in adata.obs["Organ"].unique():
        idx = np.where(adata.obs['Organ'] == organ)[0]
        if len(idx) > 20000:
            idx = np.random.choice(idx, 9000, replace=False)
        keep.extend(idx)
    adata = adata[keep].copy()
    sc.pp.normalize_total(adata, target_sum=1e4, layer='raw_counts')
    sc.pp.log1p(adata, layer='raw_counts')

    
# 2. Split into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42)

# 3. Instantiate the classifier
    clf = LogisticRegression(
    penalty='l1',
    solver='saga',
    multi_class='multinomial',
    C=0.1,
    tol=1e-4,
    max_iter=2000,
    n_jobs=-1,
    random_state=42
)

    # 4. Fit on the training data
    clf.fit(X_train, y_train)

    # 5. Predict “hard” labels on the test set
    y_pred = clf.predict(X_test)

    # 6. Predict class probabilities
    y_proba = clf.predict_proba(X_test)
    print(classification_report(y_test, y_pred))


    genes = np.array(adata.var_names)

    # 3) Method A: univariate one-vs-rest ANOVA
    mv_dir = os.path.join(output_dir, "univariate_ovo")
    for organ in np.unique(y_train):
        y_bin = (y_train == organ).astype(int)
        sel = SelectKBest(f_classif, k=n_genes).fit(X_train, y_bin)
        top = genes[sel.get_support()]
        od = os.path.join(mv_dir, organ.replace(" ", "_"))
        os.makedirs(od, exist_ok=True)
        pd.Series(top, name='gene')\
          .to_csv(os.path.join(od, f"top_{n_genes}_genes_{organ}.csv"),
                  index=False)
        print(f"[Univariate] {organ}: saved {len(top)} genes → {od}")

    # 4) Method B: multinomial L1-penalized Logistic Regression
    lg_dir = os.path.join(output_dir, "logistic_l1")
    clf = LogisticRegression(
        penalty='l1',
        solver='saga',
        multi_class='multinomial',
        C=0.1,
        tol=1e-4,
        max_iter=2000,
        n_jobs=-1,
        random_state=42
    )
    print("Fitting LogisticRegression…")
    clf.fit(X, y)

    for i, organ in enumerate(clf.classes_):
        coefs = clf.coef_[i]
        idx = np.argsort(np.abs(coefs))[::-1][:n_genes]
        top = genes[idx]
        od = os.path.join(lg_dir, organ.replace(" ", "_"))
        os.makedirs(od, exist_ok=True)
        pd.Series(top, name='gene')\
          .to_csv(os.path.join(od, f"top_{n_genes}_genes_{organ}.csv"),
                  index=False)
        print(f"[Logistic L1] {organ}: saved {len(top)} genes → {od}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Per-organ gene module discovery: univariate & L1-logistic"
    )
    p.add_argument("--input_h5ad",  type=str, required=True,
                   help="Path to HVG-filtered .h5ad file")
    p.add_argument("--organ_label", type=str, default="Organ",
                   help="Column in adata.obs with organ labels")
    p.add_argument("--output_dir",   type=str, required=True,
                   help="Directory to write organ-specific CSVs")
    p.add_argument("--n_genes",      type=int, default=100,
                   help="Number of genes per organ per method")
    args = p.parse_args()

    main(
        input_h5ad=args.input_h5ad,
        organ_label=args.organ_label,
        output_dir=args.output_dir,
        n_genes=args.n_genes
    )

