#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import numpy as np
from scipy.io import mmread
from pathlib import Path
import argparse

def read_velocity_sample(sample_dir, sample_id):
    vdir = Path(sample_dir) / "output" / "Velocyto" / "filtered"

    spliced   = mmread(vdir / "spliced.mtx.gz").tocsr().T
    unspliced = mmread(vdir / "unspliced.mtx.gz").tocsr().T

    features = pd.read_csv(vdir / "features.tsv.gz", header=None, sep="\t")
    barcodes = pd.read_csv(vdir / "barcodes.tsv.gz", header=None, sep="\t")[0].values

    gene_names = features[1].values
    gene_ids   = features[0].values

    adata_v = ad.AnnData(X=spliced.copy())
    adata_v.var_names = gene_names
    adata_v.var["gene_ids"] = gene_ids

    adata_v.obs_names = [f"{bc}-{sample_id}" for bc in barcodes]
    adata_v.obs["Sample_ID"] = sample_id

    adata_v.layers["spliced"]   = spliced
    adata_v.layers["unspliced"] = unspliced

    adata_v.var_names_make_unique()
    return adata_v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True,
                    help="Directory containing Sanger/ and GSE143753/")
    ap.add_argument("--output_dir", required=True,
                    help="Directory to write outputs")
    ap.add_argument("--merged_name", default="velocity_merged.h5ad")
    args = ap.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------
    # FIX: explicitly gather samples from Sanger + GSE143753
    # -------------------------------------------------------
    subdirs = ["Sanger", "GSE143753"]
    samples = []

    for sd in subdirs:
        subpath = input_dir / sd
        if not subpath.exists():
            print(f"Warning: folder {sd} not found, skipping.")
            continue

        # sample directories inside each dataset
        for d in subpath.iterdir():
            if d.is_dir():
                samples.append(d)
    # samples is now list of Path objects, not names

    print(f"Found {len(samples)} sample folders across Sanger + GSE143753.")
    # -------------------------------------------------------

    vel_adatas = []

    for sample_path in samples:
        sid = sample_path.name
        vdir = sample_path / "output" / "Velocyto" / "filtered"

        if not vdir.exists():
            print(f"Skipping {sid}: no velocity folder.")
            continue

        try:
            adata_v = read_velocity_sample(sample_path, sid)
        except Exception as e:
            print(f"Skipping {sid} due to error: {e}")
            continue

        # save per-sample
        out_sdir = output_dir / sid
        out_sdir.mkdir(exist_ok=True)
        sample_h5 = out_sdir / f"{sid}_velocity.h5ad"
        adata_v.write_h5ad(sample_h5)
        print(f"Saved {sample_h5}")

        vel_adatas.append(adata_v)

    if len(vel_adatas) == 0:
        raise ValueError("No samples processed.")

    # merge
    merged = ad.concat(vel_adatas, join="outer", axis=0, index_unique=None)
    merged_path = output_dir / args.merged_name
    merged.write_h5ad(merged_path)
    print(f"Merged velocity AnnData saved to {merged_path}")


if __name__ == "__main__":
    main()
