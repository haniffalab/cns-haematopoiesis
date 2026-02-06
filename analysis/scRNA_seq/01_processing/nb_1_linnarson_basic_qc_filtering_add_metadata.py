#!/usr/bin/env python3
# coding: utf-8

import argparse
import scanpy as sc
import pandas as pd
import helper_functions
import numpy as np
from scipy.stats import median_abs_deviation


def add_metadata_to_obs(
    adata,
    df: pd.DataFrame,
    df_key: str,
    obs_key: str = None,
):
    """
    Merge EVERY column of `df` into `adata.obs` by matching
      df[df_key]  ↔  adata.obs[obs_key] (or adata.obs_names if obs_key is None).

    Keeps duplicates in df_key by dropping all but the first.
    Leaves adata.obs index untouched.
    """
    # 1. Sanity checks
    assert df_key in df.columns, f"{df_key!r} not in metadata DataFrame"
    if obs_key is not None:
        assert obs_key in adata.obs.columns, f"{obs_key!r} not in adata.obs"

    # 2. Deduplicate metadata on the key, keep only first
    df_unique = df.drop_duplicates(subset=[df_key], keep="first").copy()
    # 3. Build a lookup Series indexed by the key
    df_indexed = df_unique.set_index(df_key)
    # 4. Pick the “left” key from adata.obs
    if obs_key is None:
        # map from obs_names
        left_index = pd.Series(adata.obs_names, index=adata.obs_names, name=df_key)
    else:
        # map from an existing obs column
        left_index = adata.obs[obs_key].rename(df_key)
    # 5. For every column in df_indexed (including df_key itself), map into adata.obs
    for col in df_indexed.columns:
        adata.obs[col] = left_index.map(df_indexed[col])
    # 6. Report if anything expected went missing
    missing = adata.obs[df_indexed.columns].isna().any(axis=1).sum()
    if missing:
        print(f"⚠️  {missing} cells did not get full metadata from `{df_key}` join")
    return adata


def upper_thresh(x, nmads=3):
    return x.median() + nmads * median_abs_deviation(x)


def main():
    parser = argparse.ArgumentParser(
        description="QC/filtering on AnnData with your pre-computed metrics"
    )
    parser.add_argument(
        "--adata-file", required=True,
        help="Path to the input unfiltered .h5ad"
    )
    parser.add_argument(
        "--metadata-file", required=True,
        help="Path to your metadata CSV (index_col=0)"
    )
    parser.add_argument(
        "--output-file", required=True,
        help="Where to write the filtered .h5ad"
    )
    args = parser.parse_args()

    # load AnnData
    adata = sc.read_h5ad(args.adata_file)

    # load metadata
    metadata = pd.read_csv(args.metadata_file, index_col=0)

    # add metadata
    add_metadata_to_obs(
        adata,
        metadata,
        df_key='Sample_ID',
        obs_key='Sanger_ID'
    )

    # now run your already-written filtering/outlier code:
    samples = list(adata.obs["Sanger_ID"].unique())

    # flag mt outliers
    metric = "log1p_total_counts_mt"
    th = (
        adata
        .obs
        .groupby("Sanger_ID")[metric]
        .transform(lambda x: upper_thresh(x, nmads=3))
    )
    adata.obs["mt_high_outlier"] = adata.obs[metric] > th

    # assign QC categories
    conditions = [
        (adata.obs['total_counts_outlier'] == True),
        (adata.obs['gene_counts_outlier']  == True),
        (adata.obs['mt_outlier']          == True),
        (
            (adata.obs['gene_counts_outlier']== False)
            & (adata.obs['total_counts_outlier'] == False)
            & (adata.obs['mt_high_outlier']==False)
        ),
    ]
    values = ['Low_total_count', 'Low_nFeature', 'High_Mito', 'Pass']
    adata.obs['QC'] = np.select(conditions, values, default="Pass")

    # apply combined mask of QC=='Pass' AND any other hard filters you already have
    # (assuming you did those flag columns elsewhere)
    
    
    
    linnarson_paper_metadata= pd.read_csv('/nfs/team298/sm54/BoneAtlasProject/metadata/sample_metadata/Linnarson_Metadata_from_Paper.csv')
    # get the list of Sample_IDs that failed
    failed_ids = linnarson_paper_metadata.loc[linnarson_paper_metadata["QC"] == "Failed","Sample_ID"].unique().tolist()

    print(f"Found {len(failed_ids)} failed samples:", failed_ids)

    # if you want to subset your AnnData to just the failed samples:
    adata_failed = adata[adata.obs["Sanger_ID"].isin(failed_ids)].copy()
    print(f"adata_failed has {adata_failed.n_obs} cells across {len(failed_ids)} samples")
    adata=adata[(adata.obs['QC']=="Pass") & (adata.obs['n_genes_by_counts']>200) &(adata.obs['total_counts']<60000)& (adata.obs['total_counts']>500) ]
    # failed_ids = list of Sample_IDs with QC == "Fail"
# build mask of all cells NOT in a failed sample
#     keep_mask = ~adata.obs["Sanger_ID"].isin(failed_ids)

# # subset
#     adata_pass = adata[keep_mask]

    # write out
    adata.write_h5ad(args.output_file)
    print(f"Filtered AnnData written to {args.output_file}")
    


if __name__ == "__main__":
    main()
