#!/usr/bin/env python
# coding: utf-8

"""
Stage 1: Read metadata, select samples, generate per-sample .h5ad and merged AnnData.
Logs to <output_dir>/logs/generate_h5ad.log.
"""

import argparse
import logging
import os
import pandas as pd
import scanpy as sc
from helper_functions import generate_h5ads

def parse_args():
    p = argparse.ArgumentParser(description="Stage 1: Generate and merge h5ads")
    p.add_argument("--metadata",   required=True, help="CSV metadata (index_col=0)")
    p.add_argument("--dataset",    default="Bone_Atlas", help="Dataset to select")
    p.add_argument("--raw_dir",    required=True, help="Dir of raw .h5 outputs")
    p.add_argument("--output_dir", required=True, help="Where to write h5ads")
    p.add_argument("--merged_name", default="bone_atlas_unfiltered.h5ad",
                   help="Merged AnnData filename")
    p.add_argument("--log_level",  default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    p.add_argument("--log_dir",    help="Where to write logs (default: output_dir/logs)")
    return p.parse_args()

def main():
    args = parse_args()
    logging.basicConfig(level=getattr(logging,args.log_level))
    sc.settings.verbosity = args.log_level.lower()

    os.makedirs(args.output_dir, exist_ok=True)
    logs = args.log_dir or os.path.join(args.output_dir, 'logs')
    os.makedirs(logs, exist_ok=True)
    fh = logging.FileHandler(os.path.join(logs,'generate_h5ad.log'))
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(fh)
    logging.info(f"Writing logs to {logs}/generate_h5ad.log")

    meta = pd.read_csv(args.metadata, index_col=0)
    #samples = meta.loc[meta.Dataset==args.dataset,'Sample_ID'].tolist()
     # pick only the rows in our dataset, but keep both Sample_ID and alias
    samp_df = meta.loc[meta.Dataset == args.dataset, ['Sample_ID','alias']]
    samples = []
    for sample_id, alias in samp_df.itertuples(index=False):
        dir1 = os.path.join(args.raw_dir, sample_id)
        dir2 = os.path.join(args.raw_dir, alias) if pd.notna(alias) else None

        if os.path.isdir(dir1):
            samples.append(sample_id)
        elif dir2 and os.path.isdir(dir2):
            samples.append(alias)
            logging.info(f"Using alias for {sample_id}: {alias}/ exists")
        else:
            logging.warning(f"Skipping {sample_id}: neither {sample_id}/ nor {alias}/ exist")

    logging.info(f"After checking raw_dir, {len(samples)} samples will be processed")

    from collections import Counter
    counts = Counter(samples)
    samples = [s for s in samples if counts[s]==1]
    logging.info(f"Found {len(samples)} unique samples")

    generate_h5ads(
        sample_list=samples,
        raw_dir=args.raw_dir,
        out_dir=args.output_dir,
        merged_fname=args.merged_name
    )
    logging.info(f"Finished processing {len(samples)} samples")

if __name__=="__main__":
    main()
