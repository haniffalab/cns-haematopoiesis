#!/usr/bin/env python
# coding: utf-8

"""
Stage 2 Bulk: Annotate all samples in merged .h5ad with CellTypist in one job.
Logs to <output_dir>/logs/celltypist_bulk.log.
"""

import argparse
import logging
import os
import scanpy as sc
from helper_functions import annotate_celltypist

def parse_args():
    p = argparse.ArgumentParser(description="Bulk CellTypist annotation")
    p.add_argument("--merged_h5ad", required=True, help="Path to merged .h5ad")
    p.add_argument("--model_path", required=True, help="CellTypist model .pkl")
    p.add_argument("--output_dir", required=True, help="Where to write results")
    p.add_argument("--adata_status",  default="raw", help=["is the data raw or filtered? bone atlas or cross organ or brain"])
    p.add_argument("--log_level",  default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    p.add_argument("--log_dir",    help="Where to write logs")
    return p.parse_args()

def main():
    args = parse_args()
    logging.basicConfig(level=getattr(logging,args.log_level))
    sc.settings.verbosity = args.log_level.lower()

    os.makedirs(args.output_dir, exist_ok=True)
    logs = args.log_dir or os.path.join(args.output_dir,'logs')
    os.makedirs(logs, exist_ok=True)
    fh = logging.FileHandler(os.path.join(logs,'celltypist_bulk.log'))
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(fh)

    logging.info(f"Bulk CellTypist on {args.merged_h5ad}")
    annotate_celltypist(
        merged_h5ad_path=args.merged_h5ad,
        model_path=args.model_path,
        out_dir=args.output_dir,
        adata_status=args.adata_status
    )
    logging.info("Bulk CellTypist complete")

if __name__=="__main__":
    main()