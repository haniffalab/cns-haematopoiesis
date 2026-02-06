#!/usr/bin/env python
# coding: utf-8

"""
Stage 3: Run Scrublet on all samples in merged .h5ad in one job.
Logs to <output_dir>/logs/doublet_bulk.log.
"""

import argparse
import logging
import os
import scanpy as sc
from helper_functions import calculate_doublets

def parse_args():
    p = argparse.ArgumentParser(description="Bulk Scrublet doublet detection")
    p.add_argument("--merged_h5ad", required=True, help="Path to annadata .h5ad")
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
    fh = logging.FileHandler(os.path.join(logs,'doublet_bulk.log'))
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(fh)

    logging.info(f"Bulk Scrublet on {args.merged_h5ad}")
    calculate_doublets(
        merged_h5ad_path=args.merged_h5ad,
        out_dir=args.output_dir,
        adata_status=args.adata_status
    )
    logging.info("Bulk Scrublet complete")

if __name__=="__main__":
    main()