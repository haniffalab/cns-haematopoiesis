"""
utils.py

Common utility functions for BoneAtlasProject.
"""

import os
import pandas as pd
import logging

def ensure_dir(path):
    """
    Ensure that a directory exists.
    If it does not, create it (including any necessary parent directories).
    """
    os.makedirs(path, exist_ok=True)

def find_h5_file(sample_dir, key='filtered'):
    """
    Locate the single .h5 file in sample_dir containing 'key'.
    Returns the full path, or raises FileNotFoundError/ValueError.
    """
    files = [f for f in os.listdir(sample_dir) 
             if key in f and f.endswith('.h5')]
    if len(files) == 0:
        raise FileNotFoundError(f"No H5 file containing '{key}' in {sample_dir}")
    if len(files) > 1:
        raise ValueError(f"Multiple H5 files containing '{key}' in {sample_dir}: {files}")
    return os.path.join(sample_dir, files[0])


def get_logger(name=__name__):
    """
    Return a logger that prints INFO-level messages to stdout.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        fmt = '%(asctime)s %(levelname)-8s %(message)s'
        handler.setFormatter(logging.Formatter(fmt))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger