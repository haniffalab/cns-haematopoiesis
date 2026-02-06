# helper_function.py
# Minimal helpers to save/load ONLY .obs for one AnnData at a time,
# and to fully REPLACE existing obs on restore.

from __future__ import annotations
import os
import pandas as pd

def save_obs(adata, out_path: str, index_name: str = "obs_name") -> str:
    """
    Save adata.obs (+ obs_names) to a single file.
    Supports .parquet (recommended) or .csv.gz.
    """
    df = adata.obs.copy()

    # Insert obs_names and avoid writing the pandas index to Parquet
    df.insert(0, index_name, adata.obs_names.astype(str).to_numpy())
    df = df.reset_index(drop=True)

    # Workaround for pyarrow JSON metadata issue (np.bool_ etc.)
    try:
        df.attrs.clear()
    except Exception:
        df.attrs = {}

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    if out_path.endswith(".parquet"):
        df.to_parquet(out_path, compression="zstd", index=False)  # <- index=False
        return out_path
    elif out_path.endswith(".csv.gz"):
        df.to_csv(out_path, index=False)
        return out_path
    else:
        out_path = out_path + ".parquet"
        df.to_parquet(out_path, compression="zstd", index=False)
        return out_path


def load_obs(path: str, index_name: str = "obs_name") -> pd.DataFrame:
    """
    Load an obs table saved by save_obs(); returns a DataFrame indexed by obs_name.
    """
    if path.endswith(".parquet"):
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)
    if index_name not in df.columns:
        raise KeyError(f"Expected '{index_name}' column in {path}.")
    return df.set_index(index_name)

def replace_obs(adata, new_obs_df: pd.DataFrame, strict: bool = True) -> None:
    """
    Replace adata.obs entirely with new_obs_df, aligning by obs_names.
    - strict=True: require exact 1:1 index match; raise if mismatch.
    """
    idx_ad  = pd.Index(adata.obs_names.astype(str))
    idx_new = pd.Index(new_obs_df.index.astype(str))
    missing = idx_ad.difference(idx_new)
    extras  = idx_new.difference(idx_ad)
    if strict and (len(missing) or len(extras)):
        raise ValueError(
            f"Index mismatch: {len(missing)} missing and {len(extras)} extra rows. "
            f"Set strict=False to allow alignment (extras dropped; missing become NaN)."
        )
    adata.obs = new_obs_df.reindex(idx_ad)
    print(f"[replace_obs] adata.obs replaced with {adata.obs.shape[1]} cols for {adata.n_obs} rows.")

def restore_obs(adata, path: str, strict: bool = True, index_name: str = "obs_name") -> None:
    """
    Convenience: load_obs(path) + replace_obs(adata, ...).
    """
    df = load_obs(path, index_name=index_name)
    replace_obs(adata, df, strict=strict)