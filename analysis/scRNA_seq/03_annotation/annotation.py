import numpy as np
import pandas as pd
from typing import Union
from scipy.stats import fisher_exact

def _rgg_clusters(r):
    names = r["names"]
    # structured array with named fields (common)
    if hasattr(names, "dtype") and getattr(names.dtype, "names", None):
        return list(names.dtype.names)
    # dict-like
    if isinstance(names, dict):
        return list(names.keys())
    raise ValueError("Unsupported format for r['names'] in rank_genes_groups")

def _rgg_pick(r, field, cluster):
    """Pick the array for a given field & cluster, robust to multiple formats."""
    if field not in r:
        return None
    arr = r[field]
    # structured array
    if hasattr(arr, "dtype") and getattr(arr.dtype, "names", None) and cluster in arr.dtype.names:
        return arr[cluster]
    # dict-like
    if isinstance(arr, dict):
        return arr.get(cluster, None)
    # fallback: list aligned to cluster order
    names = r["names"]
    if hasattr(names, "dtype") and getattr(names.dtype, "names", None):
        clusts = list(names.dtype.names)
        if cluster in clusts:
            idx = clusts.index(cluster)
            try:
                return arr[idx]
            except Exception:
                return None
    return None



def annotate_broad_topN(
    adata,
    CATEGORY_MARKERS,
    cluster_key="leiden_res_3",
    deg_key="rank_genes_groups",
    top_n=100,
    min_lfc=0.0,
    max_padj=1.0,
    min_overlap=2,
    prefer_weighted=True,
    up_only=True,
    background="markers",  # "markers" (union of marker genes) or "var" (adata.var_names) 
):
    """
    - Extract top-N DEGs per cluster.
    - Score overlap/enrichment vs CATEGORY_MARKERS.
    - Pick best category per cluster and write adata.obs['broad_anno'].
    - Returns (labels_dict, scores_df)
    """
    deg = extract_topN_from_rgg(
        adata, deg_key=deg_key, top_n=top_n,
        up_only=up_only, min_lfc=min_lfc, max_padj=max_padj
    )

    # per-cluster gene sets + logFC maps
    cluster_sets, lfc_maps = {}, {}
    for c, sub in deg.groupby("cluster"):
        gs = set(sub["gene_norm"])
        cluster_sets[c] = gs
        if "logfoldchanges" in sub and sub["logfoldchanges"].notna().any():
            lfc_maps[c] = {g: float(sub.loc[sub["gene_norm"]==g, "logfoldchanges"].iloc[0]) for g in gs}
        else:
            lfc_maps[c] = {g: 0.0 for g in gs}

    # background universe for enrichment
    if background == "var":
        bg_all = set(map(str.upper, adata.var_names))
    else:
        bg_all = set().union(*CATEGORY_MARKERS.values()) if CATEGORY_MARKERS else set()

    def score_cluster(deg_genes, lfc_map, background_genes=None):
        deg_set = set(deg_genes)
        bg = set(background_genes) if background_genes is not None else (bg_all | deg_set)
        M = max(len(bg), 1)

        out = {}
        for cat, markers in CATEGORY_MARKERS.items():
            marks = set(markers) & bg
            a = len(deg_set & marks); b = len(deg_set) - a
            c = len(marks)   - a;     d = M - a - b - c
            _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
            w = sum(max(0.0, lfc_map.get(g, 0.0)) for g in (deg_set & marks))
            out[cat] = {"overlap": a, "weighted": w, "p_enrich": p}
        return out

    rows, chosen = [], {}
    for c, gs in cluster_sets.items():
        scores = score_cluster(gs, lfc_maps[c])
        sdf = pd.DataFrame(scores).T.reset_index().rename(columns={"index":"category"})
        sdf = sdf.sort_values(
            by=["weighted","overlap","p_enrich"] if prefer_weighted else ["overlap","weighted","p_enrich"],
            ascending=[False, False, True]
        )
        top = sdf.iloc[0]
        label = str(top["category"]) if int(top["overlap"]) >= min_overlap else "Stroma"
        chosen[c] = label
        sdf["cluster"] = c
        rows.append(sdf)

    score_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["category","cluster"])

    # Attach to obs
    adata.obs["broad_anno"] = (
        adata.obs[cluster_key].astype(str).map(chosen).fillna("Undefined").astype("category")
    )

    # Keep extras
    adata.uns["broad_anno_scores"] = score_df
    adata.uns["broad_anno_by_cluster_auto"] = {str(k): str(v) for k, v in chosen.items()}

    return chosen, score_df


def extract_topN_from_rgg(
    adata,
    deg_key="rank_genes_groups",
    top_n=100,
    up_only=True,
    min_lfc=0.0,
    max_padj=1.0,
    sort_preference=("scores","logfoldchanges","pvals_adj")
) -> pd.DataFrame:
    """
    Returns a tidy DF with columns:
      cluster, gene, logfoldchanges, scores, pvals_adj, gene_norm
    keeping the top-N per cluster by the first available metric in sort_preference.
    """
    if deg_key not in adata.uns:
        raise KeyError(f"{deg_key} not in adata.uns")

    r = adata.uns[deg_key]
    clusters = _rgg_clusters(r)
    out = []

    for c in clusters:
        names = _rgg_pick(r, "names", c)
        if names is None:
            continue
        lfc   = _rgg_pick(r, "logfoldchanges", c)
        scr   = _rgg_pick(r, "scores", c)
        padj  = _rgg_pick(r, "pvals_adj", c)

        n = len(names)
        def _as_list(x, fill=np.nan):
            if x is None:
                return [fill]*n
            try:
                xx = list(x)
            except Exception:
                xx = [x]*n
            if len(xx) < n:
                xx += [fill]*(n - len(xx))
            return xx

        df = pd.DataFrame({
            "cluster": [c]*n,
            "gene": [str(g) for g in names],
            "logfoldchanges": _as_list(lfc, np.nan),
            "scores": _as_list(scr, np.nan),
            "pvals_adj": _as_list(padj, np.nan),
        })

        # Significance + up-reg filters (optional but recommended)
        if max_padj < 1.0 and df["pvals_adj"].notna().any():
            df = df[df["pvals_adj"].fillna(1.0) <= max_padj]
        if up_only and df["logfoldchanges"].notna().any():
            df = df[df["logfoldchanges"].fillna(0) > min_lfc]

        # Ranking
        chosen = None
        for field in sort_preference:
            if field in df and df[field].notna().any():
                df = df.sort_values(field, ascending=(field == "pvals_adj"))
                chosen = field
                break

        out.append(df.head(top_n))

    res = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["cluster","gene"])
    res["gene_norm"] = res["gene"].str.upper().str.strip()
    return res




def top_genes_for_cluster(
    adata,
    cluster,
    deg_key: str = "rank_genes_groups",
    top_n: int = 100,
    up_only: bool = True,
    min_lfc: float = 0.0,
    max_padj: float = 1.0,
    sort_preference = ("scores", "logfoldchanges", "pvals_adj"),
    out_csv: str = None,
    return_genes_only: bool = False,
) -> Union[pd.DataFrame, list]: 
    """
    Return the top-N differentially expressed genes for one cluster.

    Parameters
    ----------
    adata : AnnData
    cluster : str|int
        Cluster name as used in adata.uns[deg_key]['names'] (e.g. "12", "B", 0).
    deg_key : str
        Key in adata.uns with rank_genes_groups result.
    top_n : int
        Number of genes to return (after filtering).
    up_only : bool
        If True and logfoldchanges exists, keep genes with logFC > min_lfc.
    min_lfc : float
        Minimum log fold-change if up_only=True.
    max_padj : float
        Keep genes with adjusted p-value <= max_padj (if pvals_adj exists).
    sort_preference : tuple
        Order of metrics to sort by. First available column is used.
        'scores' and 'logfoldchanges' are sorted DESC; 'pvals_adj' ASC.
    out_csv : str|None
        If provided, write the table to this path.
    return_genes_only : bool
        If True, return just the list of gene names.

    Returns
    -------
    pd.DataFrame or list
        DataFrame with columns: cluster, gene, logfoldchanges, scores, pvals_adj, gene_norm
        or list of gene names if return_genes_only=True.
    """
    if deg_key not in adata.uns:
        raise KeyError(f"{deg_key!r} not found in adata.uns")

    r = adata.uns[deg_key]

    # --- helpers to handle different Scanpy storage formats ---
    def _clusters(rr):
        names = rr["names"]
        if hasattr(names, "dtype") and getattr(names.dtype, "names", None):
            return list(names.dtype.names)                 # structured array
        if isinstance(names, dict):
            return list(names.keys())                      # dict-like
        raise ValueError("Unsupported format for r['names'] in rank_genes_groups")

    def _pick(rr, field, c):
        if field not in rr:
            return None
        arr = rr[field]
        # structured array with named fields
        if hasattr(arr, "dtype") and getattr(arr.dtype, "names", None):
            if c in arr.dtype.names:
                return arr[c]
        # dict-like
        if isinstance(arr, dict):
            return arr.get(c, None)
        # fallback: list aligned to cluster order
        names = rr["names"]
        if hasattr(names, "dtype") and getattr(names.dtype, "names", None):
            cl = list(names.dtype.names)
            if c in cl:
                idx = cl.index(c)
                try:
                    return arr[idx]
                except Exception:
                    return None
        return None

    # --- find the cluster exactly as stored in rgg ---
    clist = _clusters(r)
    c_str = str(cluster)
    if c_str not in clist:
        # sometimes cluster labels are not strings; try exact match first
        if cluster in clist:
            c_key = cluster
        else:
            raise KeyError(f"Cluster {cluster!r} not found. Available: {clist[:10]}{'...' if len(clist)>10 else ''}")
    else:
        c_key = c_str

    names = _pick(r, "names", c_key)
    if names is None:
        raise KeyError(f"No 'names' found for cluster {c_key!r}")

    lfc  = _pick(r, "logfoldchanges", c_key)
    scr  = _pick(r, "scores", c_key)
    padj = _pick(r, "pvals_adj", c_key)

    n = len(names)
    def _as_list(x, fill=np.nan):
        if x is None: return [fill]*n
        try: xx = list(x)
        except Exception: xx = [x]*n
        if len(xx) < n: xx += [fill]*(n - len(xx))
        return xx

    df = pd.DataFrame({
        "cluster": [str(c_key)] * n,
        "gene": [str(g) for g in names],
        "logfoldchanges": _as_list(lfc, np.nan),
        "scores": _as_list(scr, np.nan),
        "pvals_adj": _as_list(padj, np.nan),
    })

    # --- filters ---
    if max_padj < 1.0 and "pvals_adj" in df and df["pvals_adj"].notna().any():
        df = df[df["pvals_adj"].fillna(1.0) <= max_padj]
    if up_only and "logfoldchanges" in df and df["logfoldchanges"].notna().any():
        df = df[df["logfoldchanges"].fillna(0) > min_lfc]

    # --- ranking by first usable metric ---
    for field in sort_preference:
        if field in df and df[field].notna().any():
            ascending = (field == "pvals_adj")
            df = df.sort_values(field, ascending=ascending)
            break

    df = df.head(top_n).copy()
    df["gene_norm"] = df["gene"].str.upper().str.strip()

    if out_csv:
        df.to_csv(out_csv, index=False)

    if return_genes_only:
        return df["gene"].tolist()

    return df


def attach_broad_annotations(adata, mapping, cluster_key="leiden_res_3", out_key="broad_anno"):
    lab_map = {str(k): str(v) for k, v in mapping.items()}
    adata.obs[out_key] = (
        adata.obs[cluster_key].astype(str).map(lab_map).fillna("Unassigned").astype("category")
    )
    adata.uns[f"{out_key}_by_cluster_manual"] = lab_map
    # Optional: compact table of cluster→label
    tbl = (adata.obs[[cluster_key, out_key]]
           .drop_duplicates()
           .sort_values(by=cluster_key, key=lambda s: pd.to_numeric(s, errors="ignore")))


