import os
import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation
import pandas as pd
import celltypist
import scrublet as scr
import scipy.stats
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
from utils import ensure_dir
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as smm
# Utility Functions

def is_outlier(adata, metric: str, nmads: int):
    """
    Identify outliers in an AnnData.obs metric using Median Absolute Deviation (MAD).
    """
    M = adata.obs[metric]
    mad = median_abs_deviation(M)
    return (M < np.median(M) - nmads * mad) | (M > np.median(M) + nmads * mad)
    
def is_upper_outlier(adata, metric: str, nmads: int):
    """
    Flag cells whose obs[metric] is more than nmads*MAD above the median.
    """
    M   = adata.obs[metric]
    mad = median_abs_deviation(M)
    thresh = np.median(M) + nmads * mad
    return M > thresh


def bh_correct_using_statsmodels(pvalues):
    """
    Perform Benjamini–Hochberg FDR correction on a list/array of p-values.
    Returns the array of adjusted p-values (q-values).
    """
    _, qvals, _, _ = smm.multipletests(pvalues, method='fdr_bh')
    return qvals

# 1) Data preparation and QC




def generate_h5ads(
    sample_list,
    raw_dir,
    out_dir,
    sample_key: str = 'filtered',
    compression: str = 'gzip',
    merged_fname: str = 'merged.h5ad',
):
    """
    Reads 10x HDF5s, annotates QC metrics, writes per-sample and merged AnnData.
    Skips any sample whose directory or .h5 file isn’t found.

    Returns
    -------
    per_sample_files : dict
        Mapping Sanger_ID -> path to its .h5ad.
    merged_path : str
        Path to the merged .h5ad file.
    """
    ensure_dir(out_dir)
    per_sample_files = {}
    adatas = []
    processed = []

    for sanger in sample_list:
        samp_raw = os.path.join(raw_dir, sanger)
        if not os.path.isdir(samp_raw):
            print(f"Warning: directory not found for sample {sanger!r}; skipping.")
            continue

        files = [f for f in os.listdir(samp_raw)
                 if sample_key in f and f.endswith('.h5')]
        if len(files) != 1:
            print(f"Warning: {sanger!r} expected one .h5, found {len(files)} ({files}); skipping.")
            continue

        path = os.path.join(samp_raw, files[0])
        try:
            ad = sc.read_10x_h5(path)
        except FileNotFoundError:
            print(f"Warning: file {path!r} not found; skipping sample {sanger}.")
            continue

        # Annotate barcodes and IDs
        ad.obs['barcode_raw'] = ad.obs_names.astype(str)
        ad.obs['Sanger_ID']   = sanger
        ad.obs['Cell_ID']     = ad.obs['barcode_raw'] + '-' + sanger
        ad.obs_names = ad.obs['Cell_ID']

        # Gene metadata
        ad.var['SYMBOL'] = ad.var_names
        ad.var_names_make_unique()
        ad.var['mt'] = ad.var_names.str.startswith('MT-')
        ad.var['rb'] = ad.var_names.str.startswith(('RPS','RPL'))

        # QC metrics and outliers
        sc.pp.calculate_qc_metrics(ad, qc_vars=['mt','rb'], inplace=True, log1p=True)
        ad.obs['total_counts_outlier'] = is_outlier(ad, 'log1p_total_counts', 5)
        ad.obs['gene_counts_outlier']  = is_outlier(ad, 'log1p_n_genes_by_counts', 5)
        ad.obs['mt_outlier']           = is_upper_outlier(ad, 'log1p_total_counts_mt', 3)

        # Write per-sample
        sample_out = os.path.join(out_dir, sanger)
        ensure_dir(sample_out)
        sample_path = os.path.join(sample_out, f"{sanger}.h5ad")
        ad.write_h5ad(sample_path, compression=compression)

        per_sample_files[sanger] = sample_path
        adatas.append(ad)
        processed.append(sanger)

    if not adatas:
        raise RuntimeError("No samples were processed — check your sample_list and directories.")

    # Merge all processed and write
    merged = sc.concat(
        adatas,
        join='outer',
        label='Sanger_ID',
        keys=processed,
        index_unique=None,
    )
    merged_path = os.path.join(out_dir, merged_fname)
    merged.write_h5ad(merged_path, compression=compression)

    return per_sample_files, merged_path


# 2) CellTypist annotation

import os
import logging
import pandas as pd
import scanpy as sc
import celltypist
import scrublet as scr
import scipy.stats
import statsmodels.stats.multitest as smm
from utils import ensure_dir
# … other imports …

def annotate_celltypist(
    merged_h5ad_path: str,
    model_path: str,
    out_dir: str,
    adata_status: str,
     target_sum: float = 1e4,
):
    """
    Runs CellTypist per sample, saves per-sample and global predictions.
    merged_h5ad_path: path to the anndata,
    model_path: path to celltypist model,
    out_dir: path to directory where the output should be saved,
    target_sum: float = 1e4,
    adata_status: is the adata filtered or not, options: raw/filtered 
    """
    ensure_dir(out_dir)
    adata = sc.read_h5ad(merged_h5ad_path)

    # Load model
    model_name = os.path.splitext(os.path.basename(model_path))[0]
    model = celltypist.models.Model.load(model=model_path)

    all_prob_dfs = []
    all_pred_dfs = []
    processed, skipped = [], []

    # Precompute cell counts per sample for logging
    counts = adata.obs['Sanger_ID'].value_counts().to_dict()

    for sanger in adata.obs['Sanger_ID'].unique():
        samp_out = os.path.join(out_dir, sanger)
        ensure_dir(samp_out)

        pred_csv = os.path.join(samp_out, f"predicted_labels_{model_name}_{sanger}_{adata_status}.csv")
        prob_csv = os.path.join(samp_out, f"probability_matrix_{model_name}_{sanger}_{adata_status}.csv")

        # PROCESS  sample
        logging.info(f"[CellTypist] Processing {sanger} ({counts[sanger]} cells)")
        subset = adata[adata.obs['Sanger_ID'] == sanger].copy()
        sc.pp.normalize_total(subset, target_sum=target_sum)
        sc.pp.log1p(subset)

        result = celltypist.annotate(
            subset,
            model=model,
            majority_voting=True,
        )

        # DataFrames
        pred_df = result.predicted_labels.copy()
        prob_df = result.probability_matrix.copy()

        # Save per-sample outputs
        pred_df.to_csv(pred_csv)
        result.decision_matrix.to_csv(os.path.join(
            samp_out, f"decision_matrix_{model_name}_{sanger}_{adata_status}.csv"))
        prob_df.to_csv(prob_csv)
        result.adata.write_h5ad(os.path.join(
            samp_out, f"{model_name}_{sanger}_{adata_status}_celltypist_annotated.h5ad"))

        # Collect for global
        prob_df.index.name = 'Cell_ID'
        all_prob_dfs.append(prob_df)

        if 'Cell_ID' not in pred_df.columns:
            tmp = pred_df.rename_axis('Cell_ID').reset_index()
        else:
            tmp = pred_df.copy()
        tmp['Sanger_ID'] = sanger
        all_pred_dfs.append(tmp)

        processed.append(sanger)

    # Global concatenation & save
    big_probs = pd.concat(all_prob_dfs, axis=0)
    big_probs.to_csv(os.path.join(out_dir, 'all_samples_probability_matrix.csv'))

    big_preds = pd.concat(all_pred_dfs, axis=0)
    cols = ['Cell_ID','Sanger_ID'] + [c for c in big_preds.columns
                                     if c not in ('Cell_ID','Sanger_ID')]
    big_preds[cols].to_csv(os.path.join(out_dir, f"all_samples_predictions_{adata_status}.csv"),
                           index=False)

    # Summary logging
    total = len(adata.obs['Sanger_ID'].unique())
    logging.info(f"[CellTypist] Done: {len(processed)} processed, {len(skipped)} skipped, out of {total} samples.")
    for s, n in counts.items():
        logging.info(f"  Sample {s}: {n} cells")

    return big_probs, big_preds



def calculate_doublets(
    merged_h5ad_path,
    out_dir: str,
    adata_status: str,
    doublet_rate_factor: float = 0.008,
    sim_doublet_ratio: float = 2,
    q_threshold: float = 0.1,
):
    """
    1) Load merged AnnData (either a file path or an AnnData object).
    2) For each sample in adata.obs['Sanger_ID']:
       - run Scrublet → per-cell scrublet_score & default predicted_doublets
       - save histogram & UMAP plots
       - run Scanpy preprocessing, PCA, neighbors, Leiden clustering
       - compute cluster-level median scrublet_score
       - fit a median-centred, MAD-variance normal null on cluster medians
         (MAD computed only on values above the median to avoid zero-truncation)
       - compute raw p-value per cluster, BH-correct those *cluster* p-values
       - map BH-adjusted q back to cells as bh_pval
       - add a second call `predicted_doublets_pval = (bh_pval < q_threshold)`
       - save per-sample CSV with columns:
           ['scrublet_score','predicted_doublets',
            'scrublet_cluster_score','bh_pval','predicted_doublets_pval']
    3) Save combined CSV and failed_runs.log in out_dir.

    Reference
    ---------
    Supplementary Section S6 of “A step-by-step workflow for low-level analysis
    of single-cell RNA-seq data with Bioconductor” (F1000Res 2016;5:2122),
    PMCID:PMC6522369.

    Parameters
    ----------
    merged_h5ad_path : str or AnnData
        Path to merged .h5ad file, or an AnnData object.
    out_dir : str
        Top-level output directory.
    adata_status : str
        Tag for the AnnData state (e.g. "raw" or "filtered").
    doublet_rate_factor : float
        Multiplier for expected doublet rate = (n_cells/1000)*factor.
    sim_doublet_ratio : float
        Ratio of simulated doublets per observed cell for Scrublet.
    q_threshold : float
        BH-adjusted p-value cutoff for calling a cluster doublet-rich.

    Returns
    -------
    combined : pd.DataFrame
        Concatenated per-cell tables for all samples.
    failed : list
        Sample IDs that failed processing.
    """
    import matplotlib.pyplot as plt

    # 1) load AnnData
    if isinstance(merged_h5ad_path, str):
        adata = sc.read_h5ad(merged_h5ad_path)
    else:
        adata = merged_h5ad_path.copy()

    sc.settings.verbosity = 1
    os.makedirs(out_dir, exist_ok=True)

    scorenames = [
        'scrublet_score',
        'predicted_doublets',
        'scrublet_cluster_score',
        'bh_pval',
        'predicted_doublets_pval'
    ]
    counts = adata.obs['Sanger_ID'].value_counts().to_dict()
    processed, failed, all_results = [], [], []

    for sample in adata.obs['Sanger_ID'].unique():
        sample_dir = os.path.join(out_dir, sample)
        os.makedirs(sample_dir, exist_ok=True)
        plots_dir = os.path.join(sample_dir, f"scrublet-plots_{adata_status}")
        os.makedirs(plots_dir, exist_ok=True)
        scrub_csv = os.path.join(
            sample_dir,
            f"doublet_score_{sample}_{adata_status}.csv"
        )

        try:
            logging.info(f"[Doublet] Processing {sample} ({counts[sample]} cells)")
            adata_sample = adata[adata.obs['Sanger_ID']==sample].copy()
            adata_sample.layers['raw_counts'] = adata_sample.X.copy()

            # Run Scrublet
            dbl_rate = adata_sample.n_obs/1000 * doublet_rate_factor
            scrub = scr.Scrublet(
                adata_sample.X,
                expected_doublet_rate=dbl_rate,
                sim_doublet_ratio=sim_doublet_ratio
            )
            doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
            adata_sample.obs['scrublet_score']     = doublet_scores
            adata_sample.obs['predicted_doublets'] = predicted_doublets

            # Plot histogram
            try:
                scrub.plot_histogram()
                plt.savefig(os.path.join(
                    plots_dir, f'doublet_score_histogram_{adata_status}.png'))
                plt.close()
            except Exception as e:
                print(f"Warning: Could not plot histogram for {sample}: {e}")

            # UMAP embedding
            try:
                emb = scr.get_umap(scrub.manifold_obs_, n_neighbors=10, min_dist=0.3)
                scrub.set_embedding('UMAP', emb)
                scrub.plot_embedding('UMAP', order_points=True)
                plt.savefig(os.path.join(
                    plots_dir, f'UMAP_{adata_status}.png'))
                plt.close()
            except Exception as e:
                print(f"Warning: Could not plot UMAP for {sample}: {e}")

            # Scanpy preprocessing & Leiden
            sc.pp.filter_genes(adata_sample, min_cells=3)
            adata_sample.X = adata_sample.X.astype('float32')
            sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
            sc.pp.log1p(adata_sample)
            sc.pp.highly_variable_genes(
                adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5
            )
            adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
            adata_sample.X = adata_sample.layers['raw_counts'].copy()
            sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
            sc.pp.log1p(adata_sample)
            sc.tl.pca(adata_sample, svd_solver='arpack')
            sc.pp.neighbors(adata_sample)
            sc.tl.leiden(adata_sample)

            # Over-cluster
            for clus in np.unique(adata_sample.obs['leiden']):
                sc.tl.leiden(adata_sample, restrict_to=('leiden',[clus]))
                adata_sample.obs['leiden'] = adata_sample.obs['leiden_R']

            # Compute cluster medians
            clusters = np.unique(adata_sample.obs['leiden'])
            cluster_meds = {
                cl: np.median(
                    adata_sample.obs.loc[
                        adata_sample.obs['leiden']==cl,
                        'scrublet_score'
                    ]
                ) for cl in clusters
            }
            adata_sample.obs['scrublet_cluster_score'] = \
                adata_sample.obs['leiden'].map(cluster_meds)

            # Null on >median cluster medians
            all_meds = np.array(list(cluster_meds.values()))
            med = np.median(all_meds)
            mad = np.median(all_meds[all_meds>med] - med)
            scale = 1.4826 * mad
            raw_p = {
                cl: 1 - scipy.stats.norm.cdf(cluster_meds[cl], loc=med, scale=scale)
                for cl in clusters
            }

            # BH-correct cluster p-values
            clusts, pvals = zip(*raw_p.items())
            qvals = bh_correct_using_statsmodels(np.array(pvals))
            bh_dict = dict(zip(clusts, qvals))
            adata_sample.obs['bh_pval'] = \
                adata_sample.obs['leiden'].map(bh_dict)

            # New p-val–based call
            adata_sample.obs['predicted_doublets_pval'] = \
                adata_sample.obs['bh_pval'] < q_threshold

            # Save per-sample CSV
            scrublet_df = pd.DataFrame(
                index=adata_sample.obs_names,
                columns=scorenames
            )
            for col in scorenames:
                scrublet_df[col] = adata_sample.obs[col]
            scrublet_df.to_csv(scrub_csv)

            all_results.append(scrublet_df)
            processed.append(sample)
            logging.info(f"[Doublet] Finished {sample}")

        except Exception as e:
            logging.error(f"[Doublet] Error processing {sample}: {e}")
            failed.append(sample)

    # Combine all results
    combined = pd.concat(all_results) if all_results else pd.DataFrame()
    combined.to_csv(os.path.join(
        out_dir, f'combined_scrublet_scores_{adata_status}.csv'
    ))

    if failed:
        with open(os.path.join(out_dir,'failed_runs.log'),'w') as f:
            f.write("\n".join(failed))
        logging.info(f"[Doublet] Logged {len(failed)} failures")

    return combined, failed
