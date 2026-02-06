import scanpy as sc
import pandas as pd
import numpy as np
import time
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix

# Load data
# adata= sc.read_h5ad("/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/haem_compartment/bone_atlas_scvi_runs/processed_adata/adata_nlatent30_nlayers2_nhidden128_nhvg3000_leiden.h5ad")

adata= sc.read_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_best_embedding_with_hvg_label_propagation.h5ad') 

#adata = adata[~adata.obs['majority_voting_Level2'].isna()].copy()

                     

# Ensure 'published_anno' exists
if 'published_anno' not in adata.obs.columns:
    raise KeyError("Column 'published_anno' not found in adata.obs.")

# Count NA entries
na_count = adata.obs['published_anno'].isna().sum()
print(f"Total NA entries in 'published_anno': {na_count}")

if na_count == 0:
    raise ValueError("No NA entries found in 'published_anno'. Check data validity.")

# Define 'Condition' column
adata.obs['Condition'] = np.where(adata.obs['published_anno'].isna(), 'UNASSIGNED', 'REFERENCE')
adata.obs['Condition'] = adata.obs['Condition'].astype(str)


def predict_label2(merged_adata, anno_col='annotation_reference', k=50, min_score=0.1):
    """
    Predict annotation labels for query cells (Condition == UNASSIGNED) 
    using k-nearest neighbors from reference cells (Condition == REFERENCE).

    Parameters:
    -----------
    merged_adata : AnnData
        Annotated data matrix with scVI embeddings and published annotations.

    anno_col : str
        Column in .obs with reference annotations to propagate.

    k : int
        Number of neighbors to consider for each query cell.

    min_score : float
        Minimum confidence score threshold for assigning predicted labels. 
        Below this threshold, cells are labeled 'low_confidence'.

    Returns:
    --------
    pd.DataFrame
        DataFrame with index as query cell names and a single column 'propagated_anno'
        containing predicted cell type labels (or 'low_confidence').
    """
    start_time = time.time()

    if anno_col not in merged_adata.obs.columns:
        raise KeyError(f"Column '{anno_col}' not found in merged_adata.obs")

    # Extract scVI embedding
    X_emb = merged_adata.obsm["X_scVI"].astype(np.float32)

    # Identify query and reference cells
    is_query = merged_adata.obs["Condition"] == "UNASSIGNED"
    is_reference = merged_adata.obs["Condition"] == "REFERENCE"

    X_emb_ref = X_emb[is_reference, :]
    X_emb_que = X_emb[is_query, :]

    print(f"Number of reference cells: {X_emb_ref.shape[0]}")
    print(f"Number of query cells: {X_emb_que.shape[0]}")

    if X_emb_ref.shape[0] == 0 or X_emb_que.shape[0] == 0:
        raise ValueError("Either reference or query cells are empty. Check 'Condition' labels.")

    # Build KDTree and query for neighbors
    tree = cKDTree(X_emb_ref)
    k_index_ref = tree.query(X_emb_que, k=k)[1]  # neighbor indices

    print(f"k-NN search completed in {time.time() - start_time:.2f} seconds.")

    # Build sparse knn matrix
    num_query, num_ref = sum(is_query), sum(is_reference)
    knn_mat = csr_matrix((np.ones(num_query * k), k_index_ref.ravel(), np.arange(0, num_query * k + 1, k)),
                          shape=(num_query, num_ref))

    # Subset ref cells used in neighbors
    keep_ref_ixs = np.unique(k_index_ref.ravel())
    keep_ref_ixs.sort()
    small_knn_mat = knn_mat[:, keep_ref_ixs]

    # Reference annotations for dummy matrix
    annos = merged_adata[is_reference].obs[anno_col].iloc[keep_ref_ixs].copy()
    dummy_df = pd.get_dummies(annos)
    dummy_mat = dummy_df.values.astype(np.float32)

    # Probabilities of annotations
    new_anno = small_knn_mat.dot(dummy_mat)
    n_neighbors = np.array(small_knn_mat.sum(axis=1)).flatten()
    n_neighbors[n_neighbors == 0] = 1
    new_anno_prob = (new_anno.T / n_neighbors).T

    best_label = dummy_df.columns[new_anno_prob.argmax(axis=1)].values
    best_label_score = new_anno_prob.max(axis=1)

    # Apply minimum confidence threshold
    best_label = best_label.astype(str)
    best_label[best_label_score <= min_score] = "low_confidence"

    print("Label propagation completed.")

    # Create DataFrame of propagated labels for query cells only
    propagated_labels_df = pd.DataFrame({
        'propagated_anno': best_label
    }, index=merged_adata.obs_names[is_query])

    return propagated_labels_df


# Run label propagation
propagated_labels_df = predict_label2(adata, anno_col='published_anno', k=35, min_score=0.3)

# Save the propagated annotations only for query cells
propagated_labels_df.to_csv(
    "/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/haem_compartment/propagated_labels_query_cells_version3_cpu_knn_longer_time.csv"
)

print("Saved propagated labels for query cells.")
