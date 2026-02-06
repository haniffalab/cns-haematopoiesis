import scanpy as sc
import celltypist
import os
import sys
import time
import numpy as np

def train_celltypist_model(adata_path, output_dir, label_column='Organ_Celltype_Age'):
    """
    Train a CellTypist model for a given AnnData file and save the model.

    Parameters:
    -----------
    adata_path : str
        Path to the input AnnData file (.h5ad).
    output_dir : str
        Directory where the trained model will be saved.
    label_column : str
        The column name in adata.obs containing cell type labels.
    """

    basename = os.path.basename(adata_path).replace('.h5ad', '')
    print(f"\n==== Processing {basename} ====")

    # Load data
    adata = sc.read_h5ad(adata_path)
    adata.X= adata.layers["raw_counts"].copy()

    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Rough training for feature selection
    print("Training rough model for feature selection...")
    t0 = time.time()
    model_fs = celltypist.train(adata, label_column, n_jobs=2, max_iter=6, use_SGD=True)
    print(f"Feature selection model trained in {time.time() - t0:.2f} seconds.")

    # Select genes
    gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -600, axis=1)[:, -600:]
    gene_index = np.unique(gene_index)
    print(f"Number of genes selected: {len(gene_index)}")

    # Final training
    print("Training final model on selected genes...")
    t0 = time.time()
    model = celltypist.train(adata[:, gene_index], label_column, check_expression=False, n_jobs=2, max_iter=100)
    print(f"Final model trained in {(time.time() - t0)/60:.2f} minutes.")

    # Save model
    os.makedirs(output_dir, exist_ok=True)
    output_model_path = os.path.join(output_dir, f"{basename}_celltypist_model.pkl")
    model.write(output_model_path)

    print(f"Model saved to {output_model_path}\n")


if __name__ == '__main__':
    # Example usage:
    # python train_celltypist_model_batch.py /path/to/adata.h5ad /path/to/output_dir

    if len(sys.argv) != 3:
        print("Usage: python train_celltypist_model_batch.py <input_adata_path> <output_dir>")
        sys.exit(1)

    adata_path = sys.argv[1]
    output_dir = sys.argv[2]

    train_celltypist_model(adata_path, output_dir)

