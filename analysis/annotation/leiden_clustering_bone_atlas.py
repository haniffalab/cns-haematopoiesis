import scanpy as sc

adata= sc.read_h5ad("/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/best_integration/processed_adata/Best_Integration_bone_atlas_nlatent30_nlayers3_nhidden128_nhvg3000.h5ad")

adata.X= adata.layers['raw_counts'].copy()

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.tl.leiden(adata, resolution= 3, key_added='leiden_res_3',neighbors_key='X_scVI')

sc.tl.rank_genes_groups(adata, groupby='leiden_res_3', method='wilcoxon', corr_method='benjamini-hochberg', groups='all', reference='rest', n_genes=500, use_raw=False, log_transformed=True)


adata.write_h5ad("/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/best_integration/processed_adata/Best_Integration_bone_atlas_nlatent30_nlayers3_nhidden128_nhvg3000_with_leiden_clustering.h5ad")