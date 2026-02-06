import scanpy as sc 
import pandas as pd
pd.set_option('display.max_rows', 70)
pd.set_option('display.max_columns', 70)  # show all columns
pd.set_option('display.width', None)        # don’t wrap columns
adata= sc.read_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata_scvi.h5ad')
adata.layers['raw_counts']= adata.X.copy()
# Preprocessing
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, batch_key="10X_Chemistry")  # Use n_hvg argument

from scanpy.tl import score_genes_cell_cycle
# 1) Define your S‐phase and G2M‐phase marker lists
s_genes = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP','RAD51AP1',
    'GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2',
    'CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1',
    'CHAF1B','BRIP1','E2F8'
]
g2m_genes = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2',
    'CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2',
    'AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDC20',
    'TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA3','HMMR','TPX2','CDCA8'
]

sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=s_genes,
    g2m_genes=g2m_genes,
    copy=False
)

adata.X= adata.layers['raw_counts'].copy()
del adata.layers['raw_counts']
adata.write_h5ad('/lustre/scratch124/cellgen/haniffa/users/sm54/data/Bone_Atlas/bone_atlas_filtered_with_metadata_scvi.h5ad')

