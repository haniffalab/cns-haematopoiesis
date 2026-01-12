#!/usr/bin/env python3

import scanpy as sc
import liana as li
import pandas as pd
import numpy as np
import plotnine as p9
import os
from liana.method import cellphonedb

# ------------------------------------------------
# MAKE PLOTS DIRECTORY
# ------------------------------------------------
os.makedirs("plots", exist_ok=True)

# ------------------------------------------------
# LOAD DATA WITH RAW COUNTS
# ------------------------------------------------
adata = sc.read_h5ad(
    "/nfs/team298/sm54/BoneAtlasProject/data/visium/adata_spine_9_10_for_cell2location_with_raw_counts.h5ad"
)
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

adata.raw = adata.copy()

# Ensure raw exists
if adata.raw is None:
    adata.raw = adata.copy()


###############################################################
# 1) RANK AGGREGATE — HIGH RESOLUTION COLUMN
###############################################################
print("Running Rank Aggregate for HIGH mapping...")

li.mt.rank_aggregate(
    adata,
    groupby="cell2location_mapping_high_V1",
    resource_name="consensus",
    expr_prop=0.1,
    verbose=True,
    key_added="liana_res",
    use_raw=True
)

# Save results internally and externally
adata.uns["liana_res_highV1"] = adata.uns["liana_res"].copy()
adata.uns["liana_res_highV1"].to_csv("liana_spine_highV1_results.csv", index=False)

# DOTPLOT — HIGH
p = li.pl.dotplot(
    adata=adata,
    colour="magnitude_rank",
    size="specificity_rank",
    inverse_size=True,
    inverse_colour=True,
    top_n=15,
    orderby="magnitude_rank",
    orderby_ascending=True,
    uns_key="liana_res",
    figure_size=(8, 7)
)
p.save("plots/liana_spine_highV1_dotplot.pdf")


###############################################################
# 2) RANK AGGREGATE — FINE RESOLUTION COLUMN
###############################################################
print("Running Rank Aggregate for FINE mapping...")

li.mt.rank_aggregate(
    adata,
    groupby="cell2location_mapping_fine_V1",
    resource_name="consensus",
    expr_prop=0.1,
    verbose=True,
    key_added="liana_res",
    use_raw=True
)

adata.uns["liana_res_fineV1"] = adata.uns["liana_res"].copy()
adata.uns["liana_res_fineV1"].to_csv("liana_spine_fineV1_results.csv", index=False)

# DOTPLOT — FINE
p = li.pl.dotplot(
    adata=adata,
    colour="magnitude_rank",
    size="specificity_rank",
    inverse_size=True,
    inverse_colour=True,
    top_n=15,
    orderby="magnitude_rank",
    orderby_ascending=True,
    uns_key="liana_res",
    figure_size=(8, 7)
)
p.save("plots/liana_spine_fineV1_dotplot.pdf")


###############################################################
# 3) CELLPHONEDB — HIGH RESOLUTION (FOR TILEPLOTS)
###############################################################
print("Running CellPhoneDB for HIGH mapping (slow)...")

cellphonedb(
    adata,
    groupby="cell2location_mapping_high_V1",
    resource_name="consensus",
    expr_prop=0.1,
    key_added="cpdb_highV1",
    verbose=True
)

adata.uns["cpdb_highV1"].to_csv("liana_spine_highV1_cpdb_results.csv")

# TILEPLOT — HIGH
tp = li.pl.tileplot(
    adata=adata,
    fill="means",
    label="props",
    label_fun=lambda x: f"{x:.2f}",
    top_n=10,
    orderby="cellphone_pvals",
    orderby_ascending=True,
    uns_key="cpdb_highV1",
    source_labels=["FIBRO_1", "FIBRO_2", "MACROPHAGE (LYVE1_HIGH)"],
    target_labels=["CHON_RESTING", "CHON_PROGENITOR"],
    figure_size=(10, 8)
)
tp.save("plots/liana_spine_highV1_tileplot.pdf")


###############################################################
# 4) CELLPHONEDB — FINE RESOLUTION (FOR TILEPLOTS)
###############################################################
print("Running CellPhoneDB for FINE mapping (slow)...")

cellphonedb(
    adata,
    groupby="cell2location_mapping_fine_V1",
    resource_name="consensus",
    expr_prop=0.1,
    key_added="cpdb_fineV1",
    verbose=True
)

adata.uns["cpdb_fineV1"].to_csv("liana_spine_fineV1_cpdb_results.csv")

# TILEPLOT — FINE
tp = li.pl.tileplot(
    adata=adata,
    fill="means",
    label="props",
    label_fun=lambda x: f"{x:.2f}",
    top_n=10,
    orderby="cellphone_pvals",
    orderby_ascending=True,
    uns_key="cpdb_fineV1",
    source_labels=["Chondrocyte", "Fibroblast", "MACROPHAGE (LYVE1_HIGH)"],
    target_labels=["Osteogenic", "Fibroblast"],
    figure_size=(10, 8)
)
tp.save("plots/liana_spine_fineV1_tileplot.pdf")


###############################################################
# SAVE UPDATED ANNDATA
###############################################################
adata.write_h5ad("adata_spine_liana_processed.h5ad")

print("All LIANA analyses (RankAggregate + CPDB + Dotplots + Tileplots) completed successfully.")
