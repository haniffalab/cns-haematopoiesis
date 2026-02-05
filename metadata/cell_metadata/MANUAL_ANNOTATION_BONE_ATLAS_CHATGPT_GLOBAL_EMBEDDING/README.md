# Manual Annotation Results

This folder contains the outputs of manual/broad annotation for your Bone Atlas clusters, computed **using your own marker CSV (exact category names)** and your **uploaded DEGs**.

## How labels were assigned
- **DEGs input:** `DEGs_boneatlas_sig.csv.gz`
  - filtered to FDR ≤ **0.05** and logFC > **0.25** (up‑regulated only)
- Compared each cluster’s DE list to **all marker gene sets** from `FineCelltypes_compact_genesets.csv`.
- Scored each category per cluster by:
  1) **weighted** = sum of positive logFC among overlapping genes,
  2) **overlap** = number of overlapping genes,
  3) **p_enrich** = Fisher’s exact test (greater), background = union of all markers ∪ DEGs.
- Ranking priority: **weighted**, then **overlap**, then **p_enrich**.
- Assigned the top category **if overlap ≥ 2**; otherwise **Stroma**.

Parameters used are recorded in `annotation_params.json`.

## Files

- **annotation_labels.from_DEGs_sig.json**  
  Dict of `"cluster" -> "broad_label"` from your marker sets.

- **annotation_labels.table.csv**  
  Flat two-column table of the same mapping (easy to join into `adata.obs`).

- **annotation_label_counts.csv**  
  Counts of how many clusters received each label.

- **annotation_scores.from_DEGs_sig.csv.gz**  
  Full scoring matrix: every `cluster × category` with columns  
  `overlap, p_enrich, weighted`. Use this to audit or re-rank.

- **annotation_top3.from_DEGs_sig.csv**  
  The top 3 categories per cluster with their scores.

- **annotation_drivers.from_DEGs_sig.csv**  
  For the chosen label in each cluster, lists the overlapping **driver genes** (top 20).

- **BroadCategory_to_Genes.FROM_CSV.EXACT.json**  
  The **full** marker map built from your CSV (deduped per category; exact labels).

- **marker_dict_5per_broad.FROM_CSV.EXACT.json**  
  A compact panel (4–5 markers per category) chosen by frequency within your CSV.

- **annotation_params.json**  
  Provenance: inputs, thresholds, and options used to generate these outputs.

## Reuse tips

- To map labels into your `AnnData`:
  ```python
  import json, pandas as pd
  labels = json.load(open("annotation_labels.from_DEGs_sig.json"))
  df = pd.read_csv("annotation_labels.table.csv")  # if you prefer a table join
  ```

- For plotting with your compact marker panel:
  ```python
  import json
  marker_dic = json.load(open("marker_dict_5per_broad.FROM_CSV.EXACT.json"))
  # sc.pl.dotplot(adata, marker_dic, groupby="leiden_res_3")
  ```

---
Generated on: 2025-08-27 14:06:18 UTC
