Annotation set: annotation_v2 (cross-chat bundle)
-------------------------------------------------
This folder aggregates all clusters shared so far in this conversation series.
Each `cluster_XX.json` contains:
  - cluster: numeric cluster id
  - provisional_label: short human-readable label
  - top_genes: ~10 representative genes
  - notes: brief rationale/context (including user-provided hints)

Coverage (cluster → label):
  54 → Synovial fibroblast (user-specified) / chondrocyte-like ECM-high
  55 → Cycling neuroectodermal progenitor (Notch+/cell-cycle+)
  56 → Interzone-like synovial fibroblast (CREB5+/PRG4+)
  57 → PHOX2B+ sensory/autonomic neuron-like (immature)
  58 → DCX+ immature neuron / neuroblast
  59 → Neutrophil / granulocyte
  60 → Tissue macrophage (MRC1+/CD163+; LYVE1+ subset)
  61 → Schwann cell / peripheral glia
  62 → Erythroid (reticulocyte/late)
  63 → Doublets/multiplets — needs splitting
  64 → Chondrocyte (immature, COL2A1+/ACAN+/HAPLN1+)
  65 → Arachnoid fibroblast / barrier-like meningeal cell
  66 → Classical/inflammatory monocyte
  67 → Erythroid (early–mid maturation)
  69 → Meningeal/arachnoid fibroblast (axon-guidance rich)
  70 → Pre-/naive B cell (IGHM+/IGLL1+, cycling fraction)
  73 → Meningeal/arachnoid fibroblast (MEIS2+/SFRP1+)

Notes:
- Labels are provisional and reflect transcriptional programs plus user context.
- Cluster 63 was identified as a likely doublet population and split into successors (e.g., 64).
- If you share more clusters, re-run packaging to append them here.

Files:
- cluster_XX.json (per cluster)
- annotation_v2_combined.json (single-file rollup of all clusters)
