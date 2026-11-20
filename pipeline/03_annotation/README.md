# Embedding, clustering and labels

Clusters and labels are stored in the cohort object (`obsm["X_umap"]`, `obsm["X_umap_tnk"]`, `ct_*_final`,
`is_ct_final_resolved`); later steps read them from there.

1. 465 features: the 461 panel genes measured in every section, three HLA allele-specific probes and an ORAI2 exon-junction probe; mutation-allele and wild-type locus probes, negative controls and protein channels are left out. Log-normalized; PCA (200 components, 82.8% of the variance), Harmony on section, cosine
   15-nearest-neighbor graph, Leiden (resolution 1.0).
2. Batch-restricted islands and small artifact clusters are removed (`02_counts_qv/artifact_cells.csv`); graph,
   UMAP (`min_dist` 0.1, `spread` 1.5) and Leiden are recomputed.
3. Cell types per cluster from Wilcoxon marker genes (`sc.tl.rank_genes_groups`). The early myeloid-progenitor
   cluster is subclustered on the QV-filtered counts (PCA 200, Harmony on section, cosine 15-NN, Leiden 0.4). Clusters
   without a coherent lineage signal over the section-level background profile are Indeterminate (635,222 cells).
4. T/NK cells are screened for other-lineage signal: Leiden (resolution 1.0) on the cohort embedding, and a cluster
   with at least 2 of its top-5 Wilcoxon markers in the erythroid, myeloid, stromal, endothelial, B-cell, macrophage
   or megakaryocyte marker sets is T/NK mixed (33,270 cells). The other T/NK cells (52,542) are subclustered into
   T and NK states (PCA 200, Harmony on section, 15-NN with Euclidean distance, Leiden resolution 1.25). Their
   display UMAP (`obsm["X_umap_tnk"]`) is computed separately (PCA 200, cosine 15-NN).
5. 5 lineages, 23 tier-2 types, 29 named tier-3 types (4,001,365 resolved cells) and the broad labels of Fig. 3a;
   see `src/mfx/taxonomy.py`.
