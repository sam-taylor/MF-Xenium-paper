# MF Xenium spatial paper: code

Code that builds the cohort AnnData object and computes the data behind the figure panels and reported statistics.

## Layout

```
config/paths.yaml   data locations (copy to config/paths.local.yaml)
src/mfx/            shared code: paths, cohort object access, cell-type taxonomy, niche groups
pipeline/           builds the cohort object (steps 01-06, see pipeline/README.md)
analysis/           one script per analysis; results are written to <work>/results/<topic>/
slurm/run.sbatch    batch wrapper
```

Setup: `conda env create -f environment.yml`, set the paths in `config/paths.local.yaml`, run `pipeline/` in step order, then
`analysis/`. Seeds are fixed. Proseg and the embedding and clustering (`pipeline/03_annotation`) are stochastic, and their results are
stored in the cohort object.

## Source data

| Panel | Table(s) in `<work>/results/` | Script |
|---|---|---|
| 1f | chimerism/y_bias_per_section.csv | analysis/chimerism/chimerism_tables.py |
| 2f, ED2h, ED2i | vaf/ (per-section Xenium and clinical VAFs are listed in the script) | analysis/mutants/xenium_vs_clinical_vaf_correlation.py |
| ED2d, ED2e | chimerism/residual_host_by_type.csv, mesenchymal_host_donor_by_patient.csv | analysis/chimerism/chimerism_tables.py |
| 3c, 3d, 3e | mutants/mutant_density.csv | analysis/mutants/mutant_density.py |
| 3g | deg/cell_level_wilcoxon.csv | analysis/deg/cell_level_wilcoxon.py |
| 4a, 4e | niches/umap_subsample.npy, umap_projected_jak2.npy; niche_clone_enrichment/pooled_by_group.csv | analysis/niches/niche_umap.py, niche_umap_project.py, clone_enrichment_pooled.py |
| 4b | niche_composition/niche_*_celltype_log2_enrichment.csv | analysis/niches/niche_composition_table.py |
| 4d | niche_composition/dm_contrasts.csv, dm_stage_variation.csv | analysis/niches/niche_group_counts.py, niche_composition_dm.R |
| 4f | niche_clone_enrichment/section_gee.csv | analysis/niches/clone_enrichment.py |
| 5b, ED6b, ED6c | imputation/per_type.csv, per_sample.csv, spearman.csv, host_fraction.csv | analysis/imputation/imputation_tables.py |
| 5d, 5e | deg/global_pooled/deg_global.csv | analysis/deg/pseudobulk_deg.py (`--design global_pooled`) |
| 5f, ED7b | deg/by_type/ | analysis/deg/pseudobulk_deg.py (`--design by_type`) |
| ED7a | deg/global_ct/ | analysis/deg/pseudobulk_deg.py (`--design global_ct`) |
| 5g | proximity/effector_distance/ | analysis/proximity/effector_distance.py |
| 5h | proximity/es/ | analysis/proximity/proximity_es.py |
| 5i | proximity/immune_near_far/effector_genes.csv | analysis/proximity/immune_near_far_deg.py |
| ED5a | niches/radius_neighbor_curve.csv | analysis/niches/radius_neighbor_curve.py |
| ED5e | niches/k20_ward_linkage.csv; niche_composition/niche_k20_celltype_log2_enrichment.csv; niches/niche_neighbor_density.csv | pipeline/06_spatial/niche_fit.py, analysis/niches/niche_composition_table.py, niche_neighbor_density.py |
| ED6a | imputation/burden_vs_clinical_vaf.csv, burden_vs_clinical_vaf_summary.csv | analysis/imputation/burden_vs_clinical_vaf.py |
| ED6d-g | niche_host_clone/ | analysis/niches/host_clone_enrichment.py |
| ED8a, ED8b | proximity/nk/ | analysis/proximity/nk_proximity.py |
| ED9 | checkpoint/ | analysis/checkpoint/checkpoint_burden.py |

Numbers quoted in the Methods:

| Statement | Table(s) in `<work>/results/` | Script |
|---|---|---|
| negative-control share of retained calls at QV ≥ 15 | qc/negative_control_fraction.csv | analysis/qc/negative_control_fraction.py |
| variance of the 200 principal components; transcripts per Indeterminate and resolved cell | celltyping/pca_variance.csv, transcripts_by_status.csv | analysis/celltyping/annotation_stats.py |
| clone-zone neighborhood size (50-200 neighbors) and distance to the 100th neighbor | imputation/neighborhood_size.csv, neighborhood_distance.csv | analysis/imputation/neighborhood_size.py |
| Fig. 2f regression and patient mixed model | vaf/fig2f_regression.csv | analysis/mutants/xenium_vs_clinical_vaf_correlation.py |

Panels that show or tabulate columns of the cohort object directly need no code: 3a, 3b (`obsm["X_umap"]`, `["X_umap_tnk"]`,
`ct_*_final`), 3f (cell-type and chimerism columns by section), 4c, 4g (`obsm["spatial"]`, niche columns), 5c (probe and clone-zone
columns), ED1a (section table), ED3a (lineage counts per section), ED3b, ED3c (`obsm["protein_mfi"]`), ED4a, ED4b (marker expression
by cell type), ED5b-d, ED5f (neighbor counts of one cell, niche groups of one section), Supplementary Table 5 (cells per section).
