# MF Xenium spatial paper: code

Code that builds the cohort AnnData object and computes the data behind the figure panels and reported statistics.

## Repository structure

```
config/paths.yaml   data locations (copy to config/paths.local.yaml)
src/mfx/            shared code: paths, cohort object access, cell-type taxonomy, niche groups, clone-zone rule
pipeline/           builds the cohort object (steps 01-06, see pipeline/README.md)
analysis/           one script per analysis; results are written to <work>/results/<topic>/
slurm/run.sbatch    batch wrapper
environment.yml     Python environment and software versions
```

## Contents

### Figure panels

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

### Numbers quoted in the Methods

| Statement | Table(s) in `<work>/results/` | Script |
|---|---|---|
| negative-control share of retained calls at QV ≥ 15 | qc/negative_control_fraction.csv | analysis/qc/negative_control_fraction.py |
| variance of the 200 principal components; transcripts per Indeterminate and resolved cell | celltyping/pca_variance.csv, transcripts_by_status.csv | analysis/celltyping/annotation_stats.py |
| clone-zone neighborhood size (50-200 neighbors) and distance to the 100th neighbor | imputation/neighborhood_size.csv, neighborhood_distance.csv | analysis/imputation/neighborhood_size.py |
| Fig. 2f regression and patient mixed model | vaf/fig2f_regression.csv | analysis/mutants/xenium_vs_clinical_vaf_correlation.py |

### Panels without code

Panels that show or tabulate columns of the cohort object directly need no code: 3a, 3b (`obsm["X_umap"]`, `["X_umap_tnk"]`,
`ct_*_final`), 3f (cell-type and chimerism columns by section), 4c, 4g (`obsm["spatial"]`, niche columns), 5c (probe and clone-zone
columns), ED1a (section table), ED3a (lineage counts per section), ED3b, ED3c (`obsm["protein_mfi"]`), ED4a, ED4b (marker expression
by cell type), ED5b-d, ED5f (neighbor counts of one cell, niche groups of one section), Supplementary Table 5 (cells per section).


## Installation

```bash
conda env create -f environment.yml    # Python 3.14 environment; installs src/mfx as the package `mfx`
conda activate mfx
```

`analysis/niches/niche_composition_dm.R` additionally needs R with brms 2.23, rstan 2.32, posterior 1.7 and data.table.

## Quickstart

1. Copy `config/paths.yaml` to `config/paths.local.yaml` and set the locations: `work` (intermediates and results),
   `cohort_object` (the cohort AnnData, zarr), and for the pipeline `sections_csv` and `protein_dir`. `MFX_CONFIG` can point to a
   config file elsewhere.
2. Run a script, for example the pseudobulk differential expression of Fig. 5d, 5e:

   ```bash
   python analysis/deg/pseudobulk_deg.py --design global_pooled
   ```

   Tables are written to `<work>/results/<topic>/`, the folders named in the tables above.
3. On a Slurm cluster: `sbatch slurm/run.sbatch analysis/deg/pseudobulk_deg.py --design global_pooled`. Per-section pipeline steps
   run as arrays (`pipeline/README.md`).

Scripts read the cohort object and change it only when run with `--write-object`.

## Data

The analysis scripts start from the cohort object, one AnnData (zarr) holding the cells of all sections with counts, cell-type
labels, probe calls, chimerism calls, the spatial graph and niche labels. Its availability is given in the Data availability
statement of the paper. `pipeline/` rebuilds it from the Xenium output of each section (Proseg segmentation inputs, transcript
tables and morphology images), listed per section in `sections_csv`.

## Things to know before running

- **Run order.** The pipeline steps run in numerical order. Some analyses read tables written by an earlier script:
  `pipeline/06_spatial/niche_fit.py` before `analysis/niches/niche_umap.py` and `niche_umap_project.py`;
  `niche_group_counts.py` before `niche_composition_dm.R`; the per-section recounts of `pipeline/02_counts_qv` before
  `analysis/qc/negative_control_fraction.py`.
- **Panels are drawn from the tables.** The scripts compute the values shown in each panel and write them as tables; the plotting
  code is not part of this repository, except for Fig. 2f, ED2h and ED2i, which the VAF script also draws.
- **Source data in a script.** The per-section Xenium and clinical VAFs of Fig. 2f, ED2h and ED2i are entered in
  `analysis/mutants/xenium_vs_clinical_vaf_correlation.py` (they are listed in Supplementary Table 1), and `analysis/imputation/burden_vs_clinical_vaf.py`
  reads the clinical VAFs from that script.
- **Stochastic steps.** Seeds are fixed. Proseg and the embedding and clustering (`pipeline/03_annotation`) are stochastic, and their
  results are stored in the cohort object rather than recomputed.
- **Resources.** The full-cohort scripts were run with 8 CPUs and up to 250 GB of memory (`slurm/run.sbatch`; adjust the partition
  and resources to your site).

## License

MIT, see [`LICENSE`](LICENSE).
