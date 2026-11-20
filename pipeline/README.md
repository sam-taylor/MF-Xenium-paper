# Pipeline

Paths are set in `config/paths.yaml`. Per-section steps (marked "array") run once per section, as a Slurm array over the rows of `sections_csv` or over a list of
section directories (`sbatch --array=1-70 slurm/run.sbatch <script> ...`; create `logs/` first).

| Step | Script | Produces |
|---|---|---|
| 01 Segmentation | `01_segmentation/proseg_section.sh` | Proseg cells, polygons, transcript assignments |
| 02 Counts | `02_counts_qv/qv_recount_section.py` (array), `assemble_counts.py` | QV ≥ 15 counts, cell filter, outlier flag, normalized matrix |
| 03 Annotation | `03_annotation/README.md` | embedding, clusters, cell-type labels |
| 04 Mutations | `04_mutations/probe_recovery_section.py` (array), `probe_counts.py`, `clone_zone_imputation.py` | re-assigned probe counts, JAK2 call, clone zone |
| 05 Chimerism, protein | `05_chimerism_protein/chimerism_calls.py`, `quantify_protein_section.py` (array), `assemble_protein.py` | donor/host calls, protein and stain intensities |
| 06 Spatial | `06_spatial/radius_graph.py`, `niche_fit.py` | 50 µm radius graph; neighborhood composition, 20 niches, 8 niche groups |

Counts: `assemble_counts.py` writes `<work>/cohort_counts.zarr`. Step 03 adds the embedding and labels to it, and later steps
read and extend the result as `cohort_object` (`config/paths.yaml`). The per-section recount keeps features detected in at least
5 cells of the section (`scanpy.pp.filter_genes`, mutation probes and sex genes included).

`02_counts_qv/artifact_cells.csv` lists the cells of the artifact clusters of the cohort embedding; they are excluded from the counts.
Scripts that change the cohort object do so only with `--write-object`.
