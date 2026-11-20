#!/usr/bin/env python3
"""Fig. 3c, 3d, 3e: density of JAK2 V617F mutant cells (>= 1 ALT probe transcript) on the cohort UMAP, per MF4 and MF6 biopsy.
Input: cohort object (obs sample_id, probe call; obsm X_umap)
Output: <work>/results/mutants/mutant_density.csv (mutant cells with their density, scaled 0-1 within the biopsy)"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from mfx import cohort, paths

SAMPLES = ["MF4-T1", "MF4-T2", "MF4-T3", "MF4-T4", "MF6-T2", "MF6-T3", "MF6-T4", "MF6-T5", "MF6-T6"]

obs = cohort.read_obs(["sample_id", cohort.PROBE_COL])
keep = (obs["sample_id"].isin(SAMPLES) & cohort.read_bool(cohort.PROBE_COL)).to_numpy()
mut = ad.AnnData(obs=pd.DataFrame({"sample_id": pd.Categorical(obs["sample_id"][keep])}, index=cohort.obs_names()[keep]),
                 obsm={"X_umap": np.asarray(cohort.read_obsm("X_umap"))[keep]})
sc.tl.embedding_density(mut, basis="umap", groupby="sample_id")
out = mut.obs.rename(columns={"umap_density_sample_id": "mutant_density"})[["sample_id", "mutant_density"]]
out.rename_axis("cell_id").to_csv(paths.results("mutants") / "mutant_density.csv")
print(out.groupby("sample_id", observed=True).size().to_dict())
