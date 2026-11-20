#!/usr/bin/env python3
"""Project the other cells of the JAK2 sections onto the niche UMAP (Fig. 4e) with sc.tl.ingest.

Reference: the 200,000 fitted cells (analysis/niches/niche_umap.py) and the neighbor graph of their CLR vectors;
query: every other cell with >= 10 neighbors in a section of a JAK2-mutant patient.

Input:  <work>/results/niches/{nhood_composition_clr,imputed,umap_subsample}.npy
Output: <work>/results/niches/umap_projected_jak2.npy (fitted cells as fitted, query cells projected, NaN elsewhere)
"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from umap.umap_ import find_ab_params

from mfx import cohort, paths

N_NEIGHBORS = 150
MIN_DIST = 0.2
SEED = 42
N_CHUNKS = 12
D = paths.results("niches")

clr = np.load(D / "nhood_composition_clr.npy", mmap_mode="r")
umap = np.load(D / "umap_subsample.npy")
imputed = np.load(D / "imputed.npy")
in_jak2 = cohort.read_obs(["sample_id"]).sample_id.isin(cohort.sections().query("is_jak2_patient").sample_id).to_numpy()
fitted = np.isfinite(umap[:, 0])
query = np.where(in_jak2 & ~fitted & ~imputed)[0]

ref_rows = np.where(fitted)[0]
ref = ad.AnnData(obs=pd.DataFrame(index=ref_rows.astype(str)),
                 obsm={"X_clr": np.asarray(clr[ref_rows], dtype=np.float32), "X_umap": umap[ref_rows]})
sc.pp.neighbors(ref, use_rep="X_clr", n_neighbors=N_NEIGHBORS, random_state=SEED)
a, b = find_ab_params(1.0, MIN_DIST)
ref.uns["umap"] = {"params": {"a": a, "b": b, "random_state": SEED}}

out = umap.copy()
for rows in np.array_split(query, N_CHUNKS):
    q = ad.AnnData(obs=pd.DataFrame(index=rows.astype(str)), obsm={"X_clr": np.asarray(clr[rows], dtype=np.float32)})
    sc.tl.ingest(q, ref, embedding_method="umap")
    out[rows] = q.obsm["X_umap"]
np.save(D / "umap_projected_jak2.npy", out)
print(f"projected {len(query):,} cells onto {len(ref_rows):,} fitted cells")
