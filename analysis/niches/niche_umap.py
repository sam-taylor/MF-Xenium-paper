#!/usr/bin/env python3
"""UMAP of neighborhood composition (Fig. 4a, 4e): 200,000 cells with >= 10 neighbors, 150-nearest-neighbor graph on
the CLR vectors, PAGA over the niche groups as the initial layout.

Input:  <work>/results/niches/{nhood_composition_clr,imputed}.npy, obs niche group
Output: <work>/results/niches/umap_subsample.npy (NaN outside the subsample); --write-object: obsm X_umap_niche_clr50_nn150_md02_final
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from mfx import cohort, paths
from mfx import niches as nc

N_CELLS = 200_000
N_NEIGHBORS = 150
MIN_DIST = 0.2
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()
D = paths.results("niches")

imputed = np.load(D / "imputed.npy")
clr = np.load(D / "nhood_composition_clr.npy", mmap_mode="r")
group = cohort.read_str(cohort.NICHE_GROUP_COL)
valid = np.where(~imputed)[0]
sub = np.sort(np.random.RandomState(SEED).choice(valid, size=N_CELLS, replace=False))

a = ad.AnnData(obs=pd.DataFrame({"group": pd.Categorical(group[sub], nc.GROUP_ORDER)}, index=sub.astype(str)),
               obsm={"X_clr": np.asarray(clr[sub])})
sc.pp.neighbors(a, use_rep="X_clr", n_neighbors=N_NEIGHBORS, random_state=SEED)
sc.tl.paga(a, groups="group")
sc.pl.paga(a, plot=False)  # computes the node positions used as initial layout
sc.tl.umap(a, init_pos="paga", min_dist=MIN_DIST, random_state=SEED)

emb = np.full((len(group), 2), np.nan, dtype=np.float32)
emb[sub] = a.obsm["X_umap"]
np.save(D / "umap_subsample.npy", emb)
if args.write_object:
    cohort.write_obsm(cohort.NICHE_UMAP_KEY, emb)
print(f"embedded {N_CELLS:,} cells")
