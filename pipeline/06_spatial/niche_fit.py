#!/usr/bin/env python3
"""Spatial niches: CLR composition of the 50 um neighborhood over the 31 tier-3 labels, k-means (k = 20) on cells
with >= 10 neighbors (the others take the majority niche of their 15 nearest fitted cells), Ward grouping into 8 groups.

Input:  cohort object (obsp spatial_connectivities, obs ct_tier3_final)
Output: <work>/results/niches/{nhood_composition_clr,nhood_neighbor_counts,labels_k20,imputed}.npy, k20_ward_linkage.csv
        (scipy linkage matrix of the 20 niche centroids), k20_ward8_map.csv;
        --write-object: the niche obs columns of the cohort object
"""
import argparse

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import OneHotEncoder

from mfx import cohort, paths
from mfx import niches as nc

PSEUDOCOUNT = 1e-5
MIN_NEIGHBORS = 10
KNN_IMPUTE = 15
SEED = 42
K = 20
N_GROUPS = 8

ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true", help="write the labels into the cohort object")
args = ap.parse_args()
OUT = paths.results("niches")

conn = cohort.read_obsp("spatial_connectivities").tocsr()
conn.data[:] = 1.0
conn.setdiag(0)
conn.eliminate_zeros()
onehot = OneHotEncoder(dtype=np.float32).fit_transform(cohort.read_str("ct_tier3_final")[:, None])
counts = (conn @ onehot).toarray().astype(np.float32)
n_nb = counts.sum(1)
logc = np.log(counts / np.maximum(n_nb, 1)[:, None] + PSEUDOCOUNT)
clr = (logc - logc.mean(1, keepdims=True)).astype(np.float32)
np.save(OUT / "nhood_composition_clr.npy", clr)
np.save(OUT / "nhood_neighbor_counts.npy", n_nb.astype(np.int32))
valid = n_nb >= MIN_NEIGHBORS
print(f"{valid.sum():,} fitted cells, {(~valid).sum():,} edge cells; median neighbors {np.median(n_nb):.0f}", flush=True)

fitted = clr[valid]
lab_v = MiniBatchKMeans(n_clusters=K, random_state=SEED, batch_size=10_000, n_init=10, max_iter=300).fit_predict(fitted)
labels = np.empty(len(clr), np.int32)
labels[valid] = lab_v
labels[~valid] = KNeighborsClassifier(n_neighbors=KNN_IMPUTE, n_jobs=-1).fit(fitted, lab_v).predict(clr[~valid])
np.save(OUT / "labels_k20.npy", labels)
np.save(OUT / "imputed.npy", ~valid)

centroids = np.vstack([fitted[lab_v == n].mean(0) for n in range(K)])
tree = linkage(centroids, method="ward")
pd.DataFrame(tree, columns=["left", "right", "distance", "n_niches"]).to_csv(OUT / "k20_ward_linkage.csv", index=False)
ward = fcluster(tree, t=N_GROUPS, criterion="maxclust") - 1
pd.Series(ward, name="ward_cut").rename_axis("k20").to_csv(OUT / "k20_ward8_map.csv")

if args.write_object:
    # group names in mfx.niches refer to this Ward partition of the niche ids
    assert {frozenset(ks) for ks in nc.GROUPS_K20.values()} == {frozenset(np.where(ward == w)[0]) for w in np.unique(ward)}
    cohort.write_obs_columns({
        "niche_clr50_km_k20_final": pd.Categorical(labels.astype(str), categories=[str(i) for i in range(K)]),
        "is_niche_clr50_km_k20_final_imputed": ~valid,
        "niche_group8_clr50_km_k20_final_label": pd.Categorical([nc.GROUP_OF_K20[k] for k in labels], categories=nc.GROUP_ORDER),
        "nhood_n_neighbors_radius50": n_nb.astype(np.int32),
    })
