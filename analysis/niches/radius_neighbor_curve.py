"""ED5a: neighbors per cell at a range of radii (choice of the 50 um niche radius).

Input:  cohort object (obsm spatial, obs sample_id); radius graph is built per section
Output: <work>/results/niches/radius_neighbor_curve.csv
"""
import numpy as np
import pandas as pd
import zarr
from scipy.spatial import cKDTree

from mfx import paths

RADII = np.array([10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 175, 200, 250])
g = zarr.open_group(str(paths.cohort_object()), mode="r", use_consolidated=False)
codes = g["obs"]["sample_id"]["codes"][:]
xy = g["obsm"]["spatial"][:]
n = len(codes)
counts = np.empty((len(RADII), n), np.int32)

for sid in np.unique(codes):
    m = np.flatnonzero(codes == sid)
    tree = cKDTree(xy[m])
    for i, r in enumerate(RADII):
        counts[i, m] = np.asarray(tree.query_ball_point(xy[m], float(r), return_length=True)) - 1
    print(f"section {sid} {len(m):,}", flush=True)

ref = g["obs"]["nhood_n_neighbors_radius50"][:]
assert np.array_equal(counts[RADII == 50][0], ref)

pd.DataFrame({
    "radius_um": RADII,
    "mean": counts.mean(1),
    "sd": counts.std(1),
    "median": np.median(counts, 1),
    "q25": np.percentile(counts, 25, axis=1),
    "q75": np.percentile(counts, 75, axis=1),
    "n_cells": n,
}).to_csv(paths.results("niches") / "radius_neighbor_curve.csv", index=False)
