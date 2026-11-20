#!/usr/bin/env python3
"""50 um radius graph within each section (no self edges): obsp spatial_connectivities (1 per edge) and spatial_distances
(um), as squidpy.gr.spatial_neighbors(coord_type="generic", radius=50) builds them per section.
Input: cohort object (obsm spatial, obs sample_id). Default: build the graph and compare it with the one stored in the
cohort object; --write-object: write it."""
import argparse

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

from mfx import cohort, paths

RADIUS = 50.0

ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()

xy = np.asarray(cohort.read_obsm("spatial"))[:, :2].astype(np.float64)
sid = cohort.read_obs(["sample_id"]).sample_id.astype(str).to_numpy()
n = len(xy)
rows, cols, dist = [], [], []
for s in np.unique(sid):
    idx = np.where(sid == s)[0]
    d = cKDTree(xy[idx]).sparse_distance_matrix(cKDTree(xy[idx]), RADIUS, output_type="coo_matrix")
    keep = d.row != d.col
    rows.append(idx[d.row[keep]]); cols.append(idx[d.col[keep]]); dist.append(d.data[keep])
r, c, dd = np.concatenate(rows), np.concatenate(cols), np.concatenate(dist)
distances = sp.csr_matrix((dd.astype(np.float32), (r, c)), shape=(n, n))
connectivities = sp.csr_matrix((np.ones(len(r), dtype=np.float32), (r, c)), shape=(n, n))
print(f"{connectivities.nnz:,} edges; mean {connectivities.nnz / n:.1f} neighbours per cell")

if args.write_object:
    cohort.write_obsp("spatial_connectivities", connectivities)
    cohort.write_obsp("spatial_distances", distances)
    cohort.write_uns("spatial_graph_canonical", {"method": "radius", "radius_um": RADIUS, "per_sample": True,
                                                "include_self": False, "note": "radius graph within each section, 50 um, no self edges"})
else:
    stored = cohort.read_obsp("spatial_connectivities").tocsr()
    same = (stored != 0).astype(np.int8) != (connectivities != 0).astype(np.int8)
    deg_new, deg_old = np.diff(connectivities.indptr), np.diff(stored.indptr)
    out = paths.results("niches")
    with open(out / "radius_graph_check.txt", "w") as f:
        f.write(f"edges built {connectivities.nnz} stored {stored.nnz} differing {same.nnz}\n"
                f"cells with a different neighbour count {int((deg_new != deg_old).sum())}\n")
    print(open(out / "radius_graph_check.txt").read())
