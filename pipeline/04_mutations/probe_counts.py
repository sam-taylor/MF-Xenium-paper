#!/usr/bin/env python3
"""Per-cell mutation-probe counts from the section recoveries (probe_recovery_section.py), matched on
(xenium_id, cell_id): probe_alt_wt_xy (Proseg-assigned transcripts + recovered) and the
JAK2 V617F call jak2_alt_pos_xy (>= 1 ALT transcript in probe_alt_wt_xy). With --write-object they are added to the
cohort object (obsm probe_alt_wt_xy; obs jak2_alt_pos_xy; uns probe_xy_recovery).

Input:  <work>/probe_recovery/<xenium_id>_counts.parquet, cohort object
Output: <work>/probe_counts.npz
"""
import argparse
import glob

import numpy as np
import pandas as pd

from mfx import cohort, paths

ALT_GENE = "JAK2_ALT_1848:T"
ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()

counts = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(paths.path("work") / "probe_recovery/*_counts.parquet")))],
                   ignore_index=True)
obs = cohort.read_obs(["xenium_id", "cell_id"])
cell_index = pd.MultiIndex.from_arrays([obs["xenium_id"].astype(str), obs["cell_id"].astype(int)])
rows = cell_index.get_indexer(pd.MultiIndex.from_arrays([counts["xenium_id"], counts["cell"]]))
keep = rows >= 0
print(f"{int((~keep).sum()):,} cell x probe rows fall outside the cohort cells and are ignored")
probe = pd.Categorical(counts["gene"])
probes = list(probe.categories)
assigned = np.zeros((len(obs), len(probes)), dtype=np.int32)
recovered = np.zeros_like(assigned)
np.add.at(assigned, (rows[keep], probe.codes[keep]), counts.loc[keep, "n_assigned"].to_numpy())
np.add.at(recovered, (rows[keep], probe.codes[keep]), counts.loc[keep, "n_recovered"].to_numpy())
xy = assigned + recovered
alt = probes.index(ALT_GENE)
pos_xy = xy[:, alt] >= 1
print(f"{ALT_GENE}: positive cells {int((assigned[:, alt] >= 1).sum()):,} assigned only, {int(pos_xy.sum()):,} with recovery")
np.savez_compressed(paths.path("work") / "probe_counts.npz", probes=np.array(probes), assigned=assigned, xy=xy, pos_xy=pos_xy)

if args.write_object:
    names = cohort.obs_names()
    cohort.write_obsm("probe_alt_wt_xy", pd.DataFrame(xy, columns=probes, index=names))
    cohort.write_obs_columns({"jak2_alt_pos_xy": pos_xy})
    cohort.write_uns("probe_xy_recovery", {"polygons": "Proseg cell-polygons.geojson.gz (non-overlapping 2D)",
                                          "tie_break": "nearest cell centroid", "threshold": "alt >= 1"})
