#!/usr/bin/env python3
"""Neighbors within 50 µm and share of Indeterminate cells in each niche and niche group (ED5e, right).

Input:  cohort object (obs niche, niche group, nhood_n_neighbors_radius50, ct_tier3_final)
Output: <work>/results/niches/niche_neighbor_density.csv (rows: the 20 niches, the 8 groups and the cohort)
"""
import pandas as pd

from mfx import cohort, paths
from mfx import taxonomy as tx

MIN_NEIGHBORS_FIT = 10

obs = cohort.read_obs([cohort.NICHE_COL, cohort.NICHE_GROUP_COL, "nhood_n_neighbors_radius50", "ct_tier3_final"])
obs["niche"] = "K" + obs.pop(cohort.NICHE_COL).astype(str)
obs["group"] = obs.pop(cohort.NICHE_GROUP_COL).astype(str)
obs["n"] = obs.pop("nhood_n_neighbors_radius50")
obs["below_fit"] = obs.n < MIN_NEIGHBORS_FIT
obs["indeterminate"] = obs.pop("ct_tier3_final").astype(str) == tx.INDETERMINATE

tables = []
for level, keys in (("niche", ["group", "niche"]), ("group", ["group"]), ("cohort", [lambda i: "all cells"])):
    g = obs.groupby(keys)
    t = g.n.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).unstack().set_axis(["q10", "q25", "median", "q75", "q90"], axis=1)
    t = t.join(g.agg(n_cells=("n", "size"), mean=("n", "mean"), frac_below_fit=("below_fit", "mean"),
                     frac_indeterminate=("indeterminate", "mean")))
    t = t.reset_index().rename(columns={"index": "group"})
    t.insert(0, "level", level)
    tables.append(t)
pd.concat(tables)[["level", "group", "niche", "n_cells", "q10", "q25", "median", "q75", "q90", "mean", "frac_below_fit",
                   "frac_indeterminate"]].to_csv(paths.results("niches") / "niche_neighbor_density.csv", index=False)
