#!/usr/bin/env python3
"""Mutant-cell enrichment per niche group, pooled over the cells of all JAK2 sections (Fig. 4e).

Input:  cohort object (obs sample_id, niche group, jak2_alt_pos_xy, clone-zone column)
Output: <work>/results/niche_clone_enrichment/pooled_by_group.csv
"""
import numpy as np
import pandas as pd

from mfx import cohort, paths

COLS = [cohort.PROBE_COL, cohort.CLONE_COL]
OUT = paths.results("niche_clone_enrichment")

obs = cohort.read_obs(["sample_id", cohort.NICHE_GROUP_COL])
obs["group"] = obs.pop(cohort.NICHE_GROUP_COL).astype(str)
for col in COLS:
    obs[col] = cohort.read_bool(col)
sub = obs[obs.sample_id.isin(cohort.sections().query("is_jak2_patient").sample_id)]
print(f"{len(sub):,} cells in {sub.sample_id.nunique()} JAK2 sections")
rows = []
for col in COLS:
    p0 = sub[col].mean()
    for g, d in sub.groupby("group"):
        rows.append({"label": col, "group": g, "n_cells": len(d), "n_clone": int(d[col].sum()),
                     "clone_frac": d[col].mean(),
                     "log2_enrichment": np.log2((d[col].sum() + 0.5) / (len(d) * p0 + 0.5))})
res = pd.DataFrame(rows)
res.to_csv(OUT / "pooled_by_group.csv", index=False)
print(res.round(3).to_string(index=False))
