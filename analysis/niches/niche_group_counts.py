#!/usr/bin/env python3
"""Cells per section and niche group, the input of the Dirichlet-multinomial model (Fig. 4d).

Input:  cohort object (obs sample_id, niche group); sections of MFUT, pre-transplant, remission and relapse
Output: <work>/results/niche_composition/{niche_group_counts,niche_group_columns}.csv
"""
import pandas as pd

from mfx import cohort, paths
from mfx import niches as nc

OUT = paths.results("niche_composition")
secs = cohort.sections()
secs = secs[secs.group.isin(cohort.GROUPS)].set_index("sample_id")
obs = cohort.read_obs(["sample_id", cohort.NICHE_GROUP_COL])
ORDER = sorted(nc.GROUP_ORDER)
counts = pd.crosstab(obs.sample_id, obs[cohort.NICHE_GROUP_COL].astype(str))[ORDER]
counts = counts.loc[counts.index.isin(secs.index)]
cols = {g: f"n{i}" for i, g in enumerate(ORDER)}
df = counts.rename(columns=cols).join(secs[["patient_id", "group"]])
df["total"] = counts.sum(1)
df.index.name = "sample_id"
df.to_csv(OUT / "niche_group_counts.csv")
pd.DataFrame({"niche": list(cols), "response": list(cols.values())}).to_csv(OUT / "niche_group_columns.csv", index=False)
print(df.shape, df.group.value_counts().to_dict())
