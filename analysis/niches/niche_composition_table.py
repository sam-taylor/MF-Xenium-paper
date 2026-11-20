#!/usr/bin/env python3
"""Cell-type composition of the k20 niches and the log2 enrichment of each cell type per niche and niche group (Fig. 4b).

Input:  cohort object (obs niche, ct_tier3_final)
Output: <work>/results/niche_composition/{niche_celltype_composition,niche_table,niche_k20_celltype_log2_enrichment,
        niche_group_celltype_log2_enrichment}.csv; enrichment = log2 of the cell-type fraction in the niche over the
        cohort fraction (pseudocount 1e-4), resolved cell types only
"""
import numpy as np
import pandas as pd

from mfx import cohort, paths
from mfx import niches as nc
from mfx import taxonomy as tx

EPS = 1e-4
OUT = paths.results("niche_composition")
obs = cohort.read_obs([cohort.NICHE_COL, "ct_tier3_final"]).astype(str)
comp = pd.crosstab(obs[cohort.NICHE_COL].astype(int).rename("niche"), obs.ct_tier3_final.rename("ct"))
comp.to_csv(OUT / "niche_celltype_composition.csv")

types = [t for t in tx.TYPES.tier3 if t not in tx.UNRESOLVED]
cohort_frac = comp.sum() / comp.values.sum()


def enrichment(counts):
    frac = counts.div(counts.sum(1), axis=0)
    return np.log2((frac[types] + EPS) / (cohort_frac[types] + EPS)).rename(columns=tx.LABEL3)


# short niche names: the two most enriched cell types, leaving out the rare osteoclasts and pDC
naming = enrichment(comp).drop(columns=[tx.label3("Osteoclast_MMP9_ITGB3"), tx.label3("pDC_IRF8_MPEG1")])
names = {k: "–".join(naming.loc[k].sort_values(ascending=False).index[:2]) for k in comp.index}
table = pd.DataFrame([{"niche_id": k, "niche": f"k{k} {names[k]}", "group": g}
                      for g, ks in nc.GROUPS_K20.items() for k in ks])
table.to_csv(OUT / "niche_table.csv", index=False)

by_niche = enrichment(comp.loc[nc.K20_ORDER])
by_niche.index = table.niche
by_niche.insert(0, "group", table.group.to_numpy())
by_niche.to_csv(OUT / "niche_k20_celltype_log2_enrichment.csv")
enrichment(comp.loc[nc.K20_ORDER].groupby(table.group.to_numpy()).sum()).loc[nc.GROUP_ORDER] \
    .rename_axis("niche_group").to_csv(OUT / "niche_group_celltype_log2_enrichment.csv")
print(table[["niche_id", "niche", "group"]].to_string(index=False))
