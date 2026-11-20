#!/usr/bin/env python3
"""Niche-group enrichment of clone-zone cells (vs all cells) and of host cells (vs donor cells) (ED6d-g).

Input:  cohort object (obs sample_id, niche group, chimerism_status, disease_cat_trans, clone-zone column)
Output: <work>/results/niche_host_clone/{clone_vs_all,host_vs_donor,host_vs_clone_pooled,scatter_stats}.csv; the pooled table uses
        post-transplant JAK2 sections with >= 100 host and >= 100 donor calls; scatter_stats holds the Pearson and
        Spearman correlation across the niche groups of the pooled host and clone-zone enrichments
"""
import numpy as np
import pandas as pd
from scipy import stats

from mfx import cohort, paths
from mfx import niches as nc

OUT = paths.results("niche_host_clone")
MIN_CALLED = 100
EPS = 1e-6

obs = cohort.read_obs(["sample_id", cohort.NICHE_GROUP_COL, "chimerism_status", "disease_cat_trans"])
obs["group"] = obs.pop(cohort.NICHE_GROUP_COL).astype(str)
obs["mut"] = cohort.read_bool(cohort.CLONE_COL)
obs["host"] = obs.chimerism_status.astype(str).eq("host")
obs["donor"] = obs.chimerism_status.astype(str).eq("donor")
secs = cohort.sections()
jak2 = obs.sample_id.isin(secs[secs.is_jak2_patient].sample_id)

# clone vs all cells, JAK2 sections
j = obs[jak2]
n_all, n_mut = len(j), int(j.mut.sum())
rows = []
for g in nc.GROUP_ORDER:
    d = j[j.group == g]
    fa, fm = len(d) / n_all, int(d.mut.sum()) / n_mut
    rows.append({"label": g, "n_cells": len(d), "n_mut": int(d.mut.sum()), "frac_all": fa, "frac_mut": fm,
                 "log2_enrichment": np.log2((fm + EPS) / (fa + EPS))})
pd.DataFrame(rows).to_csv(OUT / "clone_vs_all.csv", index=False)

# host vs donor, all cells with a call in post-transplant sections
post = obs[~obs.disease_cat_trans.astype(str).eq("Disease_0")]
n_host, n_donor = int(post.host.sum()), int(post.donor.sum())
rows = []
for g in nc.GROUP_ORDER:
    d = post[post.group == g]
    fh, fd = d.host.sum() / n_host, d.donor.sum() / n_donor
    rows.append({"label": g, "frac_host": fh, "frac_donor": fd, "log2_host_vs_donor": np.log2((fh + EPS) / (fd + EPS)),
                 "n_host_all": n_host, "n_donor_all": n_donor})
pd.DataFrame(rows).to_csv(OUT / "host_vs_donor.csv", index=False)

# pooled over post-transplant JAK2 sections with both host and donor calls
called = obs.groupby("sample_id")[["host", "donor"]].sum()
keep = called.index[(called.host >= MIN_CALLED) & (called.donor >= MIN_CALLED)]
p = obs[jak2 & obs.sample_id.isin(keep) & ~obs.disease_cat_trans.astype(str).eq("Disease_0")]
secs_used = sorted(p.sample_id.unique())
pooled = p.groupby("group").agg(n=("mut", "size"), n_mut=("mut", "sum"), n_host=("host", "sum"),
                                n_donor=("donor", "sum")).reindex(nc.GROUP_ORDER)
frac = pooled / pooled.sum()
pooled["log2_mut"] = np.log2(frac.n_mut / frac.n)
pooled["log2_host"] = np.log2(frac.n_host / frac.n_donor)
pooled["n_sections"] = len(secs_used)
pooled["n_patients"] = len({s.split("-")[0] for s in secs_used})
pooled.to_csv(OUT / "host_vs_clone_pooled.csv")
pearson = stats.pearsonr(pooled.log2_mut, pooled.log2_host)
spearman = stats.spearmanr(pooled.log2_mut, pooled.log2_host)
pd.DataFrame([{"n_groups": len(pooled), "n_sections": len(secs_used), "n_patients": pooled.n_patients.iloc[0],
               "pearson_r": pearson.statistic, "pearson_p": pearson.pvalue,
               "spearman_rho": spearman.statistic, "spearman_p": spearman.pvalue}]).to_csv(OUT / "scatter_stats.csv", index=False)
print(f"{len(secs_used)} sections: {', '.join(secs_used)}")
print(pooled.round(3).to_string())
