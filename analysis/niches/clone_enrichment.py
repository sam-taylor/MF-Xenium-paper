#!/usr/bin/env python3
"""Clone enrichment per niche group (Fig. 4f): in each disease-active JAK2 section, observed JAK2 V617F probe-positive
cells among MK, erythroid and myeloid cells against the number expected from the section's rate. Per group, observed ~
Poisson with offset log(expected), intercept only, GEE with sections clustered by patient (exchangeable correlation);
BH across groups.

Input:  cohort object (obs sample_id, patient_id, ct_tier3_final, niche group, jak2_alt_pos_xy)
Output: <work>/results/niche_clone_enrichment/{per_section,section_gee}.csv
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from mfx import cohort, paths
from mfx import niches as nc
from mfx import taxonomy as tx

LABELS = {"JAK2 V617F probe+": cohort.PROBE_COL}
MIN_CELLS = 50
OUT = paths.results("niche_clone_enrichment")

obs = cohort.read_obs(["sample_id", "ct_tier3_final", cohort.NICHE_GROUP_COL])
for lab, col in LABELS.items():
    obs[lab] = cohort.read_bool(col)
secs = cohort.disease_active_sections()
cells = obs[obs.sample_id.isin(secs.sample_id) & obs.ct_tier3_final.astype(str).isin(tx.CLONE_ELIGIBLE_TIER3)].copy()
cells["group"] = cells[cohort.NICHE_GROUP_COL].astype(str)
print(f"{len(cells):,} eligible cells in {cells.sample_id.nunique()} sections", flush=True)

rows = []
for sid, d in cells.groupby("sample_id"):
    for lab in LABELS:
        x = d[lab].to_numpy()
        rate = x.mean()
        for g in nc.GROUP_ORDER:
            m = (d.group == g).to_numpy()
            n = int(m.sum())
            if n < MIN_CELLS:
                continue
            o = x[m].sum()
            rows.append({"sample_id": sid, "label": lab, "niche_group": g, "n_cells": n, "observed": int(o), "expected": n * rate})
per = pd.DataFrame(rows).merge(secs[["sample_id", "patient_id", "group"]], on="sample_id", how="left")
per.to_csv(OUT / "per_section.csv", index=False)

res = []
for lab in LABELS:
    block = []
    for g in nc.GROUP_ORDER:
        d = per[(per.label == lab) & (per.niche_group == g) & (per.expected > 0)].sort_values("patient_id")
        fit = sm.GEE(d.observed, np.ones((len(d), 1)), groups=d.patient_id, family=sm.families.Poisson(),
                     offset=np.log(d.expected), cov_struct=sm.cov_struct.Exchangeable()).fit()
        lo, hi = fit.conf_int()[0][0], fit.conf_int()[1][0]
        block.append({"label": lab, "niche_group": g, "n_sections": len(d), "n_patients": d.patient_id.nunique(),
                      "log2_rate_ratio": fit.params[0] / np.log(2), "ci_lo": lo / np.log(2), "ci_hi": hi / np.log(2),
                      "pvalue": fit.pvalues[0], "observed": int(d.observed.sum()), "expected": float(d.expected.sum()), "converged": bool(fit.converged)})
    block = pd.DataFrame(block)
    block["fdr"] = multipletests(block.pvalue, method="fdr_bh")[1]
    res.append(block)
res = pd.concat(res, ignore_index=True)
res.to_csv(OUT / "section_gee.csv", index=False)
print(res.round(3).to_string(index=False))
