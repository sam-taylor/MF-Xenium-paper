#!/usr/bin/env python3
"""ED9: checkpoint-gene positivity (>= 1 count) and per-cell burden in clone-zone myeloid cells (MFUT, pre-transplant, relapse)
and non-clone-zone myeloid cells (remission); section-level burden tests; UpSet combination counts (ED9a-d, f) and burden
fractions (ED9e, h).
Input: cohort object. Output: <work>/results/checkpoint/"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import kruskal

from mfx import cohort, counts, paths

GENES = ["CD47", "VSIR", "CD274", "CD80", "TIGIT", "CD86", "IDO1", "HAVCR2", "PDCD1LG2", "TNFRSF9",
         "ICOS", "LAG3", "CTLA4", "PDCD1"]
MUT_GROUPS = ["MFUT", "Pre-Tx", "Relapse"]
TOP = 15  # UpSet combinations shown
OUT = paths.results("checkpoint")

obs = cohort.read_obs(["sample_id", "patient_id", "ct_tier1_final", "ct_tier2_final"])
secs = cohort.sections().set_index("sample_id")
obs["group"] = obs.sample_id.map(secs.group)
clone = cohort.read_bool(cohort.PROBE_COL) | cohort.read_bool(cohort.CLONE_COL)
myeloid = (obs.ct_tier1_final.astype(str) == "Myeloid").to_numpy() & cohort.read_bool("is_ct_final_resolved")
jak2 = ~cohort.read_bool("is_non_jak2_control")
obs["set"] = np.where(myeloid & jak2 & clone & obs.group.isin(MUT_GROUPS).to_numpy(), obs.group,
                      np.where(myeloid & ~clone & (obs.group == "Remission").to_numpy(), "WT-Remission", ""))
keep = np.where(obs.set != "")[0]

X, _ = counts.layer_counts(keep, GENES)
cells = obs.iloc[keep][["sample_id", "patient_id", "group", "set", "ct_tier2_final"]].copy()
cells[GENES] = (X.toarray() > 0)
cells["burden"] = cells[GENES].sum(1).to_numpy()
cells["is_mk"] = cells.ct_tier2_final.astype(str).str.startswith("Megakaryocyte")
cells["patient_id"] = cells.patient_id.astype(str)
for c in ("group", "set", "ct_tier2_final"):
    cells[c] = cells[c].astype(str)

populations = {k: cells[cells.set == k] for k in (*MUT_GROUPS, "WT-Remission")}
populations |= {f"{k}_MK": populations[k][populations[k].is_mk] for k in MUT_GROUPS}
combos, totals, burden = [], [], {}
for name, c in populations.items():
    pos = c[GENES]
    n = pos.apply(lambda r: "+".join(np.array(GENES)[r.to_numpy()]), axis=1).value_counts().drop("", errors="ignore")
    combos.append(n.head(TOP).rename_axis("combination").reset_index(name="n_cells").assign(panel=name))
    totals.append(pos.sum().rename_axis("gene").reset_index(name="n_cells").assign(panel=name))
    burden[name] = c.burden.clip(upper=4).value_counts(normalize=True).reindex(range(5), fill_value=0)
pd.concat(combos)[["panel", "combination", "n_cells"]].to_csv(OUT / "upset_combinations.csv", index=False)
pd.concat(totals)[["panel", "gene", "n_cells"]].to_csv(OUT / "upset_gene_totals.csv", index=False)
bf = pd.DataFrame(burden).T.add_prefix("burden_").rename(columns={"burden_4": "burden_4plus"}).rename_axis("population")
bf.assign(n_cells=[len(c) for c in populations.values()]).to_csv(OUT / "burden_fractions.csv")

summary = pd.DataFrame([{"set": k, "n_cells": int((cells.set == k).sum()),
                         "frac_ge1": float((cells.loc[cells.set == k, "burden"] >= 1).mean()),
                         "frac_ge3": float((cells.loc[cells.set == k, "burden"] >= 3).mean()),
                         "n_sections": cells.loc[cells.set == k, "sample_id"].nunique()}
                        for k in (*MUT_GROUPS, "WT-Remission")])
summary.to_csv(OUT / "summary.csv", index=False)

sec = cells[cells.set.isin(MUT_GROUPS)].groupby(["sample_id", "patient_id", "group"])["burden"].mean().reset_index()
sec.to_csv(OUT / "section_mean_burden.csv", index=False)
kw = kruskal(*[sec.loc[sec.group == k, "burden"] for k in MUT_GROUPS]).pvalue
fit = smf.ols("burden ~ C(group, Treatment('Pre-Tx'))", sec).fit(cov_type="cluster", cov_kwds={"groups": sec["patient_id"]})
rel = "C(group, Treatment('Pre-Tx'))[T.Relapse]"
pd.DataFrame([{"n_sections": len(sec), "kruskal_p": kw, "relapse_vs_pretx_coef": fit.params[rel],
               "relapse_vs_pretx_p": fit.pvalues[rel]}]).to_csv(OUT / "burden_test.csv", index=False)
