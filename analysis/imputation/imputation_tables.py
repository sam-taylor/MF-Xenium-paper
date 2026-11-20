#!/usr/bin/env python3
"""Clone-zone imputation summaries: per_type.csv (Fig. 5b), per_sample.csv and spearman.csv (ED6b), host_fraction.csv (ED6c).
Input: cohort object. Output: <work>/results/imputation/"""
import pandas as pd
from scipy.stats import spearmanr

from mfx import cohort, paths
from mfx import taxonomy as tx

OUT = paths.results("imputation")
ACTIVE_5A = {"MFUT", "Pre-Tx", "Relapse", "CTMS"}
ACTIVE_7A = {"MFUT", "Pre-Tx", "Relapse"}
NON_HEM_T1 = {"Stromal", "Endothelial"}

obs = cohort.read_obs(["sample_id", "ct_tier1_final", "ct_tier2_final", "chimerism_status_sex_only", "calculated_days"])
obs["probe"] = cohort.read_bool(cohort.PROBE_COL)
obs["imp"] = cohort.read_bool(cohort.CLONE_COL)
obs["resolved"] = cohort.read_bool("is_ct_final_resolved")
secs = cohort.sections().set_index("sample_id")
obs["jak2"] = obs.sample_id.map(secs.is_jak2_patient).fillna(False).astype(bool)
obs["group"] = obs.sample_id.map(secs.group).fillna("")

a = obs[obs.jak2 & obs.group.isin(ACTIVE_5A) & obs.resolved]
t = a.groupby(a.ct_tier2_final.astype(str)).agg(n=("imp", "size"), imputed=("imp", "mean"), probe=("probe", "mean"))
t.index = [tx.label2(c) for c in t.index]
t.sort_values("imputed", ascending=False).to_csv(OUT / "per_type.csv")

b = obs[obs.jak2 & obs.group.isin(ACTIVE_7A) & obs.resolved]
hem = ~b.ct_tier1_final.astype(str).isin(NON_HEM_T1)
per_sample = pd.DataFrame({"n_probe_hem": (b.probe & hem).groupby(b.sample_id).sum(),
                           "n_imputed": b.imp.groupby(b.sample_id).sum(),
                           "group": b.groupby("sample_id").group.first()})
per_sample.to_csv(OUT / "per_sample.csv")
rho, p = spearmanr(per_sample.n_probe_hem, per_sample.n_imputed)
pd.DataFrame([{"n_sections": len(per_sample), "spearman_rho": rho, "spearman_p": p}]).to_csv(OUT / "spearman.csv", index=False)

days = pd.to_numeric(obs.calculated_days, errors="coerce")
h = obs[obs.jak2 & obs.chimerism_status_sex_only.astype(str).isin(["host", "donor"]) & (days >= 0) & obs.resolved]
host = h.chimerism_status_sex_only.astype(str) == "host"
sets = {"Probe−": ~h.probe, "Probe+": h.probe, "Imputed−": ~h.probe & ~h.imp, "Imputed+": ~h.probe & h.imp}
hf = pd.DataFrame({k: {"n": int(m.sum()), "host_frac": float(host[m].mean())} for k, m in sets.items()}).T
hf.to_csv(OUT / "host_fraction.csv")
