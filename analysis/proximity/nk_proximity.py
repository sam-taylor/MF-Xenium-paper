#!/usr/bin/env python3
"""ED8a, ED8b: fraction of probe-positive cells with a CD56bright or CD56dim NK cell within 50 um per section, and the
distances to the nearest CD56bright NK cell, with group tests.
Input: cohort object. Output: <work>/results/proximity/nk/"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree
from scipy.stats import kruskal

from mfx import cohort, paths

GROUPS = ["MFUT", "Pre-Tx", "Relapse"]
BRIGHT, DIM = "NK_Cell_GNLY_CD56bright", "NK_Cell_GNLY_CD56dim"
MIN_MUT, MIN_TYPE, MIN_TYPE_KDE = 20, 10, 3
RADIUS = 50  # um, the neighborhood radius used throughout
OUT = paths.results("proximity/nk")

obs = cohort.read_obs(["sample_id", "patient_id", "ct_tier3_final"])
obs["sample_id"] = obs.sample_id.astype(str)
xy = np.asarray(cohort.read_obsm("spatial"))[:, :2]
grp = cohort.sections().query("is_jak2_patient").set_index("sample_id").group.to_dict()
mut = cohort.read_bool(cohort.PROBE_COL)
ct = obs.ct_tier3_final.astype(str).to_numpy()

rows, kde_d = [], {k: [] for k in GROUPS}
for sid, idx in obs.groupby("sample_id").indices.items():
    gname = grp.get(sid)
    if gname not in GROUPS or mut[idx].sum() < MIN_MUT:
        continue
    m = idx[mut[idx]]
    for t in (BRIGHT, DIM):
        tt = idx[ct[idx] == t]
        if len(tt) < MIN_TYPE_KDE:
            continue
        d, _ = cKDTree(xy[tt]).query(xy[m], k=1)
        if t == BRIGHT:
            kde_d[gname].extend(d.tolist())
        if len(tt) >= MIN_TYPE:
            rows.append({"sample_id": sid, "patient_id": str(obs.patient_id.iloc[idx[0]]), "group": gname, "type": t,
                         "n_type": len(tt), "n_mut": len(m), "frac_within": float(np.mean(d <= RADIUS))})
nn = pd.DataFrame(rows)
nn.to_csv(OUT / "section_metrics.csv", index=False)
pd.DataFrame([(k, v) for k, vs in kde_d.items() for v in vs], columns=["group", "distance_um"]).to_csv(OUT / "bright_distances.csv", index=False)

b, d = nn[nn.type == BRIGHT].set_index("sample_id"), nn[nn.type == DIM].set_index("sample_id")
df = b[["patient_id", "group", "frac_within", "n_type"]].join(d[["frac_within", "n_type"]], rsuffix="_dim", how="inner")
df["delta"] = df.frac_within - df.frac_within_dim
df.to_csv(OUT / "bright_minus_dim.csv")
kw_a = kruskal(*[df.loc[df.group == k, "delta"] for k in GROUPS]).pvalue
fit = smf.ols("delta ~ C(group, Treatment('Pre-Tx'))", df.reset_index()).fit(
    cov_type="cluster", cov_kwds={"groups": df["patient_id"].to_numpy()})
rel = "C(group, Treatment('Pre-Tx'))[T.Relapse]"
kw_b = kruskal(*[b.loc[b.group == k, "frac_within"] for k in GROUPS]).pvalue
pd.DataFrame([{"ed8a_kw_p": kw_a, "ed8a_ols_relapse_pp": 100 * fit.params[rel], "ed8a_ols_p": fit.pvalues[rel],
               "ed8b_kw_p": kw_b}]).to_csv(OUT / "summary.csv", index=False)
