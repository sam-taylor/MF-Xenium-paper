#!/usr/bin/env python3
"""Chimerism tables: Y-chromosome transcript bias per section (Fig. 1f), residual host cells by type in remission sections
(ED2d), host and donor share of mesenchymal cells per patient (ED2e).
Input: cohort object (raw counts, chimerism_status)
Output: <work>/results/chimerism/{y_bias_per_section,residual_host_by_type,mesenchymal_host_donor_by_patient}.csv"""
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

from mfx import cohort, paths

OUT = paths.results("chimerism")
Y_GENES, X_GENE = ["DDX3Y", "RPS4Y1", "ZFY"], "XIST"
MESENCHYMAL = ("Stromal_", "Endothelial_", "Pericyte_", "Adipocyte_")

g = cohort.open_group()
sid = cohort.read_obs(["sample_id"])["sample_id"].to_numpy()
raw_var = read_elem(g["raw"]["var"]).index
cols = [raw_var.get_loc(x) for x in Y_GENES + [X_GENE]]
C = sparse_dataset(g["raw"]["X"])
sums = [np.asarray(C[s:s + 500_000][:, cols].todense()) for s in range(0, len(sid), 500_000)]
t = pd.DataFrame(np.vstack(sums), columns=Y_GENES + [X_GENE]).groupby(sid).sum()
t["y_bias"] = t[Y_GENES].sum(1) / t.sum(1)
t.to_csv(OUT / "y_bias_per_section.csv")

obs = cohort.read_obs(["sample_id", "patient_id", "ct_tier2_final", "is_ct_final_resolved", "chimerism_status",
                       "calculated_days"])
post = (pd.to_numeric(obs["calculated_days"], errors="coerce") >= 0) & obs["is_ct_final_resolved"].astype(bool)
o = obs[post].copy()
o["call"] = o["chimerism_status"].astype(str)
o = o[o["sample_id"].isin(o.loc[o["call"].isin(["host", "donor"]), "sample_id"].unique())]
o["ct"] = o["ct_tier2_final"].astype(str)
o["mes"] = o["ct"].str.startswith(MESENCHYMAL)
group = o["sample_id"].map(cohort.sections().set_index("sample_id")["group"])
hd = o[o["call"].isin(["host", "donor"]) & (group == "Remission")]
e = hd.groupby("ct").agg(host=("call", lambda c: (c == "host").sum()), n=("call", "size"), mes=("mes", "first"))
e["pct_host"] = 100 * e["host"] / e["host"].sum()
e = e[e["pct_host"] >= 1].sort_values(["mes", "pct_host"], ascending=[False, False])
e.to_csv(OUT / "residual_host_by_type.csv")

m = o[o["mes"] & o["call"].isin(["host", "donor"])]
f = m.groupby(m["patient_id"].astype(str)).agg(host=("call", lambda c: (c == "host").sum()),
                                                 donor=("call", lambda c: (c == "donor").sum()), n=("call", "size"))
f["pct_host"], f["pct_donor"] = 100 * f["host"] / f["n"], 100 * f["donor"] / f["n"]
f.loc[sorted(f.index, key=lambda p: int(p[2:]))].to_csv(OUT / "mesenchymal_host_donor_by_patient.csv")
