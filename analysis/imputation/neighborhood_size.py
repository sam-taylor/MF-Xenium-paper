#!/usr/bin/env python3
"""Choice of the clone-zone neighborhood size (Methods, Clone-Zone Imputation).

For K = 50, 100, 150 and 200 nearest neighbors, the clone-zone rule is applied unchanged and scored against the sex-based
host and donor calls of resolved cells in post-transplant sections of JAK2 V617F patients (host-derived non-lymphoid
hematopoietic cells belong to the clone, donor-derived cells do not). Also reports the distance to the 100th neighbor and the
local probe-positive rate that a single same-type probe-positive neighbor gives at K = 100.
Input: cohort object. Output: <work>/results/imputation/neighborhood_size.csv, neighborhood_distance.csv"""
import numpy as np
import pandas as pd

from mfx import clonezone as cz
from mfx import cohort, paths

KS = (50, 100, 150, 200)
OUT = paths.results("imputation")

obs = cz.load_cells()
bg, global_bg = cz.background(obs)
bg_rate = obs.ct2.map(dict(zip(bg.ct2, bg.bg_rate_used))).fillna(global_bg).to_numpy(float)
extra = cohort.read_obs(["chimerism_status_sex_only", "calculated_days"])
status = extra.chimerism_status_sex_only.astype(str).to_numpy()
days = pd.to_numeric(extra.calculated_days, errors="coerce").to_numpy()
jak2 = ~obs.is_non_jak2.to_numpy()
called = jak2 & (days >= 0) & obs.resolved.to_numpy() & obs.is_nlh.to_numpy() & np.isin(status, ["host", "donor"])
host, donor = called & (status == "host"), called & (status == "donor")
probe = obs.is_probe.to_numpy()

rows = []
for k in KS:
    same, prb, kth = cz.neighbor_counts(obs, k, return_kth_distance=True)
    label, _, _ = cz.clone_zone(obs, bg_rate, same, prb, k)
    imp = label & ~probe
    rows.append({"k": k, "n_clone_zone": int(label.sum()),
                 "host_cells": int(host.sum()), "donor_cells": int(donor.sum()),
                 "host_labeled_frac": float(label[host].mean()), "donor_labeled_frac": float(label[donor].mean()),
                 "imputed_host_share": float(host[imp & called].mean())})
    if k == cz.K:
        elig = jak2 & obs.is_nlh.to_numpy() & obs.resolved.to_numpy()
        d = pd.DataFrame([{"cells": name, "n": int(np.isfinite(kth[m]).sum()),
                           "median_um": float(np.nanmedian(kth[m])), "q25_um": float(np.nanpercentile(kth[m], 25)),
                           "q75_um": float(np.nanpercentile(kth[m], 75)), "frac_within_120um": float(np.nanmean(kth[m] <= 120))}
                          for name, m in (("all cells", np.ones(len(obs), bool)), ("eligible cells, JAK2 patients", elig))])
        d["one_probe_neighbor_rate"] = 1 / cz.K
        d["max_background_x3"] = cz.X_BG * float(bg.bg_rate_used.max())
        d.to_csv(OUT / "neighborhood_distance.csv", index=False)
pd.DataFrame(rows).to_csv(OUT / "neighborhood_size.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
