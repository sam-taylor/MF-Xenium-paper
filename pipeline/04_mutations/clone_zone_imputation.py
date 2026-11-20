#!/usr/bin/env python3
"""Clone zone: JAK2 V617F probe-positive cells plus eligible neighbours whose 100 nearest same-section cells are enriched for
same-type probe-positive cells (>= 3x the type's background rate in non-JAK2 sections).
Input: cohort object. Output: <work>/results/clone_zone/; with --write-object the obs columns of the cohort object."""
import argparse

import numpy as np
import pandas as pd

from mfx import clonezone as cz
from mfx import cohort, paths

SCORE_COL = "score_mut_imp_local_spatial_knn_cttier2final_selected"

ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()
OUT = paths.results("clone_zone")

obs = cz.load_cells()
bg, global_bg = cz.background(obs)
bg.to_csv(OUT / "background_by_type.csv", index=False)
obs["bg_rate"] = obs.ct2.map(dict(zip(bg.ct2, bg.bg_rate_used))).fillna(global_bg)
same_count, probe_count = cz.neighbour_counts(obs)
label, score, eligible = cz.clone_zone(obs, obs.bg_rate.to_numpy(float), same_count, probe_count)
probe = obs.is_probe.to_numpy()
np.save(OUT / "label.npy", label)
np.save(OUT / "score.npy", score)

ctrl_a = obs.is_non_jak2.to_numpy()
rem = obs.group.eq("Remission").to_numpy() & ~ctrl_a
n_probe_eligible = int((eligible & probe).sum())
rule = pd.DataFrame([{
    "k": cz.K, "min_same_frac": cz.MIN_SAME_FRAC, "x_bg": cz.X_BG, "min_probe": cz.MIN_PROBE,
    "n_clone_zone": int(label.sum()), "n_probe_positive_eligible": n_probe_eligible,
    "probe_capture": float((label & probe & eligible).sum()) / max(n_probe_eligible, 1),
    "non_jak2_rate": float(label[ctrl_a].sum()) / max(int(ctrl_a.sum()), 1),
    "remission_rate": float(label[rem].sum()) / max(int(rem.sum()), 1),
}])
rule.to_csv(OUT / "rule_summary.csv", index=False)
pd.DataFrame({"sample_id": obs.sample_id, "probe": probe, "clone_zone": label}).groupby("sample_id").sum().to_csv(OUT / "counts_by_section.csv")

if args.write_object:
    cohort.write_obs_columns({cohort.CLONE_COL: label, SCORE_COL: score})
