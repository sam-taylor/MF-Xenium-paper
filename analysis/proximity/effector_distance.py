#!/usr/bin/env python3
"""Fig. 5g: composite GZMB/GNLY/PRF1 score of T and NK cells against distance to the nearest JAK2-mutant cell;
running mean and 90th percentile over distance-sorted cells, per-section distance-bin fold changes (summarized per bin)
and a 0-50 vs 200-500 um Wilcoxon signed-rank test.
Input: cohort object. Output: <work>/results/proximity/effector_distance/"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon

from mfx import cohort, counts, paths

GROUPS = {"MFUT", "Pre-Tx", "Relapse"}
MIN_ANCHORS = 20
IMMUNE = ["T_Cell_CD8_Memory", "T_Cell_CD8_Naive", "T_Cell_CD4_Memory", "T_Cell_CD4_Naive", "T_Cell_Reg",
          "T_Cell_Prolif", "T_Cell_ISG", "NK_Cell_GNLY_CD56bright", "NK_Cell_GNLY_CD56dim"]
EFFECTORS = ["GZMB", "GNLY", "PRF1"]
BINS = [0, 25, 50, 75, 100, 150, 200, 300, 500]
CENTERS = [(BINS[i] + BINS[i + 1]) / 2 for i in range(len(BINS) - 1)]
EPS = 1e-4
WINDOW = 1500  # cells per running window
OUT = paths.results("proximity/effector_distance")

obs = cohort.read_obs(["sample_id", "ct_tier3_final"])
xy = np.asarray(cohort.read_obsm("spatial"))[:, :2]
probe = cohort.read_bool(cohort.PROBE_COL)
secs = cohort.sections().set_index("sample_id")
ok = obs.sample_id.map(secs.is_jak2_patient & secs.group.isin(GROUPS)).fillna(False).astype(bool).to_numpy()
is_imm = obs.ct_tier3_final.astype(str).isin(IMMUNE).to_numpy()
sample_idx = pd.Series(np.arange(len(obs)))[ok].groupby(obs.sample_id[ok].values).apply(np.array).to_dict()
imm_rows = np.where(is_imm & ok)[0]
X, _ = counts.log_norm(imm_rows, EFFECTORS)
score = pd.Series(X.toarray().mean(1), index=imm_rows)

parts = []
for sid, idx in sample_idx.items():
    a = idx[probe[idx]]
    if len(a) < MIN_ANCHORS:
        continue
    im = idx[is_imm[idx] & ~probe[idx]]
    if len(im) < 10:
        continue
    d, _ = cKDTree(xy[a]).query(xy[im], k=1)
    keep = d <= BINS[-1]
    parts.append(pd.DataFrame({"sample_id": sid, "d": d[keep], "score": score[im].to_numpy()[keep]}))
cells = pd.concat(parts)
cells.to_csv(OUT / "immune_cells.csv.gz", index=False)
cells["bin"] = pd.cut(cells.d, BINS, labels=False, include_lowest=True)
sb = cells.groupby(["sample_id", "bin"]).agg(n=("score", "size"), score=("score", "mean")).reset_index()
sb = sb[sb.n >= 5]
base = sb[sb["bin"].isin([6, 7])].groupby("sample_id")["score"].mean()
sb["lfc"] = np.log2((sb.score + EPS) / (sb.sample_id.map(base) + EPS))
sb = sb.dropna(subset=["lfc"])
sb["center"] = sb["bin"].map(dict(enumerate(CENTERS)))
sb.to_csv(OUT / "section_bins.csv", index=False)
sb.groupby("center")["lfc"].agg(["mean", "sem"]).to_csv(OUT / "bin_summary.csv")

cells = cells.sort_values("d").reset_index(drop=True)
roll = cells["score"].rolling(WINDOW, center=True, min_periods=WINDOW // 2)
pd.DataFrame({"d": cells["d"].rolling(WINDOW, center=True, min_periods=WINDOW // 2).median(),
              "mean": roll.mean(), "q90": roll.quantile(0.90)}).dropna().to_csv(OUT / "running_mean.csv", index=False)
prox = cells[cells.d <= 50].groupby("sample_id")["score"].mean()
dist = cells[(cells.d > 200) & (cells.d <= 500)].groupby("sample_id")["score"].agg(["mean", "size"])
dist = dist[dist["size"] >= 5]["mean"]
common = prox.index.intersection(dist.index)
lfc = np.log2((prox[common] + EPS) / (dist[common] + EPS))
w = wilcoxon(lfc)
rho, rp = spearmanr(sb.center, sb.lfc)
stats = pd.DataFrame([{"n_sections": len(common), "median_lfc_0_50_vs_200_500": lfc.median(),
                       "n_sections_lower": int((lfc < 0).sum()), "wilcoxon_p": w.pvalue, "spearman_rho_section_bins": rho,
                       "spearman_p": rp, "n_cells": len(cells), "frac_cells_score0": float((cells.score == 0).mean())}])
stats.to_csv(OUT / "stats.csv", index=False)
