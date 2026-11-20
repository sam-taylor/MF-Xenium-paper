#!/usr/bin/env python3
"""Fig. 5h: proximity enrichment score of every gene along the distance to the nearest JAK2-mutant cell, per section
(signed peak of the normalised running sum), tested across sections with a Wilcoxon signed-rank test (BH).
Median and interquartile range across sections of the running sum of GZMB, GNLY and PRF1 are written as effector_curves.csv.
Input: cohort object. Output: <work>/results/proximity/es/{es_per_section.csv, es_gene_table.csv, effector_curves.csv}"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

from mfx import cohort, counts, paths

GROUPS = {"MFUT", "Pre-Tx", "Relapse"}
MIN_MUT_CELLS = 50
MIN_NONMUT_CELLS = 200
DIST_MAX = 600.0
DIST_BINS = np.linspace(0, DIST_MAX, 201)
DIST_CENTERS = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])
CURVE_GENES = ["GZMB", "GNLY", "PRF1"]
OUT = paths.results("proximity/es")


def running_sum(expr_ranked):
    """Normalised running sum of every column of a (cells x genes) matrix sorted by distance (row 0 = closest)."""
    total = expr_ranked.sum(axis=0)
    rs = np.cumsum((expr_ranked - total / expr_ranked.shape[0]).astype(np.float64), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, rs / total, np.nan).astype(np.float32)


rows = np.where(cohort.read_bool("is_ct_final_resolved"))[0]
obs = cohort.read_obs(["sample_id"]).iloc[rows].reset_index(drop=True)
xy = np.asarray(cohort.read_obsm("spatial"))[rows, :2]
probe = cohort.read_bool(cohort.PROBE_COL)[rows]
imp = cohort.read_bool(cohort.CLONE_COL)[rows]
X, genes = counts.log_norm(rows)
X = np.asarray(X.todense(), dtype=np.float32)

secs = cohort.sections().query("is_jak2_patient").set_index("sample_id")
totals = pd.Series(probe).groupby(obs.sample_id.values).sum()
valid = sorted(s for s in secs.index if secs.group[s] in GROUPS and int(totals.get(s, 0)) >= MIN_MUT_CELLS)
indices = obs.groupby("sample_id").indices
curve_idx = [genes.get_loc(g) for g in CURVE_GENES]

per, curves, density = [], {g: [] for g in CURVE_GENES}, []
for sid in valid:
    idx = indices[sid]
    a, nb = probe[idx], ~probe[idx] & ~imp[idx]
    if a.sum() < MIN_MUT_CELLS or nb.sum() < MIN_NONMUT_CELLS:
        continue
    d, _ = cKDTree(xy[idx][a]).query(xy[idx][nb], k=1)
    o = np.argsort(d)
    d_sorted, rws = d[o], idx[nb][o]
    rs = running_sum(X[rws])
    # signed peak: the extreme running-sum value of larger magnitude
    es = np.where(rs.max(axis=0) >= -rs.min(axis=0), rs.max(axis=0), rs.min(axis=0))
    ok = ~np.isnan(es)
    per.append(pd.DataFrame({"sample_id": sid, "gene": np.asarray(genes)[ok], "es": es[ok].astype(float),
                             "n_cells_used": int(nb.sum()), "n_anchors": int(a.sum())}))
    for g, k in zip(CURVE_GENES, curve_idx):
        at = np.interp(DIST_CENTERS, d_sorted, running_sum(X[rws, k][:, None])[:, 0])
        curves[g].append(np.where(DIST_CENTERS <= d_sorted.max(), at, np.nan))
    density.append(np.histogram(d_sorted[d_sorted <= DIST_MAX], bins=DIST_BINS)[0] / len(d_sorted))
per = pd.concat(per, ignore_index=True)
per.to_csv(OUT / "es_per_section.csv", index=False)

rows_out = []
for gene, d in per.groupby("gene"):
    v = d.es
    if len(v) < 5:
        continue
    med = v.median()
    rows_out.append({"gene": gene, "median_es": float(med), "mean_es": float(v.mean()), "n_samples": int(len(v)),
                     "n_concordant": int((np.sign(v) == np.sign(med)).sum()), "wilcoxon_p": float(wilcoxon(v).pvalue)})
tab = pd.DataFrame(rows_out)
tab["fdr"] = multipletests(tab.wilcoxon_p, method="fdr_bh")[1]
tab.sort_values("wilcoxon_p").to_csv(OUT / "es_gene_table.csv", index=False)
curve_tab = {"distance_um": DIST_CENTERS, "neighbour_density": np.mean(density, 0)}
for g, c in curves.items():
    c = np.array(c)
    curve_tab |= {f"{g}_median": np.nanmedian(c, 0), f"{g}_q1": np.nanpercentile(c, 25, 0), f"{g}_q3": np.nanpercentile(c, 75, 0)}
pd.DataFrame(curve_tab).to_csv(OUT / "effector_curves.csv", index=False)
