#!/usr/bin/env python3
"""Numbers quoted in the Methods (Xenium Cell Typing): variance explained by the leading principal components of the 465
log-normalized features, and gene-expression transcripts per cell of Indeterminate, T/NK mixed and resolved cells.
Input: cohort object (X, obs). Output: <work>/results/celltyping/{pca_variance,transcripts_by_status}.csv"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import sparse_dataset

from mfx import cohort, paths
from mfx import taxonomy as tx

OUT = paths.results("celltyping")

X = sp.csr_matrix(sparse_dataset(cohort.open_group()["X"])[:], dtype=np.float64)
n = X.shape[0]
mu = np.asarray(X.mean(axis=0)).ravel()
cov = (X.T @ X).toarray() / n - np.outer(mu, mu)
ev = np.sort(np.linalg.eigvalsh(cov))[::-1]
frac = np.cumsum(ev) / ev.sum()
pd.DataFrame({"n_components": [10, 20, 50, 100, 200], "variance_explained": [frac[k - 1] for k in (10, 20, 50, 100, 200)]}
             ).to_csv(OUT / "pca_variance.csv", index=False)

obs = cohort.read_obs(["ct_tier3_final", "total_counts"])
ct = obs.ct_tier3_final.astype(str)
status = np.select([ct.eq(tx.INDETERMINATE), ct.eq(tx.TNK_MIXED)], ["Indeterminate", "T/NK mixed"], "Resolved")
d = obs.total_counts.groupby(status).describe(percentiles=[0.25, 0.5, 0.75])
d["frac_below_20"] = obs.total_counts.lt(20).groupby(status).mean()
d.to_csv(OUT / "transcripts_by_status.csv")
print(open(OUT / "pca_variance.csv").read(), d.round(2).to_string())
