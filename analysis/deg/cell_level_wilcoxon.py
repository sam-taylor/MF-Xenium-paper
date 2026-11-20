#!/usr/bin/env python3
"""Fig. 3g: JAK2 mutant (probe-positive) vs wild-type cells, Wilcoxon rank-sum on the log-normalised matrix, pooled over cell
types, in patients with >= 50 mutant cells at any timepoint.
Input: cohort object. Output: <work>/results/deg/cell_level_wilcoxon.csv"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata.io import sparse_dataset

from mfx import cohort, paths

MIN_MUT_PER_TIMEPOINT = 50

obs = cohort.read_obs(["patient_id", "timepoint", cohort.PROBE_COL, "is_ct_final_resolved"])
resolved = obs["is_ct_final_resolved"].fillna(False).to_numpy(bool)
mut = obs[cohort.PROBE_COL].fillna(False).to_numpy(bool) & resolved
patient = obs["patient_id"].astype(str)
n_mut = pd.Series(mut).groupby([patient, obs["timepoint"].astype(str)]).sum()
keep_pt = n_mut[n_mut >= MIN_MUT_PER_TIMEPOINT].index.get_level_values(0).unique()
idx = np.where(resolved & patient.isin(keep_pt).to_numpy())[0]

X = sparse_dataset(cohort.open_group()["X"])
a = ad.AnnData(X=sp.csr_matrix(X[idx]), obs=pd.DataFrame(
    {"JAK2_status": pd.Categorical(np.where(mut[idx], "Mutant", "WT"))}, index=cohort.obs_names()[idx]),
    var=pd.DataFrame(index=cohort.var_names()))
sc.tl.rank_genes_groups(a, "JAK2_status", groups=["Mutant"], reference="WT", method="wilcoxon", use_raw=False)
df = sc.get.rank_genes_groups_df(a, group="Mutant")
df.to_csv(paths.results("deg") / "cell_level_wilcoxon.csv", index=False)
print(f"{len(keep_pt)} patients, {a.n_obs:,} cells, mutant {int(mut[idx].sum()):,}")
