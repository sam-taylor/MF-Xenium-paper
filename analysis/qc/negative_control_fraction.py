#!/usr/bin/env python3
"""Share of negative-control codewords among the transcripts kept after the QV >= 15 filter, per section (Methods,
Xenium Pre-processing). Input: <work>/qv_recount/<xenium_id>.h5ad from pipeline/02_counts_qv/qv_recount_section.py.
Output: <work>/results/qc/negative_control_fraction.csv"""
import anndata as ad
import numpy as np
import pandas as pd

from mfx import paths

rows = []
for f in sorted((paths.path("work") / "qv_recount").glob("*.h5ad")):
    a = ad.read_h5ad(f)
    ft = a.var.feature_type.to_numpy()
    tot = np.asarray(a.layers["counts"].sum(0)).ravel()
    gex, neg = tot[ft == "Gene Expression"].sum(), tot[ft == "Negative Control"].sum()
    rows.append({"xenium_id": f.stem, "gene_expression": int(gex), "negative_control": int(neg), "neg_fraction": neg / (gex + neg)})
t = pd.DataFrame(rows)
t.to_csv(paths.results("qc") / "negative_control_fraction.csv", index=False)
print(t.neg_fraction.describe())
