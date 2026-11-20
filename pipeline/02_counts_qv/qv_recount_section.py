"""Counts of one section from Proseg-assigned transcripts, QV >= 15 (row of the section map as argument)."""
import os
import re
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc
import scipy.sparse as sp

from mfx import paths
from mfx.qc import outlier_flag

MIN_QV = 15
PROBE = re.compile(r"(?:_ALT_|_WT_|_WT$|_MUT$|_alt_)")
SEX_GENES = {"XIST", "DDX3Y", "EIF1AY", "KDM5D", "RPS4Y1", "USP9Y", "UTY", "ZFY"}
HLA_ALLELES = ["HLA_DPA1*02:01:01:01_b1+b2", "HLA_DPB1*02:01:02:01_b1+b2", "HLA_DPB1*13:01:01:01_b1+b2"]
OUT = paths.path("work") / "qv_recount"
OUT.mkdir(parents=True, exist_ok=True)

row = pd.read_csv(paths.path("sections_csv")).iloc[int(os.environ.get("SLURM_ARRAY_TASK_ID") or sys.argv[1])]
x = row.xenium_id

t = pq.read_table(row.proseg_metadata, columns=["transcript_id", "gene", "assignment"]).to_pandas()
t = t.merge(pq.read_table(row.transcripts, columns=["transcript_id", "qv"]).to_pandas(), on="transcript_id", how="left")
t = t[t.assignment.notna()].copy()
t["cell"] = t.assignment.astype(np.int64)
t["gene"] = t.gene.astype(str)

artifact = pd.read_csv(os.path.join(os.path.dirname(__file__), "artifact_cells.csv"))
t = t[~t.cell.isin(artifact.loc[artifact.xenium_id == x, "cell_id"])]

exempt = t.gene.str.contains(PROBE) | t.gene.isin(SEX_GENES)
tq = t[exempt | (t.qv >= MIN_QV)]
cells, feats = np.sort(tq.cell.unique()), np.array(sorted(tq.gene.unique()))
m = sp.csr_matrix((np.ones(len(tq), dtype=np.float32), (np.searchsorted(cells, tq.cell), np.searchsorted(feats, tq.gene))),
                  shape=(len(cells), len(feats)))
m.sum_duplicates()

feature_type = np.where([f.startswith("NegControl") for f in feats], "Negative Control",
                        np.where([bool(PROBE.search(f)) for f in feats], "Mutation Probe", "Gene Expression"))
a = ad.AnnData(m, obs=pd.DataFrame({"cell_id": cells.astype(str), "xenium_id": x}, index=[f"{x}_{c}" for c in cells]),
               var=pd.DataFrame({"feature_type": feature_type}, index=feats))

gex = a[:, a.var.feature_type == "Gene Expression"].copy()
sc.pp.filter_cells(gex, min_counts=10)
a = a[gex.obs_names].copy()
sc.pp.filter_genes(a, min_cells=5)
a.obs["total_counts_all_features"] = np.asarray(a.X.sum(1)).ravel()

gex = a[:, a.var.feature_type == "Gene Expression"].copy()
sc.pp.calculate_qc_metrics(gex, percent_top=(10, 20, 50, 150), inplace=True)
gex.obs["outlier"] = outlier_flag(gex.obs)
a.obs = gex.obs
a.layers["counts"] = a.X.copy()

kept = a.obs.cell_id.astype(np.int64).to_numpy()
hla = t[t.gene.isin(HLA_ALLELES) & t.cell.isin(kept)]
hla_m = np.zeros((a.n_obs, len(HLA_ALLELES)), dtype=np.int32)
pos = pd.Series(np.arange(a.n_obs), index=kept)
np.add.at(hla_m, (pos.loc[hla.cell].to_numpy(), pd.Index(HLA_ALLELES).get_indexer(hla.gene)), 1)
a.obsm["hla_alleles_any_qv"] = hla_m
a.write_h5ad(OUT / f"{x}.h5ad")
