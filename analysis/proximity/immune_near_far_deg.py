#!/usr/bin/env python3
"""Fig. 5i: T and NK cells <= 50 um vs > 150 um from the nearest JAK2-mutant cell, pseudobulk DESeq2 per section and zone.
Input: cohort object. Output: <work>/results/proximity/immune_near_far/{deg.csv, effector_genes.csv, pseudobulk_meta.csv}"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.spatial import cKDTree

from mfx import cohort, counts, paths

GROUPS = {"MFUT", "Pre-Tx", "Relapse"}
MIN_ANCHORS = 20
MIN_PB_CELLS = 5
PROX, DIST = 50, 150
IMMUNE = ["T_Cell_CD8_Memory", "T_Cell_CD8_Naive", "T_Cell_CD4_Memory", "T_Cell_CD4_Naive", "T_Cell_Reg",
          "T_Cell_Prolif", "T_Cell_ISG", "NK_Cell_GNLY_CD56bright", "NK_Cell_GNLY_CD56dim"]
GENES_5D = ["IFNG", "TNF", "GZMB", "PRF1", "GNLY", "KLRD1", "KLRK1", "FCGR3A", "PDCD1", "HAVCR2", "TIGIT",
            "LAG3", "CTLA4", "TOX", "IL2", "CXCL10", "IL6", "CD274"]
OUT = paths.results("proximity/immune_near_far")

obs = cohort.read_obs(["sample_id", "ct_tier3_final"])
xy = np.asarray(cohort.read_obsm("spatial"))[:, :2]
probe = cohort.read_bool(cohort.PROBE_COL)
secs = cohort.sections().set_index("sample_id")
ok = obs.sample_id.map(secs.is_jak2_patient & secs.group.isin(GROUPS)).fillna(False).astype(bool).to_numpy()
is_imm = obs.ct_tier3_final.astype(str).isin(IMMUNE).to_numpy()
sample_idx = pd.Series(np.arange(len(obs)))[ok].groupby(obs.sample_id[ok].values).apply(np.array).to_dict()

imm_rows = np.where(is_imm & ok)[0]
C, genes = counts.raw_gex(imm_rows)
order = [genes.index(v) for v in cohort.var_names() if v in genes]  # panel order of the object
C, genes = C[:, order], [genes[i] for i in order]
row_of = pd.Series(np.arange(len(imm_rows)), index=imm_rows)

parts = []
for sid, idx in sample_idx.items():
    a = idx[probe[idx]]
    if len(a) < MIN_ANCHORS:
        continue
    im = idx[is_imm[idx] & ~probe[idx]]
    if len(im) < 10:
        continue
    d, _ = cKDTree(xy[a]).query(xy[im], k=1)
    for zone, m in (("proximal", d <= PROX), ("distal", d > DIST)):
        parts.append(pd.DataFrame({"sample_id": sid, "zone": zone, "row": im[m]}))
cells = pd.concat(parts, ignore_index=True)
cells["pb"] = cells.sample_id + "|" + cells.zone
n = cells.pb.value_counts()
use = cells[cells.pb.isin(n[n >= MIN_PB_CELLS].index)]
both = use.groupby("sample_id").zone.nunique()
use = use[use.sample_id.isin(both[both == 2].index)]
pbs = use.pb.unique()
M = sp.csr_matrix((np.ones(len(use)), (pd.Categorical(use.pb, pbs).codes, row_of[use.row].to_numpy())),
                  shape=(len(pbs), len(imm_rows)))
count_df = pd.DataFrame(np.asarray((M @ C).todense()).round().astype(int), index=pbs, columns=genes)
meta = use.drop_duplicates("pb").set_index("pb").loc[pbs, ["sample_id", "zone"]].copy()
meta["n_cells"] = n.loc[pbs].to_numpy()
meta["zone"] = pd.Categorical(meta.zone, ["distal", "proximal"])
meta.to_csv(OUT / "pseudobulk_meta.csv")

keep = (count_df >= 5).sum(axis=0) >= (len(count_df) * 0.3)
count_df = count_df.loc[:, keep]
dds = DeseqDataSet(counts=count_df, metadata=meta, design="~ sample_id + zone", quiet=True)
dds.deseq2()
st = DeseqStats(dds, contrast=["zone", "proximal", "distal"], quiet=True)
st.summary()
res = st.results_df.sort_values("padj")
res.to_csv(OUT / "deg.csv")
sel = res.reindex(GENES_5D)[["log2FoldChange", "lfcSE", "stat", "padj"]]
sel["passed_count_filter"] = [bool(keep.get(g, False)) for g in GENES_5D]
sel.to_csv(OUT / "effector_genes.csv")
