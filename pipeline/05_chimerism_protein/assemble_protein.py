#!/usr/bin/env python3
"""Join the per-section intensity tables (quantify_protein_section.py) onto the cohort cells, matched on
(xenium_id, cell_id): obsm protein_mfi (8 proteins; NaN outside the 12 sections with the protein panel) and obs
seg_18s_mfi, seg_alphasma_vimentin_mfi (segmentation stains). With --write-object they are added to the cohort object.

Input:  <protein_dir>/<xenium_id>_protein_mfi.parquet, cohort object
Output: <work>/protein_mfi.parquet
"""
import argparse

import numpy as np
import pandas as pd

from mfx import cohort, paths

NAMES = {"prot_cd45_mfi": "CD45", "prot_cd3e_mfi": "CD3E", "prot_cd4_mfi": "CD4", "prot_cd8a_mfi": "CD8a",
         "prot_cd20_mfi": "CD20", "prot_alphasma_mfi": "alpha-SMA", "prot_vimentin_mfi": "Vimentin",
         "prot_e_cadherin_mfi": "E-cadherin"}
SEG = ["seg_18s_mfi", "seg_alphasma_vimentin_mfi"]
ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()

obs = cohort.read_obs(["xenium_id", "cell_id"])
xid, cid = obs["xenium_id"].astype(str).to_numpy(), obs["cell_id"].astype(str).to_numpy()
out = pd.DataFrame(np.nan, index=cohort.obs_names(), columns=SEG + list(NAMES.values()), dtype="float32")
for x in pd.unique(xid):
    t = pd.read_parquet(paths.path("protein_dir") / f"{x}_protein_mfi.parquet").rename(columns=NAMES)
    m = xid == x
    cols = [c for c in out.columns if c in t.columns]
    out.loc[m, cols] = t.reindex(cid[m])[cols].to_numpy(dtype="float32")
out.to_parquet(paths.path("work") / "protein_mfi.parquet")
print(out.notna().sum())

if args.write_object:
    cohort.write_obsm("protein_mfi", out[list(NAMES.values())])
    cohort.write_obs_columns({c: out[c].to_numpy() for c in SEG})
