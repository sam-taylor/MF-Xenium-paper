#!/usr/bin/env python3
"""Donor / host calls per cell for sex-mismatched transplants.

host = host-sex transcripts (DDX3Y/RPS4Y1/ZFY or XIST, by pt_donor_sex) or a JAK2 ALT probe assigned by Proseg or, in MF4, donor-absent HLA alleles; donor = donor-sex transcripts; both/neither = unassigned.
chimerism_status_sex_only uses the sex genes alone.

Input:  cohort object, <work>/probe_counts.npz, <work>/qv_recount/<xenium_id>.h5ad
Output: <work>/chimerism_calls.csv.gz; with --write-object also obs chimerism_status, chimerism_status_sex_only,
        chimer_host, chimer_donor, chimer_unassigned
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

from mfx import cohort, paths

Y_GENES, XIST = ["DDX3Y", "RPS4Y1", "ZFY"], "XIST"
ap = argparse.ArgumentParser()
ap.add_argument("--write-object", action="store_true")
args = ap.parse_args()

pds = cohort.read_str("pt_donor_sex")
pid = cohort.read_str("patient_id")
g = cohort.open_group()
raw_var = list(read_elem(g["raw"]["var"]).index)
pos = sparse_dataset(g["raw"]["X"])[:][:, [raw_var.index(v) for v in Y_GENES + [XIST]]].toarray() > 0
has_y, has_xist = pos[:, :len(Y_GENES)].any(1), pos[:, len(Y_GENES)]

# HLA alleles absent from the MF4 donor (DPA1*02, DPB1*13), counted at any QV in the section recounts
obs = cohort.read_obs(["xenium_id", "cell_id"])
xid, cid = obs["xenium_id"].astype(str).to_numpy(), obs["cell_id"].astype(str).to_numpy()
has_hla = np.zeros(len(pid), dtype=bool)
for x in pd.unique(xid):
    sec = ad.read_h5ad(paths.path("work") / "qv_recount" / f"{x}.h5ad")
    m = np.where(xid == x)[0]
    j = pd.Series(np.arange(sec.n_obs), index=sec.obs["cell_id"].astype(str).to_numpy()).reindex(cid[m])
    ok = j.notna().to_numpy()
    has_hla[m[ok]] = (sec.obsm["hla_alleles_any_qv"][j[ok].astype(int).to_numpy()][:, [0, 2]] > 0).any(1)
has_hla &= pid == "MF4"
pc = np.load(paths.path("work") / "probe_counts.npz", allow_pickle=True)
jak2_assigned = pc["assigned"][:, list(pc["probes"]).index("JAK2_ALT_1848:T")] > 0

male_host = pds == "XY/XX"  # pt_donor_sex lists host first
mismatch = male_host | (pds == "XX/XY")
host_sex = np.where(male_host, has_y, has_xist)
donor_sex = np.where(male_host, has_xist, has_y)


def call(host, donor):
    return np.select([mismatch & host & ~donor, mismatch & donor & ~host], ["host", "donor"], "unassigned")


sex_only = call(host_sex, donor_sex)
final = call(host_sex | jak2_assigned | has_hla, donor_sex)
print("final   ", {k: int((final == k).sum()) for k in ("host", "donor", "unassigned")})
print("sex only", {k: int((sex_only == k).sum()) for k in ("host", "donor", "unassigned")})
pd.DataFrame({"chimerism_status": final, "chimerism_status_sex_only": sex_only},
             index=cohort.obs_names()).to_csv(paths.path("work") / "chimerism_calls.csv.gz")

if args.write_object:
    cohort.write_obs_columns({"chimerism_status": pd.Categorical(final, categories=["donor", "host", "unassigned"]),
                              "chimerism_status_sex_only": pd.Categorical(sex_only, categories=["donor", "host", "unassigned"]),
                              "chimer_host": final == "host", "chimer_donor": final == "donor",
                              "chimer_unassigned": final == "unassigned"})
