#!/usr/bin/env python3
"""Clone-zone vs other cells of the same section: pseudobulk DESeq2 (pydeseq2), design global_ct (ED7a), global_pooled
(Fig. 5d, 5e) or by_type (Fig. 5f, ED7b).
Input: cohort object. Output: <work>/results/deg/<design>/"""
import argparse

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from mfx import cohort, counts, paths

MIN_CELLS = 10
MIN_SAMPLES = 5

ap = argparse.ArgumentParser()
ap.add_argument("--design", choices=["global_ct", "global_pooled", "by_type"], required=True)
args = ap.parse_args()
OUT = paths.results(f"deg/{args.design}")

obs = cohort.read_obs(["sample_id", "ct_tier1_final", "ct_tier3_final"])
obs["mut"] = cohort.read_bool(cohort.CLONE_COL)
obs["jak2"] = obs.sample_id.isin(cohort.sections().query("is_jak2_patient").sample_id)
rows = np.where(cohort.read_bool("is_ct_final_resolved"))[0]
obs = obs.iloc[rows].reset_index(drop=True)
X, genes = counts.raw_gex(rows)

obs["ct"] = obs.ct_tier3_final.astype(str)
if args.design == "global_pooled":
    eligible = ~obs.ct_tier1_final.astype(str).isin(["Stromal", "Endothelial", "Lymphoid"])
    obs["ct"] = np.where(eligible, "eligible_lineages", "other")
sample_ids = sorted(obs.sample_id[obs.jak2].unique())
is_mut = (obs.mut & obs.jak2).to_numpy()
is_wt = (~obs.mut & obs.jak2).to_numpy()
mut_by_ct = obs.ct[is_mut].value_counts()
eligible_cts = sorted(mut_by_ct[mut_by_ct >= MIN_CELLS * MIN_SAMPLES].index)
if args.design == "global_pooled":
    eligible_cts = ["eligible_lineages"]
sample_arr, ct_arr = obs.sample_id.to_numpy(), obs.ct.to_numpy()


def pseudobulks(cts):
    """List of (label, count vector, metadata) for every section x type x condition with enough cells."""
    out = []
    for sid in sample_ids:
        smask = sample_arr == sid
        for ct in cts:
            cmask = ct_arr == ct
            for cond, m in (("mutant", smask & cmask & is_mut), ("wt", smask & cmask & is_wt)):
                n = int(m.sum())
                if n >= MIN_CELLS:
                    out.append((f"{sid}__{ct}__{cond}", np.asarray(X[m].sum(axis=0)).ravel(),
                                {"sample_id": sid, "ct": ct, "condition": cond, "n_cells": n}))
    return out


def filter_genes(count_df):
    return count_df.loc[:, (count_df >= 5).sum(axis=0) >= len(count_df) * 0.3]


def run_deseq(count_df, meta, design):
    dds = DeseqDataSet(counts=count_df, metadata=meta, design=design, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "mutant", "wt"], quiet=True)
    stat.summary()
    return stat.results_df


if args.design != "by_type":
    pb = pseudobulks(eligible_cts)
    count_df = pd.DataFrame([p[1] for p in pb], index=[p[0] for p in pb], columns=genes)
    meta = pd.DataFrame([p[2] for p in pb], index=[p[0] for p in pb])
    per = meta.groupby(["ct", "condition"]).size().unstack(fill_value=0)
    valid = per.index[(per.mutant >= MIN_SAMPLES) & (per.wt >= MIN_SAMPLES)]
    keep = meta.ct.isin(valid)
    count_df, meta = count_df[keep], meta[keep].copy()
    meta["condition"] = pd.Categorical(meta.condition, categories=["wt", "mutant"])
    meta["ct"] = pd.Categorical(meta.ct)
    meta["sample_id"] = pd.Categorical(meta.sample_id)
    meta.to_csv(OUT / "pseudobulk_metadata.csv")
    design = "~ sample_id + condition" if args.design == "global_pooled" else "~ sample_id + ct + condition"
    res = run_deseq(filter_genes(count_df), meta, design).sort_values("padj")
    res.to_csv(OUT / "deg_global.csv")
    print(f"{design}: {int((res.padj < 0.05).sum())} genes FDR < 0.05")
else:
    allres = {}
    for ct in eligible_cts:
        pb = pseudobulks([ct])
        meta = pd.DataFrame([p[2] | {"sample": p[0]} for p in pb]).set_index("sample")
        if (meta.condition == "mutant").sum() < MIN_SAMPLES or (meta.condition == "wt").sum() < MIN_SAMPLES:
            continue
        count_df = filter_genes(pd.DataFrame(np.array([p[1] for p in pb], dtype=int), index=meta.index, columns=genes))
        res = run_deseq(count_df, meta, "~ sample_id + condition")
        res.index.name = "gene"
        res.to_csv(OUT / f"deg_{ct[:50].replace('/', '_')}.csv")
        allres[ct] = res
        print(f"{ct}: {int((res.padj < 0.05).sum())} genes FDR < 0.05")
    pd.concat([r[r.padj < 0.05].assign(cell_type=ct) for ct, r in allres.items()]).to_csv(OUT / "deg_all_significant.csv")
