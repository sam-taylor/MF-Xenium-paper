#!/usr/bin/env python3
"""Probe-positive and clone-zone fractions of haematopoietic cells per section beside the clinical JAK2 V617F VAF (ED6a).
Input: cohort object; clinical VAFs and groups of Fig. 2f (analysis/xenium_vs_clinical_vaf_correlation.py), replaced by the
VAF measured on the FFPE biopsy itself where the cohort object has one (jak2_v617f_vaf_ffpe).
Output: <work>/results/imputation/burden_vs_clinical_vaf.csv, burden_vs_clinical_vaf_summary.csv"""
import ast

import pandas as pd

from mfx import cohort, paths

OUT = paths.results("imputation")
NON_HEM_T1 = {"Stromal", "Endothelial"}
ACTIVE = ["MFUT", "Pre-Tx", "Relapse"]

src = ast.parse((paths.ROOT / "analysis/xenium_vs_clinical_vaf_correlation.py").read_text())
rows = [r for node in src.body if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") in ("DATA", "DATA_MF8")
        for r in ast.literal_eval(node.value)]
clin = pd.DataFrame(rows, columns=["sample_id", "xenium_vaf", "clinical_vaf", "group"]).set_index("sample_id")

obs = cohort.read_obs(["sample_id", "ct_tier1_final", "jak2_v617f_vaf_ffpe"])
ffpe = obs.drop_duplicates("sample_id").set_index("sample_id").jak2_v617f_vaf_ffpe
clin["clinical_vaf"] = ffpe.reindex(clin.index).fillna(clin.clinical_vaf)
obs["probe"] = cohort.read_bool(cohort.PROBE_COL)
obs["clone"] = cohort.read_bool(cohort.CLONE_COL)
obs = obs[cohort.read_bool("is_ct_final_resolved") & ~obs.ct_tier1_final.astype(str).isin(NON_HEM_T1)]

t = obs.groupby("sample_id").agg(n_hem=("probe", "size"), n_probe=("probe", "sum"), n_clone_zone=("clone", "sum"))
t["probe_frac"] = t.n_probe / t.n_hem
t["clone_zone_frac"] = t.n_clone_zone / t.n_hem
jak2 = cohort.sections().query("is_jak2_patient").sample_id
t = t.loc[t.index.isin(jak2)].join(clin[["clinical_vaf", "group"]], how="inner")
t = t[t.group.isin(ACTIVE) & (t.clinical_vaf > 0)].sort_values("clinical_vaf")
t["clone_zone_over_vaf"] = t.clone_zone_frac / t.clinical_vaf
t.to_csv(OUT / "burden_vs_clinical_vaf.csv")

pd.DataFrame([{"n_sections": len(t),
               "median_clinical_vaf": t.clinical_vaf.median(),
               "median_probe_frac": t.probe_frac.median(),
               "max_probe_frac": t.probe_frac.max(),
               "median_clone_zone_frac": t.clone_zone_frac.median(),
               "median_clone_zone_over_vaf": t.clone_zone_over_vaf.median(),
               "n_clone_zone_below_vaf": int((t.clone_zone_frac < t.clinical_vaf).sum()),
               "pooled_clone_zone_over_probe": t.n_clone_zone.sum() / t.n_probe.sum()}]).to_csv(
    OUT / "burden_vs_clinical_vaf_summary.csv", index=False)
