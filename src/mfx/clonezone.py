"""Clone-zone rule shared by the pipeline step (pipeline/04_mutations/clone_zone_imputation.py) and the neighborhood-size
comparison (analysis/imputation/neighborhood_size.py)."""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from mfx import cohort

K = 100
MIN_SAME_FRAC = 0.15
X_BG = 3.0
MIN_PROBE = 1
MIN_N_CT = 10_000
NON_HEM_TIER1 = {"Stromal", "Endothelial", "Indeterminate"}


def load_cells() -> pd.DataFrame:
    """Per-cell table used by the rule: section, centroid, tier-1/2 type, probe call, control and eligibility flags."""
    obs = cohort.read_obs(["sample_id", "x_centroid", "y_centroid", "ct_tier1_final", "ct_tier2_final"]).reset_index(drop=True)
    obs["ct1"] = obs.pop("ct_tier1_final").astype(str)
    obs["ct2"] = obs.pop("ct_tier2_final").astype(str)
    obs["is_probe"] = cohort.read_bool(cohort.PROBE_COL)
    obs["is_non_jak2"] = cohort.read_bool("is_non_jak2_control")
    obs["resolved"] = cohort.read_bool("is_ct_final_resolved")
    obs["is_hem"] = ~obs.ct1.isin(NON_HEM_TIER1)
    obs["is_nlh"] = obs.is_hem & ~obs.ct1.eq("Lymphoid")  # non-lymphoid hematopoietic
    secs = cohort.sections().set_index("sample_id")
    obs["group"] = obs.sample_id.map(secs.group).fillna("")
    return obs


def background(obs: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Probe-positive rate per tier-2 type in patients without the mutation (pooled rate below MIN_N_CT control cells)."""
    ctrl = obs.is_non_jak2 & obs.is_hem & obs.resolved
    bg = obs[ctrl].groupby("ct2").agg(n=("is_probe", "size"), n_probe=("is_probe", "sum")).reset_index()
    global_bg = float(obs.loc[ctrl, "is_probe"].mean())
    bg["bg_rate"] = bg.n_probe / bg.n
    bg["bg_rate_used"] = np.where(bg.n >= MIN_N_CT, bg.bg_rate, global_bg)
    return bg, global_bg


def neighbor_counts(obs: pd.DataFrame, k: int = K, return_kth_distance: bool = False):
    """For each cell, the number of its k nearest same-section cells of its own tier-2 type, and how many of those are
    probe-positive; optionally the distance to its k-th neighbor (NaN in sections with <= k cells)."""
    n = len(obs)
    same_count = np.zeros(n, dtype=np.uint16)
    probe_count = np.zeros(n, dtype=np.uint16)
    kth = np.full(n, np.nan, dtype=np.float32)
    codes = pd.Categorical(obs.ct2).codes.astype(np.int16)
    probe = obs.is_probe.to_numpy()
    xy = obs[["x_centroid", "y_centroid"]].to_numpy(np.float32)
    for idx in obs.groupby("sample_id").indices.values():
        idx = np.asarray(idx, dtype=np.int64)
        if len(idx) <= k:
            continue
        dist, nbr = cKDTree(xy[idx]).query(xy[idx], k=k + 1, workers=-1)
        nbr_g = idx[nbr[:, 1:]]
        same = codes[nbr_g] == codes[idx, None]
        same_count[idx] = same.sum(1)
        probe_count[idx] = (same & probe[nbr_g]).sum(1)
        kth[idx] = dist[:, -1]
    return (same_count, probe_count, kth) if return_kth_distance else (same_count, probe_count)


def clone_zone(obs: pd.DataFrame, bg_rate: np.ndarray, same_count, probe_count, k: int = K):
    """Clone-zone label and local same-type probe-positive rate for every cell."""
    local_rate = probe_count / np.maximum(same_count, 1)
    enriched = (((same_count / k) >= MIN_SAME_FRAC) & (probe_count >= MIN_PROBE) & (local_rate >= X_BG * bg_rate))
    eligible = ~obs.is_non_jak2.to_numpy() & obs.is_nlh.to_numpy() & obs.resolved.to_numpy()
    label = eligible & (obs.is_probe.to_numpy() | enriched)
    return label, local_rate.astype(np.float32), eligible
