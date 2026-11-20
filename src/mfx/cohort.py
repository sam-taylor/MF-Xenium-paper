"""Access to the cohort object and its section table."""
import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem

from . import paths

GROUPS = ["MFUT", "Pre-Tx", "Remission", "Relapse"]  # CTMS (JAK-inhibitor treated) sections are not compared
PROBE_COL = "jak2_alt_pos_xy"                                         # >= 1 JAK2 V617F mutant-probe transcript
CLONE_COL = "is_mut_imp_local_spatial_knn_cttier2final_selected"      # probe+ or imputed clone-zone cell
NICHE_COL = "niche_clr50_km_k20_final"
NICHE_GROUP_COL = "niche_group8_clr50_km_k20_final_label"
NICHE_UMAP_KEY = "X_umap_niche_clr50_nn150_md02_final"


def open_group(mode: str = "r"):
    return zarr.open_group(str(paths.cohort_object()), mode=mode, use_consolidated=False)


def read_obs(cols) -> pd.DataFrame:
    """Selected obs columns (sample_id stripped of whitespace)."""
    g = open_group()
    df = pd.DataFrame({c: read_elem(g["obs"][c]) for c in cols})
    if "sample_id" in df:
        df["sample_id"] = df["sample_id"].astype(str).str.strip()
    return df


def read_str(col: str) -> np.ndarray:
    return np.asarray(read_elem(open_group()["obs"][col])).astype(str)


def read_bool(col: str) -> np.ndarray:
    return pd.Series(np.asarray(read_elem(open_group()["obs"][col]))).fillna(False).astype(bool).to_numpy()


def read_obsm(key: str):
    return read_elem(open_group()["obsm"][key])


def read_obsp(key: str):
    return read_elem(open_group()["obsp"][key])


def var_names() -> pd.Index:
    return pd.Index([str(v) for v in open_group()["var"]["_index"][:]])


def sections() -> pd.DataFrame:
    """One row per section: sample_id, patient_id, group (fig2_group), days from transplant,
    disease_cat_trans and is_jak2_patient (False for patients without a JAK2 mutation, the background
    controls of the clone-zone imputation)."""
    o = read_obs(["sample_id", "patient_id", "fig2_group", "calculated_days", "disease_cat_trans", "is_non_jak2_control"])
    t = o.drop_duplicates("sample_id").rename(columns={"fig2_group": "group", "calculated_days": "days"})
    t["patient_id"] = t["patient_id"].astype(str)
    t["group"] = t["group"].astype(str)
    t["is_jak2_patient"] = ~t.pop("is_non_jak2_control").astype(str).eq("True")
    return t.reset_index(drop=True)


def disease_active_sections() -> pd.DataFrame:
    """Sections of JAK2-mutant patients with active disease: untreated MF and the pre-transplant and
    relapse biopsies (disease_cat_trans 'Disease_*'); remission and JAK-inhibitor sections are left out."""
    s = sections()
    active = s.disease_cat_trans.astype(str).str.startswith("Disease") | s.sample_id.str.startswith("MFUT")
    return s[s.is_jak2_patient & active & (s.group != "CTMS")].reset_index(drop=True)


# Writing results into the cohort object (pipeline steps)
def write_obs_columns(columns: dict) -> None:
    """Add or replace obs columns of the cohort object (dict name -> array/Series/Categorical)."""
    import anndata as ad
    from anndata.io import write_elem

    ad.settings.zarr_write_format = 3
    g = open_group("r+")
    obs = read_elem(g["obs"])
    for name, values in columns.items():
        obs[name] = pd.Series(values.array if isinstance(values, pd.Series) else values, index=obs.index)
    del g["obs"]
    write_elem(g, "obs", obs)
    zarr.consolidate_metadata(str(paths.cohort_object()))


def write_obsm(key: str, value) -> None:
    import anndata as ad
    from anndata.io import write_elem

    ad.settings.zarr_write_format = 3
    g = open_group("r+")
    if key in g["obsm"]:
        del g["obsm"][key]
    write_elem(g["obsm"], key, value)
    zarr.consolidate_metadata(str(paths.cohort_object()))


def write_obsp(key: str, value) -> None:
    import anndata as ad
    from anndata.io import write_elem

    ad.settings.zarr_write_format = 3
    g = open_group("r+")
    if key in g["obsp"]:
        del g["obsp"][key]
    write_elem(g["obsp"], key, value)
    zarr.consolidate_metadata(str(paths.cohort_object()))


def write_uns(key: str, value) -> None:
    import anndata as ad
    from anndata.io import write_elem

    ad.settings.zarr_write_format = 3
    g = open_group("r+")
    if key in g["uns"]:
        del g["uns"][key]
    write_elem(g["uns"], key, value)
    zarr.consolidate_metadata(str(paths.cohort_object()))


def obs_names() -> np.ndarray:
    """Cell ids (obs index) of the cohort object, in row order."""
    return np.array([str(v) for v in open_group()["obs"]["_index"][:]])
