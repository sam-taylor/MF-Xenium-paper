"""Stack the section recounts into the cohort count matrices (raw, layers['counts'], normalized X)."""
import anndata as ad
import pandas as pd
import scanpy as sc

from mfx import paths

ad.settings.zarr_write_format = 3
TARGET_SUM = 50

sections = pd.read_csv(paths.path("sections_csv"))
parts = [ad.read_h5ad(paths.path("work") / "qv_recount" / f"{x}.h5ad") for x in sections.xenium_id]
feature_type = pd.concat([p.var["feature_type"] for p in parts]).groupby(level=0).first()
new = ad.concat(parts, join="outer", fill_value=0, merge="first")
new.var["feature_type"] = feature_type.reindex(new.var_names).to_numpy()

panel = new[:, new.var.feature_type == "Gene Expression"].copy()
sc.pp.normalize_total(panel, target_sum=TARGET_SUM)
sc.pp.log1p(panel)
out = ad.AnnData(panel.X, obs=new.obs, var=pd.DataFrame(index=panel.var_names))
out.layers["counts"] = new[:, panel.var_names].layers["counts"]
out.raw = ad.AnnData(new.layers["counts"], obs=pd.DataFrame(index=new.obs_names), var=new.var[["feature_type"]])
out.write_zarr(paths.path("work") / "cohort_counts.zarr")
