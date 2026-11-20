#!/usr/bin/env python3
"""Recover mutation-probe transcripts that Proseg left unassigned or called background, by containment in its
non-overlapping 2D cell polygons (ties: nearest cell centroid). One section per call (argument or $SLURM_ARRAY_TASK_ID).

Input:  section map (config `sections_csv`); Proseg transcript-metadata.parquet and cell-polygons.geojson.gz
Output: <work>/probe_recovery/<xenium_id>_counts.parquet   per cell x probe: n_assigned, n_recovered
        <work>/probe_recovery/<xenium_id>_summary.csv      per probe: transcripts, dropped, recovered
"""
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from shapely.geometry import shape

from mfx import paths

OUT = paths.path("work") / "probe_recovery"
OUT.mkdir(parents=True, exist_ok=True)

task = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 0))
row = pd.read_csv(paths.path("sections_csv")).iloc[task]
x = row.xenium_id
proseg_dir = Path(row.proseg_metadata).parent
print(x, row.sample_id, flush=True)

tx = pd.read_parquet(proseg_dir / "transcript-metadata.parquet",
                     columns=["x", "y", "gene", "assignment", "background"])
probe = tx[tx["gene"].str.contains(r"_ALT_|_WT_", case=False, na=False)]
bg = probe["assignment"].isna() | probe["background"].astype(bool)
assigned, dropped = probe[~bg], probe[bg]

with gzip.open(proseg_dir / "cell-polygons.geojson.gz", "rt") as fh:
    features = json.load(fh)["features"]
geoms = [shape(ft["geometry"]) for ft in features]
cells = np.array([ft["properties"]["cell"] for ft in features], dtype=np.int64)

# polygons containing each dropped transcript; the nearest centroid wins where more than one does
point, poly = STRtree(geoms).query(shapely.points(dropped["x"].to_numpy(), dropped["y"].to_numpy()), predicate="within")
centroid = shapely.centroid(geoms)
dist = np.hypot(shapely.get_x(centroid[poly]) - dropped["x"].to_numpy()[point],
                shapely.get_y(centroid[poly]) - dropped["y"].to_numpy()[point])
hit = pd.DataFrame({"point": point, "cell": cells[poly], "dist": dist}).sort_values(["point", "dist"]).drop_duplicates("point")
recovered = dropped.iloc[hit["point"].to_numpy()].assign(cell=hit["cell"].to_numpy())
print(f"probe transcripts {len(probe):,}, dropped {len(dropped):,}, recovered {len(recovered):,}", flush=True)

counts = pd.concat([
    assigned.assign(cell=assigned["assignment"].astype(int)).groupby(["cell", "gene"]).size().rename("n_assigned"),
    recovered.groupby(["cell", "gene"]).size().rename("n_recovered"),
], axis=1).fillna(0).astype(np.int32).reset_index()
counts.insert(0, "xenium_id", x)
counts.to_parquet(OUT / f"{x}_counts.parquet")

summary = pd.DataFrame({"n_total": probe.groupby("gene").size(), "n_assigned": assigned.groupby("gene").size(),
                        "n_dropped": dropped.groupby("gene").size(),
                        "n_recovered": recovered.groupby("gene").size()}).fillna(0).astype(int)
summary.to_csv(OUT / f"{x}_summary.csv")
