"""Per-cell protein and segmentation-stain intensity of one section: mean pixel value of the morphology-focus image
inside each Proseg polygon, round(10 * max(0, mean - offset)) with offset 100 for the 8 protein channels and 0 for the
two segmentation stains.

Usage:  python quantify_protein_section.py cell-polygons.geojson.gz <xenium_output_dir> <protein_mfi.parquet>
Output: parquet indexed by Proseg cell number, one column per channel.
"""
import gzip
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.draw import polygon
from tifffile import TiffFile

OME = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
SEG = {"18S": "seg_18s_mfi", "alphaSMA/Vimentin": "seg_alphasma_vimentin_mfi"}
PROT = {"CD4": "prot_cd4_mfi", "CD20": "prot_cd20_mfi", "CD8A": "prot_cd8a_mfi", "CD3E": "prot_cd3e_mfi",
        "Vimentin": "prot_vimentin_mfi", "alphaSMA": "prot_alphasma_mfi", "CD45": "prot_cd45_mfi",
        "E-Cadherin": "prot_e_cadherin_mfi"}
FILE_KEY = {"18s": "18S", "alphasma_vimentin": "alphaSMA/Vimentin", "cd4": "CD4", "cd20": "CD20", "cd8a": "CD8A",
            "cd3e": "CD3E", "vimentin": "Vimentin", "alphasma": "alphaSMA", "cd45": "CD45", "e-cadherin": "E-Cadherin"}
MIN_PIXELS = 4


def ome(path):
    with TiffFile(str(path)) as tf:
        return ET.fromstring(tf.ome_metadata)


def channel_files(morph_dir):
    """{channel: file}: one file per channel (ch*_<name>.ome.tif) or multi-channel files named by OME channel order."""
    per_channel = sorted(morph_dir.glob("ch*.ome.tif"))
    if per_channel:
        return {FILE_KEY[f.name[:-8].split("_", 1)[1]]: f for f in per_channel if f.name[:-8].split("_", 1)[1] in FILE_KEY}
    files = sorted(morph_dir.glob("morphology_focus_*.ome.tif"))
    names = [c.get("Name", "") for c in ome(files[0]).findall(".//ome:Channel", OME)]
    return {n: f for n, f in zip(names, files) if n in SEG}


def largest_ring(geom):
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    rings = [np.asarray(p[0], dtype=float) for p in polys if p]
    area = lambda r: 0.5 * abs(np.dot(r[:, 0], np.roll(r[:, 1], -1)) - np.dot(r[:, 1], np.roll(r[:, 0], -1)))
    return max(rings, key=area)


polygons, xenium_dir, output = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
files = channel_files(xenium_dir / "morphology_focus")
first = next(iter(files.values()))
px = ome(first).find(".//ome:Pixels", OME).get("PhysicalSizeX")
with TiffFile(str(first)) as tf:
    height, width = tf.pages[0].shape

with gzip.open(polygons, "rt") as fh:
    features = json.load(fh)["features"]
ids = [int(f["properties"]["cell"]) for f in features]
rings = [largest_ring(f["geometry"]) for f in features]
# polygon coordinates are in um (scaled to pixels) unless they already span the image
scale = 1 / float(px) if max(r[:, 0].max() for r in rings[:500]) / width < 0.5 else 1.0

out = {}
for channel, file in files.items():
    col = {**SEG, **PROT}[channel]
    offset = 0.0 if channel in SEG else 100.0
    with TiffFile(str(file)) as tf:
        img = tf.pages[0].asarray().astype(np.float32)
    mfi = np.full(len(rings), np.nan, dtype=np.float32)
    for i, ring in enumerate(rings):
        xs, ys = ring[:, 0] * scale, ring[:, 1] * scale
        x0, x1 = max(0, int(np.floor(xs.min()))), min(width, int(np.ceil(xs.max())) + 1)
        y0, y1 = max(0, int(np.floor(ys.min()))), min(height, int(np.ceil(ys.max())) + 1)
        if x1 <= x0 or y1 <= y0:
            continue
        rr, cc = polygon(ys - y0, xs - x0, shape=(y1 - y0, x1 - x0))
        if rr.size >= MIN_PIXELS:
            mfi[i] = round(10 * max(0.0, float(img[y0:y1, x0:x1][rr, cc].mean()) - offset))
    out[col] = mfi

output.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(out, index=pd.Index([str(i) for i in ids], name="cell_id")).to_parquet(output)
