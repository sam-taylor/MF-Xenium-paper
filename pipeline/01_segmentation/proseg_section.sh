#!/usr/bin/env bash
# Proseg 3.0.7 re-segmentation of one section, initialised from the onboard Xenium segmentation.
# Usage:  bash proseg_section.sh <onboard_output_dir> <output_dir>
# --min-qv 0 and --include-neg-ctrls keep all transcripts; the QV filter is applied in pipeline/02_counts_qv.
set -euo pipefail
src=$1
out=$2
mkdir -p "$out"
proseg \
  --xenium \
  --use-cell-initialization \
  --max-transcript-nucleus-distance 100 \
  --prior-seg-reassignment-prob 0.9 \
  --nuclear-reassignment-prob 0.2 \
  --ab-nihlo-bubble-prob 0.1 \
  --min-qv 0 \
  --include-neg-ctrls \
  --output-transcript-metadata "$out/transcript-metadata.parquet" \
  --output-cell-polygons "$out/cell-polygons.geojson.gz" \
  --output-cell-metadata "$out/cell-metadata.parquet" \
  --output-expected-counts "$out/expected-counts.parquet" \
  --nthreads "${SLURM_CPUS_PER_TASK:-12}" \
  "$src/transcripts.parquet"
