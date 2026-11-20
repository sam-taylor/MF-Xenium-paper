"""Data locations from config/paths.yaml.

Resolution order: $MFX_CONFIG, config/paths.local.yaml, config/paths.yaml (relative to the repository
root). Outputs of every step go under `work`.
"""
import os
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
_CFG = None


def config() -> dict:
    global _CFG
    if _CFG is None:
        env = os.environ.get("MFX_CONFIG")
        candidates = [Path(env)] if env else []
        candidates += [ROOT / "config/paths.local.yaml", ROOT / "config/paths.yaml"]
        for p in candidates:
            if p.is_file():
                _CFG = yaml.safe_load(p.read_text())
                break
        else:
            raise FileNotFoundError("no paths config found; copy config/paths.yaml to config/paths.local.yaml")
    return _CFG


def path(key: str) -> Path:
    return Path(config()[key])


def cohort_object() -> Path:
    return path("cohort_object")


def results(name: str) -> Path:
    """Output folder of one analysis step: <work>/results/<name>."""
    p = path("work") / "results" / name
    p.mkdir(parents=True, exist_ok=True)
    return p
