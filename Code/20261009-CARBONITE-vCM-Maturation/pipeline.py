"""
pipeline.py

Shared config and file discovery for the 20261009 CARBONITE-vCM calcium pipeline
(steps 1-4). Everything path-related lives here so each step finds the same
files the same way.

The discovery rules follow clark-deep-phenotyping/iPSC-CM/3-align-mi.py (Alex):
files are globbed, the well is parsed out of the filename, and exactly one file
per well is required -- if there are zero or several, stop rather than guess.
"""

import re
import sys
from pathlib import Path

import yaml

CODE_DIR = Path(__file__).resolve().parent
REPO = CODE_DIR.parents[1]                      # CARBONITE-vCMs-Analysis/
CONFIG = CODE_DIR / "experiment.yaml"
OUT_ROOT = REPO / "Calcium"

WELL_RE = re.compile(r"(?:^|_)([A-P][0-9]{2})_")   # B02_... (stain) or ..._Stream_B02_... (video)


def die(msg):
    sys.exit(f"error: {msg}")


def load_config():
    cfg = yaml.safe_load(CONFIG.read_text())
    cfg["data_root"] = Path(cfg["data_root"])
    return cfg


def plate_cfg(cfg, plate):
    if plate not in cfg["plates"]:
        die(f"unknown plate {plate!r}; choose from {list(cfg['plates'])}")
    p = dict(cfg["plates"][plate])
    p["stain_dir"] = cfg["data_root"] / p["stain_dir"]
    p["video_dir"] = cfg["data_root"] / p["video_dir"]
    for k in ("stain_dir", "video_dir"):
        if not p[k].is_dir():
            die(f"{plate} {k} not found: {p[k]}")
    return p


def out_dir(plate, *parts):
    d = OUT_ROOT / plate
    for part in parts:
        d = d / part
    d.mkdir(parents=True, exist_ok=True)
    return d


def wells_from(paths, kind):
    """Map well -> path, requiring exactly one file per well."""
    by_well = {}
    for p in paths:
        m = WELL_RE.search(p.name)
        if m:
            by_well.setdefault(m.group(1), []).append(p)
    dupes = {w: [p.name for p in v] for w, v in by_well.items() if len(v) > 1}
    if dupes:
        die(f"{len(dupes)} well(s) have multiple {kind} files, e.g. {next(iter(dupes.items()))}")
    return {w: v[0] for w, v in sorted(by_well.items())}


def stain_files(pcfg, channel):
    """well -> stain image for one channel ID (e.g. 'w1')."""
    files = wells_from(sorted(pcfg["stain_dir"].glob(f"*_s1_{channel}.tif")), f"stain {channel}")
    if not files:
        die(f"no *_s1_{channel}.tif in {pcfg['stain_dir']}")
    return files


def video_files(cfg, pcfg):
    files = wells_from(sorted(pcfg["video_dir"].glob(cfg["video_glob"])), "video")
    if not files:
        die(f"no '{cfg['video_glob']}' in {pcfg['video_dir']}")
    return files


def condition(cfg, well):
    c = cfg["conditions_by_row"].get(well[0])
    if c is None:
        die(f"experiment.yaml has no condition for row {well[0]} ({well})")
    return c


def excluded(cfg, plate):
    return dict((cfg.get("exclude_wells") or {}).get(plate) or {})


def select_wells(available, requested):
    if not requested:
        return sorted(available)
    requested = [w.upper() for w in requested]
    missing = [w for w in requested if w not in available]
    if missing:
        die(f"requested well(s) not available: {missing}")
    return requested
