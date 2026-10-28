#!/usr/bin/env python3
"""
1-cellpose-nuclei.py

Cellpose nuclear segmentation of the fixed plate's DAPI channel (w1). These
nuclei are the calcium ROIs: step 3 maps them onto the Rhod-2 movie.

Adapted from clark-deep-phenotyping/iPSC/4-cellpose-staining-run.py (Alex), with
his settings unchanged:
    cellprob_threshold = -2.0, flow_threshold = 0.8, diameter = None
His script also passes model_type="nuclei" and channels=[0, 0]; in Cellpose
v4.0.1+ both are ignored (v4 always uses its single CPSAM model on the first
channels), so they are dropped here rather than passed and warned about. The
result is the same model and the same call.

Reads (read-only):  <stain_dir>/<well>_..._s1_w1.tif        2048 x 2048 DAPI
Writes:             Calcium/<plate>/cellpose/<well>_nuc_mask.tif     uint16 labels (gitignored)
                    Calcium/<plate>/cellpose/<well>_nuc_overlay.png  sanity check
                    Calcium/<plate>/cellpose/segmentation_summary.csv

Usage (GPU node, carbonite-vcms env):
    python 1-cellpose-nuclei.py plate1
    python 1-cellpose-nuclei.py plate2 --wells B02 C05
    python 1-cellpose-nuclei.py plate1 --overwrite
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage.segmentation import find_boundaries

import pipeline

CELLPROB_THRESHOLD = -2.0
FLOW_THRESHOLD = 0.8
DIAMETER = None
CROP = 512                 # full-res side of the zoomed overlay panel


def norm(img, lo=1, hi=99.5):
    p_lo, p_hi = np.percentile(img, [lo, hi])
    return np.clip((img.astype(np.float32) - p_lo) / (p_hi - p_lo + 1e-8), 0, 1)


def save_overlay(well, plate_name, dapi, masks, path):
    """Full field + a centre crop, DAPI grey with nucleus outlines in yellow."""
    g = norm(dapi)
    rgb = np.repeat(g[..., None], 3, axis=-1)
    rgb[find_boundaries(masks, mode="inner")] = [1, 1, 0]
    y0 = x0 = (dapi.shape[0] - CROP) // 2
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    ax[0].imshow(rgb)
    ax[0].add_patch(plt.Rectangle((x0, y0), CROP, CROP, fill=False, ec="cyan", lw=1))
    ax[0].set_title(f"{plate_name}  {well}: {masks.max()} nuclei")
    ax[1].imshow(rgb[y0:y0 + CROP, x0:x0 + CROP])
    ax[1].set_title(f"centre {CROP} x {CROP} px")
    for a in ax:
        a.axis("off")
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plate", help="plate1 or plate2")
    p.add_argument("--wells", nargs="*", default=None)
    p.add_argument("--overwrite", action="store_true", help="re-segment wells that already have a mask")
    p.add_argument("--cpu", action="store_true", help="force CPU (slow, testing only)")
    args = p.parse_args()

    cfg = pipeline.load_config()
    pcfg = pipeline.plate_cfg(cfg, args.plate)
    files = pipeline.stain_files(pcfg, cfg["channels"]["nuclei"])
    wells = pipeline.select_wells(files, args.wells)
    out = pipeline.out_dir(args.plate, "cellpose")

    from cellpose import models      # imported late so --help works without a GPU
    model = models.CellposeModel(gpu=not args.cpu)
    print(f"{args.plate}: segmenting {len(wells)} wells on {model.device} -> {out}")

    rows = []
    for i, well in enumerate(wells, 1):
        mask_path = out / f"{well}_nuc_mask.tif"
        if mask_path.exists() and not args.overwrite:
            masks = tifffile.imread(mask_path)
            print(f"  [{i}/{len(wells)}] {well}: mask exists ({masks.max()} nuclei), skipping")
        else:
            dapi = tifffile.imread(files[well])
            masks, _, _ = model.eval(dapi, cellprob_threshold=CELLPROB_THRESHOLD,
                                     flow_threshold=FLOW_THRESHOLD, diameter=DIAMETER)
            masks = masks.astype("uint16")
            tifffile.imwrite(mask_path, masks)
            save_overlay(well, pcfg["name"], dapi, masks, out / f"{well}_nuc_overlay.png")
            print(f"  [{i}/{len(wells)}] {files[well].name} -> {mask_path.name} ({masks.max()} nuclei)")
        areas = np.bincount(masks.ravel())[1:]
        areas = areas[areas > 0]
        rows.append(dict(well=well, condition=pipeline.condition(cfg, well),
                         n_nuclei=int(masks.max()),
                         median_area_px=float(np.median(areas)) if areas.size else np.nan,
                         stain_file=files[well].name))

    summary = out / "segmentation_summary.csv"
    df = pd.DataFrame(rows)
    if summary.exists() and args.wells:
        old = pd.read_csv(summary)
        df = pd.concat([old[~old.well.isin(df.well)], df]).sort_values("well")
    df.to_csv(summary, index=False)
    print(f"\nWrote {summary}")
    print(df.drop(columns="stain_file").to_string(index=False))


if __name__ == "__main__":
    main()
