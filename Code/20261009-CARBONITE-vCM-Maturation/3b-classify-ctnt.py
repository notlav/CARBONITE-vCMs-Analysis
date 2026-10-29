#!/usr/bin/env python3
"""
3b-classify-ctnt.py

Call each step-1 DAPI nucleus cardiomyocyte (cTnT+) or not (cTnT-) from the
fixed plate's cTnT stain (w2), so steps 5/6 can drop the flat non-CM nuclei.

New for this repo (no Alex counterpart). Readout: cTnT is cytoplasmic, so the
nucleus interior is a poor readout. Each nucleus gets a thin perinuclear ring,
RING_PX full-res px wide (RING_PX * 0.678 um), grown with expand_labels so
neighbouring rings never overlap (each ring pixel goes to its nearest nucleus).
The score is

    ctnt_ratio = median cTnT in the ring / the well's background

Background = mode of the well's cTnT pixel histogram (camera offset plus
autofluorescence, ~530-590 counts in every well). The median rather than the
mean, so a ring that only touches a neighbouring CM's edge on one side is not
called positive. A nucleus boxed in by neighbours gets no pixels from
expand_labels (2-186 per well, mostly inside dense CM clusters); for those the
ring is instead its own RING_PX dilation minus all nuclei, shared with the
neighbours' rings (ring_shared = True). Only if that is empty too (a nucleus
wedged between others; 444 of 96,473, both plates) is ctnt_ratio NaN, and the
nucleus is called cTnT-.

Threshold (chosen 2026-10-29 from all 31 wells, both plates): log10(ctnt_ratio)
is cleanly bimodal, with a narrow background peak at ratio ~1.0 and a broad
cTnT+ mode at ~3-8. The valley between them is wide and flat: the pooled KDE
minimum is at 1.49 (Plate 1 alone 2.14, Plate 2 alone 1.45). Moving the cut
anywhere from 1.4 to 2.2 changes the CM fraction by 1.6 percentage points on
Plate 1 (26.5 -> 24.9%) and 4.2 on Plate 2 (68.6 -> 64.4%; its valley is
shallower because dense wells have more nuclei touching CM cytoplasm).
CTNT_RATIO = 2.0 sits inside the valley on both plates and is the
conservative end of it. One threshold for both plates: they
were stained and imaged together (2026-10-09) and their background peaks match.
The script re-derives the valleys and prints the sensitivity table every run,
so a change in the data shows up.

Known limit: at 10X, a cTnT- nucleus lying on top of or tightly against a CM
(common in the dense Plate 2 wells) can still read positive. The QC overlays
are there to judge how often.

Reads (read-only):  <stain_dir>/<well>_..._s1_w2.tif                cTnT, 2048 x 2048
                    Calcium/<plate>/cellpose/<well>_nuc_mask.tif     step-1 labels
Writes:             Calcium/<plate>/ctnt_classification.csv  one row per nucleus;
                        joins calcium_features.csv on Well + ObjectNumber
                    Calcium/<plate>/ctnt_qc/<well>.png       nuclei coloured by class on cTnT
                    Calcium/<plate>/ctnt_qc/ctnt_summary.csv one row per well
                    Calcium/qc/ctnt_threshold.png            pooled distributions + threshold

Usage (no GPU needed; runs on the Mac as well as the cluster):
    python 3b-classify-ctnt.py plate1
    python 3b-classify-ctnt.py plate2 --wells B03 C05
    python 3b-classify-ctnt.py plate2 --ratio 1.5      # try another threshold
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation, find_objects, gaussian_filter1d
from scipy.signal import argrelmin
from scipy.stats import gaussian_kde
from skimage.morphology import disk
from skimage.segmentation import expand_labels, find_boundaries

import pipeline

RING_PX = 6            # perinuclear ring width, full-res px (~4 um)
CTNT_RATIO = 2.0       # cTnT+ if ring median >= this x the well background
POS_COLOR = (0.1, 1.0, 0.1)      # cTnT+ outlines
NEG_COLOR = (1.0, 0.2, 1.0)      # cTnT- outlines
CROP = 400             # full-res side of each zoom panel
SENSITIVITY = [1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0]


def data_root(cfg):
    """experiment.yaml points at the cluster mount. On the Mac the same folder is
    the repo's parent (/Volumes/SaucermanLab/...), so fall back to that."""
    if not cfg["data_root"].is_dir() and (pipeline.REPO.parent / "20261009-Lavi-plate1").is_dir():
        cfg["data_root"] = pipeline.REPO.parent
    return cfg


# --- Measurement -------------------------------------------------------------
def background(img):
    """Mode of the pixel histogram (1-count bins, lightly smoothed)."""
    lo, hi = np.percentile(img, [0.1, 50])
    h, e = np.histogram(img, bins=np.arange(np.floor(lo), np.ceil(hi) + 2))
    return float(e[np.argmax(gaussian_filter1d(h.astype(float), 5))])


def label_median(lab, vals, n):
    """Median of vals per label 0..n (NaN where a label has no pixels), plus counts."""
    order = np.lexsort((vals, lab))
    lab, vals = lab[order], vals[order]
    idx = np.arange(n + 1)
    start = np.searchsorted(lab, idx)
    count = np.searchsorted(lab, idx, side="right") - start
    med = np.full(n + 1, np.nan)
    ok = count > 0
    med[ok] = 0.5 * (vals[start[ok] + (count[ok] - 1) // 2] + vals[start[ok] + count[ok] // 2])
    return med, count


def shared_ring_median(ctnt, nuc, labels):
    """Median cTnT in each label's own RING_PX dilation minus every nucleus (rings may overlap)."""
    out = np.full(len(labels), np.nan)
    boxes = find_objects(nuc)
    fp = disk(RING_PX)
    for i, lab in enumerate(labels):
        sl = boxes[lab - 1]
        y0, x0 = max(sl[0].start - RING_PX, 0), max(sl[1].start - RING_PX, 0)
        y1, x1 = sl[0].stop + RING_PX, sl[1].stop + RING_PX
        sub = nuc[y0:y1, x0:x1]
        ring = binary_dilation(sub == lab, structure=fp) & (sub == 0)
        if ring.any():
            out[i] = np.median(ctnt[y0:y1, x0:x1][ring])
    return out


def measure(ctnt, nuc):
    """Per-nucleus ring and nucleus cTnT; ring label image for the overlay."""
    n = int(nuc.max())
    ring = np.where(nuc > 0, 0, expand_labels(nuc, RING_PX))
    sel = ring.ravel() > 0
    ring_med, ring_px = label_median(ring.ravel()[sel], ctnt.ravel()[sel], n)
    ring_sum = np.bincount(ring.ravel(), weights=ctnt.ravel(), minlength=n + 1)
    area = np.bincount(nuc.ravel(), minlength=n + 1)
    nuc_sum = np.bincount(nuc.ravel(), weights=ctnt.ravel(), minlength=n + 1)
    labels = np.flatnonzero(area)
    labels = labels[labels > 0]
    bg = background(ctnt)
    df = pd.DataFrame({
        "ObjectNumber": labels,                                  # = cellpose label in the step-1 mask
        "area_full_px": area[labels],
        "ring_px": ring_px[labels],
        "ctnt_ring_median": ring_med[labels],
        "ctnt_ring_mean": ring_sum[labels] / np.maximum(ring_px[labels], 1),
        "ctnt_nuc_mean": nuc_sum[labels] / area[labels],
        "ctnt_bg": bg,
    })
    df["ring_shared"] = df.ring_px == 0
    boxed = df.ObjectNumber[df.ring_shared].to_numpy()
    df.loc[df.ring_shared, "ctnt_ring_median"] = shared_ring_median(ctnt, nuc, boxed)
    df.loc[df.ring_shared, "ctnt_ring_mean"] = np.nan
    df["ctnt_ratio"] = df.ctnt_ring_median / bg
    return df, ring


# --- QC ----------------------------------------------------------------------
def norm(img, lo=1, hi=99.5):
    p_lo, p_hi = np.percentile(img, [lo, hi])
    return np.clip((img - p_lo) / (p_hi - p_lo + 1e-8), 0, 1)


def save_overlay(title, ctnt, nuc, pos_labels, path):
    """cTnT grey, nucleus outlines green (cTnT+) / magenta (cTnT-): full field +
    two zooms (centre, and the densest-CM corner of a 4x4 grid)."""
    rgb = np.repeat(norm(ctnt)[..., None], 3, axis=-1)
    edge = find_boundaries(nuc, mode="inner")
    pos = np.isin(nuc, pos_labels)
    rgb[edge & pos] = POS_COLOR
    rgb[edge & ~pos & (nuc > 0)] = NEG_COLOR

    N = ctnt.shape[0]
    centre = ((N - CROP) // 2, (N - CROP) // 2)
    # second zoom: the grid cell with the most mixed classes, to show the boundary cases
    g = N // 4
    pos_blocks = pos.reshape(4, g, 4, g).mean((1, 3))
    nuc_blocks = (nuc > 0).reshape(4, g, 4, g).mean((1, 3))
    mix = np.minimum(pos_blocks, nuc_blocks - pos_blocks)
    by, bx = np.unravel_index(np.argmax(mix), mix.shape)
    mixed = (by * g + (g - CROP) // 2, bx * g + (g - CROP) // 2)

    fig, ax = plt.subplots(1, 3, figsize=(18, 6.3))
    ax[0].imshow(rgb)
    for (y0, x0), c in [(centre, "cyan"), (mixed, "yellow")]:
        ax[0].add_patch(plt.Rectangle((x0, y0), CROP, CROP, fill=False, ec=c, lw=1))
    n_pos, n_all = int(len(pos_labels)), int(len(np.unique(nuc)) - 1)
    ax[0].set_title(f"{title}: {n_pos}/{n_all} cTnT+ ({100 * n_pos / max(n_all, 1):.0f}%)")
    for a, (y0, x0), name in [(ax[1], centre, "centre"), (ax[2], mixed, "most mixed")]:
        a.imshow(rgb[y0:y0 + CROP, x0:x0 + CROP], interpolation="nearest")
        a.set_title(f"{name} {CROP} x {CROP} px  (green = cTnT+, magenta = cTnT-)")
    for a in ax:
        a.axis("off")
    plt.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def valleys(x, grid, between=(1.0, 5.0)):
    """KDE local minima of log10 ratio, kept only between the two modes (tail wiggles dropped)."""
    k = gaussian_kde(x, bw_method=0.05)(grid)
    v = [float(10 ** grid[i]) for i in argrelmin(k, order=20)[0]]
    return [x for x in v if between[0] < x < between[1]], k


def threshold_report(ratio):
    """Pool every ctnt_classification.csv that exists, print valleys + sensitivity,
    and plot the distributions with the threshold."""
    tables = [pd.read_csv(p) for p in sorted(pipeline.OUT_ROOT.glob("plate*/ctnt_classification.csv"))]
    df = pd.concat(tables, ignore_index=True).dropna(subset=["ctnt_ratio"])
    grid = np.linspace(-0.3, 1.3, 801)

    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    print("\nlog10(ring / background) valleys (KDE local minima):")
    colors = {"pooled": "k", "plate1": "C0", "plate2": "C1"}
    for name, sub in [("pooled", df)] + list(df.groupby("Plate")):
        v, k = valleys(np.log10(sub.ctnt_ratio.clip(lower=0.3)), grid)
        print(f"  {name:>7}: {[round(x, 2) for x in v]}  (n = {len(sub)})")
        ax[0].plot(10 ** grid, k, color=colors.get(name, "0.5"), label=f"{name} (n={len(sub)})",
                   lw=1.5 if name == "pooled" else 1)
    print(f"\n% cTnT+ by threshold (this run uses {ratio}):")
    print(pd.DataFrame({t: df.groupby("Plate").ctnt_ratio.apply(lambda r: 100 * (r >= t).mean())
                        for t in SENSITIVITY}).round(1).to_string())

    ax[0].set_xscale("log")
    ax[0].axvline(ratio, color="0.4", ls="--", lw=1)
    ax[0].set_xlabel("cTnT ring median / well background")
    ax[0].set_ylabel("Density (KDE of log10)")
    ax[0].set_title("Pooled and per plate")
    ax[0].legend(fontsize=8)
    bins = np.logspace(-0.3, 1.3, 60)
    for (plate, well), sub in df.groupby(["Plate", "Well"]):
        h, _ = np.histogram(sub.ctnt_ratio, bins=bins, density=True)
        ax[1].plot(bins[:-1], h, lw=0.6, color="C0" if plate == "plate1" else "C1", alpha=0.7)
    ax[1].set_xscale("log")
    ax[1].axvline(ratio, color="k", ls="--", lw=1)
    ax[1].set_xlabel("cTnT ring median / well background")
    ax[1].set_title("Per well (blue = Plate 1, orange = Plate 2)")
    fig.suptitle(f"cTnT+ threshold: ring median >= {ratio} x background (ring = {RING_PX} px)", fontsize=10)
    fig.tight_layout()
    path = pipeline.OUT_ROOT / "qc" / "ctnt_threshold.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\nWrote {path}")


# --- Driver ------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plate", help="plate1 or plate2")
    p.add_argument("--wells", nargs="*", default=None)
    p.add_argument("--ratio", type=float, default=CTNT_RATIO, help="cTnT+ threshold (ring / background)")
    p.add_argument("--no-overlays", action="store_true")
    args = p.parse_args()

    cfg = data_root(pipeline.load_config())
    pcfg = pipeline.plate_cfg(cfg, args.plate)
    files = pipeline.stain_files(pcfg, cfg["channels"]["align"])         # w2 = cTnT
    mask_dir = pipeline.OUT_ROOT / args.plate / "cellpose"
    exclude = pipeline.excluded(cfg, args.plate)
    for w, why in exclude.items():
        print(f"  excluding {w} per experiment.yaml: {why}")
    wells = pipeline.select_wells([w for w in files if w not in exclude], args.wells)
    qc_dir = pipeline.out_dir(args.plate, "ctnt_qc")
    print(f"{args.plate}: classifying {len(wells)} wells, ring {RING_PX} px, cTnT+ if ratio >= {args.ratio}")

    tables, summary = [], []
    for well in wells:
        mask_path = mask_dir / f"{well}_nuc_mask.tif"
        if not mask_path.exists():
            pipeline.die(f"{mask_path} not found -- run 1-cellpose-nuclei.py first")
        nuc = tifffile.imread(mask_path).astype(np.int64)
        ctnt = tifffile.imread(files[well]).astype(np.float32)
        df, _ = measure(ctnt, nuc)
        df["cTnT_pos"] = df.ctnt_ratio >= args.ratio
        df.insert(0, "Plate", args.plate)
        df.insert(1, "Well", well)
        df.insert(3, "Condition", pipeline.condition(cfg, well))
        tables.append(df)
        if not args.no_overlays:
            save_overlay(f"{args.plate} {well}", ctnt, nuc,
                         df.ObjectNumber[df.cTnT_pos].to_numpy(), qc_dir / f"{well}.png")
        summary.append(dict(well=well, condition=pipeline.condition(cfg, well), n_nuclei=len(df),
                            n_ctnt_pos=int(df.cTnT_pos.sum()),
                            n_ring_shared=int(df.ring_shared.sum()),
                            pct_ctnt_pos=round(100 * df.cTnT_pos.mean(), 1),
                            background=df.ctnt_bg.iloc[0],
                            median_ratio=round(float(df.ctnt_ratio.median()), 2)))
        print(f"  {well}: {len(df)} nuclei, bg {df.ctnt_bg.iloc[0]:.0f}, "
              f"cTnT+ {df.cTnT_pos.sum()} ({100 * df.cTnT_pos.mean():.1f}%)")

    out = pd.concat(tables, ignore_index=True)
    csv = pipeline.OUT_ROOT / args.plate / "ctnt_classification.csv"
    sm = qc_dir / "ctnt_summary.csv"
    summary = pd.DataFrame(summary)
    if args.wells and csv.exists():
        old = pd.read_csv(csv)
        out = pd.concat([old[~old.Well.isin(wells)], out], ignore_index=True)
        olds = pd.read_csv(sm)
        summary = pd.concat([olds[~olds.well.isin(wells)], summary])
    out.sort_values(["Well", "ObjectNumber"]).to_csv(csv, index=False)
    summary.sort_values("well").to_csv(sm, index=False)
    print(f"\nWrote {csv}  ({len(out)} nuclei)")
    print(f"Wrote {sm}")
    print(summary.sort_values("well").to_string(index=False))
    threshold_report(args.ratio)


if __name__ == "__main__":
    main()
