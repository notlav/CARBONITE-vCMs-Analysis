#!/usr/bin/env python3
"""
2-align-mi.py

Align the fixed cTnT stain (w2) to each well's Rhod-2 calcium movie with mutual
information (MI). Step 3 uses the shift to move the DAPI nuclei masks from
step 1 into calcium coordinates.

Adapted from clark-deep-phenotyping/iPSC-CM/3-align-mi.py (Alex). The numerics
are his, unchanged: robust 1-99% normalisation, 32-bin MI, coarse-to-fine
integer search (coarse x4 over +/-100 px, then full-res +/-6 px), frame 0 as the
calcium reference, and the Theil-Sen plate-plane QC. What changed is only file
discovery (pipeline.py / experiment.yaml) and where outputs go.

Why cTnT and not DAPI: the live movie has no nuclear channel. Rhod-2 is
cytosolic and fills the same cardiomyocytes that cTnT outlines, so cTnT is the
stain that shares texture with the movie. cTnT is used here for registration
only -- the calcium ROIs are the DAPI nuclei.

Scale: stain 2048 px at 0.678 um/px, calcium 682 px at 2.034 um/px (3x3 binning,
same 10X objective). The stain is resized to 682 x 682, then one integer
translation is fitted. Sign convention, shared with step 3:
    positive X_shift -> move stain right (calcium coordinates)
    positive Y_shift -> move stain down  (calcium coordinates)

Plate QC: the shift comes from a stage offset, so across a plate it lies on a
plane (X varies with column, Y with row). A well off that plane did not find the
true optimum. The script flags (qc_pass = False); a human inspects the overlay
and records any exclusion in experiment.yaml.

Writes: Calcium/<plate>/alignment/alignment_mi.csv
        Calcium/<plate>/alignment/overlays/<well>.png   (before / after)

Usage:
    python 2-align-mi.py plate1
    python 2-align-mi.py plate2 --wells B02 C05
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import shift as ndi_shift
from scipy.stats import theilslopes
from skimage.transform import resize, rescale

import pipeline

TARGET_SIZE = 682          # calcium frame side length
CAL_FRAME = 0              # which calcium frame to register against

# Coarse-to-fine search parameters. Alex's, except COARSE_RANGE: his 25 (+/-100 px)
# is too small here -- on plate 2 B02/C05 both fits sat on the window edge at
# Y=106 while the MI landscape over +/-240 px peaks cleanly at Y~116. 50 gives
# +/-200 px; fits that still land on the edge are flagged (at_search_edge).
COARSE_FACTOR = 4          # downsample factor for the coarse pass
COARSE_RANGE = 50          # coarse search radius (== +/- 200 px at full res)
FINE_RADIUS = COARSE_FACTOR + 2   # full-res refinement radius around coarse best
N_BINS = 32                # joint-histogram bins for MI

# MetaMorph burns a timestamp into every Rhod-2 frame: rows 0-15, cols 0-327,
# pixel values 0 / 65535, identical position in all 32 movies. Replaced with the
# frame median before registration so the text cannot drive the fit.
TEXT_BOX = (slice(0, 16), slice(0, 328))


# --- Core numerics (unchanged from iPSC-CM/3-align-mi.py) --------------------
def robust_normalize(img, pmin=1, pmax=99):
    """Percentile-clip to [0, 1] so a few bright pixels don't dominate the
    histogram bins. Returns float image in [0, 1]."""
    img = img.astype(np.float64)
    lo, hi = np.percentile(img, [pmin, pmax])
    if hi <= lo:
        return np.zeros_like(img)
    return np.clip((img - lo) / (hi - lo), 0.0, 1.0)


def mutual_information(a, b, bins=N_BINS):
    """Mutual information (nats) between two equal-length 1-D samples in [0, 1]."""
    if a.size == 0:
        return -np.inf
    hist, _, _ = np.histogram2d(a, b, bins=bins, range=[[0, 1], [0, 1]])
    pxy = hist / hist.sum()
    px = pxy.sum(axis=1)
    py = pxy.sum(axis=0)
    nz = pxy > 0
    px_py = px[:, None] * py[None, :]
    return float(np.sum(pxy[nz] * np.log(pxy[nz] / px_py[nz])))


def overlap_crops(ref, mov, dx, dy):
    """Overlapping pixels of `ref` and `mov` when `mov` is translated by
    (dx right, dy down). No wrap-around."""
    h, w = ref.shape
    y0, y1 = max(0, dy), min(h, h + dy)
    x0, x1 = max(0, dx), min(w, w + dx)
    if y1 <= y0 or x1 <= x0:
        return np.empty(0), np.empty(0)
    ref_c = ref[y0:y1, x0:x1]
    mov_c = mov[y0 - dy:y1 - dy, x0 - dx:x1 - dx]
    return ref_c.ravel(), mov_c.ravel()


def mi_at(ref, mov, dx, dy):
    a, b = overlap_crops(ref, mov, dx, dy)
    # Guard against tiny overlaps giving spuriously high MI.
    if a.size < 0.25 * ref.size:
        return -np.inf
    return mutual_information(a, b)


def coarse_to_fine_align(ref, mov):
    """Integer (dx, dy) that maximizes MI of `mov` shifted onto `ref`."""
    ref_c = rescale(ref, 1.0 / COARSE_FACTOR, order=1, anti_aliasing=True)
    mov_c = rescale(mov, 1.0 / COARSE_FACTOR, order=1, anti_aliasing=True)

    best = (-np.inf, 0, 0)
    for dy in range(-COARSE_RANGE, COARSE_RANGE + 1):
        for dx in range(-COARSE_RANGE, COARSE_RANGE + 1):
            score = mi_at(ref_c, mov_c, dx, dy)
            if score > best[0]:
                best = (score, dx, dy)
    cx, cy = best[1] * COARSE_FACTOR, best[2] * COARSE_FACTOR

    best = (-np.inf, cx, cy)
    for dy in range(cy - FINE_RADIUS, cy + FINE_RADIUS + 1):
        for dx in range(cx - FINE_RADIUS, cx + FINE_RADIUS + 1):
            score = mi_at(ref, mov, dx, dy)
            if score > best[0]:
                best = (score, dx, dy)
    mi, dx, dy = best
    return dx, dy, mi


# --- Visualization -----------------------------------------------------------
def overlay_rgb(cal, stain):
    """cal -> red, stain -> green."""
    rgb = np.zeros((*cal.shape, 3), dtype=float)
    rgb[..., 0] = robust_normalize(cal)
    rgb[..., 1] = robust_normalize(stain) * 0.6
    return rgb


def save_overlay(title, cal, stain, dx, dy, mi, path):
    stain_shifted = ndi_shift(stain.astype(float), shift=(dy, dx),
                              order=0, mode="constant", cval=0)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].imshow(overlay_rgb(cal, stain))
    axes[0].set_title(f"{title}  before")
    axes[1].imshow(overlay_rgb(cal, stain_shifted))
    axes[1].set_title(f"{title}  after  (X={dx}, Y={dy}, MI={mi:.3f})")
    for a in axes:
        a.axis("off")
    fig.suptitle("Rhod-2 frame 0 = red, cTnT stain = green (yellow = overlap)")
    plt.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# --- Loading -----------------------------------------------------------------
def load_resized_stain(path):
    stain = tifffile.imread(path)
    return resize(stain, (TARGET_SIZE, TARGET_SIZE), order=1,
                  preserve_range=True).astype(np.float64)


def load_cal_frame(path):
    # key= decodes one frame instead of the whole ~935 MB stack.
    cal = tifffile.imread(path, key=CAL_FRAME).astype(np.float64)
    cal[TEXT_BOX] = np.median(cal)
    return cal


# --- Quality control (unchanged) ---------------------------------------------
# Theil-Sen (median of pairwise slopes), not least squares with outlier
# trimming: Alex found trimming fit the plane to a failed row and flagged the
# good wells instead.
QC_MIN_WELLS = 6       # below this a plane fit is meaningless
QC_FLOOR_PX = 3.0      # never flag anything closer than this to the plane


def plate_qc(df):
    """Add plane_residual_px and qc_pass columns."""
    df = df.copy()
    if len(df) < QC_MIN_WELLS:
        print(f"\nQC skipped: {len(df)} wells is too few to fit a plate plane")
        df["plane_residual_px"] = np.nan
        df["qc_pass"] = True
        return df

    row = np.array([ord(w[0]) - ord("A") for w in df["well"]], dtype=float)
    col = np.array([int(w[1:]) for w in df["well"]], dtype=float)

    sx, ix, *_ = theilslopes(df["X_shift"].to_numpy(float), col)
    sy, iy, *_ = theilslopes(df["Y_shift"].to_numpy(float), row)
    pred_x, pred_y = sx * col + ix, sy * row + iy
    dist = np.hypot(df["X_shift"] - pred_x, df["Y_shift"] - pred_y)

    mad = np.median(np.abs(dist - np.median(dist))) * 1.4826
    thresh = max(QC_FLOOR_PX, np.median(dist) + 4 * mad)

    df["plane_residual_px"] = dist.round(2)
    df["qc_pass"] = dist <= thresh

    print(f"\nPlate QC (Theil-Sen): X = {sx:+.2f}*col {ix:+.1f}   "
          f"Y = {sy:+.2f}*row {iy:+.1f}")
    print(f"  median residual {np.median(dist):.2f}px, flagging above {thresh:.1f}px")
    failed = df.loc[~df["qc_pass"]].sort_values("plane_residual_px", ascending=False)
    for _, r in failed.iterrows():
        i = df.index[df["well"] == r["well"]][0]
        print(f"  FAIL {r['well']}: got({r['X_shift']:>4.0f},{r['Y_shift']:>4.0f}) "
              f"expected({pred_x[i]:>6.1f},{pred_y[i]:>6.1f})  "
              f"off {r['plane_residual_px']:>5.1f}px  MI={r['MI']:.3f}")
    if failed.empty:
        print(f"  all {len(df)} wells on the plane")
    else:
        print(f"  {len(failed)} of {len(df)} wells failed -- inspect their overlays; "
              f"consider excluding them in experiment.yaml")
    return df


# --- Driver ------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plate", help="plate1 or plate2")
    p.add_argument("--wells", nargs="*", default=None)
    args = p.parse_args()

    cfg = pipeline.load_config()
    pcfg = pipeline.plate_cfg(cfg, args.plate)
    stain = pipeline.stain_files(pcfg, cfg["channels"]["align"])
    videos = pipeline.video_files(cfg, pcfg)
    both = sorted(set(stain) & set(videos))
    for kind, extra in (("stain", set(stain) - set(videos)), ("video", set(videos) - set(stain))):
        if extra:
            print(f"  note: {len(extra)} well(s) have {kind} only: {sorted(extra)}")
    wells = pipeline.select_wells(both, args.wells)

    out = pipeline.out_dir(args.plate, "alignment")
    overlay_dir = pipeline.out_dir(args.plate, "alignment", "overlays")
    print(f"{args.plate}: aligning {len(wells)} wells")

    rows = []
    for well in wells:
        cal = load_cal_frame(videos[well])
        st = load_resized_stain(stain[well])
        dx, dy, mi = coarse_to_fine_align(robust_normalize(cal), robust_normalize(st))
        edge = max(abs(dx), abs(dy)) >= COARSE_RANGE * COARSE_FACTOR
        print(f"  {well}: X_shift={dx:>4d}  Y_shift={dy:>4d}  MI={mi:.4f}"
              + ("  WARNING: at search edge -- widen COARSE_RANGE" if edge else ""))
        save_overlay(f"{args.plate} {well}", cal, st, dx, dy, mi, overlay_dir / f"{well}.png")
        rows.append(dict(well=well, condition=pipeline.condition(cfg, well),
                         X_shift=dx, Y_shift=dy, MI=mi, at_search_edge=edge,
                         stain_file=stain[well].name, video_file=videos[well].name))

    df = pd.DataFrame(rows)
    csv_path = out / "alignment_mi.csv"
    if args.wells and csv_path.exists():
        old = pd.read_csv(csv_path)
        df = pd.concat([old[~old.well.isin(df.well)], df]).sort_values("well").reset_index(drop=True)
        df = df.drop(columns=[c for c in ("plane_residual_px", "qc_pass") if c in df])
    df = plate_qc(df)
    df.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path}")
    print(df.drop(columns=["stain_file", "video_file"]).to_string(index=False))


if __name__ == "__main__":
    main()
