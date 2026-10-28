#!/usr/bin/env python3
"""
3-extract-nuclear-transients.py

Raw Rhod-2 trace for every DAPI nucleus: map the step-1 cellpose nuclei into
calcium coordinates with the step-2 shift, then take each nucleus's mean pixel
value in every frame.

Adapted from two of Alex's scripts in clark-deep-phenotyping:
  - iPSC-CM/4-extract-calcium-transients.py: per-label mean via bincount,
    frame timing from the per-frame DateTime tags, QC/exclusion handling.
  - iPSC/3-2-stain-anchored-calcium.py: the nucleus "core" rule for a 3x3
    binned live camera, used instead of iPSC-CM's nearest-neighbour resize.

Why the core rule: iPSC-CM's step 4 shrinks the 2048 px mask to 682 px by
nearest neighbour, which is fine for whole cells. A nucleus here is only ~6
calcium pixels across (median ~340 full-res px^2), so most of its calcium pixels
straddle the nucleus edge and mix in background or a neighbour. A calcium pixel
covers exactly a 3x3 block of stain pixels (0.678 vs 2.034 um/px, same
objective); it is kept for a nucleus only if >= 7 of those 9 carry that label.
Nuclei whose core is under 27 full-res px^2 (3 calcium px) are dropped.

Geometry (same convention as step 2): the step-2 shift (X_shift, Y_shift) is in
calcium pixels; x3 gives the full-res shift. Calcium pixel (r, c) covers stain
pixels rows 3r-3*Y_shift .. +2, cols 3c-3*X_shift .. +2.

Other rules:
  - Nuclei with any core pixel inside the burned-in timestamp box (rows 0-15,
    cols 0-327 of every frame) are dropped: the text changes frame to frame.
  - Every nucleus gets a trace, cardiomyocyte or not; CM filtering on cTnT
    comes later and is cheap, re-reading 28 GB of video is not.
  - Wells failing step-2 QC (qc_pass False) or listed in experiment.yaml
    exclude_wells are skipped.
  - Frame interval = span of the DateTime tags / (n - 1). The tags have 1 ms
    resolution, so a median of frame-to-frame gaps reads a false 10.00 ms
    (Alex's finding); these movies read 9.830 ms = 101.73 Hz.

Writes: Calcium/<plate>/nuclear_traces.parquet   one row per nucleus, raw CalciumTrace (gitignored)
        Calcium/<plate>/traces_summary.csv       one row per well
        Calcium/<plate>/roi_qc/<well>.png        nucleus cores on the Rhod-2 movie

Usage:
    python 3-extract-nuclear-transients.py plate1
    python 3-extract-nuclear-transients.py plate2 --wells B02 C05
"""

import argparse
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage.segmentation import find_boundaries

import pipeline

BIN = 3               # calcium px = BIN x BIN stain px
N_LIVE = 682          # calcium frame side length
WIN = 3               # purity window, full-res px per side (Alex's)
PURITY = 7            # of WIN*WIN full-res pixels that must share the label (Alex's)
MIN_AREA = 27         # full-res px^2 of core; 3 calcium px at BIN=3 (Alex's)
TEXT_BOX = (slice(0, 16), slice(0, 328))   # burned-in timestamp, calcium px
QC_CROP = 200         # calcium px side of the QC zoom


# --- Nucleus cores on the calcium grid (from iPSC/3-2-stain-anchored-calcium.py)
def shifted_labels(stain_lab, shift_full, N):
    """L[r, c] = stain_lab[r - dy, c - dx] on an N x N full-res canvas."""
    dy, dx = shift_full
    L = np.zeros((N, N), stain_lab.dtype)
    r0, r1 = max(0, dy), min(N, stain_lab.shape[0] + dy)
    c0, c1 = max(0, dx), min(N, stain_lab.shape[1] + dx)
    L[r0:r1, c0:c1] = stain_lab[r0 - dy:r1 - dy, c0 - dx:c1 - dx]
    return L


def core_labels(stain_lab, shift_full, n=N_LIVE, b=BIN):
    """Calcium-grid label image of nucleus cores (>= PURITY of each b x b block agree)."""
    L = shifted_labels(stain_lab, shift_full, n * b)
    blocks = L.reshape(n, b, n, b).transpose(0, 2, 1, 3).reshape(n, n, b * b)
    centre = blocks[..., (b * b) // 2]
    agree = (blocks == centre[..., None]).sum(-1)
    return np.where((centre > 0) & (agree >= PURITY), centre, 0)


def build_rois(stain_lab, X_shift, Y_shift):
    """Core label image + per-nucleus table (geometry, why dropped)."""
    shift_full = (int(Y_shift) * BIN, int(X_shift) * BIN)
    core = core_labels(stain_lab, shift_full)
    n_stain = int(stain_lab.max())

    cnt = np.bincount(core.ravel(), minlength=n_stain + 1)
    in_text = np.unique(core[TEXT_BOX])
    keep = np.flatnonzero(cnt >= int(np.ceil(MIN_AREA / BIN ** 2)))
    keep = keep[(keep > 0) & ~np.isin(keep, in_text)]
    core = np.where(np.isin(core, keep), core, 0)

    lab = stain_lab.ravel()
    rows, cols = np.divmod(np.arange(lab.size), stain_lab.shape[1])
    area = np.bincount(lab, minlength=n_stain + 1)
    sy = np.bincount(lab, weights=rows, minlength=n_stain + 1)[keep] / area[keep]
    sx = np.bincount(lab, weights=cols, minlength=n_stain + 1)[keep] / area[keep]
    rois = pd.DataFrame({
        "ObjectNumber": keep,                          # = cellpose label in the step-1 mask
        "stain_x": sx, "stain_y": sy,                  # full-res stain px (for cTnT/MYL2/ANP later)
        "live_x": (sx + shift_full[1] - (BIN - 1) / 2) / BIN,
        "live_y": (sy + shift_full[0] - (BIN - 1) / 2) / BIN,
        "area_full_px": area[keep],
        "n_core_px": cnt[keep],
        "core_frac": cnt[keep] * BIN ** 2 / area[keep],
    })
    in_field = np.unique(shifted_labels(stain_lab, shift_full, N_LIVE * BIN))
    stats = dict(n_stain=n_stain, n_in_field=int((in_field > 0).sum()),
                 n_text_dropped=int((in_text > 0).sum()), n_kept=len(keep))
    return core, rois, stats


# --- Video -------------------------------------------------------------------
def traces_and_times(path, core, labels):
    """Mean pixel per label per frame (bincount: one pass per frame for all
    nuclei), plus seconds-from-frame-0 from the DateTime tags. Frames are read
    one page at a time so the ~935 MB stack is never held in memory."""
    flat = core.ravel()
    counts = np.bincount(flat, minlength=int(labels.max()) + 1)
    with tifffile.TiffFile(path) as tif:
        n = len(tif.pages)
        out = np.empty((labels.size, n), dtype=np.float32)
        ts, mean_img = [], np.zeros(core.shape, np.float64)
        for t, page in enumerate(tif.pages):
            frame = page.asarray()
            sums = np.bincount(flat, weights=frame.ravel(), minlength=counts.size)
            out[:, t] = sums[labels] / counts[labels]
            mean_img += frame
            ts.append(datetime.strptime(page.tags["DateTime"].value, "%Y%m%d %H:%M:%S.%f"))
    t_s = np.array([(x - ts[0]).total_seconds() for x in ts])
    dt = (t_s[-1] - t_s[0]) / (len(t_s) - 1)
    return out, dt, n, mean_img / n


# --- QC ----------------------------------------------------------------------
def norm(img, lo=1, hi=99.5):
    p_lo, p_hi = np.percentile(img, [lo, hi])
    return np.clip((img - p_lo) / (p_hi - p_lo + 1e-8), 0, 1)


def save_qc(title, mean_img, core, path):
    """Mean Rhod-2 projection, full field and a centre zoom, core ROIs in yellow."""
    img = mean_img.copy()
    img[TEXT_BOX] = np.median(img)
    rgb = np.repeat(norm(img)[..., None], 3, axis=-1)
    rgb[find_boundaries(core, mode="inner")] = [1, 1, 0]
    y0 = x0 = (N_LIVE - QC_CROP) // 2
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    ax[0].imshow(rgb)
    ax[0].add_patch(plt.Rectangle((x0, y0), QC_CROP, QC_CROP, fill=False, ec="cyan", lw=1))
    ax[0].add_patch(plt.Rectangle((0, 0), TEXT_BOX[1].stop, TEXT_BOX[0].stop,
                                  fill=False, ec="red", lw=1, ls="--"))
    ax[0].set_title(f"{title}: mean Rhod-2 + nucleus cores")
    ax[1].imshow(rgb[y0:y0 + QC_CROP, x0:x0 + QC_CROP], interpolation="nearest")
    ax[1].set_title(f"centre {QC_CROP} x {QC_CROP} calcium px")
    for a in ax:
        a.axis("off")
    plt.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --- Driver ------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plate", help="plate1 or plate2")
    p.add_argument("--wells", nargs="*", default=None)
    args = p.parse_args()

    cfg = pipeline.load_config()
    pcfg = pipeline.plate_cfg(cfg, args.plate)
    videos = pipeline.video_files(cfg, pcfg)
    mask_dir = pipeline.OUT_ROOT / args.plate / "cellpose"
    align_csv = pipeline.OUT_ROOT / args.plate / "alignment" / "alignment_mi.csv"
    if not align_csv.exists():
        pipeline.die(f"{align_csv} not found -- run 2-align-mi.py first")
    align = pd.read_csv(align_csv).set_index("well")

    failed = align.index[~align["qc_pass"].astype(bool)].tolist()
    if failed:
        print(f"  skipping {len(failed)} well(s) that failed step-2 QC: {failed}")
    exclude = pipeline.excluded(cfg, args.plate)
    for w, why in exclude.items():
        print(f"  excluding {w} per experiment.yaml: {why}")
    usable = [w for w in align.index if w not in failed and w not in exclude]
    wells = pipeline.select_wells(usable, args.wells)
    qc_dir = pipeline.out_dir(args.plate, "roi_qc")
    print(f"{args.plate}: extracting {len(wells)} wells")

    tables, summary = [], []
    for well in wells:
        mask_path = mask_dir / f"{well}_nuc_mask.tif"
        if not mask_path.exists():
            pipeline.die(f"{mask_path} not found -- run 1-cellpose-nuclei.py first")
        a = align.loc[well]
        core, rois, st = build_rois(tifffile.imread(mask_path), a.X_shift, a.Y_shift)
        tr, dt, n_frames, mean_img = traces_and_times(videos[well], core,
                                                      rois.ObjectNumber.to_numpy())
        save_qc(f"{args.plate} {well}", mean_img, core, qc_dir / f"{well}.png")

        rois.insert(0, "Plate", args.plate)
        rois.insert(1, "Well", well)
        rois.insert(2, "Condition", pipeline.condition(cfg, well))
        rois["FrameInterval_s"] = dt
        rois["CalciumTrace"] = list(tr)
        tables.append(rois)
        summary.append(dict(well=well, condition=pipeline.condition(cfg, well),
                            X_shift=a.X_shift, Y_shift=a.Y_shift, **st,
                            pct_kept=round(100 * st["n_kept"] / max(st["n_stain"], 1), 1),
                            median_core_px=float(rois.n_core_px.median()),
                            core_frac_med=round(float(rois.core_frac.median()), 3),
                            n_frames=n_frames, frame_interval_ms=round(dt * 1e3, 4)))
        print(f"  {well}: shift=({a.X_shift:>4},{a.Y_shift:>4})  stain={st['n_stain']}  "
              f"in field={st['n_in_field']}  kept={st['n_kept']}  "
              f"(text box dropped {st['n_text_dropped']})  "
              f"core px median={rois.n_core_px.median():.0f}  "
              f"{n_frames} frames @ {1/dt:.2f} Hz")

    out = pd.concat(tables, ignore_index=True)
    pq = pipeline.OUT_ROOT / args.plate / "nuclear_traces.parquet"
    sm = pipeline.OUT_ROOT / args.plate / "traces_summary.csv"
    if args.wells and pq.exists():
        old = pd.read_parquet(pq)
        out = pd.concat([old[~old.Well.isin(wells)], out], ignore_index=True)
        olds = pd.read_csv(sm)
        summary = pd.concat([olds[~olds.well.isin(wells)], pd.DataFrame(summary)])
    summary = pd.DataFrame(summary).sort_values("well")
    out.sort_values(["Well", "ObjectNumber"]).to_parquet(pq, engine="pyarrow", index=False)
    summary.to_csv(sm, index=False)
    v = summary.frame_interval_ms
    if v.max() - v.min() > 0.05:
        print(f"  WARNING: wells disagree on frame interval: {v.min():.3f}-{v.max():.3f} ms")
    print(f"\nWrote {pq}  ({len(out)} nuclei)")
    print(f"Wrote {sm}")
    print(summary.drop(columns=["X_shift", "Y_shift"]).to_string(index=False))


if __name__ == "__main__":
    main()
