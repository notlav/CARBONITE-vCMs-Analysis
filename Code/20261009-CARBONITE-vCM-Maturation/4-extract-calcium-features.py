#!/usr/bin/env python3
"""
4-extract-calcium-features.py

Turn each nucleus's raw Rhod-2 trace into dF/F0 features.

Adapted from clark-deep-phenotyping/iPSC-CM/5-extract-calcium-features.py
(Alex). The signal maths lives in calcium.py (copied from his, unchanged), with
his defaults:

    trim 1 frame each end   frame 0 is artifactual in a few % of ROIs
    F0 = p10                dF/F0 = (F - F0) / F0; a beating cell sits near
                            baseline most of the time
    smooth 3 frames
    peak height 5x noise    per nucleus, against its own noise, not a fixed
                            dF/F0, so low-amplitude conditions are not reported
                            as "beating less"
    time_to_peak            from max upstroke velocity (CARBONITE1 definition)

Frame timing is each well's own FrameInterval_s from step 3 (101.73 Hz here).
Nuclei with < 2 detected transients get NaN kinetics and beating = False rather
than being dropped: "does not beat" is a phenotype. Non-CM nuclei are still in
the table; cTnT-based CM filtering comes later.

Reads:  Calcium/<plate>/nuclear_traces.parquet
Writes: Calcium/<plate>/calcium_features.csv   one row per nucleus, no traces

Usage:
    python 4-extract-calcium-features.py plate1
    python 4-extract-calcium-features.py plate2 --snr 4
"""

import argparse

import pandas as pd

import calcium
import pipeline

ROI_COLS = ["stain_x", "stain_y", "live_x", "live_y", "area_full_px", "n_core_px",
            "core_frac", "FrameInterval_s"]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plate", help="plate1 or plate2")
    p.add_argument("--f0-pct", type=float, default=calcium.F0_PCT)
    p.add_argument("--smooth", type=int, default=calcium.SMOOTH)
    p.add_argument("--trim", type=int, default=calcium.TRIM)
    p.add_argument("--snr", type=float, default=calcium.SNR_K)
    p.add_argument("--floor", type=float, default=calcium.FLOOR)
    p.add_argument("--min-sep", type=float, default=calcium.MIN_SEP)
    p.add_argument("--width", type=float, default=calcium.WIDTH)
    args = p.parse_args()

    cfg = pipeline.load_config()
    src = pipeline.OUT_ROOT / args.plate / "nuclear_traces.parquet"
    if not src.exists():
        pipeline.die(f"{src} not found -- run 3-extract-nuclear-transients.py first")
    df = pd.read_parquet(src)

    exclude = pipeline.excluded(cfg, args.plate)
    hit = sorted(set(df.Well) & set(exclude))
    for w in hit:
        print(f"  excluding {w} per experiment.yaml: {exclude[w]}")
    df = df[~df.Well.isin(exclude)].reset_index(drop=True)
    if df.empty:
        pipeline.die("no nuclei left after exclusions")
    print(f"{args.plate}: {len(df)} nuclei, {df.Well.nunique()} wells, "
          f"dt={df.FrameInterval_s.mean()*1000:.3f} ms")

    feat = calcium.features(df, args.f0_pct, args.smooth, args.snr, args.floor,
                            args.min_sep, args.width, args.trim)
    # features() returns rows in df order; join positionally (ObjectNumber is
    # only unique within a well).
    if len(feat) != len(df):
        pipeline.die(f"feature rows {len(feat)} != nucleus rows {len(df)}")
    out = feat.join(df[ROI_COLS])

    dst = pipeline.OUT_ROOT / args.plate / "calcium_features.csv"
    out.to_csv(dst, index=False)
    print(f"  beating (>=2 transients): {out.beating.sum()}/{len(out)} "
          f"({100*out.beating.mean():.0f}%)")
    print(f"\nWrote {dst}")
    print(out.groupby("Condition").agg(
        nuclei=("beating", "size"), beating=("beating", "mean"),
        beat_rate_hz=("beat_rate", "median"), amplitude=("amplitude", "median"),
        ttp_ms=("time_to_peak", lambda s: 1000 * s.median()),
        d50_ms=("duration_50", lambda s: 1000 * s.median()),
    ).round(3).to_string())


if __name__ == "__main__":
    main()
