"""
calcium.py

Calcium trace -> dF/F0 -> per-transient features, for the 20261009
CARBONITE-vCM per-nucleus pipeline. Shared by 4-extract-calcium-features.py
(the pipeline) and 5-plot-calcium.ipynb (plots), so a parameter chosen while
looking at plots is the same one the pipeline applies.

Copied from clark-deep-phenotyping/iPSC-CM/calcium.py (Alex) on 2026-10-28.
The signal maths is unchanged, so features stay comparable with his data.
Changes: this docstring, and features() carries a Plate column. Alex's notes below refer to his datasets
(icell-Dec25, PURA, NR2F2), where each default was settled.

Input is the RAW per-nucleus mean Rhod-2 intensity that
3-extract-nuclear-transients.py stores. dF/F0 = (F - F0) / F0 with F0 = the
10th percentile of the (trimmed) trace; no background subtraction (Alex's
convention).

Every parameter here was settled by looking at plots, not assumed. See the
docstrings for what each one is defending against.
"""

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d
from scipy import stats
from scipy.signal import find_peaks

# Defaults, settled in explore-calcium.py against icell-Dec25.
F0_PCT = 10          # percentile for F0
SMOOTH = 3           # moving-average frames
TRIM = 1             # frames dropped each end (frame 0 is artifactual in ~4% of cells)
SNR_K = 5.0          # peak height as a multiple of each cell's own noise
FLOOR = 0.02         # absolute dF/F floor
MIN_SEP = 0.5        # min peak separation (s)
WIDTH = 0.07         # min peak width (s)


# --- Signal ------------------------------------------------------------------
def dff(trace, f0_pct, smooth, trim=1):
    """dF/F against a percentile baseline, optionally smoothed.

    F0 as a low percentile rather than the mean: a beating cell spends most of
    its time near baseline, but the mean is dragged up by the transients it is
    supposed to be measured against.

    trim drops that many frames from each end first. The first frame carries an
    artifact in ~4% of cells -- for dim cells a small blip over a tiny F0
    explodes into dF/F of 100+. Dropping the two edge frames takes the maximum
    peak dF/F across the plate from 110.6 to 0.8 while moving the median from
    0.200 to 0.198: the entire extreme tail was frame-0, and nothing real is
    lost. The retired step 4 handled this by overwriting frame 0 with the
    baseline percentile; trimming is preferred because it discards a bad
    measurement rather than fabricating a good one.
    """
    t = np.asarray(trace, dtype=float)
    if trim:
        t = t[trim:len(t) - trim]
    f0 = np.percentile(t, f0_pct)
    d = (t - f0) / max(f0, 1e-9)
    if smooth and smooth > 1:
        d = uniform_filter1d(d, size=smooth)
    return d


def noise_floor(d, smooth=SMOOTH):
    """Noise sd of the trace, robust to the transient.

    THIRD ATTEMPT, and the first two were both biased by the transient they were
    supposed to ignore:

      sd of the sub-median half   inflated for SLOW cells, whose lower half
                                  contains the rise and fall
      sd of frame-to-frame diffs  inflated for FAST cells, because a fast
                                  upstroke IS a large frame-to-frame change

    The second is worse than it looks. Simulated at the real regime (noise sd
    0.005, transient 0.5, dt 9.8 ms) sd(diff) returns 0.0141 for a 50 ms upstroke
    and 0.0065 for a 350 ms one -- it varies 2.1x with the very feature it must be
    independent of, and reads 1.3-2.8x the true noise throughout. In the data it
    put the fast-vs-slow SNR delta at -0.891; with MAD the same comparison is
    +0.046, i.e. there is no real SNR difference at all.

    MAD of the diffs fixes it. The upstroke contributes a handful of large diffs
    out of ~800; a median absolute deviation ignores them by construction, while
    an sd squares them. Same simulation: MAD returns 0.0036-0.0039 across every
    upstroke speed, an 11% spread against sd's 210%.

    scale="normal" converts MAD to an sd-equivalent for Gaussian noise; /sqrt(2)
    because differencing doubles the variance; sqrt(smooth) undoes dff()'s
    smoothing (for a uniform filter of width w, sd(diff) = s*sqrt(2)/w while
    sd(smoothed) = s/sqrt(w), so the ratio is exactly sqrt(w)).
    """
    if len(d) < 3:
        return float("nan")
    mad = stats.median_abs_deviation(np.diff(d), scale="normal")
    return float(mad / np.sqrt(2) * np.sqrt(max(smooth, 1)))


def peak_height(d, snr_k, abs_floor, smooth=SMOOTH):
    """Detection threshold for THIS cell: k x its own noise, with a floor.

    Deliberately not a fixed dF/F. Amplitude varies by condition -- PURA's
    median peak is 0.148 vs 0.239 for NR2F2 -- so a fixed bar culls the
    small-amplitude conditions and reports them as "beating less" when they are
    beating fine with smaller transients. Measured against each cell's own
    noise, PURA beats as often as Untreated (96% vs 96%); at a fixed 0.075 it
    appears to beat less (86% vs 92%). The absolute floor stops a dead-flat
    trace from detecting its own noise.

    This replaces the retired per-condition thresholds (PURA 0.01, SP1 0.03,
    else 0.075), which were the same problem patched by hand -- and PURA's 0.01
    sat at 2x its noise, so it fired on noise.
    """
    return max(snr_k * noise_floor(d, smooth), abs_floor)


def detect(d, dt, height, min_sep_s, width_s):
    """Peak detection in seconds, converted to frames using this well's dt."""
    return find_peaks(d, height=height,
                      distance=max(1, int(round(min_sep_s / dt))),
                      width=max(1, int(round(width_s / dt))))[0]


# time_to_peak follows CARBONITE1-archive/process-images.py: the transient starts
# at the point of maximum upstroke VELOCITY, not at a fraction of amplitude.
#
#     diff_peaks, _ = find_peaks(np.diff(ca[peaks[0]-200:peaks[0]]),
#                                distance=200, height=0.02)
#     ttp = t[peaks[0]] - t[diff_peaks[0]]
#
# That measures the upper part of the rise, so it is systematically shorter than
# a 10%-of-amplitude definition. It is kept because the fast_rise (<0.2 s) and
# dome (>0.33 s) cutoffs in fig-cal-groups.py are calibrated to it: change the
# definition and those numbers silently stop meaning what they meant.
#
# The one translation: the original ran on a min-max normalised 0-1 trace, where
# height=0.02 is 2% of the cell's full range per frame. Here the trace is dF/F,
# so the equivalent threshold scales with this transient's own amplitude. dF/F is
# kept over min-max because min-max normalises amplitude away, and amplitude is
# what exposed NR2F2.
TTP_WINDOW_S = 2.0     # the original's 200 frames at its assumed 100 Hz
TTP_DV_FRAC = 0.02     # the original's height=0.02 on a 0-1 trace


def upstroke_start(d, dt, pk, lo, amp):
    """Index of maximum upstroke velocity before peak `pk`, or None."""
    w = int(round(TTP_WINDOW_S / dt))
    # Bounded by the previous peak as well as the window. The original used a
    # fixed 200-frame lookback, which can reach back into the PREVIOUS
    # transient's upstroke -- and with distance=200 only the tallest derivative
    # peak survives, so a bigger preceding transient would win and the start
    # would be taken from the wrong beat.
    a = max(lo, pk - w)
    if pk - a < 3:
        return None
    dv = np.diff(d[a:pk + 1])
    hits, _ = find_peaks(dv, height=TTP_DV_FRAC * amp, distance=max(1, w))
    return (a + hits[0]) if len(hits) else None


def kinetics(d, dt, peaks):
    """Per-transient shape. Returns medians across a cell's transients.

    Reported separately from whether the cell beats, so amplitude and rhythm are
    never conflated -- that distinction is the whole point of the SNR rule above.

    For each peak, the baseline is the minimum between the previous peak and
    this one, so a transient riding on a neighbour's decay is measured against
    where it actually started.

      amplitude     transient rise above its OWN local baseline, median across
                    the cell's transients. Distinct from peak_dff, which is a
                    single max over the whole trace against the global F0:
                    amplitude is robust (a median over transients) and immune to
                    a lone artifactual frame.
      time_to_peak  from the point of maximum upstroke velocity to the peak,
                    matching CARBONITE1-archive/process-images.py so its
                    fast_rise (<0.2 s) / dome (>0.33 s) cutoffs still apply.
                    See upstroke_start(). NOTE this is a different instant from
                    the 10%-of-amplitude crossing that the transient overlays
                    align on -- the overlays show the whole rise, this measures
                    its upper part.
      duration_X    width at (100-X)% of amplitude above baseline: duration_50
                    is width at half maximum, duration_90 is 90% recovered
      tau_decay     exponential fit to the decay, log-linear
    """
    amps, ttp, d20, d50, d90, taus = [], [], [], [], [], []
    for i, pk in enumerate(peaks):
        lo = peaks[i - 1] if i > 0 else 0
        hi = peaks[i + 1] if i + 1 < len(peaks) else len(d)
        pre = d[lo:pk + 1]
        if len(pre) < 3:
            continue
        base = float(pre.min())
        amp = float(d[pk]) - base
        if amp <= 0:
            continue
        amps.append(amp)

        start = upstroke_start(d, dt, pk, lo, amp)
        if start is not None:
            ttp.append((pk - start) * dt)

        # All three levels or none: taking medians over different subsets lets a
        # transient contribute its d50 but not its d90 (when the next peak
        # arrives before it decays to 10%), which produced the impossible
        # d90 < d50. Durations must describe the same transients to be compared.
        post = d[pk:hi]
        widths = {}
        for frac in (0.2, 0.5, 0.9):
            lvl = base + (1 - frac) * amp
            l = np.where(pre <= lvl)[0]
            r = np.where(post <= lvl)[0]
            if len(l) and len(r):
                widths[frac] = ((pk + r[0]) - (lo + l[-1])) * dt
        if len(widths) == 3:
            d20.append(widths[0.2]); d50.append(widths[0.5]); d90.append(widths[0.9])

        # decay: peak -> 10% of amplitude, fit log(y-base) vs t
        end = np.where(post <= base + 0.1 * amp)[0]
        seg = post[:end[0] + 1] if len(end) else post
        y = seg - base
        ok = y > 1e-6
        if ok.sum() >= 5:
            t = np.arange(len(seg))[ok] * dt
            slope = np.polyfit(t, np.log(y[ok]), 1)[0]
            if slope < 0:
                taus.append(-1.0 / slope)

    med = lambda a: float(np.median(a)) if a else np.nan
    return dict(amplitude=med(amps),
                time_to_peak=med(ttp), duration_20=med(d20), duration_50=med(d50),
                duration_90=med(d90), tau_decay=med(taus),
                d50_over_d90=(med(d50) / med(d90)) if (d50 and d90 and med(d90)) else np.nan)


def features(df, f0_pct, smooth, snr_k, abs_floor, min_sep_s, width_s, trim=1):
    rows = []
    for r in df.itertuples():
        dt = float(r.FrameInterval_s)
        d = dff(r.CalciumTrace, f0_pct, smooth, trim)
        nf = noise_floor(d, smooth)
        h = peak_height(d, snr_k, abs_floor, smooth)
        pk = detect(d, dt, h, min_sep_s, width_s)
        dur = (len(d) - 1) * dt
        rec = dict(
            Plate=getattr(r, "Plate", None),
            Well=r.Well, ObjectNumber=r.ObjectNumber, Condition=r.Condition,
            row=r.Well[0], col=int(r.Well[1:]),
            n_peaks=len(pk),
            beat_rate=len(pk) / dur,
            peak_dff=float(d.max()),
            baseline_noise=nf,
            snr=float(d.max() / max(nf, 1e-9)),
            threshold_used=h,
        )
        rec.update(kinetics(d, dt, pk))
        rows.append(rec)
    f = pd.DataFrame(rows)
    f["beating"] = f.n_peaks >= 2
    return f


