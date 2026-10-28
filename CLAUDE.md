# CARBONITE-vCMs Analysis

Analysis code for a project studying how **maturation media** affects the development of
**ventricular hiPSC-derived cardiomyocytes (vCMs)**. Owner: Lavie Ngo (Saucerman Lab).

## Background references

- **Maturation media:** Feyen et al. 2020, Cell Reports — https://pubmed.ncbi.nlm.nih.gov/32697997/
  ("Feyen media" in this repo means this metabolic maturation medium).
- **Assay (CARBONITE):** https://pmc.ncbi.nlm.nih.gov/articles/PMC12918877/ — this is the source of the
  assay design and readouts. In short: live Rhod-2 calcium movies give single-cell
  transients, and the same wells are then fixed and immunostained.
  - Stains: **MYL2** (ventricular marker), **ANP** (NPPA; immaturity/stress marker),
    **cTnT** (CM identity/sarcomeres), **DAPI** (nuclei).

## Pilot experiment design (20261009 experiments)

- Differentiation D0: 2026-08-31. On D18 (2026-09-18), cells were replated into two 96-well plates
  and kept in maintenance media until D25 (2026-09-25). On D25, wells were moved to
  their assigned media.
- Both plates use the same 2×8 layout (columns 02–09):
  - **Row B (B02–B09):** maintenance media (control)
  - **Row C (C02–C09):** Feyen maturation media (treatment)
- **Plate 1:** 1 week in the assigned media. Calcium imaging (CARBONITE) on 2026-10-02, then fixed.
- **Plate 2:** 2 weeks in the assigned media. Calcium imaging on 2026-10-09, then fixed.
- Both plates were stained and imaged on 2026-10-09.
- Biological unit for now: one well = one replicate (n = 8 per condition per plate).

## Data locations (outside this repo, read-only)

Parent folder: `/Volumes/SaucermanLab/Lavie/Experiments/CARBONITE-vCMs/`

| What | Plate 1 | Plate 2 |
|---|---|---|
| Fixed/stained images, per channel (`<well>_..._s1_w1..w4.tif`) | `20261009-Lavi-plate1/all/` | `20261009-Lavi-plate2/all/` |
| All-channel merged stack | `20261009-Lavi-plate1/20261009-CARBONITE-vCM-Plate1-all-channels-merged.tif` | `20261009-Lavi-plate2/20261009-CARBONITE-vCM-Plate2-all-channels-merged.tif` |
| Rhod-2 calcium movies (TIFF stacks, `*_Stream_<well>_s1_Texas Red_Low.tif`) | `20261009-Lavi-plate1/Lavi-CARBONITE-plate1-10-02-26/3188/` | `20261009-Lavi-plate2/LN-20261009-CARBONITE-Plate2/3243/` |
| Calcium movies (AVI exports, `*_Stream_<well>_s1_Texas Red_Low.avi`) | `20261009-Lavi-plate1/calcium movies/` | `20261009-Lavi-plate2/calcium movies/` |
| Raw MetaXpress output (`<date>/<plate id>/TimePoint_1`) | `20261009-Lavi-plate1/2026-10-10/3251/` | `20261009-Lavi-plate2/2026-10-10/3252/` |

- All paths in this table were confirmed by the user on 2026-10-28. Each plate folder also
  has a `thumb/` folder of thumbnails, which the analysis does not use.
- **Fixed-image channels:** `w1` = DAPI, `w2` = cTnT, `w3` = ANP, `w4` = MYL2.
- **Calcium movies:** Rhod-2 imaged in the Texas Red channel, one stack per well,
  **1000 frames at 102 fps** (about 9.8 s per movie). Use the frame rate to put the time
  axis in seconds.
  - Select stacks with the glob `*_Stream_<well>_s1_Texas Red_Low.tif`, which gives exactly 16 per plate
    (B02–B09, C02–C09). The naming differs by plate: Plate 1 files start with
    `Lavi-CARBONITE-plate1-10-02-26_` and Plate 2 files start with `Plate2_`.
  - Plate 1's `3188/` folder also holds two files from other experiments,
    `Plate1-9426-t55_D03_s1_w1.tif` and `Plate1-10226-t5min_D05_s1_w1.tif`. Exclude them.
- Never modify or overwrite raw data. Write outputs to a separate folder.

## Code

- `Code/20261009-CARBONITE-vCM-Maturation/` — notebooks for this experiment.
  - `...-montages.ipynb`: RGB channel-merge montages per condition. Adapted from the
    montage notebooks in the IRI experiments
    (`/Volumes/SaucermanLab/Lavie/LRP1/Experiments/IRI-experiments/`). Keep the same style:
    a reusable `create_montage()` function, `CHANNEL_CONFIG` to map channels to colors,
    an explicit well→condition dict, and `downsample` for display.
- Stack: Python in Jupyter notebooks, with `numpy`, `matplotlib`, `tifffile`, and `pathlib`.
- Notebook style: start with markdown cells for the title, author, and experiment info. Then
  reusable functions, then short calling cells. Use retina inline figures.

## Analysis goals (current round)

1. **Montages:** merged RGB montages of the fixed stains for each plate and condition, in the
   same style as the IRI montage notebooks.
2. **Calcium transients:** automated version of ImageJ's "Plot Z-axis Profile" (mean
   intensity per frame) for every well and replicate, plotted in a notebook. Single-cell
   or ROI-based extraction may follow.
   - Likely metrics, still to be decided: F/F0 normalization, beat rate, amplitude,
     time to peak, decay (e.g., time to 50%/90% decay or tau), and transient duration.
3. **Main comparison:** Plate 1 (1 week) vs. Plate 2 (2 weeks), maintenance vs. Feyen, using
   simple plots and interpretation. The results will guide future maturation
   experiments (e.g., maturation before ischemia-reperfusion injury assays).
