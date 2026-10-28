# CARBONITE-vCMs Analysis

Analysis code for a project studying how **maturation media** affects the development of
**ventricular hiPSC-derived cardiomyocytes (vCMs)**. Owner: Lavie Ngo (Saucerman Lab).

## Background references

- **Maturation media:** Feyen et al. 2020, Cell Reports — https://pubmed.ncbi.nlm.nih.gov/32697997/
  ("Feyen media" in this repo means this metabolic maturation medium).
- **Assay (CARBONITE):** https://pmc.ncbi.nlm.nih.gov/articles/PMC12918877/ — this is the source of the
  assay design and readouts. In short: live Rhod-2 calcium movies give single-cell
  transients, and the same wells are then fixed and immunostained.
  - Stains: **MYL2** (ventricular marker), **ANP** (NPPA; atrial marker),
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

Parent folder:
- **Cluster (Rivanna/Tardis; the calcium pipeline runs here):** `/sfs/gpfs/tardis/project/SaucermanLab/Lavie/Experiments/CARBONITE-vCMs/`
- **Mac (mounted share; montages were run here):** `/Volumes/SaucermanLab/Lavie/Experiments/CARBONITE-vCMs/`

Both point to the same files; the relative paths below apply to both. This repo is
`CARBONITE-vCMs-Analysis/` inside that parent folder.

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
- **Calcium movies:** Rhod-2 imaged in the Texas Red channel, one stack per well, 682 × 682 px,
  uint16, **1000 frames at 101.73 Hz** (9.830 ms/frame, 9.82 s per movie). Nominal rate is
  102 fps; the true rate comes from the per-frame `DateTime` tags (span / (n − 1); the tags
  have 1 ms resolution, so a median of frame gaps falsely reads 10.00 ms).
  - The AVI exports are the same 1000 frames at 102 fps (9.8 s). They look ~30 s long only
    because movie players cannot display 102 fps and play them in slow motion.
  - **Geometry:** stain images are 2048 px at 0.678 µm/px; calcium movies are 682 px at
    2.034 µm/px (3×3 camera binning), same 10X Plan Apo objective, so one calcium pixel =
    exactly a 3×3 block of stain pixels.
  - Select stacks with the glob `*_Stream_<well>_s1_Texas Red_Low.tif`, which gives exactly 16 per plate
    (B02–B09, C02–C09). The naming differs by plate: Plate 1 files start with
    `Lavi-CARBONITE-plate1-10-02-26_` and Plate 2 files start with `Plate2_`.
  - Plate 1's `3188/` folder also holds two files from other experiments,
    `Plate1-9426-t55_D03_s1_w1.tif` and `Plate1-10226-t5min_D05_s1_w1.tif`. Exclude them.
- Never modify or overwrite raw data. Write outputs to a separate folder.

### Known data issues

- **Plate 2, wells C05 and C06 (Feyen, 2 weeks):** hair or debris is visible in the fixed/stained
  images (C06: bright streak; C05: thin fiber). The C05 fiber is also visible in the Rhod-2 movie
  (mean projection, `Calcium/plate2/roi_qc/C05.png`). Both wells aligned normally and are kept,
  but they are flagged (hollow points) in every calcium plot. C05 is the low-amplitude Feyen outlier.
- **Plate 1, well C07 (Feyen, 1 week): EXCLUDED from calcium analysis** (recorded in
  `experiment.yaml`). Its Rhod-2 field is a uniform confluent sheet with nothing matching its
  sparse cTnT patches, so alignment failed: MI ~0.005–0.01 at every shift, for both frame 0 and
  the movie SD projection, and the fit was 295 px off the plate plane. Nuclei cannot be placed on
  the movie. Plate 1 Feyen is n = 7.
- **Burned-in timestamp text** in every Rhod-2 frame: rows 0–15, cols 0–327 (calcium px), pixel
  values 0/65535, changing every frame. Masked for alignment; nuclei under it are dropped.
- **Stain↔calcium stage offset:** the fixed image sits about (X −20…0, Y +113…+120) calcium px
  from the movie. Shifts lie on a plate plane (X with column, Y with row) within ~1 px.
  Consequence: the top ~117 rows (~17%) of every movie have no nuclei, and 70–90% of each
  well's segmented nuclei land inside the movie.
- **MYL2 and cTnT display ranges:** the ImageJ ranges saved in the merged stacks set the minimum too
  high (MYL2 15424-20855 shows only ~0.1-3% of pixels). For the single-channel montages, these two
  channels use one range pooled over all 32 wells (min = p50, max = p99.9), the same for both plates.
  ANP and DAPI keep the ImageJ ranges. The 4-channel merged montages use the same ranges with
  DAPI blue, cTnT green, ANP red, MYL2 magenta (not the cyan/yellow LUTs saved in the stacks).
  This is for display only; no intensity quantification has been done.

## Code

- `Code/20261009-CARBONITE-vCM-Maturation/` — notebooks for this experiment.
  - `...-montages.ipynb`: RGB channel-merge montages per condition. Adapted from the
    montage notebooks in the IRI experiments
    (`/Volumes/SaucermanLab/Lavie/LRP1/Experiments/IRI-experiments/`). Keep the same style:
    a reusable `create_montage()` function, `CHANNEL_CONFIG` to map channels to colors,
    an explicit well→condition dict, and `downsample` for display.
  - **Per-nucleus calcium pipeline** (built 2026-10-28). Steps 1–4 are scripts, step 5 is
    a notebook. Adapted from Alex's pipeline in
    `/sfs/gpfs/tardis/project/SaucermanLab/Alex/clark-deep-phenotyping/iPSC-CM/` (plus the
    nuclei Cellpose settings and the "core" rule from its `iPSC/` sibling). Each file's header
    names its source. **Never run, import, or write anything in Alex's directory.** Copy or adapt
    into this repo only.

    | Step | File | What it does |
    |---|---|---|
    | config | `experiment.yaml` | Paths, row B/C → condition, `exclude_wells` with reasons |
    | shared | `pipeline.py` | Config loading and file discovery (one file per well or stop) |
    | 1 | `1-cellpose-nuclei.py` | Cellpose on DAPI (w1), Alex's settings: cellprob −2.0, flow 0.8, diameter None; Cellpose 4.2.1 (CPSAM model). GPU. |
    | 2 | `2-align-mi.py` | cTnT (w2) resized to 682 → integer shift onto Rhod-2 frame 0 by mutual information; Theil-Sen plate-plane QC. Search widened to ±200 px (Alex's ±106 was too small). |
    | 3 | `3-extract-nuclear-transients.py` | Nuclei → calcium grid by the core rule (a calcium px belongs to a nucleus only if ≥7 of its 9 stain px share the label; min 3 calcium px); mean Rhod-2 per nucleus per frame; dt from timestamps |
    | 4 | `4-extract-calcium-features.py` + `calcium.py` | ΔF/F0 = (F − F0)/F0, F0 = p10, trim 1 frame each end, smooth 3; peaks at 5× each nucleus's own noise (floor 0.02); metrics per nucleus. `calcium.py` is Alex's, maths unchanged. No camera-offset subtraction (Alex's convention). |
    | 5 | `5-plot-calcium.ipynb` | Example traces, plate-layout mean traces, heatmaps, spatial maps, peak ΔF/F0 distributions, per-well summaries |

    Run order: `python 1-cellpose-nuclei.py plate1` (likewise for steps 2, 3, 4), for each plate,
    from this folder. Each step takes `--wells` for a subset. Whole run is about 15 min on a GPU node.
  - **Environment (cluster):** conda env `carbonite-vcms` (`/home/lln5zuj/miniconda3/envs/carbonite-vcms`,
    Python 3.11, Cellpose 4.2.1, PyTorch with CUDA). Jupyter kernel "carbonite-vcms (cellpose)".
    Always set `PYTHONNOUSERSITE=1` (the env and kernel already do), because
    `~/.local/lib/python3.11` otherwise leaks 155 packages into it.
  - **Outputs:** `Calcium/<plate>/`:
    - `cellpose/<well>_nuc_overlay.png`: segmentation check
    - `alignment/overlays/<well>.png` and `alignment_mi.csv`: alignment check
    - `roi_qc/<well>.png`: nucleus cores on the mean Rhod-2 image
    - `traces_summary.csv` and `calcium_features.csv` (one row per nucleus)
    - `nuclear_traces.parquet`: raw traces
    - `cellpose/*_nuc_mask.tif`: label masks

    Figures from step 5 go in `Calcium/figures/`, with per-well tables in `Calcium/well_summary_*.csv`
    and logs in `Calcium/logs/`. Masks and parquet files are gitignored.
- Stack: Python in Jupyter notebooks, with `numpy`, `matplotlib`, `tifffile`, and `pathlib`.
- Notebook style: start with markdown cells for the title, author, and experiment info. Then
  reusable functions, then short calling cells. Use retina inline figures.

## Analysis goals (current round)

1. **Montages: done.** Merged and single-channel montages plus selected-well panels are in `Montages/`.
2. **Calcium transients: per-nucleus pipeline done for 31/32 wells** (Plate 1 C07 excluded).
   This replaced the original whole-field "Plot Z-axis Profile" plan, at the user's request.
   - **Next:** a cTnT-based cardiomyocyte filter (use stain cTnT around each nucleus). All DAPI
     nuclei are currently included, and most are flat non-CMs. Then review the beat criterion,
     since ≥2 peaks at 5× noise passes many flat nuclei. `MIN_PEAK_DFF = 0.2` in the notebook is a
     placeholder for the filter.
3. **Main comparison:** Plate 1 (1 week) vs. Plate 2 (2 weeks), maintenance vs. Feyen, using
   simple plots and interpretation. The results will guide future maturation
   experiments (e.g., maturation before ischemia-reperfusion injury assays).

## Current calcium results (exploratory, 2026-10-28)

Not CM-filtered, no statistics. "Active" means ≥2 peaks and peak ΔF/F0 ≥ 0.2. The numbers are
medians of per-well medians (from `Calcium/well_summary_active-peak0.2.csv`).

| | P1 Maint (n=8) | P1 Feyen (n=7) | P2 Maint (n=8) | P2 Feyen (n=8) |
|---|---|---|---|---|
| Active nuclei per well | 180 | **7** | 936 | 711 |
| % active | 8.6 | **0.3** | 38.6 | 35.1 |
| Beat rate (Hz) | 0.31 | 0.26 | 0.36 | 0.31 |
| Amplitude (ΔF/F0) | 0.19 | 0.10 | 0.10 | **0.25** |
| Duration at 50% / 90% (ms) | 219 / 637 | 118 / 546 | 182 / 467 | **452 / 1175** |

How to read these (and the caveats to repeat when interpreting):
- **Plate 1 Feyen is almost silent at 1 week.** Only C03 shows small transients; the other wells are
  flat with a slow upward drift. Its kinetics rest on 1–67 nuclei per well, so they are not
  trustworthy. Its 595 ms time to peak and the 14 s decay outlier are likely the drift. One
  hypothesis to test, not a conclusion: a transient suppression of beating after switching to the
  metabolic (low-glucose, fatty-acid) medium, recovered by 2 weeks.
- **At 2 weeks, Feyen wells beat about as often as maintenance** (similar % active), with **larger
  and longer transients** (amplitude ~2.5×, durations ~2.5×). Do not read "longer" as more
  mature without checking: slower kinetics are not the usual maturation signature, and beat rate
  and the CM filter can change this.
- **Amplitude vs. peak ΔF/F0 disagree on Plate 2 maintenance** (peak 0.38 vs. Feyen 0.31, but
  amplitude 0.10 vs. 0.25). Three maintenance wells (B04, B05, B07) beat at ~0.7–0.8 Hz. There,
  transients ride on incomplete decay, so the per-transient amplitude (rise above a local
  baseline) is small while the global peak sums up. Beat rate confounds amplitude and duration.
- **Beat rate is quantized** (2 beats in 9.82 s = 0.20 Hz, 3 = 0.31 Hz, …). Each nucleus has
  only 1–3 transients, and some are already underway at frame 0.
- Example nuclei with high signal can be artifacts (e.g., slow rises without decay in Plate 1 C04,
  possibly dying cells or motion). Check the traces before trusting extremes.
- Plate also differs in imaging day and cell age (~D32 vs. D39), so compare Feyen vs. maintenance
  within a plate first. Cross-plate differences mix time in media with session effects.
