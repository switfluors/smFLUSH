# Balanced Detection Method and MC-simulation

MATLAB code for single-molecule emission-spectrum **centroid** estimation by the
**Balanced Method**, plus a Monte-Carlo study of its precision and accuracy.
The repository has two finalized, independent scripts (in `Final/`):

| Script | Purpose |
|--------|---------|
| `BalancedDetection.m` | **Experiment.** Analyze a real ND2 movie: extract single-molecule spectra and compute their spectral centroids by the Balanced Method. |
| `MCsimulation_final.m` | **Simulation.** Monte-Carlo model of centroid precision/accuracy vs. photon budget and analysis-window choice. |
| `Alexa_Fluor_647.csv` | Reference emission spectrum used by the simulation (col 1 = wavelength nm, col 3 = intensity). |
| `bfmatlab/` | Bundled **Bio-Formats** MATLAB toolbox. Provides `bfopen`, used by the experiment script to read `.nd2` microscope stacks. Add it to the MATLAB path before running (see below). |

---

## 1. `BalancedDetection.m` — experimental analysis

Loads a multi-frame ND2 stack and processes it end to end:

1. Build a wavelength calibration from known reference lines.
2. Locate the spatial slit; estimate and subtract a global background.
3. Detect single-molecule events and extract their dispersed spectra onto a
   common 1 nm wavelength grid (400–900 nm).
4. Filter events by photon budget and refine the spectral centroid with the
   **iterative Balanced Method**: a fixed-width window centered on the previous
   round's centroid, with the window half-width set once by the 1/N method.
5. Produce a 3-panel summary figure — (a) spectra heatmap, (b) mean spectrum,
   (c) centroid histogram — and log results.

**Requirements:** the **Bio-Formats** toolbox (`bfmatlab`) on the MATLAB path.
The `bfmatlab/` folder is bundled in this repository — it is the Open Microscopy
Environment's MATLAB reader that lets MATLAB open proprietary `.nd2` files via
`bfopen`. Point `addpath` at wherever you place this folder.

**Edit before running** (top of the script):
```matlab
addpath('D:\bfmatlab')                              % Bio-Formats location (the bundled bfmatlab/ folder)
nd2FilePath    = 'D:\...\00.nd2';               % input movie
outputFilePath = 'D:\...\G1.xlsx';                  % Excel log; output folder derived from it
```
Other tunable knobs: `photonMin/photonMax` (event filter), `level_frac` and
`width_factor` (Balanced-Method window), `num`/`convTol` (iteration limit/stop).

**Outputs** (written to the output folder):
- `BalancedDetection.mat` — centroids `c`, photon counts `s`, spectra `spe`, `halfWidth`.
- `<file>.pdf` — the 3-panel figure.
- `<file>_summary.csv` — one-row per-file summary.
- `<file>.xlsx` — running log: sheet 2 = per-file summary, sheet 3 = per-round
  convergence trace.

---

## 2. `MCsimulation_final.m` — Monte-Carlo simulation

Convolves a reference dye spectrum with the microscope PSF, disperses it onto a
16-row detector, adds Poisson (signal + background) and Gaussian (readout) noise,
subtracts the background, and computes the intensity-weighted centroid over a
sliding analysis window. For each photon budget `Psig` it simulates **one** set
of noisy spectra and applies **all** Width × Center windows to it, so differences
across the result maps reflect the window choice alone.

- **Part 1 — Homogeneous emitters:** identical molecules; sweeps Width × Center ×
  `Psig`. → `Heatmaps_5Psig.pdf`.
- **Part 2 — Spectral heterogeneity:** adds per-molecule peak-position jitter
  `delta ~ N(0, sigma)`; repeats for several `sigma`.
  → `Figure4_Contrast_Plot.pdf`, `Figure5_Heatmap_Std.pdf`,
  `Figure5_Heatmap_Error.pdf`, and a sweet-spot precision plot.

Set `ENABLE_HETERO = false` to run Part 1 only (homogeneous; equivalent to the
`sigma = 0` case) and skip Part 2.

**Requirements:** **Image Processing** and **Statistics and Machine Learning**
toolboxes (`fspecial`, `poissrnd`).

**Key parameters** (top of the script):
```matlab
Psig_File_Vals     = 100:500:5000;    % signal photons
Center_Wavelengths = 550:5:800;       % analysis window center (nm)
Width_Vals         = 50:50:300;       % analysis window width (nm)
Sigma_Shift_Vals   = [0 0.5 1.5 2.5]; % peak-shift std (nm), Part 2
numSpectra         = 1000;            % molecules per condition
```

**Outputs:** `Std_3D` / `Std_4D` (centroid std, nm) and matching `PercentError`
arrays stay in the workspace; PDFs are written to the working directory. No
intermediate `.mat` files — everything runs in memory.

---

## Usage

```matlab
% Run from the repository root (so the .csv spectra are on the path):
run MCsimulation_final.m     % simulation
run BalancedDetection.m       % experimental analysis (set the paths first)
```

Tested on MATLAB R2024a

## Notes
- Simulation results are stochastic (`rng('shuffle')`). For reproducible runs,
  replace it with `rng(0)`.
