# Repo-specific Copilot instructions

This repository is primarily a MATLAB research codebase (signal processing/TDOA/up-sampling) with a small Python helper. The guidance below helps an AI coding assistant be immediately productive.

- Purpose: implement and analyse TDOA/upsampling pipelines for lightning/echo signals (scripts under `2023/`, `2019/`, and `evaluate/`). See example pipeline: [2023/find_peak_upsample_corr.m](2023/find_peak_upsample_corr.m).
- Primary languages: MATLAB (.m). Minimal Python tools and notebooks exist (e.g. `x_y_z_to_angle.py`, `other/Analysis Signal Interval.ipynb`).

- Project layout quick map:
  - `2023/`, `2019/` — experiment scripts and data processing entry scripts (run these in MATLAB). Example: [2023/find_peak_upsample_corr.m](2023/find_peak_upsample_corr.m).
  - `evaluate/` — figure generation and evaluation utilities (paper-format plotting). Example: [evaluate/evalute_upsample.m](evaluate/evalute_upsample.m).
  - `other/` — utility helpers: TDMS readers, filters, plotting helpers. Example: [other/filter_bp.m](other/filter_bp.m), [other/TDMS_readTDMSFile.m](other/TDMS_readTDMSFile.m).
  - top-level small tools: `x_y_z_to_angle.py` (coordinate helper).

- How scripts are structured (important patterns):
  - Most MATLAB files are scripts (executable top-level code) that include helper function definitions at the bottom of the same `.m` file. Edit within the same file when modifying algorithm flow. See [2023/find_peak_upsample_corr.m](2023/find_peak_upsample_corr.m) for the pattern.
  - Many scripts use hard-coded numeric constants and relative Windows-style paths (e.g. `'..\\20230718175104.9180CH1.dat'`). Before running, ensure MATLAB's current folder is set so relative paths resolve, or update paths in the script.
  - Processing flow is typically: read signal -> bandpass (see `filter_bp`) -> smoothing/windowing -> upsample (`interp1`, spline) -> cross-correlation/gccphat (`xcorr`) -> TDOA estimation -> lsq optimization (`lsqnonlin`) -> write result text files.

- Required MATLAB toolboxes/assumptions:
  - Signal Processing functions (`butter`, `filtfilt`, `xcorr`, `designfilt`) are used.
  - Optimization (`lsqnonlin`) is used in some scripts — ensure Optimization Toolbox is available.
  - Scripts assume large sampling rates (Fs ~ 200e6) and may expect enough memory for long arrays.

- Typical developer workflow (to run an analysis):
  1. Open MATLAB and set Current Folder to the repository root (`f:\\LHW`) or the specific experiment folder (e.g., `2023`).
  2. Inspect and, if necessary, fix relative paths in the script.
  3. Run the main script in the editor/command window: e.g. `run find_peak_upsample_corr` (or press Run in the editor while the script is active).
  4. Outputs: scripts usually write text results (e.g. `result_4-6.txt`) and figures (in `evaluate/` scripts). Check the script header/comments for plot/save behavior.

- Notes for editing and testing:
  - Keep helper functions inside the same `.m` file unless factoring out is deliberate — many scripts rely on local subfunctions.
  - Prefer small, local edits and run on a single experiment (files under `2023/`) to verify behavior before broad refactors.
  - Use MATLAB's profiler and `plot`/`disp` for quick verification; scripts are not covered by unit tests.

- Useful files to inspect when changing core logic:
  - `[2023/find_peak_upsample_corr.m](2023/find_peak_upsample_corr.m)` — canonical processing pipeline.
  - `[other/filter_bp.m](other/filter_bp.m)` — bandpass filter utility used across scripts.
  - `[evaluate/evalute_upsample.m](evaluate/evalute_upsample.m)` — example figure-generation and plotting conventions.

If any part of this is unclear or you want more details (toolbox versions, frequently used constants, or to add runnable helper scripts), tell me which area to expand and I will iterate.
