CBAM -- Dynare models repository

Overview

This repository contains Dynare model files and helper scripts for a two-country NK CBAM model.

Files and purpose

- `dynare/call_file.m and ss_residuals.m`: used to adjust steady-state residuals to match GDP targets.
- `dynare/cbam_estimation.mod` and `dynare/ss_pf.m`: files used for Bayesian estimation and steady-state computation.
- `dynare/pf_cbam.mod`: perfect-foresight model used to simulate the CBAM transition; relies on `ss_pf.m` for steady-state solving.
- `dynare/compare_pf.m`: script to compare trajectories (sticky vs flexible) and run simulations for different implementation dates. Features: default horizon = 40, optional figure saving, multi-date comparison.
- `dynare/old/`: historical and experimental versions. These are excluded from recent commits.

Quick start

- Compare three implementation dates of CBAM:
  ```matlab
  compare_pf([1,13,30])
  ```
  Compare flexible and sticky and one date: 
  ```matlab
  compare_pf
  ```
- To use estimated parameters when running `compare_pf`, answer `y` to the interactive prompt; the script will patch `pf_cbam.mod` with estimated values if available.

Data

The `data/` folder contains raw and cleaned datasets and the scripts used to prepare them.

- `data/clean/`: cleaned, ready-to-use CSVs and small derived tables used by the Dynare scripts. Examples include `dynare_data_eu_row.csv`, `panel_all_indicators.csv` and several OECD series (GDP, CPI, consumption, trade, etc.). These files are safe to commit and share.
- `data/raw/`: original raw sources (large, not committed here).
- `data/cleaning.R` and `data/prepare.R`: R scripts that implement the cleaning pipeline. Run these locally to regenerate the files in `data/clean/` from `data/raw/`.

Commit policy for data

- Only cleaned, compact datasets live in version control (`data/clean/`). Raw sources are intentionally excluded because they are large.
- When you want the README to include metadata for each data file (source, frequency, transformations), I can auto-generate a table from the `data/clean/` folder or you can provide the descriptions and I'll add them.

