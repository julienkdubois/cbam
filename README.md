CBAM -- Dynare models repository

Overview

This repository contains Dynare model files and helper scripts for a two-country DSGE model with a Carbon Border Adjustment Mechanism (CBAM).

Files and purpose

- `dynare/call_file/ss_residuals.m`: helper used to adjust steady-state residuals to match GDP targets.
- `dynare/cbam_estimation.mod` and `dynare/ss_pf.m`: files used for Bayesian estimation and steady-state computation.
- `dynare/pf_cbam.mod`: perfect-foresight model used to simulate the CBAM transition; relies on `ss_pf.m` for steady-state solving.
- `dynare/compare_pf.m`: script to compare trajectories (sticky vs flexible) and run simulations for different implementation dates. Features: default horizon = 40, optional figure saving, multi-date comparison.
- `dynare/old/`: historical and experimental versions. These are excluded from recent commits.

Quick start

- Compare three implementation dates:
  ```matlab
  compare_pf([1,13,30])
  ```
- To use estimated parameters when running `compare_pf`, answer `y` to the interactive prompt; the script will patch `pf_cbam.mod` with estimated values if available.

Notes and recent changes

- On 30/04/2026 the handling of `theta1` was harmonized between the estimation code and the perfect-foresight steady-state (fix applied in `ss_pf.m`). `compare_pf.m` was updated for multi-date comparisons and figure saving.

Data section (placeholder)

- The README will include a dedicated data description and file list. Tell me which data files to document and I will add that section.
