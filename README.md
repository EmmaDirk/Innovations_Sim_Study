# Simulation Study

This repo contains **7 scripts** that together form a simulation engine:

- `00_packages.R`: Loads all dependencies (except the `here` package).

- `01_population_data_generation.R`: Generates population data according to the following data-generating mechanism (DGM):

  - **Z** ~ Normal(0, 1)
  - **U** ~ Normal(0, 1)
  - **Z** independent of **U**
  - **X** = b2 * Z + b5 * U
  - Centered interaction: **XZc** = (X * Z) - b2
  - **Y** = b1 * X + b6 * XZc + b3 * Z + b4 * U

  Residual variances are chosen such that each variable has variance 1.

- `02_sampling.R`: Contains two sampling functions:
  - **SRS**: draws a simple random sample from the population and drops **U**.
  - **NPS**: draws a non-probability sample from the population (keeps all variables) using:

    - Selection index (unit i):
      - eta_i = aX * X_i + aY * Y_i + aZ * Z_i + aU * U_i

    - Weights:
      - w_i = exp(eta_i)

    - Sampling rule: draw **exactly n units without replacement** with probabilities proportional to w_i:
      - P(i selected) proportional to w_i

    Units with larger eta_i (depending on aX, aY, aZ, aU) are more likely to be selected, creating selection bias when any a != 0.

- `03_analysis.R`: Contains three analysis functions (all return the estimated effect of **X** on **Y** plus an SE).
  - **SRS analysis**: OLS regression `Y ~ X + Z` (note: **U** is not available in the SRS by design).
  - **NPS analysis**: OLS regression `Y ~ X + Z + U`.
  - **Combined (two-step IPW)**:
    1) Fit a sample-membership model (logistic regression) on stacked SRS+NPS data: `source_nps ~ X + Z + Y`.
    2) Compute inverse-odds weights for the NPS: `w = (1 - p_hat) / p_hat` (optionally stabilized + trimmed).
    3) Fit a weighted outcome model on the NPS: `Y ~ X + Z + U` with weights `w`.
    SEs are computed via bootstrap over the two-step procedure.

- `04_simulator.R`: Runs the Monte Carlo simulation over a grid of **aU** values.
  - For each replication and each aU, it draws:
    - an **SRS** sample,
    - an **NPS** sample (optionally two sizes: “small” and “big data”),
    - and fits the corresponding estimators (SRS, NPS, Combined IPW).
  - Stores replication-level results (`beta_hat`, `se_hat`, `bias = beta_hat - b1`) and returns:
    - `draws`: all replications
    - `summary`: averages by `(aU, method)` for plotting
  - Supports parallel execution + a progress bar.

- `05_master.R`: Orchestrates the full pipeline.
  - Generates 4 population datasets (N = 1,000,000) with different **U-confounding strengths** by varying `b4` and `b5` (0.05, 0.10, 0.20, 0.40).
  - Runs `simulation_wrapper()` on each population across `aU_vals = seq(0, 2, by = 0.1)` with:
    - `n_srs = 1000`
    - `n_nps = c(200, 5000)` (small NPS + big-data NPS)
  - Saves populations to `02_data/` and simulation outputs to `03_output/`.

- `06_checks.R`: Quick sanity checks / smoke tests.
  - Sources each script and runs small examples to confirm:
    - packages load,
    - population moments look reasonable,
    - sampling functions work,
    - expected biases show up in SRS/NPS,
    - the two-step IPW estimator runs and reduces bias in a toy setup.

- `07_plots.R`: Produces the main figure from the saved simulation outputs.
  - Loads the 4 saved `out_U_*_int.RData` objects.
  - Builds:
    - **Bias vs aU** curves from `out$summary`
    - **RMSE vs aU** curves computed from replication draws (`rmse = sqrt(mean(bias^2))`)
  - Combines into a single **2x4** plot (top row: Bias, bottom row: RMSE) and saves:
    - `03_output/big_2x4_bias_rmse_int.png`