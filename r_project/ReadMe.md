# Cash-Flow-Based SDF Estimation for Private Equity Funds

This R project implements an improved cash-flow-based Stochastic Discount Factor (SDF) estimation method for private equity funds using public market factors.

## Overview

The project develops and validates estimation methods for the factor model:

```
R_PE = α + β × R_public + ε
```

Where:
- **R_PE**: Private equity returns
- **α (alpha)**: Fund-specific excess return (abnormal performance)
- **β (beta)**: Systematic risk exposure to public market factors
- **R_public**: Public market factor returns (e.g., MSCI World, Q-Factors)
- **ε (epsilon)**: Idiosyncratic residual return

### Research Papers

This codebase supports two research papers:

1. **Beta Estimation (Work in Progress)**: Focuses on improving the estimation of systematic risk exposure (β) using cash flow data and various SDF specifications.

2. **Residual/Epsilon Estimation (Published)**: [Measuring Private Equity Fund Idiosyncratic Risk](https://doi.org/10.1016/j.jfds.2024.100141) - Focuses on identifying and quantifying the idiosyncratic component (ε) of private equity returns.

---

## Project Structure

```
r_project/
├── estim_model_optimized.R        # Core SDF estimation algorithm
├── estim_model_empirical_runner.R # Runner script for empirical data estimation
├── estim_model_simulation_runner.R# Runner script for simulation studies
├── estim_residuals.R              # Residual/epsilon estimation methodology
│
├── helper/                        # Helper functions
│   ├── getNPVs.R                  # NPV calculation (with Rcpp optimization)
│   ├── hac_estimator.R            # HAC covariance matrix estimators
│   ├── load_data.R                # Data loading utilities
│   └── plot_npvs.R                # NPV visualization functions
│
├── empirical/                     # Empirical data processing
│   ├── data_in/                   # Raw input data
│   ├── data_prepared/             # Processed/cleaned data
│   ├── data_prepared_cutoff_*/    # Data with different cutoff dates
│   ├── data_private_public/       # Merged private-public datasets
│   ├── nav_returns/               # NAV return data
│   ├── prep_pitchbook.R           # PitchBook data preparation
│   ├── prep_preqin.R              # Preqin data preparation
│   └── prep_public.R              # Public market data preparation
│
├── simulation/                    # Monte Carlo simulation framework
│   ├── simulation_improved.R      # Main simulation script
│   ├── simulation_simple.R        # Simplified simulation
│   ├── data_prepared_sim/         # Simulated cash flow data
│   ├── data_results_sim/          # Simulation results
│   └── analyze_simulation_*.R     # Result analysis scripts
│
├── results/                       # Output results
│   ├── data_out/                  # Estimation output
│   ├── data_idi_published_2024/   # Published residual estimates
│   ├── analyze_result.R           # Result analysis utilities
│   └── combine_summary.R          # Summary aggregation
│
├── chart/                         # Output charts/figures
└── charts_error/                  # Error analysis visualizations
```

---

## Methodology

### SDF Estimation Approach

The estimation minimizes the squared Net Present Value (NPV) of fund cash flows when discounted using the SDF:

1. **Linear SDF Model**: `SDF_t = 1 + RF_t + β × MKT_t + ...`
2. **Exponential-Affine SDF Model**: `SDF_t = exp(RF_t + β × log(1 + MKT_t) + ...)`

### Key Features

- **Cross-Validation Support**: Optional k-fold cross-validation for hyperparameter tuning
- **Multiple Weighting Schemes**:
  - `EW` (Equal Weighting): All funds weighted equally
  - `FW` (Fund Weighting): Weighted by fund characteristics
  - `VYP` (Vintage Year Portfolios): Aggregation by vintage year
- **Regularization Options**:
  - L1 (Lasso) regularization
  - L2 (Ridge) regularization
- **HAC Standard Errors**: Heteroskedasticity and Autocorrelation Consistent inference using Bartlett/Parzen kernels
- **Rcpp Optimization**: Performance-critical NPV calculations implemented in C++

### Data Sources

- **Private Fund Data**:
  - PitchBook fund cash flows
  - Preqin fund cash flows
- **Public Market Factors**:
  - MSCI World market factors (MKT, RF, SMB, HML, etc.)
  - Hou-Xue-Zhang Q-Factors (MKT, ME, IA, ROE, EG)
  - iBoxx bond indices (for private debt analysis)

---

## Usage

### Prerequisites

```r
# Required R packages
install.packages(c(
  "data.table",    # Fast data manipulation
  "parallel",      # Parallel processing
  "doParallel",    # Parallel foreach
  "foreach",       # Loop constructs
  "Rcpp",          # C++ integration
  "fftwtools",     # FFT for HAC estimators
  "readxl",        # Excel file reading
  "xtable"         # LaTeX table export
))
```

### Preqin Data Preparation (2026 revision)

Run from `r_project/` (an empty `.here` file makes `here()` resolve to this folder):

```r
source("empirical/prep_preqin.R")   # needs the raw export in empirical/data_in/
```

**Input.** `empirical/data_in/Preqin_Cashflow_export-14_Apr_22272f54e5-...xlsx` (licensed, gitignored; data end 2022-03-31). Cash flows per USD 10mn of commitment are taken from the quarterly `Value` rows (cumulative net cash flow differenced; the final NAV is booked as a liquidating distribution on the last date). The export's dated Capital Call/Distribution rows reconcile exactly to these cumulatives; the quarterly booking dates lie about 50 days after the true flow dates on average.

**Vintage-year portfolios (constant composition).** `helper/pool_vintage.R` sums each date's flows over *all* funds of a vintage and divides by a time-invariant denominator: the fund count N_v (EW) or the total vintage commitment (FW). Funds without a transaction on a date contribute zero, so the pooled NPV equals the weighted average of fund NPVs for every discount curve. Both EW and FW pools are per unit of committed capital. Pools carry `N.funds` and `Commitment` columns (for vintage weights q_v = n_v). This replaces the earlier reporter-mean pooling (`aggregate(CF ~ Date, df, mean)`) used for `data_prepared_2026_b/`.

**Sample rule.** A fund has a *complete history* if its last record lies more than `min.history.years` (5) after its vintage year starts, or if it is still reporting at the data end (last record within `live.window.months` (12) of the export's last NAV date). The April 2022 export truncates many histories that the February 2020 export still has (vendor artifact); booking an early NAV as a liquidation is not defensible, so these funds are excluded from the main sample.

| Folder | Content |
|---|---|
| `empirical/data_prepared_2026_c/` | Main sample: complete-history funds. VYP files (tracked) and fund-level files (gitignored, regenerated) |
| `empirical/data_prepared_2026_c_allfunds/` | Robustness sample: all funds, truncated ones flagged `Complete = FALSE` |
| `empirical/data_prepared_2026_c/preqin_sample_funnel_<region>.csv` | Fund counts per strategy and vintage: raw, late start, missing NAV, truncated, complete |
| `empirical/data_prepared_2026_b/` | Previous reporter-mean pools (kept for before/after comparison) |

North America, vintages <= 2010: 383 Buyout and 432 Venture Capital complete-history funds (495 and 538 before the sample rule). Fund-level files carry `Size` (USD mn), `Last.Date`, `Age.Last` and `Complete`.

**Tests.** `LANG=en_US.UTF-8 Rscript helper/test_code_repairs.R` checks the pooling identity, the FW normalization, the sample rule, the funnel counts, the horizon convention and the U criterion on synthetic data and on all prepared files, and that every Preqin scenario resolves to a file. It fails if a prepared file is missing.

`LANG=en_US.UTF-8 Rscript helper/test_review_fixes.R` runs isolated integration regressions for cache invalidation, simulation analysis and partition lookup, CV provenance, and all three table export paths. It uses synthetic files in a temporary project; its numerical engine is replaced only for cache/runner tests. The numerical repair suite above exercises the valuation and inference calculations.

**Provenance and caching (B9).** Each result row carries `git.commit`, `git.dirty`, `config.md5`, `scenario.id`, `data.prepared.folder` and `datetime`. Generated LaTeX tables start with comment lines naming the commits, scenarios and data folders of their input results. Cache files are named `<time with ms>_cv<fold>_m<max.month>_l<lambda>_cached_res.csv`; result folders written before this change may be missing folds or horizons, because same-second files overwrote each other. One-parameter models use Brent's method on beta in [-1, 4]. `at.bound = TRUE` flags estimates on a bound. Cross-validation folds report `validation.error` (estimation horizon) and `validation.error.H0` (common H = 0 loss, comparable across horizons). The table column SR_CV is the full-sample estimate divided by the fold SD.

**Cache identity (r5 pipeline).** `helper/result_cache.R` records the resolved estimator settings and project-relative paths/hashes of cash flows, public factors, and any simulated-factor override (the YAML file path/hash is stored as provenance only, so editing comments or unrelated scenarios does not invalidate caches). Resume compares that full identity. Analysis and plot readers recheck the recorded files; changed or missing inputs are rejected, and multiple current estimator settings in one folder require an explicit identity or separate output folders. Criterion, origin and vintage-weight suffixes use the same folder-name builder (`result.cache.tag()`) for estimation, analysis-only runs and every simulation partition. The current pipeline version is `2026-11-17-r5` (admissible domain Assumption A7 `GROSS.FLOOR = 0.08`, constrained multi-start `fit.admissible()`, parameter scaling); older rows and summaries must be regenerated and are retained on disk rather than deleted.

**Reporting.** Asymptotic and CV summaries retain their input identities and contributing commits/scenarios. CV `at.bound` means that at least one fold hit a bound (`boundary.folds` gives the count). Horizon tables, vintage-cutoff tables and direct `empirical_tables.tex` exports mark full-sample boundary estimates with a dagger; invalid asymptotic SEs and t-ratios are suppressed. Fold SDs remain descriptive stability diagnostics.

**Standard errors and fixed-$b$ inference** (`helper/inference.R`, `helper/fixedb.R`). Sandwich covariance $H^{-1}\Psi H^{-1}/n$ from numerical derivatives ($\delta = 10^{-4}$) of the unit loss (squared horizon error for L2; $w_v U_{v,H}$ of the vintage for `U_pair`), with a Bartlett kernel in vintage-year distance (`kernel_bandwidth: "llsw"`, default $S = \lceil 1.3\sqrt{V}\rceil$ per Lazarus, Lewis, Stock, and Watson 2018; or a fixed numeric bandwidth) and Kiefer–Vogelsang fixed-$b$ $p$-values/critical values at $b = S/V$ (`SE.*.indep` gives the independent-unit reference). The Wald test of $\beta_{\mathrm{MKT}} = 1$ (other parameters $0$) uses $V^{-1}$. SEs are evaluated at the estimated horizon and are not valid when a parameter sits on an active boundary of $\Theta$.

**Speed.** The objective works on prepared units (`prepare.units()` and `prepare.fast()` in `helper/objective.R`): each unit's cash flows and factor columns are extracted once per data set. `do_parallel: true` (wired through `helper/config_loader.R` and `run_estimation()`) estimates strategies or simulation samples on forked workers when `n.parallel.workers(length(types)) > 1L` (falling back to serial `%do%` automatically when only 1 worker is needed or allowed). Results are bit-identical to the serial path (equivalence tests in `helper/test_code_repairs.R`).

**Memory.** The number of workers (estimation and simulation) is capped at the smallest of `n.tasks`, `max(1, detectCores() - 2)` (falling back to `1L` if `detectCores()` is unknown/non-finite), `10`, `SDF_MAX_WORKERS`, and free-memory capacity (`helper/parallel_workers.R`, keeping 6 GB reserved at 1 GB per worker by default; override via `SDF_MAX_WORKERS` or `options(sdf.gb.per.worker = ..., sdf.reserve.gb = ...)`).

**Horizon convention and U criterion** (`helper/objective.R`). Every unit of a vintage is placed on one monthly grid from the vintage origin $t_0$ to the common data end $T_{\max}$. `max.month = m` averages the pricing error over the $m$ dates $t_0, \dots, t_0 + m - 1$; `m = 1` (and its alias `m = 0`) is the inception NPV at $t_0$, with $t_0$ not discounted. `origin_convention` sets $t_0$: `"first.cf"` (first cash-flow date of the vintage, default) or `"predetermined"` (Dec 31 of vintage year $- 1$). `error_function: "U_pair"` minimizes the cross-fund pair-product criterion on fund-level files (`use_vintage_year_pfs: false`), with vintage weights `vintage_weights: "equal"` ($q_v = 1$) or `"n.funds"` ($q_v = n_v$); see the `u_*_north_america` scenarios.

**Licensed-data boundary and rebuild order (2026 revision).**
Run all scripts from `r_project/` with `LANG=en_US.UTF-8`:
1. **Prepare licensed Preqin data** (requires raw Excel export in `empirical/data_in/`, gitignored):
   - `LANG=en_US.UTF-8 Rscript empirical/prep_preqin.R` (writes `data_prepared_2026_c/` and `data_prepared_2026_c_allfunds/`, plus checked-in vintage pools and `preqin_sample_funnel_<region>.csv`)
   - `LANG=en_US.UTF-8 Rscript empirical/prep_nav_threshold.R` (writes `data_prepared_2026_c_navtrim20/`)
2. **Empirical estimation (YAML runner):**
   - `LANG=en_US.UTF-8 Rscript run_empirical_study.R` (reads `config/empirical_scenarios.yaml`, writes to `results/data_out_2026-emp-c*`; fund-level `U_pair` scenarios require the gitignored fund-level CSVs from Step 1)
   - Empirical tables and figures: `results/main_exhibits_u.R`, `results/empirical_iv.R`, `results/robustness_samples.R`, `results/plot_empirical_estimates.R`
3. **Step 4 Monte Carlo (no licensed fund-level files required — reads counts from checked-in `preqin_sample_funnel_North America.csv`):**
   - `LANG=en_US.UTF-8 Rscript simulation/mc_production.R [design ...]` (reads `config/mc_designs.yaml`, writes checkpointed batches to `simulation/data_out_mc/`)
   - `LANG=en_US.UTF-8 Rscript simulation/mc_targets.R [--out=targets_NAME.rds design ...]`
   - `LANG=en_US.UTF-8 Rscript simulation/mc_d0.R` and `LANG=en_US.UTF-8 Rscript simulation/mc_appendix_d.R`
   - `LANG=en_US.UTF-8 Rscript simulation/mc_report.R` (validates full provenance/completeness against `config/mc_designs.yaml` and writes `results/mc_step4/summary.csv`)
   - `LANG=en_US.UTF-8 Rscript simulation/mc_paper.R` (generates Step 4 LaTeX tables and figures in `results/mc_step4/`)
4. **Tests:**
   - `LANG=en_US.UTF-8 Rscript helper/test_review_fixes.R` (synthetic integration regressions)
   - `LANG=en_US.UTF-8 Rscript helper/test_mc.R` (Step 4 Monte Carlo & Appendix D acceptance suite)
   - `LANG=en_US.UTF-8 Rscript helper/test_code_repairs.R` (real-data pooling, sample-rule, horizon, U-criterion, and inference repairs)

### Deprecated Legacy Scripts

The following older scripts are retained for historical reference only and are **deprecated** for the 2026 revision:
- `estim_model_empirical_runner.R` — replaced by `run_empirical_study.R` + `config/empirical_scenarios.yaml`
- `estim_model_simulation_runner.R` — replaced by `run_simulation_study.R` + `config/simulation_scenarios.yaml` (for G1/Appendix D) and `simulation/mc_production.R` + `config/mc_designs.yaml` (for Step 4 MC)
- `empirical/prep_preqin.Rmd` — replaced by `empirical/prep_preqin.R`
- `simulation/simulation_improved.R` — replaced by `simulation/simulation_multifactor.R` (G1/B1) and `simulation/mc_generate.R` (Step 4 MC)

### Running Legacy Empirical / Simulation Runners (Deprecated)

1. **Prepare Data**: Ensure cash flow data is in `empirical/data_prepared/`

2. **Configure Parameters** in `estim_model_empirical_runner.R`:
   ```r
   private.source <- "pitchbook"  # or "preqin"
   weighting <- "FW"              # "EW" or "FW"
   use.vintage.year.pfs <- TRUE   # Use vintage year portfolios
   do.cross.validation <- FALSE   # Enable CV for hyperparameter tuning
   sdf.model <- "linear"          # or "exp.aff"
   include.alpha.term <- FALSE    # Include alpha in estimation
   ```

3. **Execute**:
   ```r
   source("estim_model_empirical_runner.R")
   ```

### Running Simulation Studies

> **2026 revision (B1).** `run_simulation_study.R` generates data with `simulation/simulation_multifactor.R` (`create.simulation.multifactor`). Deal payoffs compound over the months after the investment month, matching the estimator's discounting, and the default `shock_model: "multiplicative"` uses mean-one lognormal shocks G_u (1 + eps_u), so the DGP prices exactly at the true parameters. `shock_model: "additive"` keeps the former additive shocks with default truncation as a stress case. Data and results generated before this change must be regenerated. The older `simulation_improved.R` described below is legacy and not used by the runner.

1. **Generate Simulated Data** (`simulation/simulation_improved.R`, legacy):
   ```r
   create.simulation(
     no.deals = 15,              # Deals per fund
     investment.period = 5,      # Investment period (years)
     max.holding.period = 10,    # Maximum holding period (years)
     alpha = 0,                  # True alpha (monthly)
     beta = 1,                   # True beta
     no.samples = 1000,          # Monte Carlo samples
     no.funds = 20,              # Funds per vintage
     min.vin = 1986,             # First vintage year
     max.vin = 2005,             # Last vintage year
     stdvs = 0.2,                # Idiosyncratic volatility
     exp.aff.sdf = FALSE         # Use exponential-affine SDF
   )
   ```

2. **Run Estimation on Simulated Data**:
   ```r
   source("estim_model_simulation_runner.R")
   ```

### Estimating Residuals

For idiosyncratic return estimation (epsilon), use:

```r
source("estim_residuals.R")
```

This implements the methodology from the published paper, using Component-wise Boosting (CWB) to estimate time-varying residual returns.

---

## Key Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `private.source` | Data source | `"pitchbook"`, `"preqin"` |
| `weighting` | Cash flow weighting | `"EW"`, `"FW"` |
| `use.vintage.year.pfs` | Aggregate to vintage year portfolios | `TRUE`, `FALSE` |
| `do.cross.validation` | Enable cross-validation | `TRUE`, `FALSE` |
| `sdf.model` | SDF specification | `"linear"`, `"exp.aff"` |
| `max.months` | Cash flow horizon (months) | `120`, `150`, `180` |
| `lambda` | Regularization strength | `0` (no regularization) |
| `kernel.bandwidth` | HAC bandwidth | `"llsw"` ($S = \lceil 1.3\sqrt{V}\rceil$ vintage years, default) or numeric vintage years |

---

## Output Files

### Estimation Results

Results are saved to `results/data_out*/` as CSV files containing:
- Estimated factor loadings (betas)
- Alpha estimates (if enabled)
- Standard errors (HAC-robust)
- Model diagnostics (R², objective function values)

### Residual Estimates

Time-series of idiosyncratic returns are saved to `results/data_idi*/`:
- Date-indexed residual returns
- Decomposition of total vs. factor-driven returns

---

## Technical Notes

### Performance Optimization

- **Rcpp Integration**: The `getNPVs()` function is implemented in C++ for ~10x speedup
- **Parallel Processing**: Multi-core support via `parallel`/`doParallel` (Unix-based systems)
- **Fast HAC Estimation**: Uses FFT-based algorithm from Heberle & Sattarhoff (2017)

### Numerical Considerations

- Cash flows are discounted using cumulative SDF products
- Monthly data frequency assumed throughout
- Missing cash flow months are filled with zeros for continuous time series

---

## Citation

If using the residual estimation methodology, please cite:

```bibtex
@article{tausch2024measuring,
  title={Measuring Private Equity Fund Idiosyncratic Risk},
  author={Tausch, Christoph},
  journal={Journal of Financial Data Science},
  year={2024},
  doi={10.1016/j.jfds.2024.100141}
}
```

---

## License

*[Add license information as appropriate]*

---

## Contact

*[Add contact information or links to related resources]*
