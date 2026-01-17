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

### Running Empirical Analysis

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

1. **Generate Simulated Data** (`simulation/simulation_improved.R`):
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
| `kernel.bandwidth` | HAC bandwidth | `12` (months) |

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
