# Replication Archive: Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference

This self-contained, anonymized replication archive accompanies the manuscript *"Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference"* and its Online Supplementary Material.

---

## 1. Software & Environment Requirements

All results in the paper were produced on macOS using:
- **R 4.3.1** (`2023-06-16`)
- **Required CRAN packages:**
  - `data.table` (`1.14.8`)
  - `yaml` (`2.3.8`)
  - `ggplot2` (`3.5.2`)
  - `here` (`1.0.1`)
  - Base/recommended packages: `parallel`, `stats`, `utils`
- **Optional (Proposition 2 analytical lookup & figure verification):** Python 3.9+ with `numpy`, `scipy`, `matplotlib` (`simulation/optimal_H_theory/verify_optimal_H_regimes.py`).

Run all R commands from the root directory of this archive (`jbes_replication/`).

---

## 2. Data Availability & Proprietary Preqin Exemption

1. **Public U.S. $q^5$ Factor Data (included):**
   - `empirical/data_prepared_2026/q_factors.csv` contains monthly U.S. $q^5$ factors (`RF`, `MKT`, `ME`, `IA`, `ROE`, `EG`) downloaded from `http://global-q.org/factors.html` (Hou, Xue, and Zhang 2015; Hou, Mo, Xue, and Zhang 2021).
   - `empirical/data_prepared_2026_c/preqin_sample_funnel_North America.csv` contains **aggregate-only** fund counts per strategy and vintage year (no fund names, identifiers, sizes, or cash flows) used solely to define vintage calendar shapes in the Monte Carlo designs.
2. **Proprietary Preqin Fund-Level Cash Flows (excluded per license):**
   - Commercial Preqin fund-level cash-flow histories (April 14, 2022 export; North American USD Buyout and Venture Capital funds, vintages through 2010) cannot be redistributed under Preqin's licensing agreement.
   - Licensed Preqin users can reproduce the exact empirical cleaning and estimation pipeline by placing the raw Preqin export in `empirical/data_in/` and running `empirical/prep_preqin.R`, `run_empirical_study.R`, `results/empirical_iv.R`, `results/main_exhibits_u.R`, `results/calendar_vs_index_shac.R`, and `results/step5_diagnostics.R`.
3. **End-to-End Verification Without a Preqin License (`run_synthetic_verification.R`):**
   - Run `LANG=en_US.UTF-8 Rscript run_synthetic_verification.R` (~15 seconds) to generate a synthetic Preqin-calibrated Buyout & Venture Capital cash-flow panel (matching exact vintage counts, deal lifetimes, lognormal fund sizes, and calibrated deal volatilities on the historical $q^5$ factor path), verify core algebraic identities (`fast.unit.errors`, `u.pair`, `fixedb.crit.t`), and execute the entire single-factor, joint $(\alpha,\beta)$, calendar-vs.-index Bartlett fixed-$b$ inference, $N_{\mathrm{eff}}$, and centered Anderson--Rubin confidence-set pipeline end-to-end.

---

## 3. Directory Structure & Script-to-Exhibit Mapping

```text
jbes_replication/
├── README.md                                  # This file
├── run_synthetic_verification.R               # Quick (~15s) end-to-end synthetic verification
├── config/mc_designs.yaml                     # Complete 72-design Monte Carlo specification
├── helper/                                    # Core estimator, domain, & fixed-b SHAC modules
│   ├── objective.R                            # Horizon-averaged error, L2, U_{v,H}, U^{cross}
│   ├── admissible.R                           # Admissible domain Theta & constrained optimizer
│   ├── inference.R                            # Bartlett SHAC sandwich covariance & N_eff
│   ├── fixedb.R (+ .rds tables)               # LLSW bandwidth & Kiefer-Vogelsang fixed-b p-values
│   ├── getNPVs.R, pool_vintage.R              # Cash-flow valuation & constant-composition pooling
│   └── estim_model_core.R, result_cache.R     # Core pipeline & checkpointing
├── empirical/                                 # Public q^5 factors + licensed-Preqin prep scripts
│   ├── data_prepared_2026/q_factors.csv       # Public U.S. q^5 factors
│   ├── data_prepared_2026_c/preqin_sample_funnel_North America.csv  # Aggregate counts only
│   └── prep_preqin.R, prep_public.R, prep_nav_threshold.R
├── simulation/                                # Full Monte Carlo & Proposition 2 theory suite
│   ├── mc_generate.R, mc_driver.R, mc_analyze.R, mc_targets.R, mc_production.R
│   ├── mc_inference_checks.R                  # AR moment replay & V=400 tail check
│   ├── mc_v400_diagnostic.R, mc_v400_table.R  # V=25..1600 studentization & tail ablation
│   ├── mc_paper.R, mc_se_study.R, mc_d0.R     # Exhibit generators & D0 / SE studies
│   └── optimal_H_theory/                      # Prop. 2 lookup tables & Figure 1 generator
└── results/                                   # Empirical analysis scripts & precomputed MC CSVs
    ├── mc_step4/summary.csv, inference_checks.csv, v400_studentization_diagnostic.csv
    ├── empirical_iv.R, main_exhibits_u.R, calendar_vs_index_shac.R
    └── step5_diagnostics.R, fw_imputation_sensitivity.R, robustness_samples.R, ...
```

### Exhibit Mapping & Approximate Runtimes (8-core workstation)

| Manuscript / Supplement Exhibit | Script | Data Needed | Approx. Runtime |
| :--- | :--- | :--- | :---: |
| **Main Fig. 1** & **Supp. Tables C.1–C.2** (Prop. 2 duration matching) | `simulation/optimal_H_theory/verify_optimal_H_regimes.py` | Public `q_factors.csv` | ~10 sec |
| **Main Tab. 2–3**, **Main Fig. 2**, **Supp. Tables C.3, D.2–D.8, E.1**, **Supp. Fig. D.1–D.2** | `simulation/mc_production.R` $\to$ `simulation/mc_paper.R` | Public `q_factors.csv` + funnel counts | ~45 min (full 72 designs); ~5 sec from included `results/mc_step4/summary.csv` |
| **Supp. Table D.1** ($V=25\text{--}1{,}600$ studentization ablation) | `simulation/mc_v400_diagnostic.R` $\to$ `simulation/mc_v400_table.R` | Public `q_factors.csv` | ~15 min (500 reps); <2 sec from included CSV |
| **AR replay & $V=400$ tail check** | `simulation/mc_inference_checks.R` | Public `q_factors.csv` + MC batches | ~2 min |
| **Main Fig. 3**, **Main Tab. 4–5**, **Supp. Tables F.1–F.10**, **Supp. Fig. F.1–F.6** | `run_empirical_study.R` $\to$ `results/empirical_iv.R` $\to$ `results/main_exhibits_u.R` $\to$ `results/step5_diagnostics.R` | Licensed Preqin export (or test pipeline via `run_synthetic_verification.R`) | ~15 min (full empirical); ~15 sec (synthetic verification) |
| **Supp. Remark D.1** (calendar-lag vs. index-lag Bartlett) | `results/calendar_vs_index_shac.R` | Licensed Preqin export (or synthetic verification) | ~1 min |
