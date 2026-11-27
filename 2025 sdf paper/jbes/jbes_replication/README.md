# Replication Archive: Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference

This anonymized archive accompanies the manuscript *"Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference"* and its Online Supplementary Material. Run every command from the root of the archive (`jbes_replication/`); the empty file `.here` marks that root for the `here` package.

Runtimes below were measured on the machine that produced the paper (macOS, 12 cores, R 4.3.1), except where marked.

---

## 1. Software

- **R 4.3.1** (2023-06-16) with base/recommended packages `parallel`, `stats`, `tools`, `utils`, `grid`.
- **CRAN packages used by the Monte Carlo study and the synthetic verification:** `here` 1.0.1, `yaml` 2.3.8, `data.table` 1.14.8, `ggplot2` 3.5.2.
- **Additionally used by the empirical (licensed-data) pipeline:** `dplyr` 1.1.4, `tidyr` 1.3.1, `scales` 1.4.0, `patchwork` 1.3.1, `readxl` 1.4.3, `xtable` 1.8.4, `optimx` 2022.4.30, `foreach` 1.5.2, `doParallel` 1.0.17. (`rstudioapi` is called only when a script is sourced interactively inside RStudio.)
- **Optional (Figure 1 and Supplement Tables C.1–C.2):** Python 3.9+ with `numpy`, `scipy`, `matplotlib` for `simulation/optimal_H_theory/verify_optimal_H_regimes.py`, and a LaTeX installation with `pgfplots` for `figure_optimal_H_theory.tex`.

---

## 2. Data

1. **Public, included.**
   - `empirical/data_prepared_2026/q_factors.csv`: monthly U.S. $q^5$ factors (`RF`, `MKT`, `ME`, `IA`, `ROE`, `EG`) from `http://global-q.org/factors.html` (Hou, Xue, and Zhang 2015; Hou, Mo, Xue, and Zhang 2021).
   - `empirical/data_prepared_2026_c/preqin_sample_funnel_North America.csv`: **aggregate** fund counts per strategy and vintage year only (no fund names, identifiers, sizes, or cash flows). The Monte Carlo designs use it to reproduce the vintage calendars.
2. **Proprietary, excluded.** Preqin fund-level cash flows (April 14, 2022 export; North American USD Buyout and Venture Capital funds, vintages through 2010) cannot be redistributed under the Preqin license. Licensed users place the raw export in `empirical/data_in/` and follow Section 4.3.
3. **Precomputed simulation outputs, included** (simulated data only, no fund-level data): `results/mc_step4/summary.csv`, `inference_checks.csv`, `v400_studentization_diagnostic.csv`, `v400_error_moments.csv`, `appendix_d_summary.csv`; `simulation/data_out_mc/summary.rds`, `d0/d0.rds`, `se_study/se_study.rds`; `simulation/optimal_H_theory/data/` (Figure 1 data).

---

## 3. Quick checks (no Preqin license; seconds)

| Command | What it does | Runtime |
| :--- | :--- | :---: |
| `LANG=en_US.UTF-8 Rscript run_synthetic_verification.R` | Builds synthetic Buyout and VC panels with the paper's vintage counts on the historical $q^5$ path. Checks the exact horizon rescaling identity (T1), the $U_{v,H}$ identity (the $O(n_v)$ formula against the direct sum over distinct fund pairs, EW and FW), and the fixed-$b$ critical values 2.790 / 2.867 / 2.722. Then runs the single-factor and $(\alpha,\beta)$ estimation, Bartlett fixed-$b$ inference, $N_{\mathrm{eff}}$, and the centered AR set end to end. Writes `output/synthetic_verification_summary.csv` (deterministic). | ~11 s |
| `LANG=en_US.UTF-8 Rscript simulation/mc_paper.R` | Rebuilds all Monte Carlo exhibits from the included simulation outputs: main Tables 2–3 and Figure 2; Supplement Tables C.3, D.2–D.8 and E.1; Figures D.1–D.2. Writes `results/mc_step4/paper/` and the Figure 1 Panel C data. | ~2 s |
| `LANG=en_US.UTF-8 Rscript simulation/mc_v400_table.R` | Supplement Table D.1 from `results/mc_step4/v400_studentization_diagnostic.csv`. | <2 s |

---

## 4. Full reproduction

### 4.1 Monte Carlo study (public data only)

Run in this order:

1. `simulation/mc_production.R`: all 72 designs of `config/mc_designs.yaml` (replication counts are in the file). Writes `simulation/data_out_mc/<design>/batch_*.rds`. **About 9.6 h wall-clock** on 12 cores.
2. `simulation/mc_targets.R`, then `simulation/mc_report.R`: computes the criterion targets and writes `simulation/data_out_mc/summary.rds` and `results/mc_step4/summary.csv`.
3. `simulation/mc_d0.R` and `simulation/mc_se_study.R` (~9 min): the Gaussian size benchmark (D0) and the SE study. `simulation/mc_inference_checks.R` (~2 min): the AR replay and the $V=400$ tail check.
4. `simulation/mc_v400_diagnostic.R` (~19 min; Supplement Table D.1) and `simulation/mc_v400_error_moments.R` (moments of $\hat\beta-1$ in the growing-$V$ designs; reads the stored batches).
5. `simulation/mc_paper.R` and `simulation/mc_v400_table.R`: the exhibits listed in Section 3.

**Figure 1 (Proposition 2).** Run `python3 simulation/optimal_H_theory/verify_optimal_H_regimes.py` (prints the lookup values of Supplement Tables C.1–C.2 and writes the Panel B data; not re-timed). Run `simulation/mc_paper.R` **after** it: the Python script also writes a Python-simulated Panel C file, and `mc_paper.R` overwrites it with the R Monte Carlo numbers cited in the figure caption and Table C.3. Then compile `simulation/optimal_H_theory/figure_optimal_H_theory.tex` with `pdflatex`.

### 4.2 Comparison to Driessen, Lin, and Phalippou (2012) (Supplement Table E.1)

The 19 rows of Table E.1 are the scenarios listed in `AD.BASE` at the top of `simulation/mc_appendix_d.R` (S1 = `base_case_zero_alpha`, ..., S10 = `big_n_v_50funds_alpha_stdv30_shifted_mkt8`, plus their `_additive` twins), defined in `config/simulation_scenarios.yaml`. Generate their simulated samples with `source("run_simulation_study.R"); run_simulation_study(scenario_ids = <those ids>, generate_data = TRUE)`, which calls `simulation/simulation_multifactor.R` and writes to `simulation/data_prepared_sim/` (about 7.3 GB). `simulation/mc_appendix_d.R` then re-estimates them with the paper's estimator (~38 min) and writes `results/mc_step4/appendix_d_summary.csv`, which is included. `simulation/mc_paper.R` builds Table E.1 from that file.

### 4.3 Empirical application (licensed Preqin data)

Run in this order:

1. `empirical/prep_preqin.R` (raw export in `empirical/data_in/`): writes the main sample `empirical/data_prepared_2026_c/` and the robustness samples `_allfunds`, `_fwmean`, `_fwdrop`. Then run `empirical/prep_nav_threshold.R` (default threshold 0.20), which writes `empirical/data_prepared_2026_c_navtrim20/`.
2. `run_empirical_study.R`: in R, `source("run_empirical_study.R"); run_empirical_study()` runs all 48 active scenarios of `config/empirical_scenarios.yaml`. From the shell, `Rscript run_empirical_study.R <scenario_id> ...` runs selected ones. Writes result caches to `results/data_out_2026-emp-c[-allfunds|-navtrim20]/`. **About 5 h wall-clock** on 12 cores.
3. `results/analyze_empirical_results.R` (~1 min): runs `results/analyze_result.R` on every cache folder and writes the inference and cross-validation summaries that steps 4–7 read. A folder without current results is reported and skipped.
4. `results/export_g1_g2_tables.R`: writes the L2 and $U_{v,H}$ data files and tables in `results/figures_bovc_us_g1/` (uses `results/plot_empirical_estimates.R`). It skips the $L_1$ tables, which are not exhibits of the paper, if the `*_l1` scenarios have no current results.
5. `results/empirical_iv.R` (IV roots, AR sets, Supplement Table F.2), `results/vc_alpha_diagnostic.R`, and `results/calendar_vs_index_shac.R` (~1 min; Supplement Remark D.1).
6. `results/main_exhibits_u.R`: main Tables 4–5 and Figure 3, plus Supplement Tables F.1, F.3 and F.4.
7. `results/step5_diagnostics.R`, `results/fw_imputation_sensitivity.R`, `results/robustness_samples.R`, `results/appendix_e_cv_folds.R`, and `results/table_published_estimates.R` (run after step 6): Supplement Tables F.5–F.12. The Appendix F figures come from steps 4 and 7.

Without a license, `run_synthetic_verification.R` exercises the same estimation and inference code on a synthetic panel (Section 3).

---

## 5. Directory structure

```text
jbes_replication/
├── README.md, .here
├── run_synthetic_verification.R        # quick end-to-end check (Section 3)
├── run_empirical_study.R               # empirical estimation driver (licensed data)
├── run_simulation_study.R              # simulated samples for Supplement Table E.1
├── estim_model_optimized.R             # estimation engine of the empirical and E.1 pipelines
├── config/                             # mc_designs.yaml (72 MC designs), empirical_scenarios.yaml, simulation_scenarios.yaml
├── helper/                             # estimators, admissible domain, Bartlett sandwich, fixed-b tables, caching
├── empirical/                          # public q5 factors, aggregate funnel counts, Preqin preparation scripts
├── simulation/                         # Monte Carlo generator, driver, analysis, exhibits; optimal_H_theory/ (Figure 1)
│   └── data_out_mc/                    # included summary.rds, d0/d0.rds, se_study/se_study.rds
├── results/                            # empirical exhibit scripts; mc_step4/ precomputed simulation CSVs
└── output/                             # synthetic verification output
```
