# Revision Fixes & Open Deliverables Checklist (Reconciled Audit, Nov 2026)

Master action list reconciling the November 20 audit (`4795963..e53326c`, `revision-2026`), the follow-up double-check audit, and the November 22 cross-audit (`e3501c6`). Ordered by dependency so no estimation batch has to be run twice.

---

## Phase 1 — Settle FW Weighting & Harden Rebuild Infrastructure (Before Reruns)

### 1.1 Fund-Size Imputation & Data Preparation
- [ ] **Decide and implement the FW fund-size imputation rule (`r_project/empirical/prep_preqin.R`, lines 181–212)**
  - Current rule imputes missing fund sizes with the pooled cross-vintage mean (`$1,842.24M` Buyout over 1985–2021; `$379.43M` VC over 1982–2021).
  - In the main pre-2011 sample, **10 of 383 Buyout** and **15 of 432 VC** funds have missing sizes, concentrated in the earliest, smallest vintages:
    - Buyout single-fund shares: **82.1%** in 1985 (`1/2`), **94.9%** in 1986 (`1/2`), **77.8%** in 1990 (`1/2`); two-fund combined shares: **51.5%** in 1988 (`2/6`), **92.8%** in 1992 (`2/4`).
    - VC two-fund combined shares: **98.8%** in 1985 (`2/3`), **87.3%** in 1986 (`2/4`); single-fund shares: **34.0%** in 1992 (`1/7`), **39.1%** in 1994 (`1/5`). (Singleton 1982 VC has 100% weight at any positive imputed size under L2.)
  - Replacing pooled means with **within-vintage medians** (with nearby-vintage fallback when needed) roughly halves the Buyout FW L2 0→30 yr decline (`0.634 → 0.336`; 30 yr `0.558 → 0.889`; inception $\text{L2}-U$ gap `-0.013 → +0.043`) while holding the fund sample fixed.
  - **Dropping missing-size funds** moves Buyout FW L2 to `1.381` ($H=0$) and `0.962` ($H=30$) (`L2-U = +0.155`), but drops Buyout $U_{v,H}$ from 25 to **22** vintages while L2 stays at 25.
  - **Deliverable:** Settle the baseline FW rule, rebuild `empirical/data_prepared_2026_c*` (main, `allfunds`, `navtrim20`), re-run `helper/test_code_repairs.R`, and report **all three rules (pooled mean, within-vintage median, and deletion)** as sensitivities with explicit vintage counts and matched-support comparisons.
  - **Downstream check:** If the baseline FW weighting rule changes, revisit the FW Monte Carlo calibration (`sdlog` in `s_core_BO` / `s_core_VC`) and predicted $\text{L2}-U$ gap.

### 1.2 Runner, Parallelism, Checkpoint, and Analysis Code Fixes
- [ ] **Fix unknown CPU count & parallel cluster handling**
  - Files: `r_project/helper/parallel_workers.R` (lines 33–48), `r_project/helper/config_loader.R`, `r_project/estim_model_optimized.R` (line 747).
  - Fix: Fall back to `1L` worker when `parallel::detectCores()` returns `NA`/non-finite before the `by.memory` comparison (add unit test for unknown core count); wire `do_parallel` through the YAML adapters in `config_loader.R`; skip `makeForkCluster()` when `n.parallel.workers(...) <= 1L`.
- [ ] **Harden full MC checkpoint, report, and Appendix D validation**
  - Files: `r_project/simulation/mc_driver.R` (lines 190–256), `r_project/simulation/mc_production.R` (lines 20–25), `r_project/simulation/mc_report.R` (lines 14–23), `r_project/simulation/mc_appendix_d.R` (lines 132–149).
  - Fix: Validate full checkpoint provenance (`spec.md5`, `public.md5`, `settings`, and exact `reps` IDs) on batch reuse; fail if requested reps are missing or duplicated; validate batches and targets in `mc_report.R` against current YAML specs (with an explicit partial-run mode if needed); validate `mc_appendix_d.R` caches beyond filename presence; exit non-zero in `mc_production.R` on any failed/incomplete design.
- [ ] **Use checked-in sample funnel counts in `mc.load.designs()`**
  - File: `r_project/simulation/mc_driver.R` (lines 258–283).
  - Fix: Read BO/VC vintage counts and calendars from the checked-in `preqin_sample_funnel_North America.csv` (verified to reproduce both production and D0 counts/calendars exactly) instead of requiring the git-ignored fund-level CSV.
- [ ] **Fix `mc_analyze.R` summary metrics & add multi-start safeguard in `mc.targets()`**
  - Files: `r_project/simulation/mc_analyze.R` (lines 85–145), `r_project/simulation/mc_driver.R` (`mc.fit()`, line 166).
  - Fix: Store `neff.Coef` in `mc.fit()` and report coefficient-specific `neff`; keep all-draw `sd(est)` **and** report a clearly labelled same-sample (`tested`-only) `se.ratio` separately (note: for 30 yr core $U$ joint $\beta$, switching denominators moves the ratio from `0.915/0.933` to `1.342/1.330`); upgrade `mc.targets()` line 98 to `admissible.starts()` as a safeguard (verified not to change saved targets—max relative objective gain $7\times 10^{-13}$).
- [ ] **Fix legacy runner, CV scope, IV LOO domain, and legacy alpha conversion bugs**
  - `r_project/estim_model_optimized.R`: guard `config.path` when sourced with `source.internally = FALSE` (lines 102–125, 275–282); decide CV fold scope explicitly (line 639: currently `estimate <- 1980:2020`, latest validation block ends 2017; does not affect full-sample G1 fits) and invalidate affected CV caches if changed.
  - `r_project/results/empirical_iv.R` (lines 79–102): hold full-sample admissible domain fixed across LOO refits; stop silently folding zero-root `min.square` fallbacks into jackknife SEs.
  - `r_project/results/analyze_result.R` (lines 511–565): fix legacy monthly-vs-quarterly alpha conversions (`^4` → `^12`, remove `^(1/3)` in legacy diagnostics outside E2).
  - `r_project/run_simulation_study.R` (line 163): fix `next` inside error handler.
  - `r_project/ReadMe.md`: point readers to current YAML runners, `r5` cache identity, fixed-$b$ inference, rebuild order, and licensed-data dependencies; mark legacy `estim_model_empirical_runner.R`, `prep_preqin.Rmd`, and `simulation_improved.R` deprecated.

---

## Phase 2 — Empirical & Monte Carlo Reruns and Exhibit Updates

### 2.1 Empirical Reruns & Comparison Support
- [ ] **Rerun FW empirical caches on finalized FW pools**
  - Re-estimate FW L2 and FW $U_{v,H}$ for main (`data_out_2026-emp-c`), `allfunds`, `navtrim20`, `max_vintage`, and `nav_disc`; regenerate all affected tables/figures.
- [ ] **Add matched-vintage support for both L2 and IV alongside full-sample estimates**
  - Files: `r_project/results/empirical_iv.R`, `r_project/results/main_exhibits_u.R`.
  - Derive eligible ($n_v \ge 2$) vintages dynamically per sample rather than hardcoding 23.
  - Keep full-sample 27-vintage VC L2/IV estimates and report matched-vintage L2 **and** IV wherever $\text{L2}-U_{v,H}$ or $U_{v,H}-\text{IV}$ gaps are interpreted (on current pools, matched 23-vintage VC L2 moves the 15 yr $\text{L2}-U$ gap `0.019 → 0.158` EW and `0.035 → 0.180` FW).
- [ ] **Extend `robustness_samples.R` to joint $(\alpha,\beta)$ $U_{v,H}$**
  - Report both the VC $H=0$ `+2%` alpha cap and the $\hat\beta$ shift across samples (survives under both weighting rules: `0.258 → 0.552` under pooled mean, `0.327 → 0.543` under vintage median).
- [ ] **Vintage-cutoff ($2010 \to 2021$) scope decision**
  - Note: Approved plan kept cutoff and NAV-haircut figures L2-only (`revision_plan.tex`, line 382). At inception, 2010→2021 moves are small for EW/FW L2 and $U$, but at 30 yr VC FW $U$ jumps `0.121 → 1.002`—so retract the large **inception** cutoff inflation without promising a flat profile across horizons; explicitly approve scope if adding $U$ sweeps.

### 2.2 IV Moment Comparison Upgrade (`r_project/results/empirical_iv.R`)
- [ ] **Replace smallest-root / jackknife IV agreement claim with identification-robust moment-test inversion**
  - Present the $H=0$ mean pricing moment as valid under (P), transparently listing all detected roots and delete-one zero-root failures (`11/25` BO EW, `2/25` BO FW on current pools).
  - Report nominal 95% confidence sets from inverting the **centered Bartlett fixed-$b$ moment test** (with refined endpoints; on current pools: BO EW `[0.949, 3.994]`, BO FW `[0.623, 3.994]`, VC EW `[-0.244, 3.994]`, VC FW `[-0.360, 3.994]`; all cover $\beta=1$ and reach the ceiling).
  - Record centering convention, critical values, and all acceptance regions for reproducibility.
  - **Qualify calibration honestly:** Root robustness and finite-$V$ coverage under persistence are separate—at true $\beta=1$, the centered fixed-$b$ moment test rejects `4.0% / 3.4%` in independent BO/VC cores, but **`9.4% / 11.1%`** under persistent M4 ($c / 2c$, $\rho_W=0.6$) and **`21.8% / 9.8%`** in stronger Design A' ($V=25/100$). Present as a nominal, calibrated-by-design diagnostic.

### 2.3 U-First Main Exhibits (`r_project/results/main_exhibits_u.R`)
- [ ] **Update Figure E1**
  - Either keep the `±2 SE` band explicitly labelled as descriptive (~87% Gaussian D0 coverage for independent vintages, ~77% at $\rho=0.6$) or replace it with pointwise fixed-$b$ 95% bands (**`±2.790 SE`** BO, **`±2.867 SE`** VC).
  - Relabel triangles as the *expected estimator gap in the core designs* (e.g., VC EW inception `0.194` = `+0.264` target shift `− 0.070` finite-sample difference).
- [ ] **Update Table E1**
  - Report tests against $\beta=1$ alongside or instead of $\beta=0$ (all four inception $U$ fits have $p > 0.15$ against $\beta=1$).
  - Qualify FW single-factor stars (`8.1% / 8.5%` rejection at $H=0$ in `s_core_*`).
  - State horizon precision claims **by estimator and coefficient**: $H=0$ minimizes RMSE for single-factor $\beta$ and $(\alpha,\beta)$, whereas in two-factor Table S2, $H=15$ improves $U$'s *second-loading* RMSE for ME, IA, and ROE while *market-beta* RMSE rises in ROE (`0.083 → 0.152`).
- [ ] **Update Table E2**
  - Address inception alpha over-rejection (`12.1%` BO / `10.9%` VC in EW cores), absent FW joint calibration, and 30 yr distortion (alpha conditional rejection `0.4%–0.5%`, beta conditional rejection `1.07% / 0.33%`, `6.4% / 8.0%` on boundary, large finite-sample biases).
  - Either suppress stars in Table E2 (showing estimates + SEs with qualified inference) or calibrate the joint specification (including FW $(\alpha,\beta)$) before assigning stars.

### 2.4 Monte Carlo Exhibits (`r_project/simulation/mc_paper.R`) & Appendix D
- [ ] **Regenerate Step 4 paper exhibits with full reporting**
  - Add tables/panels for Nov 20 hardening designs: Design B (`s_core_BO`, `s_core_VC` FW—re-run if FW weights change) and Design A' (`pd_*`, DGP (c): inception IV Wald rejects `9.0% / 10.6%` at $V=100$; 15 yr rejects structural value `32.2% / 33.8%`).
  - Report **structural rejection**, **target rejection**, and **test availability (`tested / reps`)** separately; reserve *"size"* for a valid null; clarify that assigning $\beta_0$ as the IV target without path shocks reflects moment validity rather than root identification.
  - Qualify `table_S1_core.tex` line 54 MCSE (*"0.7 points at a 5% rate with 1,000 tested draws; observed-rate MCSE ranges roughly 0.14–1.55 points"*).
  - Document the **$V=400$ studentization finding** for equal-vintage $U$: rejection `6.8%–10.8%`, median $N_{\mathrm{eff}}$ `13–27`, left/right rejection counts `22/12`, `38/16`, `27/14`, `36/12`; holding the denominator at the MC SD removes/reverses the left-tail excess (strong evidence of studentization distortion; note Martin 2012 relevance to compounded tails without claiming it uniquely caused the asymmetry).
  - Fix `table_AD_dlp12.tex` line 3 provenance header to cite `appendix_d_summary.csv` (`aa67c7a`), not `summary.rds`.
- [ ] **Appendix D (DLP12) framing / bridge**
  - Use the $n_v \times V$ grid for the within-simulator infill finding; reserve attribution of the cross-paper S3 vs DLP12 gap (`0.76 / 0.48%` vs `0.98 / 0.05%`) for a matched bridge accounting for DLP12's simulated market and return-linked shifted additive lognormal innovation.

### 2.5 Open Step 5 Empirical Deliverables
- [ ] **Vintage-level residual figure:** $PV_v(\hat\beta_0)$ with within-vintage $\pm 1.96\,\text{SE}_{\text{within},v}$ bars (final sample).
- [ ] **Fix-$A$ vs Scale-Jacobian decomposition exhibit** (final sample).
- [ ] **$(\alpha,\beta)$ identification diagnostics:** 2D L2 & $U_{v,H}$ contours ($H=0, 15$), **NPV-gradient collinearity** $\mathrm{Corr}(\partial_\alpha PV_v, \partial_\beta PV_v)$, **Hessian cross-term** $[\nabla^2 Q_0]_{\alpha\beta}$, $U_{v,H}$ sampling correlation $\mathrm{Corr}(\hat\alpha,\hat\beta)$ (kept distinct from gradient collinearity), and **Merton (1980) drift-uncertainty benchmark** (R2-2c).
- [ ] **Four-year overlapping vintage-block bootstrap intervals** (final sample).
- [ ] **Vintage-level terminal NAV shares** (R2-4.1) and **Appendix E fold diagnostics** (B7/B8).
- [ ] **Published-estimates comparison table** (BO & VC $\beta$ and $\alpha$/PME across the literature; R1-4).

---

## Phase 3 — Theory Tightening, Revision Plan Updates & Manuscript Rewrite

### 3.1 Theory, Data Discussion, and Plan Updates
- [ ] **Add dedicated "Fund-size imputation" discussion to the main paper's Data section** *(homework for manuscript rewrite; recorded in `revision_plan.tex` Section B and Section 5 outline)*
  - Document the exact missingness (`10/383` BO and `15/432` VC funds in $V \le 2010$), explain why cross-vintage mean imputation overweights early small vintages under FW (e.g., 94.9% for one 1986 BO fund, 92.8% for two 1992 BO funds, 98.8% for two 1985 VC funds), state the chosen baseline rule, and report sensitivity across pooled mean, within-vintage median, and missing-size deletion (with vintage counts and matched support).
- [ ] **Tighten Theorems T4, T7, one-deal illustration, and negative-$U$ prose (`revision_plan.tex`)**
  - **T4** (line 1171): move $n_v^{-2}$ **inside** the expectation when fund counts vary across vintages.
  - **T7** (line 1174): keep $q_v$ inside the expectation; specify conditioning on eligibility ($n_v \ge 2$), fund sizes, and realized vintage origin $t_v^0$; state weighted pair expectation explicitly if conditional fund means differ by size; require unique identification before claiming consistency; carry conditions into E1's horizon note.
  - **Negative-$U$ interpretation** (line 104): remove any claim that a negative fitted $U_{v,H}$ establishes absence of shared vintage effects.
  - **One-deal illustration** (lines 1182–1195): state the crossing conditionally on the actual sensitivities; under $G_u=1$, note the alpha zero-crossing occurs at $2T + 1$ valuation dates ($T$ in months).
- [ ] **Tighten `2025 sdf paper/tex/sdf_paper_appendix_sdf.tex`**
  - Lines 24 & 51: specify appropriate **$\mathcal{F}_u$-stopping-time** conditions for both random entry ($a$) and bounded random exit ($b$).
  - Line 62: add uniqueness qualification to the *"identifies $(\alpha,\beta)$ directly"* wording.
- [ ] **Reconcile `revision_plan.tex`**
  - Retract the large **inception** post-2010 vintage-cutoff inflation in Section H (without promising a flat profile at 30 yr), update replies to **R2-2b**, **R2-3.1**, and **R2-4.1**, and update Section K counts (`1,246 rows, 55 designs` confirmed).

### 3.2 Manuscript & Appendix Rewrite and LaTeX Verification
- [ ] **Rewrite `2025 sdf paper/tex/sdf_paper.tex`**
  - Implement Nov 20 $U_{v,H}$-first storyline, repaired Lemma 1(3) ratio, T1–T8, fund-size imputation section, updated abstract/intro/conclusion, and new Step 4 + Step 5 exhibits; remove the five missing legacy `simulation/bias_comparison*` graphics calls.
- [ ] **Update `2025 sdf paper/tex/sdf_paper_online_appendices.tex` and build wrapper `sdf_paper_appendix.tex`**
  - Wire in `sdf_paper_appendix_sdf.tex` (checking appendix lettering shifts), remove `rem:equicorr` and old proofs, insert score centering and integrated-Hessian CLT proofs, and wire in updated Appendix D + Step 4/5 exhibits.
- [ ] **Verify both LaTeX build targets cleanly**
  - Build both `sdf_paper.tex` and the standalone appendix wrapper `sdf_paper_appendix.tex` from `2025 sdf paper/tex/` with resolved cross-references, explicit paths to all generated graphics, and clean layout.

---

## Design Decision Recorded (Do Not Revisit)

- **Cross-vintage block-excluded $U$ criterion ($Q_{H,h}(\theta) = |\mathcal{P}_h|^{-1}\sum_{|v-w|>h}\bar\epsilon_{v,H}\bar\epsilon_{w,H}$):** Investigated and rejected.
  - **Algebra:** $W_h = \mathbf{1}\mathbf{1}^\top - B_h$ is indefinite. Assuming a negligible fraction of excluded pairs and suitable convergence, its limit squares the mean horizon error $(\mathbb{E}[\bar\epsilon_{v,H}(\theta)])^2$ (note: at $H>0$, $\bar\epsilon_{v,H}(\theta)$ uses a $\beta$-dependent scale $A_{v,H}(\theta)$, whereas the implemented IV comparator fixes $z_v = A_{v,H}(1)$ at a reference $\beta$).
  - **Numerics (full 1,000 draws of `m4p_shared`, AR(1) $\rho_W=0.6$, $c=0.153$):** At $H=0$, $h \ge 2$ hits the upper bound (`3.994`) in **100%** of all 1,000 draws (empirically: upper bound from $h=1$ for BO EW/FW, $h=2$ for VC EW, $h=3$ for VC FW). At $H=15$, raising $h$ from `0` to `5` increases SD from `0.200 → 1.084` and upper-bound rate from `0% → 14.6%` (`0.226 → 1.211` and `0% → 19%` on seeds `1:100`).
  - **Chosen approach:** Keep within-vintage $U_{v,H}$ as the main estimator (with explicit $c_v=0$ and (X) caveats) and use **centered fixed-$b$ moment-test inversion (AR-style confidence sets)** on the $H=0$ mean moment as an identification-robust, calibrated-by-design diagnostic (explicitly noting persistence over-rejection under finite $V$).
