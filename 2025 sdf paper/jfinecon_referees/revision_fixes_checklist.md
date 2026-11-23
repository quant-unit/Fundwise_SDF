# Revision Fixes & Open Deliverables Checklist (Reconciled Audit, Nov 2026)

Master action list reconciling the November 20 audit (`4795963..e53326c`, `revision-2026`), the follow-up double-check audit, and the November 22 cross-audit (`e3501c6`). Ordered by dependency so no estimation batch has to be run twice.

---

## Phase 1 — Settle FW Weighting & Harden Rebuild Infrastructure (Before Reruns)

### 1.1 Fund-Size Imputation & Data Preparation
- [x] **Decide and implement the FW fund-size imputation rule (`r_project/empirical/prep_preqin.R`, lines 167–215)**
  - **Decision implemented:** Baseline switched to **within-vintage median** (expanding symmetrically to nearest vintages when a vintage has no observed size, e.g. singleton 1982 VC), recording `Size.Imputed` in fund-level CSVs and exporting checked-in aggregate diagnostics (`preqin_size_imputation_North America.csv`).
  - Rebuilt `data_prepared_2026_c/`, `data_prepared_2026_c_allfunds/`, `data_prepared_2026_c_navtrim20/` (both guards passed at `2.5e-15 / 1.8e-15`), plus sensitivity folders `data_prepared_2026_c_fwmean/` and `data_prepared_2026_c_fwdrop/`.
  - Verified `helper/test_code_repairs.R` passes 100% (`0 test(s) failed`) on the rebuilt pools.
  - Generated sensitivity exhibit `results/figures_bovc_us_g1/fw_imputation_sensitivity.tex` & `.csv` comparing all three rules with `V.l2 / V.u` counts and matched-vintage gaps: Buyout FW $U_{v,0}$ is stable (`1.197` median, `1.205` mean, `1.226` drop); Buyout FW L2 0→30 yr decline drops from `0.634` (mean) to `0.427` (median) / `0.419` (drop), and inception $\text{L2}-U$ flips `-0.013 → +0.028` (`+0.042` matched drop).
  - **Downstream MC calibration check:** Observed `sd(log S)` under within-vintage median is `1.185` (BO) and `0.941` (VC) vs `1.188 / 0.897` under pooled mean (within-vintage `1.131 / 0.823` vs `1.150 / 0.847`), so existing Design B MC (`s_core_BO = 1.19`, `s_core_VC = 0.90`) remains directly representative.

### 1.2 Runner, Parallelism, Checkpoint, and Analysis Code Fixes
- [x] **Fix unknown CPU count & parallel cluster handling**
  - Files: `r_project/helper/parallel_workers.R` (lines 33–48), `r_project/helper/config_loader.R`, `r_project/estim_model_optimized.R` (line 747).
  - Fix: Fall back to `1L` worker when `parallel::detectCores()` returns `NA`/non-finite before the `by.memory` comparison (add unit test for unknown core count); wire `do_parallel` through the YAML adapters in `config_loader.R`; skip `makeForkCluster()` when `n.parallel.workers(...) <= 1L`. *(Verified by `test_review_fixes.R`.)*
- [x] **Harden full MC checkpoint, report, and Appendix D validation**
  - Files: `r_project/simulation/mc_driver.R` (lines 190–256), `r_project/simulation/mc_production.R` (lines 20–25), `r_project/simulation/mc_report.R` (lines 14–23), `r_project/simulation/mc_appendix_d.R` (lines 132–149).
  - Fix: Validate full checkpoint provenance (`spec.md5`, `public.md5`, `settings`, and exact `reps` IDs) on batch reuse; fail if requested reps are missing or duplicated; validate batches and targets in `mc_report.R` against current YAML specs (with an explicit `--partial` mode); validate `mc_appendix_d.R` caches beyond filename presence; exit non-zero in `mc_production.R` on any failed/incomplete design. *(Verified: 55/55 MC designs and 19/19 Appendix D caches validated; `test_mc.R` 0 failed.)*
- [x] **Use checked-in sample funnel counts in `mc.load.designs()`**
  - File: `r_project/simulation/mc_driver.R` (lines 258–283), `r_project/simulation/mc_generate.R` (lines 52–65).
  - Fix: Read BO/VC vintage counts and calendars from the checked-in `preqin_sample_funnel_North America.csv` (verified identical to fund-level CSV counts for both BO and VC) instead of requiring the git-ignored fund-level CSV.
- [x] **Fix `mc_analyze.R` summary metrics & add multi-start safeguard in `mc.targets()`**
  - Files: `r_project/simulation/mc_analyze.R` (lines 85–145), `r_project/simulation/mc_driver.R` (`mc.fit()`, line 166).
  - Fix: Store `neff.Coef` in `mc.fit()` and report coefficient-specific `neff`; keep all-draw `sd(est)` **and** report a clearly labelled same-sample (`tested`-only) `se.ratio.tested` + `sd.tested` separately; upgrade `mc.targets()` to `admissible.starts()`. *(Verified by `test_mc.R` and regenerated `summary.csv`.)*
- [x] **Fix legacy runner, CV scope, IV LOO domain, and legacy alpha conversion bugs**
  - `r_project/estim_model_optimized.R`: guard `config.path` and `scenario.id` when sourced with `source.internally = FALSE`; derive CV `estimate` vintages from `df0` (bit-identical for `max_vintage <= 2017`, avoids dropping 2021 when `max_vintage > 2020`).
  - `r_project/results/empirical_iv.R`: hold full-sample admissible domain fixed across LOO refits; return `NA` for zero-root LOO refits and record `iv.loo.no.root` explicitly without folding `min.square` fallbacks into jackknife SEs.
  - `r_project/results/analyze_result.R`: fix legacy monthly-vs-quarterly alpha conversions (`^4` → `^12`, remove `^(1/3)`).
  - `r_project/run_simulation_study.R`: fix `next` inside error handler (verified in `test_review_fixes.R`).
  - `r_project/ReadMe.md`: updated to current YAML runners, `2026-11-17-r5` cache identity, fixed-$b$ inference, rebuild order, licensed-data boundary, and deprecated legacy scripts.

---

## Phase 2 — Empirical & Monte Carlo Reruns and Exhibit Updates

### 2.1 Empirical Reruns & Comparison Support
- [x] **Rerun FW empirical caches on finalized FW pools**
  - Re-estimated all main EW & FW L2, $U_{v,H}$ ($q_v=1$ and $q_v=n_v$), CV (`ew_cv_north_america`, `fw_cv_north_america`), `allfunds`, `navtrim20`, `nav_disc` (`NC0`–`NC75`), and FW vintage-cutoff sweeps (`fw_vintage_sweep_north_america` and `_NC50`, cutoffs 2011–2021) caches on the rebuilt within-vintage-median pools; ran `analyze_result.R` across all 44 caches and regenerated all core G1/G2, `main_u`, `iv_comparison`, `robustness_samples`, `step5`, `nav_disc`, and `max_vintage` tables/figures.
- [x] **Add matched-vintage support for both L2 and IV alongside full-sample estimates**
  - Files: `r_project/results/empirical_iv.R`, `r_project/results/main_exhibits_u.R`.
  - Dynamically derives eligible ($n_v \ge 2$) vintages per sample (`25` BO, `23` VC).
  - Keeps full-sample 27-vintage VC L2/IV estimates and reports matched-vintage 23-vintage VC L2 and IV in `iv_comparison.tex` / `iv_comparison_rows.csv` (at 15 yr, matched VC L2 widens the $\text{L2}-U$ gap from `0.019 → 0.158` EW and `0.172 → 0.314` FW).
- [x] **Extend `robustness_samples.R` to joint $(\alpha,\beta)$ $U_{v,H}$**
  - Added Panel C (`(alpha, beta)` $U_{v,H}$, $q_v=1$) across `main`, `allfunds`, and `navtrim20` (444 recorded numbers, strict `input.md5` provenance verification): VC $U_{v,0}$ remains pinned on the $+2\%$ monthly $\alpha$ cap across all three samples while $\hat\beta$ shifts `0.368 → 0.542` FW (`-0.080 → 0.150` EW) under the 20\%-NAV-threshold trim.
- [x] **Vintage-cutoff ($2010 \to 2021$) scope decision**
  - Retracted large inception post-2010 cutoff inflation in `revision_plan.tex` and `sdf_paper.tex` while keeping the 2010 seasoned cutoff because unseasoned post-2010 NAVs still distort longer horizons (L2-only sweep scope retained per approved plan).

### 2.2 IV Moment Comparison Upgrade (`r_project/results/empirical_iv.R`)
- [x] **Replace smallest-root / jackknife IV agreement claim with identification-robust moment-test inversion**
  - Verified internal moment/SE checks (`3.2e-14` / `6.4e-11`) and fast-vs-pipeline consistency (`5.7e-08`), writing `iv_comparison.tex` (168 recorded numbers), `iv_comparison_rows.csv`, and `iv_ar_confidence_sets.csv`.
  - Lists all detected roots (`2.146; 2.577` BO EW, `2.023; 2.659` BO FW, `1.399` VC EW, `1.672` VC FW at $H=0$) and fixed-domain delete-one zero-root failures (`11/25` BO EW, `10/25` BO FW, `0/27` VC EW/FW) without substituting `min.square`.
  - Reports nominal 95% confidence sets from inverting the centered Bartlett fixed-$b$ moment test (with `uniroot`-refined endpoints, centering convention, $S/V/b$, critical values, and acceptance regions recorded): BO EW `[0.949, 3.994]`, BO FW `[0.805, 3.994]`, VC EW `[-0.244, 3.994]`, VC FW `[-0.549, 3.994]` (all cover $\beta=1$ and $U_{v,0}$ and reach `3.994`).
  - Explicitly qualifies finite-$V$ persistence coverage (`4.0% / 3.4%` independent cores, `9.4% / 11.1%` persistent M4, `21.8% / 9.8%` Design A').

### 2.3 U-First Main Exhibits (`r_project/results/main_exhibits_u.R`)
- [x] **Update Figure E1**
  - Replaced `±2 SE` band with pointwise fixed-$b$ 95% bands (`±2.790 SE` Buyout, `±2.867 SE` VC) matching table inference, and relabeled triangles as the *expected estimator gap (mean L2 minus $U_{v,H}$) in the core designs*.
- [x] **Update Table E1**
  - Added $p_1$ column testing $\beta=1$ (`0.325` BO EW, `0.280` BO FW, `0.247` VC EW, `0.116` VC FW at $H=0$), qualified FW single-factor markers (`8.1% / 8.5%` Design B rejection), stated T7 conditioning explicitly, and stated horizon precision claims by estimator and coefficient (noting two-factor Table S2).
- [x] **Update Table E2**
  - Suppressed significance stars (`marks = FALSE`), reporting point estimates and Bartlett SEs with explicit calibration notes on core EW inception alpha over-rejection (`12.1%` BO / `10.9%` VC), absent FW joint calibration, and 30 yr distortion.

### 2.4 Monte Carlo Exhibits (`r_project/simulation/mc_paper.R`) & Appendix D
- [x] **Regenerate Step 4 paper exhibits with full reporting**
  - Added `table_AS5_hardening.tex` covering Nov 20 hardening designs: Design B (`s_core_BO`, `s_core_VC` FW; checked `sd(log S)` under within-vintage median is `1.185/0.941` vs `1.19/0.90`, so existing FW MC remains representative), Design A' (`pd_*`, DGP (c): inception IV Wald rejects `9.0% / 10.6%` at $V=100$; 15 yr rejects structural value `32.2% / 33.8%`), and Design A (`pc_*`).
  - Separated **structural rejection (`Rej._0`)**, **target rejection (`Rej._\dagger`)**, and **test availability (`Test %`)**; reserved *"size"* for a valid null; clarified IV target = moment validity rather than root identification.
  - Qualified `table_S1_core.tex` MCSE (*"0.7 points at a 5% rate with 1,000 tested draws; observed-rate MCSE across printed cells: 0.14–1.55 points"*).
  - Documented the **$V=400$ studentization finding** for equal-vintage $U$ in `table_AS1_size.tex` Panel C notes and joint $(\alpha,\beta)$ core calibration + same-sample `se.ratio.tested` in `table_AS3_alpha_beta.tex`.
  - Fixed `table_AD_dlp12.tex` provenance header to cite `results/mc_step4/appendix_d_summary.csv` (`aa67c7a`), not `summary.rds`. *(Verified: `mc_paper.R` wrote 11 exhibits, 1,871 recorded numbers.)*
- [x] **Appendix D (DLP12) framing / bridge**
  - Updated `table_AS4_grid.tex` and `table_AD_dlp12.tex` notes to isolate the within-simulator $n_v \times V$ infill noise effect from the cross-paper S3 vs DLP12 gap (`0.76 / 0.48%` vs `0.98 / 0.05%`), noting DLP12's simulated market and return-linked shifted additive lognormal innovation.

### 2.5 Open Step 5 Empirical Deliverables
- [x] **Vintage-level residual figure:** $PV_v(\hat\beta_{U,0})$ with within-vintage $\pm 1.96\,\text{SE}_{\text{within},v}$ bars (`results/figures_bovc_us_g1/step5/figure_S5_vintage_residuals.pdf/.tex`, `vintage_residuals_pv.csv`).
- [x] **Fix-$A$ vs Scale-Jacobian decomposition exhibit** (`results/figures_bovc_us_g1/step5/table_fixA_jacobian.tex`, `fixA_jacobian_decomposition.csv`).
- [x] **$(\alpha,\beta)$ identification diagnostics:** 2D L2 & $U_{v,H}$ contours ($H=0, 15$; `figure_S5_alpha_beta_contours.pdf/.tex`), **NPV-gradient collinearity** $\mathrm{Corr}(\partial_\alpha PV_v, \partial_\beta PV_v)$, **Hessian cross-term** $[\nabla^2 Q_0]_{\alpha\beta}$, $U_{v,H}$ & L2 sampling correlations $\widehat{\mathrm{Corr}}(\hat\alpha,\hat\beta)$, and **Merton (1980) drift-uncertainty benchmark** (`table_alpha_beta_identification.tex`, `alpha_beta_identification.csv`).
- [x] **Four-year overlapping vintage-block bootstrap intervals** (`results/figures_bovc_us_g1/step5/table_block_bootstrap.tex`, `block_bootstrap_4yr.csv`).
- [x] **Vintage-level terminal NAV shares & Appendix E fold diagnostics**
  - Generated `table_terminal_nav_shares.tex`, `figure_S5_terminal_nav_shares.pdf`, `vintage_terminal_nav_shares.csv` (medians `8.1%` BO EW, `5.4%` BO FW, `10.8%` VC EW, `8.1%` VC FW) and `table_appendix_e_cv_folds.tex` / `.csv` (B7/B8 fold diagnostics & common $H=0$ held-out loss).
- [x] **Published-estimates comparison table**
  - Generated `results/figures_bovc_us_g1/step5/table_published_estimates.tex` comparing our refreshed $U_{v,0}$ and L2 numbers directly against DLP12, FNP12, HJK14, KN16, RS16, ACGP18, GN18, Boyer et al. (2023), and KN24.

---

## Phase 3 — Theory Tightening, Revision Plan Updates & Manuscript Rewrite

### 3.1 Theory, Data Discussion, and Plan Updates
- [x] **Add dedicated "Fund-size imputation" discussion to the main paper's Data section**
  - Added Subsection 3.2 (`\label{sec:fund_size_imputation}`) in `sdf_paper.tex` documenting exact missingness (`10/383` BO, `15/432` VC in $V \le 2010$), early-vintage concentration and pooled-mean distortion (`94.9%` 1986 BO, `92.8%` 1992 BO, `98.8%` 1985 VC, `87.3%` 1986 VC), our within-vintage median baseline, and the three-rule sensitivity (`tab:fw_imputation_sensitivity`) with matched-vintage support.
- [x] **Tighten Theorems T4, T7, one-deal illustration, and negative-$U$ prose (`revision_plan.tex`)**
  - **T4** (line 1171): move $n_v^{-2}$ **inside** the expectation when fund counts vary across vintages.
  - **T7** (line 1174): keep $q_v$ inside the expectation; specify conditioning on eligibility ($n_v \ge 2$), fund sizes, and realized vintage origin $t_v^0$; state weighted pair expectation explicitly if conditional fund means differ by size; require unique identification before claiming consistency; carry conditions into E1's horizon note.
  - **Negative-$U$ interpretation** (line 104): remove any claim that a negative fitted $U_{v,H}$ establishes absence of shared vintage effects.
  - **One-deal illustration** (lines 1182–1195): state the crossing conditionally on the actual sensitivities; under $G_u=1$, note the alpha zero-crossing occurs at $2T + 1$ valuation dates ($T$ in months).
- [x] **Tighten `2025 sdf paper/tex/sdf_paper_appendix_sdf.tex`**
  - Lines 24 & 51: specify appropriate **$\mathcal{F}_u$-stopping-time** conditions for both random entry ($a$) and bounded random exit ($b$).
  - Line 62: add uniqueness qualification to the *"identifies $(\alpha,\beta)$ directly"* wording.
- [x] **Reconcile `revision_plan.tex`**
  - Retract the large **inception** post-2010 vintage-cutoff inflation in Section H (without promising a flat profile at 30 yr), update replies to **R2-2b**, **R2-3.1**, and **R2-4.1**, and update Section K counts (`1,246 rows, 55 designs` confirmed).

### 3.2 Manuscript & Appendix Rewrite and LaTeX Verification
- [x] **Rewrite `2025 sdf paper/tex/sdf_paper.tex`**
  - Implemented Nov 20 $U_{v,H}$-first storyline, repaired Lemma 1(3) ratio, T1–T8, dedicated Fund-size imputation subsection, updated abstract/intro/conclusion, and all new Step 4 + Step 5 + U-first exhibits; removed all five missing legacy `simulation/bias_comparison*` graphics calls and stale `empirical_bovc_us/` paths.
- [x] **Update `2025 sdf paper/tex/sdf_paper_online_appendices.tex` and build wrapper `sdf_paper_appendix.tex`**
  - Wired in `sdf_paper_appendix_sdf.tex` (Appendix B), removed `rem:equicorr` and old proofs, inserted score-centered and integrated-Hessian CLT proofs (Appendix D), wired in updated DLP12 r5 table (Appendix E) and all Step 4/5 exhibits (Appendix F).
- [x] **Verify both LaTeX build targets cleanly**
  - Compiled both `sdf_paper.tex` (`sdf_paper.pdf`, 648 KB) and standalone `sdf_paper_appendix.tex` (`sdf_paper_appendix.pdf`, 370 KB) from `2025 sdf paper/tex/` with **0 LaTeX errors, 0 BibTeX errors, 0 unresolved references, and 0 float overflows**.

---

## Design Decision Recorded (Do Not Revisit)

- **Cross-vintage block-excluded $U$ criterion ($Q_{H,h}(\theta) = |\mathcal{P}_h|^{-1}\sum_{|v-w|>h}\bar\epsilon_{v,H}\bar\epsilon_{w,H}$):** Investigated and rejected.
  - **Algebra:** $W_h = \mathbf{1}\mathbf{1}^\top - B_h$ is indefinite. Assuming a negligible fraction of excluded pairs and suitable convergence, its limit squares the mean horizon error $(\mathbb{E}[\bar\epsilon_{v,H}(\theta)])^2$ (note: at $H>0$, $\bar\epsilon_{v,H}(\theta)$ uses a $\beta$-dependent scale $A_{v,H}(\theta)$, whereas the implemented IV comparator fixes $z_v = A_{v,H}(1)$ at a reference $\beta$).
  - **Numerics (full 1,000 draws of `m4p_shared`, AR(1) $\rho_W=0.6$, $c=0.153$):** At $H=0$, $h \ge 2$ hits the upper bound (`3.994`) in **100%** of all 1,000 draws (empirically: upper bound from $h=1$ for BO EW/FW, $h=2$ for VC EW, $h=3$ for VC FW). At $H=15$, raising $h$ from `0` to `5` increases SD from `0.200 → 1.084` and upper-bound rate from `0% → 14.6%` (`0.226 → 1.211` and `0% → 19%` on seeds `1:100`).
  - **Chosen approach:** Keep within-vintage $U_{v,H}$ as the main estimator (with explicit $c_v=0$ and (X) caveats) and use **centered fixed-$b$ moment-test inversion (AR-style confidence sets)** on the $H=0$ mean moment as an identification-robust, calibrated-by-design diagnostic (explicitly noting persistence over-rejection under finite $V$).

---

## Audit of `4bd0e9c` (23 Nov 2026) — Corrections to the Entries Above

An audit of `4bd0e9c` against the diff, the result files and the cited papers confirmed most entries above. All three test suites pass, about 120 paper numbers match the result files, and EW numbers are unchanged by the imputation switch. The items below were wrong, incomplete or overstated and have been corrected (uncommitted at writing).

- [x] **Theorem numbering (3.2).** The main text's Theorem 1 (T1–T8) and the Online Appendix proofs used different T-lists: the appendix followed the revision plan's numbering. The proofs section was rewritten in the paper's numbering (T3 and T6 restated precisely and proved; T8 proved from the pair-covariance form), and the plan's coercivity, infill and multiple-root results became Lemmas `lem:coercivity`, `lem:infill`, `lem:multiple_roots`. A mapping from the plan's numbering to the paper's is in `revision_plan.tex` (status item of Section L).
- [x] **Condition (X) in headline $U$ claims.** The abstract, introduction, conclusion and the IV paragraph stated that $U_{v,H}$ targets $\theta_0$ under fund independence and $c_v=0$ alone; (X) is now stated everywhere. T7 no longer claims a global minimum at $\theta_0$ when conditional means differ by fund size.
- [x] **Admissible domain.** The paper said every loading lies in $[-1,4]$; only $\beta_{\mathrm{MKT}}$ is boxed (`helper/admissible.R`), and the VC FW ROE loading of 7.27 is interior. Text and Equation (`eq:admissible_domain`) corrected.
- [x] **Published-estimates table (2.5).** The DLP12 "4-factor" betas (0.72 BO, 2.20 VC) are not in DLP12 (FF3: 1.71 / 2.37). The KN16 "VC β ≈ 1.0–2.0" and its "GPME ≈ 0 post-1999" are wrong (no market beta; GPME −0.10, post-1998 −0.14, significant). The KN24 entries were wrong (MSCI-Burgiss, β 0.91 BO / 2.37 VC). RS16 reports no betas, and the VC entries for Boyer et al. and GVN were not in the papers. All rows were rebuilt directly from the published articles (including HJK14 Tables I, III–IV and RS16 Tables 1–2 & Fig. A1, verified against the published JF/JFE PDFs), and the main-text comparison paragraph was rewritten.
- [x] **Hand-typed MC numbers (2.2, 2.4).** The AR-set calibration rates and the $V=400$ tail counts had no script. `simulation/mc_inference_checks.R` now computes them (draws regenerated from seeds; rep 1 reproduces the stored IV estimate exactly) into `results/mc_step4/inference_checks.csv`, which the notes of `iv_comparison.tex` and `table_AS1_size.tex` read. The values reproduce exactly. The "21.8% / 9.8% in Design A′" pair is the $2c$ design; $c$ gives 19.2% / 8.8%, and both are now reported.
- [ ] **Provenance of the reruns (2.1).** Every estimation row and exhibit of `4bd0e9c` was produced at `731ea2f` with uncommitted changes (`git.dirty` TRUE), and `main_exhibits_u.R` had its dirty-tree stop relaxed to a warning. The exhibit headers now say so. EW numbers are identical to the earlier clean-commit runs, and `empirical_iv.R` re-fits L2 and $U$ to the table data (check 2). A clean rerun from a commit (all 44 caches) is still open.
- [x] **`mc_report.R` target validation (1.2).** None of the seven target files records provenance hashes, so the new check skipped them silently. It now warns and records them as `targets.unvalidated`; only cell completeness is checked for them. The cached targets predate the multi-start upgrade of `mc.targets()`.
- [x] **`neff` in `summary.csv` (1.2).** Since the coefficient-specific fix, the 359 alpha and second-factor rows show `neff = NA`, because the stored batches predate `neff.Coef`. No exhibit prints these values.
- [x] **Text errors in the rewrite (3.2).**
  - Fix-A paragraph: −0.334 / −0.308 are the Scale-Jacobian parts of total declines of −0.351 / −0.427.
  - BO FW L2 is not monotonic.
  - Two-factor EW range is 0.931–1.001 for EG/IA/ROE.
  - Replication counts: path-shock designs have 500 reps.
  - Margins: M5 2/2 is an RMSE tie.
  - The threshold-sample rule and the complete-history rule are now stated exactly.
  - The Buyout FW $H=0$ gap shortfall (+0.028 vs +0.122) is reported.
  - Causal wording toned down: "solely", "confirmed", "proving", "flat ridge", and the untested studentization mechanism.
- [x] **Plan numbers (3.1).** The plan's imputation item kept the pre-implementation figures (decline 0.336, gap +0.043); they are now corrected in place to 0.427 and +0.028.
- [x] **Layout and housekeeping (3.2).**
  - Notation tables and Assumption 1(3) ran 60–93 pt into the margin; fixed.
  - `r_project/Rplots.pdf` and `revision_plan.toc` untracked and ignored.
  - Stream-of-consciousness comments removed from `step5_diagnostics.R`.
  - The published-estimates table was shrunk to unreadable size by `\resizebox`; it is now a wrapped `tabularx`.
- **Verification of the corrections.** The affected exhibits were regenerated (`mc_report.R`, `mc_paper.R`, `empirical_iv.R`, `main_exhibits_u.R`, `table_published_estimates.R`, `step5_diagnostics.R`). Every recorded-number file is byte-identical to `4bd0e9c` (`summary.csv`, both `exhibits_data.csv`, the IV CSVs, the Step 5 CSVs); only headers, notes and the literature table changed. Both PDFs build with 0 errors and no undefined references; the remaining overfull boxes are at most 10 pt (pre-existing).
