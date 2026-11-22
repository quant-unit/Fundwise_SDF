# Revision Fixes & Open Deliverables Checklist (Post-Audit, Nov 2026)

Consolidated action list combining the external review of commits `4795963..e53326c` (`revision-2026`) and the follow-up double-check audit. Ordered by dependency so no estimation batch has to be run twice.

**Legend:** `[Prior review]` = identified in the prior model's review · `[Audit]` = newly identified in the double-check audit · `[Both]` = flagged and extended in both.

---

## Phase 1 — Data & Code Infrastructure Repairs (Run Before Any Re-Estimation)

### 1.1 Data Preparation (Changes FW Empirical Results)
- [ ] **Fix cross-vintage fund-size mean imputation** `[Audit]`
  - File: `r_project/empirical/prep_preqin.R` (lines 207–211).
  - Problem: Missing fund sizes are imputed with the pooled 1982–2021 mean (`$1,842.2M` Buyout, `$379.4M` VC), giving single unknown-size funds in early vintages 51.5%–94.9% (BO) and 87.3%–98.8% (VC) of their FW vintage pool and driving the "except Buyout FW" anomaly and much of the long-horizon BO FW L2 drop.
  - Fix: Impute missing fund sizes using the **within-vintage median** (or mean) of observed same-vintage funds (falling back to nearest ±2 vintages for singleton 1982 VC).
  - Add a robustness check dropping missing-size funds entirely from FW.
  - Rebuild `r_project/empirical/data_prepared_2026_c/`, `_allfunds/`, and `_navtrim20/`, and re-run `helper/test_code_repairs.R`.

### 1.2 Runner, Checkpoint, and Inference/Analysis Code Fixes
- [ ] **Guard `parallel::detectCores()` returning `NA`** `[Prior review]`
  - File: `r_project/helper/parallel_workers.R` (lines 33–48).
  - Fix: Fall back to `1L` worker when `detectCores()` is `NA` or non-finite before comparing against `by.memory`; add a regression test covering unknown core count.
- [ ] **Harden Monte Carlo checkpointing and production exit status** `[Prior review]`
  - Files: `r_project/simulation/mc_driver.R` (lines 190–256), `r_project/simulation/mc_production.R` (lines 20–25), `r_project/simulation/mc_report.R` (lines 18–23).
  - Fix: Validate full checkpoint provenance (`spec.md5`, `public.md5`, `settings`, and exact `reps` IDs) on batch reuse; fail if requested reps are missing or duplicated; make `mc_production.R` exit non-zero on failed/incomplete designs; fail or error loudly in `mc_report.R` on empty designs.
- [ ] **Decouple Monte Carlo design loader from local licensed fund file** `[Prior review]`
  - File: `r_project/simulation/mc_driver.R` (`mc.load.designs()`, lines 258–283).
  - Fix: Read empirical vintage fund counts from checked-in aggregate counts (e.g., `preqin_sample_funnel_North America.csv` or a small checked-in CSV) so the MC suite runs without the git-ignored fund-level CSV.
- [ ] **Fix `mc_analyze.R` summary and target-search bugs** `[Audit]`
  - File: `r_project/simulation/mc_analyze.R` (lines 98, 136–140) and `r_project/simulation/mc_driver.R` (`mc.fit()`, line 166).
  - Fix: Record `neff.Coef` in `mc.fit()` and report parameter-specific `neff` (currently `neff.MKT` is copied to `Alpha` and second factors); compute `se.ratio` using `sd(est[tested])` so numerator and denominator match when `on.bound > 0`; use multi-start `admissible.starts()` instead of single-start `list(par)` in `mc.targets()` for 2-parameter models.
- [ ] **Fix legacy monthly-vs-quarterly alpha bugs in `analyze_result.R`** `[Audit]`
  - File: `r_project/results/analyze_result.R` (lines 513, 564).
  - Fix: Change `(1 + df.best$Alpha)^4 - 1` to `^12 - 1` and remove `(1 + df.f["Alpha"])^(1 / 3) - 1` (alpha is already monthly).
- [ ] **Fix `estim_model_optimized.R` legacy entry & CV year range** `[Both]`
  - File: `r_project/estim_model_optimized.R` (lines 102–125, 275–282, 649).
  - Fix: Guard `config.path` when sourced with `source.internally = FALSE`; replace hardcoded `estimate <- 1980:2020` in CV with the data's actual vintage span so `max_vintage = 2021` is never dropped.
- [ ] **Fix `empirical_iv.R` leave-one-vintage-out domain & root-failure handling** `[Both]`
  - File: `r_project/results/empirical_iv.R` (lines 79–102).
  - Fix: Hold the full-sample admissible domain fixed across LOO refits (matching `estim_model_optimized.R`); stop silently substituting `min.square` fallbacks into jackknife SEs when a refit loses every root.
- [ ] **Fix `run_simulation_study.R` error handler & update `r_project/ReadMe.md`** `[Both]`
  - Files: `r_project/run_simulation_study.R`, `r_project/ReadMe.md`.
  - Fix: Repair `next` inside error handler in `run_simulation_study.R`; update `ReadMe.md` to document YAML runners, `r5` cache identity, fixed-$b$ inference, rebuild order, and licensed-data dependencies; mark legacy `estim_model_empirical_runner.R`, `prep_preqin.Rmd`, and `simulation_improved.R` as deprecated/archived.

---

## Phase 2 — Empirical & Monte Carlo Reruns and Exhibit Updates

### 2.1 Empirical Reruns & Robustness Extensions
- [ ] **Rerun FW empirical estimation caches on repaired FW pools** `[Audit]`
  - Re-estimate FW L2 and FW $U_{v,H}$ for main (`data_out_2026-emp-c`), `allfunds`, `navtrim20`, `max_vintage`, and `nav_disc`, and regenerate dependent tables/figures.
- [ ] **Extend `robustness_samples.R` to joint $(\alpha,\beta)$ $U_{v,H}$** `[Prior review]`
  - File: `r_project/results/robustness_samples.R`.
  - Include joint-$U$ sample robustness alongside L2, explicitly reporting both the VC $H=0$ alpha cap and the $\hat\beta$ shift across samples (previously `0.258 → 0.552` for VC FW in `navtrim20` before the size-imputation fix; re-check after Phase 1.1).
- [ ] **Extend vintage-cutoff ($2010 \to 2021$) diagnostic to EW and $U_{v,H}$** `[Audit]`
  - Currently `max_vintage_*_MKT_combined_NC_data.csv` contains FW L2 only. Add EW L2 and $U_{v,H}$ (EW/FW) to verify and document the post-B4/B6 flat cutoff profile.
- [ ] **Add matched-vintage (23-vintage) VC comparison** `[Prior review]`
  - File: `r_project/results/empirical_iv.R` and `r_project/results/main_exhibits_u.R`.
  - Keep full-sample 27-vintage L2/IV estimates and add matched 23-vintage L2 ($n_v \ge 2$, excluding 1982, 1988, 1989, 1991) wherever the $\text{L2}-U_{v,H}$ gap is interpreted.

### 2.2 IV Robustness Upgrade
- [ ] **Replace fragile smallest-root / jackknife IV comparison with fixed-$b$ Anderson–Rubin (AR) confidence sets** `[Both]`
  - File: `r_project/results/empirical_iv.R`, output `r_project/results/figures_bovc_us_g1/iv_comparison.tex`.
  - Present the $H=0$ mean pricing moment as a valid restriction under (P), report all detected roots and LOO zero-root failure counts transparently (previously 11/25 BO EW and 2/25 BO FW), and report identification-robust fixed-$b$ AR 95% confidence sets (e.g., on current pools: BO EW `[0.952, 3.994]`, BO FW `[0.631, 3.994]`, VC EW `[-0.239, 3.994]`, VC FW `[-0.359, 3.994]`). Recompute on updated FW pools.

### 2.3 U-First Main Exhibits (`r_project/results/main_exhibits_u.R`)
- [ ] **Fix Figure E1 confidence bands and triangle labels** `[Both]`
  - Replace $\pm 2\,\text{SE}$ bands with fixed-$b$ 95% bands ($\pm 2.79\,\text{SE}$ BO, $\pm 2.86\,\text{SE}$ VC) or state explicitly that $\pm 2\,\text{SE}$ is not a 95% interval (~86% coverage under D0).
  - Relabel triangles as the *expected estimator gap in the core designs* (finite-sample bias differences + criterion target shift, not pure target shift).
- [ ] **Revise Table E1 stars, null hypothesis, and precision note** `[Audit]`
  - Report tests of $H_0:\beta=1$ alongside or in place of $H_0:\beta=0$.
  - Qualify FW $U_{v,H}$ inference (8.1%–8.5% size at $H=0$ in `s_core_*`).
  - Narrow the *"precision falls with $H$"* wording to single-factor and $(\alpha,\beta)$ models (explicitly noting that two-factor Table S2 favors $H=15$ over $H=0$ for ME, IA, and ROE).
- [ ] **Address Table E2 alpha calibration and long-horizon distortion** `[Both]`
  - Core EW $U_{v,H}$ rejects $\alpha=0$ in **12.1%** (BO) and **10.9%** (VC) at $H=0$; FW $(\alpha,\beta)$ has no MC calibration; $H=30$ has large finite-sample bias (`-0.232 / -0.295` in $\beta$, `+0.20% / +0.24%` in $\alpha$) and 0.4%–0.5% under-rejection.
  - Decision: either remove significance stars from Table E2 (keeping estimates and SEs with explicitly qualified inference) or run dedicated $(\alpha,\beta)$ calibration (including FW $(\alpha,\beta)$ in `s_core_*`) and calibrate/qualify markers accordingly.

### 2.4 Monte Carlo Exhibits (`r_project/simulation/mc_paper.R`) & Appendix D
- [ ] **Refresh `mc_paper.R` exhibits for Nov 20 hardening runs & size reporting** `[Both]`
  - Add tables/panels for Design B (`s_core_BO`, `s_core_VC` FW) and Design A' (`pd_*`, DGP (c)), explicitly noting Design A' IV still rejects 9.0%–10.6% at $V=100$.
  - Show **both** rejection at $\theta_0$ and rejection at $\theta_H^\dagger$ (`size` vs `size.target`), reserving "size" for a valid null.
  - Show test availability (`tested / reps`) alongside conditional rejection rates for IV (and boundary shares).
  - Document the $V=400$ over-rejection (6.8%–10.8% in Table A-S1 Panel C), low $N_{\mathrm{eff}}$ (13–27 of 400, connecting to Martin 2012 / R2-3.3), and 3:1 left-tail asymmetry driven by $\mathrm{Corr}(\hat\beta,\text{SE}) > 0$.
  - Qualify `table_S1_core.tex` line 54 MCSE (*"at a 5% rate"*) and fix `table_AD_dlp12.tex` line 3 provenance header (use `appendix_d_summary.csv`'s commit `aa67c7a`, not `summary.rds`).
- [ ] **Resolve Appendix D (DLP12) framing or add a staged attribution bridge** `[Prior review]`
  - File: `r_project/results/mc_step4/paper/table_AD_dlp12.tex` and `r_project/simulation/simulation_multifactor.R`.
  - Either frame S1–S10 explicitly as an infill-effect analysis inside this simulator (without attributing the remaining S3 vs DLP12 gap: `0.76 / 0.48%` vs `0.98 / 0.05%`), or run a small staged bridge matching DLP12's Internet Appendix (simulated market + return-linked shifted additive error) one assumption at a time.

### 2.5 Open Step 5 Empirical Deliverables (Scope & Build)
- [ ] **Vintage-level residual figure**: $PV_v(\hat\beta_0)$ with within-vintage $\pm 1.96\,\text{SE}_{\text{within},v}$ bars for BO and VC (final sample). `[Prior review]`
- [ ] **Fix-$A$ vs Scale-Jacobian decomposition exhibit** (final sample). `[Prior review]`
- [ ] **$(\alpha,\beta)$ identification diagnostics**: 2D objective contours (L2 and $U_{v,H}$ at $H=0, 15$), $U_{v,H}$ sampling correlation $\mathrm{Corr}(\hat\alpha,\hat\beta)$, and **Merton (1980) drift-uncertainty benchmark** (promised for R2-2c). `[Both]`
- [ ] **Four-year overlapping vintage-block bootstrap intervals** (final sample). `[Prior review]`
- [ ] **Published-estimates comparison table** (BO and VC $\beta$ and $\alpha$/PME across the literature; R1-4). `[Prior review]`

---

## Phase 3 — Theory Tightening, Plan Reconciliation & Manuscript Rewrite

### 3.1 Theory & Plan Corrections
- [ ] **Tighten Theorem T7** (`2025 sdf paper/jfinecon_referees/revision_plan.tex`, line 1174) `[Prior review]`
  - Keep vintage weights $q_v$ inside the population expectation; specify conditioning on eligibility ($n_v \ge 2$), fund sizes, and realized vintage origin $t_v^0$; state identifiable uniqueness explicitly; carry conditions into Table E1's horizon note.
- [ ] **Remove the negative-$U$ overclaim** (`revision_plan.tex`, line 104) `[Prior review]`
  - Do not claim that a negative fitted $U_{v,H}$ establishes the absence of shared vintage effects.
- [ ] **Fix the one-deal illustration** (`revision_plan.tex`, lines 1182–1195) `[Prior review]`
  - Correct the alpha zero-crossing under $G_u=1$ to $m = 2T + 1$ valuation dates ($H \approx 2T$, not $H \approx T$); state the $\beta$-score sign directly via $\bar S_H - S_{0,T}$ without claiming a path-independent sign.
- [ ] **Add $\mathcal{F}_u$-stopping-time assumption to SDF appendix** `[Audit]`
  - File: `2025 sdf paper/tex/sdf_paper_appendix_sdf.tex` (line 51): change *"bounded exit time $b > a$"* to *"bounded $\mathcal{F}_u$-stopping time $b > a$"*.
- [ ] **Update `revision_plan.tex` retractions, counts, and referee replies** `[Both]`
  - Add *"unseasoned post-2010 vintages materially inflate betas"* to Section H retractions and update the R2-4.1 reply.
  - Reconcile working replies to R2-2b and R2-3.1 with the Nov 20 storyline and upgraded AR-IV check.
  - Update stale Section K counts (`1,182 rows, 45 designs` → `1,246 rows, 55 designs`, or final counts after Phase 2.4).

### 3.2 Manuscript Rewrite & Build Verification
- [ ] **Rewrite `2025 sdf paper/tex/sdf_paper.tex`** `[Prior review]`
  - Replace the old $H=0$ structural-target / 12.5–15-year / bandwidth-12 narrative with the Nov 20 $U_{v,H}$-first storyline, repaired Lemma 1(3) ratio, T1–T8, updated abstract/intro/conclusion, and new Step 4 + Step 5 exhibits.
  - Remove broken legacy `simulation/bias_comparison*` graphics calls.
- [ ] **Rewrite `2025 sdf paper/tex/sdf_paper_online_appendices.tex`** `[Prior review]`
  - Wire in `sdf_paper_appendix_sdf.tex`; remove `rem:equicorr` and old consistency/CLT proofs; insert score centering, integrated-Hessian CLT proofs, and updated Appendix D + Step 4/5 appendix tables.
- [ ] **Full end-to-end build & cross-reference check** `[Prior review]`
  - Compile `sdf_paper.tex` and `sdf_paper_online_appendices.tex` cleanly from `2025 sdf paper/tex/` with zero missing figures, undefined references, or overfull-box errors.

---

## Design Decision Recorded (Do Not Revisit)

- **Cross-vintage block-excluded $U$ criterion ($Q_{H,h}(\theta) = |\mathcal{P}_h|^{-1}\sum_{|v-w|>h}\bar\epsilon_{v,H}\bar\epsilon_{w,H}$):** Investigated and rejected.
  - Algebraically, $W_h = \mathbf{1}\mathbf{1}^\top - B_h$ is indefinite and its population limit collapses to the rank-1 squared IV moment $(\mathbb{E}[\bar\epsilon_{v,H}(\theta)])^2$.
  - Numerically, at $H=0$ it hits the upper domain bound (`3.994`) in 100% of `m4p_shared` draws for $h \ge 2$ (and in both empirical Buyout for $h\ge 1$ and VC for $h\ge 2$); at $H=15$, raising $h$ from 0 to 5 increases boundary hits from 0% to 19% and SD from 0.226 to 1.211.
  - **Chosen alternative under correlated vintages:** Keep within-vintage $U_{v,H}$ as the main estimator (with explicit $c_v=0$ and (X) caveats) and use **fixed-$b$ Anderson–Rubin confidence sets** on the $H=0$ mean moment (and strictly exogenous factor-path instruments $z_v = h(\mathcal{C}_v)$ under (X)) as the robustness check.
