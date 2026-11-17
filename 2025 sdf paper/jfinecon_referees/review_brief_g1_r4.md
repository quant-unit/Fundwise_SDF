# Review brief: open decisions after the G1 rerun under the new inference (17 November 2026)

You are reviewing work on branch `revision-2026` of this repository, for the revision of the paper
"Identification and Stability of SDF Estimators for Pooled, Non-traded Cash Flows" (JFEC submission).
Your job has two parts:

1. **Check the findings** that led to three open decisions (section 3). Are they true, and is the diagnosis right?
2. **Assess the options** for each decision (section 4) and recommend one, with reasons. Say where you rely on
   literature, and whether you checked it in the source or quote it from memory.

You are not asked to fix anything. Report; do not change code, results or the paper.

## 1. Where things are

| What | Where |
|---|---|
| Revision plan; "State of evidence" section (A--J) at the front, execution log after it (newest first) | `2025 sdf paper/jfinecon_referees/revision_plan.tex` (PDF next to it) |
| Log entry to review | "17 November 2026 (G1 rerun under the new inference)" |
| Check scripts behind every number below | `2025 sdf paper/jfinecon_referees/archive/g1_r4/` (each has its run command in the header) |
| Estimator | `r_project/estim_model_optimized.R` (L2 criterion line 391, L1 line 406, start values lines 755--773, optimizer choice around line 808, jackknife `jackknife.se()` line 580) |
| Sandwich SEs, numerical derivatives | `r_project/helper/inference.R` (`num.grad` 37, `num.hess` 48, `bartlett.weight` 73, `sandwich.covariance` 82) |
| Fixed-b p-values | `r_project/helper/fixedb.R`, table `r_project/helper/fixedb_quantiles.rds` |
| Result identity / caching | `r_project/helper/result_cache.R` (`RESULT.SETTINGS` line 17) |
| Results and tables | `r_project/results/data_out_2026-emp-c/` (result rows), `r_project/results/figures_bovc_us_g1/` (tables, `*_data.csv`) |

Commits: `90e5421` (implementation of the inference decision), `af61927` (rerun, checks, documentation).
The earlier review of G1 (brief `review_brief_g1.md`) covered everything up to `6d98332`.

## 2. Context in brief

- **Model.** Cash flows of vintage-year portfolios (VYP; Buyout 25 vintages, VC 27, North America, vintages up to 2010) are discounted with the cumulative gross return of a linear benchmark, `1 + r_f + beta'f` per month. The pricing error of vintage v at horizon H is the average discounted net cash flow over the first 12H months after the vintage's first cash flow; H = 0 is the inception NPV.
- **Criteria.** L2: `sqrt(sum_v e_v^2) / V`. L1: `sum_v |e_v| / V`. Optimizers: Brent for one parameter, L-BFGS-B with bounds and parameter scaling when alpha is estimated, Nelder--Mead otherwise; start `beta_MKT = 1`, other parameters 0.
- **Identity T1** (proven, in the plan, section C): `e_{v,H}(theta) = A_{v,H}(theta) * PV_v(theta)` with `A > 0`, so `e_{v,H} = 0` iff `PV_v = 0` at every horizon; the horizon only reweights vintages.
- **Inference (decided and implemented).** Bartlett (Newey--West) sandwich over vintage years with `S = ceil(1.3 sqrt V)` = 7 (Lazarus, Lewis, Stock, Watson 2018), Kiefer--Vogelsang (2005) fixed-b p-values with `b = S/V` (0.280 Buyout, 0.259 VC) from our own simulated table; leave-one-vintage-out jackknife SEs and the effective number of vintages `N_eff` reported as robustness.
- **What the rerun showed** (section E of the state of evidence): point estimates and independent-units SEs are identical to the bandwidth-12 run in all 60,736 result rows; SEs of beta rose by a median 9% (Buyout) and 3% (VC); fewer coefficients are significant. These are not in question here unless you find an error.

## 3. Findings to check

For each, say what you did (read code, ran a script, recomputed a number) and what you found.

**F1. L1 standard errors are not valid.**
Claim: in 250 of 288 full-sample L1 table entries SE_A is below 0.01 or |t| exceeds 50 (same with bandwidth 12). In 40 of 48 single-factor L1 cells the estimate sits at a kink: one vintage's error is at most 2.3e-8 while the next smallest is 6e-4 to 0.08, and the criterion rises on both sides. The numerical Hessian of |e| (step 1e-4) then explodes and the SE collapses. In the other 8 cells the minimum is smooth and very flat (SE_A 0.5 to 12). The estimate is identical across horizons when the same vintage has zero error (Buyout EW 1.25225 at 0, 5 and 10 years, vintage 2010), consistent with T1.
Scripts: `l1_kink.R`, `table_diff.R`. Check: is the mechanism right? Are the independent-units SEs (also ~1e-4) broken for the same reason? Do any L1 SEs in the tables remain usable?

**F2. MKT+EG at H = 0 has a degenerate lower minimum.**
Claim: a multi-start search (reported point, the pipeline's start, 30 random starts; exact pipeline criterion, reproduced to 4e-15) over all 260 full-sample entries of the Buyout and VC Multi tables finds 254 minima. The other 6 are MKT+EG at H = 0 (and the H = 1 alias) for Buyout FW, VC EW, VC FW: the criterion is 11.5%, 13.4%, 16.4% lower at beta ~ 3.2, EG loading ~ 13.65, where the benchmark gross return is 0.017 in January 2021 (SD of the monthly gross return 0.26 over 1969--2025 factor data). At the reported points the lowest gross return is 0.69--0.79. The same edge point appears for NATRES (4.30, 13.87) and earlier in a NATRES CV fold (beta ~ 8.3, EG ~ 13.3).
Scripts: `stall_check.R`, `edge_minimum.R`, `edge_sdf_return.R`, `natres_eg.R`. Check: is the search adequate (box `[-1,4] x [-4,4]` for random starts)? Is the "degenerate" reading right? **Not checked by us:** which cash flows the edge solution exploits. A gross return of 0.017 in one month multiplies the discounted value of every later flow by about 59; the final NAV at the data end (2021-12-31) is the obvious candidate. Please check this.

**F3. One jackknife SE is contaminated by F2.**
Claim: VC FW MKT+EG at H = 0: leaving out vintage 2003 sends the Nelder--Mead refit (started at the estimate) to the edge (beta 2.860, EG 13.402). The jackknife SEs are therefore 2.634 and 9.780; the latter is printed in the robustness table at 0 and 0.1 years. The full-sample self-check `jk.selfcheck` (0.0011) does not detect this. Script: `jackknife_vc_fw_eg.R`.

**F4. Fixed-b critical values differ slightly from the published table.**
Claim: our table gives two-sided 5% critical values of 2.777 (b = 0.28) and 2.730 (b = 0.26); Kiefer--Vogelsang Table I gives 2.805 and 2.736. Measured Monte Carlo SD of ours: 0.009 at b = 0.28 per 1e5 draws; theirs used 5e4 draws, so roughly sqrt(2) larger. 16 of 864 full-sample L2 p-values in the tables lie in [0.045, 0.055]. Check: is a 0.028 gap at b = 0.28 really within Monte Carlo error, or does it point to a difference in the simulated statistic (`fixedb.simulate()`, partial-sum formula, `T = 1000`)?

## 4. Decisions: options and our current recommendation

**D1. L1 inference.**
- (a) Drop SE_A, t and significance markers from the L1 tables; report the estimates and the cross-validation fold dispersion only. *Our recommendation for now.*
- (b) Resampling over vintages (bootstrap, subsampling, delete-d jackknife). Its validity for a kink-located estimator with V = 25--27 would have to be shown, e.g. in the Step 4 Monte Carlo.
- (c) Analytic asymptotics for least-absolute-deviation-type estimators. These need the density of the errors at zero, hard to estimate from 25 units.
- (d) Move L1 out of the main paper, given T1: the L1 estimate is pinned by the vintage with zero error.
Questions: which do you recommend? Is the delete-one jackknife inconsistent here as it is for the sample median (Efron 1982; quoted from memory, not re-read)? Does the cross-validation dispersion say anything valid about sampling uncertainty?

**D2. Excluding degenerate SDFs (from F2 and F3).**
- (a) Constrain the estimator so the benchmark gross return stays above a floor in every sample month (for example 0.5), and apply the same constraint in the jackknife refits. Sensitivity to the floor would be reported. *Our current leaning.* The floor value is arbitrary.
- (b) A bounded parameter space for loadings (box constraints). Also arbitrary, and less directly tied to the problem.
- (c) Keep the local minimum reached from the fixed start and state this in the paper. This is not a well-defined estimator; we consider it not acceptable.
- (d) A positive-by-construction SDF (the code has an `exp.aff` model). This changes the model, not just the estimator.
- (e) A penalty (lambda > 0; the L2 criterion already has a lasso term, unused with lambda = 0).
Questions: which do you recommend? Is there a standard in the SDF or GMM literature (e.g. positivity restrictions as in Hansen--Jagannathan)? Should a jackknife refit that leaves the basin of the full-sample estimate be detected and reported instead of constrained? Would (a) change any of the 254 confirmed entries? It should not, if the floor lies below 0.69.

**D3. Precision of the fixed-b table.**
- (a) Rebuild with 1e6 draws per grid point (Monte Carlo SD about 0.003; about 2 hours). This needs a `PIPELINE.VERSION` bump: the cached Wald p-values depend on the table, which is not part of the result identity. Single-coefficient p-values in the tables are recomputed at analysis time. *Our recommendation.*
- (b) Keep the table and note the Monte Carlo error.
- (c) Use the published Kiefer--Vogelsang quantiles, interpolated. These exist only for some quantiles, so p-values would be coarse.
Question: is (a) worth it given that the size of these tests at V = 25 is itself still unchecked (planned for the Step 4 Monte Carlo)? Should the table's hash enter the result identity?

## 5. Rules

- **Do not change code, configuration, results or the paper.** Put any experiment in a folder outside the repository.
- **Do not edit `r_project/config/*.yaml`**: the config hash is part of every result's identity.
- **Licensed data**: fund-level Preqin files (`r_project/empirical/data_in/`, `r_project/empirical/data_prepared_2026_c/*`) must never be committed or sent to an external service.
- Run R from `r_project/` with `LANG=en_US.UTF-8`. Memory: 32 GB; one heavy job at a time; set `SDF_MAX_WORKERS=4` if anything else runs. Do not rerun G1 (about an hour); the scripts in `archive/g1_r4/` run in seconds to a few minutes (`stall_check.R` is the longest, several minutes in parallel).

## 6. Report format

First, for F1--F4: verdict (*Wrong*, *Overstated or vague*, *Verified by reproduction*, *Verified by reading*, *Plausible but not checked*), the evidence (command and output, or code lines, or recomputed numbers and from which file), and severity.
Then, for D1--D3: the recommended option, the reasons, what it would change in the tables or the paper, and any literature you rely on, marked as checked in the source or from memory.
Finally: anything else you noticed that bears on these decisions. Do not round off; if a number matches only approximately, give both values.
