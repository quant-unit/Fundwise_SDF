# Review brief: JFEC revision, code repairs and Gate G1 (17 November 2026)

You are reviewing work done on branch `revision-2026` of this repository for the revision of the paper
"Identification and Stability of SDF Estimators for Pooled, Non-traded Cash Flows" (JFEC submission).
Your job is to **check that the code changes are correct and that the claims in the execution log are true**.
You are not asked to improve the paper or the code. Report problems; do not fix them.

## 1. Where things are

| What | Where |
|---|---|
| Revision plan and execution log (newest entry first) | `2025 sdf paper/jfinecon_referees/revision_plan.tex` (section "Execution log"; PDF next to it) |
| Our analysis of the referee reports | `2025 sdf paper/jfinecon_referees/review_of_referee_reports.tex` |
| Referee reports | `2025 sdf paper/jfinecon_referees/referee_1.md`, `Reviewer 2 - JFEC-2026-141 .pdf` |
| Paper (current text and tables) | `2025 sdf paper/tex/sdf_paper.tex`, `sdf_paper_online_appendices.tex`, tables in `2025 sdf paper/tex/empirical_bovc_us/` |
| R code | `r_project/` (run everything from there) |

Commits to review: `4795963` through `HEAD` (16--17 November; `git log 4795963^..HEAD`). The paper's current tables were produced at `d2bdf6e` (27 March).

## 2. What was done (summary; details and numbers are in the execution log)

1. **Code repairs B1--B9** (plan table "Code and data repairs"; log entries of 16 November):
   - B1 simulator: deals compound over `(start+1):end`; mean-one multiplicative shocks (`simulation/simulation_multifactor.R`, `deal.multiplier()` line 92).
   - B2 H = 0 convention: `max.month` 0 is an alias of 1 (`helper/objective.R`, `horizon.error()` line 132).
   - B3 sandwich SEs, Hessian, Wald test (`helper/inference.R`).
   - B4/B6 constant-composition vintage pools, FW normalization, complete-history sample rule (`empirical/prep_preqin.R`, `helper/pool_vintage.R`; output `empirical/data_prepared_2026_c/`).
   - B5 common vintage grid: origin t0 per vintage, every unit padded to the common data end (`pad.all.units()` line 71); cross-fund criterion `U_pair` (`u.pair()` line 143, `u.unit.weights()` line 153).
   - B7 SR_CV column, B8 common H = 0 validation loss, B9 provenance and Brent for one-parameter models (`estim_model_optimized.R`, `results/plot_empirical_estimates.R`).
2. **Result identity and caching** (`helper/result_cache.R`): pipeline version `2026-11-16-r3`, identity hash of settings and input files, CV completion markers, row deduplication.
3. **Speed and memory, intended to leave results unchanged**: prepared units (`prepare.units()` line 85, `unit.horizon.error()` line 110); parallel loop over types on a fork cluster created once per fold (`estim_model_optimized.R`, `iter.run()` line 636, cluster line 697); each task receives only its own type's rows (line 692); worker count capped by free memory (`helper/parallel_workers.R`).
4. **Fixes found while running G1**: singular Hessian gives NA SEs; CV resume requires completion markers; optimizer parameter scaling for models with alpha (`optim.parscale()`, `helper/objective.R` line 179; used at `estim_model_optimized.R` line 787); alpha bounds read as numbers (`helper/config_loader.R`, `alpha_bound()` line 237).
5. **Gate G1 results**: 40 North America scenarios (results `r_project/results/data_out_2026-emp-c/`, tables `r_project/results/figures_bovc_us_g1/`); Appendix D simulations S1--S10 in two panels (table `r_project/simulation/driessen_simulation.tex`, inputs `r_project/simulation/data_out_2026_g1/bias_analysis/2026-11-17_024631_*`).
6. **Attribution and endpoint decision** (log, 17 November afternoon and evening): the long-horizon decline of the estimates comes from the B5 window endpoint; decision to keep the common data end; exact re-weighting identity.
7. **Audit of the log** (log, 17 November late): 17 statements were corrected; the entry lists the corrections and the open items.

## 3. What to check, in order of priority

For each item say what you did (read code, ran a command, recomputed a number) and what you found.

**A. Code repairs B1--B9 against the plan.** Does each repair implement what the plan and the 16 November log entries say? Do the tests in `helper/test_code_repairs.R` and `helper/test_review_fixes.R` test the stated property, or could they pass with a wrong implementation (tautological tests, tolerances too loose, data subsets that avoid the hard cases)? Run both:
```
cd r_project
LANG=en_US.UTF-8 Rscript helper/test_code_repairs.R
LANG=en_US.UTF-8 Rscript helper/test_review_fixes.R
```
Both currently report 0 failures (several minutes each).

**B. The G1 fixes and the speed and memory changes.** Points worth checking:
- `optim.parscale()`: is a scale of 0.001 for Alpha safe for every model that uses L-BFGS-B or nlminb (bounds, `scale = 1 / par.scale` for nlminb)?
- One cluster per fold: workers are forked once per `iter.run()` call and then serve all horizons. Is any variable the loop body depends on changed between horizons without being exported by `foreach`? (`.noexport` lists the full data objects.)
- Per-type data split (line 692): is `split(...)[types]` equivalent to the former `df[df$type == type, ]` in every case (NA types, types absent from the data, factor levels, row order)?
- CV completion markers: can a scenario be skipped as complete when a fold is missing, or rerun although complete?
- `n.parallel.workers()`: is the macOS free-memory estimate (`vm_stat`: free + inactive + speculative pages) reasonable?

**C. Result integrity.** Do the G1 tables match the current-identity rows in `results/data_out_2026-emp-c`? Does the Appendix D table match `bias_analysis/2026-11-17_024631_bias_full.csv`? Is the S-number mapping right (`simulation/diressen_simulation.R` lines 73--82)? Is S4 dropped from Panel A only because it equals S3 there?

**D. Attribution (log, 17 November afternoon).** Reproduce it (about 20 minutes; it creates and removes two temporary worktrees inside the output folder):
```
bash "2025 sdf paper/jfinecon_referees/archive/g1_attribution/reproduce.sh" /path/to/empty/folder
```
It should print that the regenerated table is identical to `archive/g1_attribution_north_america.csv`, then the stall check. Check that the one-line patch it applies (window end = each unit's own last cash flow) is the only relevant difference between the "ownB/ownC" and "newB/newC" cells, and that the comparison really holds everything else fixed.

**E. Re-weighting identity and decomposition (log, 17 November evening).** Check the claim that for every vintage, horizon and parameter the horizon error equals `PV_v(theta) * A_vH(theta)`, so the horizon changes only the weights of squared loss, not the pricing conditions. Read `horizon.error()` and `getNPVs` (`helper/getNPVs.R`) and compare with the identity as stated in our review document (`review_of_referee_reports.tex`, lines 97--98; the scale-Jacobian argument is at line 71) and with the horizon criterion of the current paper (`sdf_paper.tex`, Section 2, which frames it as the "Compounding-Horizon Wedge"; the identity itself is planned for the revised text as T1). Check the window convention (B2: `max.month = m` averages the first m monthly dates). Then run
```
cd r_project
Rscript "../2025 sdf paper/jfinecon_referees/archive/check_horizon_reweighting.R" VC 360
```
Is the slope decomposition (pricing term vs. scale term) computed correctly, and does it support what the log says, and no more?

**F. The log itself.** Read the execution log entries from 16 November (evening) to 17 November (late). For every number and causal statement: is it supported by the files? Is anything still vague, overstated or wrong? Is "verified" kept separate from "consistent with" and "not tested"? Are corrections marked?

## 4. Known open items: do not report these as new findings

They are listed in the log entry of 17 November (late). In short:
1. Long-horizon SEs may understate uncertainty when few vintages carry the weight (e.g. at beta = 1 one VC vintage, 1994, carries 42% of the 30-year EW objective).
2. Multi-factor loadings differ substantially from the paper's (24 of 104 Buyout loadings change sign); not yet traced to data or endpoint.
3. The (alpha, beta) objective is nearly flat along a ridge at 10--20 years (weak separation of alpha and beta).
4. Appendix D: the cause of the bias at sigma = 30% is untested; the Appendix D text is not yet rewritten.
5. U_{v,H} results have no table; G1 tables are not yet copied into the paper.
6. Housekeeping: the config file hash is part of result identity; committed result CSVs contain absolute local paths; the commit message of `43a5aff` says S4 where it should say S6.

If you find something that changes one of these (for example evidence on its cause), report it as a finding.

## 5. Rules

- **Do not change code, configuration, results or the paper.** Put any experiment in a folder outside the repository and pass it as the output folder (see `archive/g1_attribution/run_cell.R` for how to run a scenario into another folder with a copy of the config).
- **Do not edit `r_project/config/*.yaml`**: the config file's hash is part of every result's identity, so an edit invalidates all cached results.
- **Licensed data**: fund-level Preqin files and `r_project/empirical/data_in/` must never be committed or sent to an external service.
- Run R with `LANG=en_US.UTF-8` (the YAML files contain non-ASCII characters).
- Memory: the machine has 32 GB. Run one heavy job at a time; set `SDF_MAX_WORKERS=4` if anything else is running. Do not rerun the full G1 (about 40 minutes) or the Appendix D simulations (several hours); verify from the stored results instead.
- The old estimator at `d2bdf6e` assigns global variables (for example `out`): in your own driver scripts, use names that cannot collide (the archived drivers prefix them with a dot).
- Archived drivers of the G1 runs: `2025 sdf paper/jfinecon_referees/archive/g1_runs/` (the simulation retry driver reran the four scenarios that failed in the first pass). They write into the G1 result folders (`results/data_out_2026-emp-c`, `results/figures_bovc_us_g1`, `simulation/data_out_2026_g1`), so read them, but do not run them as they are.

## 6. Report format

One list of findings, most severe first. For each finding:
- **Claim or code** checked (log entry and sentence, or file:line).
- **Verdict**: *Wrong*, *Overstated or vague*, *Verified by reproduction*, *Verified by reading*, *Plausible but not checked*.
- **Evidence**: the command and its output, or the lines of code, or the numbers recomputed from which file.
- **Severity**: does it change a result, a table, a claim in the paper, or only the log?

Then, separately: items from section 3 you could not check, and why. Do not round off: if you only read the code, say so; if a number matches only approximately, give both values.
