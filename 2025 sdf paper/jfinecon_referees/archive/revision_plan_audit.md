# Audit of Fable's revision plan

15 November 2026. Reviewed commit `4942249896d1827f40116e3667f8a6da4e061c7b`, against its parent `14ef2fb2eda3cfb659294b87e27c6c813305e909`, the supplied task, the current referee review, and the relevant R implementation. References below point to the plan after compilation/layout repairs; its substantive wording and numbers were preserved.

**Assessment:** retain the proposed direction, but repair the mathematical claims and experimental design before implementing it. The exact rescaling identity, distinction between pricing moments and squared-loss targets, separation of fund and vintage noise, and repair-first sequence are useful. The plan addresses the requested triage and minimal/effective comparison. It does not yet provide a mathematically reliable blueprint in several important places.

**Compilation:** the original source failed because it redefined AMS-LaTeX's `\tag` with `\newcommand`. Renamed the evidence-label macro to `\evidencetag`, repaired line wrapping, compacted the contents page, and adjusted spacing/page breaks. The final PDF has 16 pages; the final build has no LaTeX warnings, unresolved references, or overfull/underfull boxes. All pages were visually inspected. No substantive plan, manuscript, data, or estimator changes were made.

## Standards: mathematics and evidence

This axis checks correctness rather than repository style; no applicable repository standards file was found. The seven findings from the independent mathematical review are retained below, with verification details.

### S1. [P1] The proposed learned instrument is not automatically valid

[Plan, instruments paragraph](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:277>) says both candidates depend only on the benchmark path and that leave-one-vintage-out prediction removes noise correlation. Candidate (ii) also depends on other vintages' realized cash flows. Under the plan's persistent common shocks, those observations can be correlated with the held-out vintage's shock.

For a simple counterexample, let vintage errors be a stationary, unit-variance AR(1) process with correlation 0.8, independent of benchmark paths. Use the mean of the other vintages' errors as the learned component of an instrument. With 25 vintages and the middle vintage held out,

\[
E[z_{-v}u_v]=\frac{1}{24}\sum_{w\ne v}0.8^{|v-w|}=0.310427,
\]

not zero. Each residual still has mean zero conditional on its benchmark path. Leaving out its own observation is insufficient.

**Repair:** keep fixed-reference horizon weights as the clean starting candidate. Present the learned gradient as a research extension, requiring a validity proof, independence or suitable conditional restrictions, and dependence-aware training splits. A useful primary-source precedent is the neighbor-excluding construction in [Semenova, Goldman, Chernozhukov and Taddy (2023), Section 3.1](https://onlinelibrary.wiley.com/doi/full/10.3982/QE1670). That paper supports the need for a dependence-aware construction; it does not establish validity for this particular proposed instrument. A reference parameter learned from the same outcomes also requires treatment as an estimated nuisance parameter.

### S2. [P1] T5 does not prove a global exploding-alpha solution

[Plan, T5](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:248>) compares the criterion at infinity with its value at the truth, then says the degenerate solution wins whenever noise exceeds the squared first call. That comparison can rule out the truth as the global minimizer, but cannot establish infinity as the minimizer.

Take an initial call of -1 and a one-period payoff of 10 with probability 0.1, zero otherwise. The true gross return is 1, and conditional pricing holds. For trial gross return \(g>0\),

\[
Q(g)=E[(-1+Y/g)^2]=1-2/g+10/g^2.
\]

Here \(Q(1)=9\), the infinite-return limit is 1, and the global minimum is \(Q(10)=0.9\). Infinity beats truth but loses to a finite parameter.

**Repair:** state the limit and the sufficient condition for the truth to be beaten, separately from any global boundary theorem. Distinguish a large-alpha diagnostic from minimization on the compact admissible set in A7.

There is also an origin problem: the positive limit in B2 requires a nonzero cash flow at the valuation origin. Under A5, the predetermined origin may precede every call. In the supplied predetermined-origin probes, **all 25 Buyout and all 27 Venture vintages have zero initial cash flow**. At such an origin, the correctly normalized inception criterion also tends to zero as alpha increases. Unifying the 0/1-month convention does not universally remove the zero-limit degeneracy. State the limit as the squared cash flow at the chosen origin, possibly zero, and reconcile this with the economic definition of the claim.

### S3. [P1] T6 mixes fixed-vintage infill with growing-vintage asymptotics

[Plan, T6](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:249>) says infill with a common component converges to the minimizer of signal plus vintage-level noise. With the plan's fixed number of vintages, common shocks remain realized random variables; infill does not average them into their variance.

For \(PV_{vi}(\theta)=\theta-1+u_v+e_{vi}\), increasing funds removes the average of \(e_{vi}\), but leaves

\[
V^{-1}\sum_v(\theta-1+u_v)^2,
\]

whose minimizer is \(1-\bar u\). With realized shocks 0.2 and 0.4, it is 0.7. The expected signal-plus-variance criterion, for centered shocks, is minimized at 1. These are different objects.

**Repair:** give the random conditional limit for fixed-V infill. Add a growing-vintage LLN before claiming convergence to the unconditional population criterion. Retain the review's identification and uniform-convergence conditions for the no-common-shock consistency result. T7 can compare population targets, but should not equate the limits of these different sampling regimes.

### S4. [P2] Beta does not have the asserted monotone effect

[Plan, opening explanation](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:68>) treats higher alpha and higher beta alike. Higher alpha increases every admissible benchmark gross return; beta multiplies signed factor returns.

For a one-period factor return -0.1, \(G=1-0.1\beta\) and the two-date weight is \(A=(1+G)/2\), so \(\partial_\beta A=-0.05\). The claimed growth of \(A^2\) with beta is false even on the positive-return domain.

**Repair:** retain the directional alpha result under its stated assumptions. Describe the beta profile and the 12.5–15-year cancellation as hypotheses about a particular corrected simulation design, requiring actual verification. Do not infer them from the alpha score.

### S5. [P2] Simulator conditional pricing remains unproved after the month repair

[Plan, evidence table](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:100>) says the off-by-one month is the only exception to conditional pricing. In [the simulator](</Users/christausch/github/Fundwise_SDF/r_project/simulation/simulation_multifactor.R:225>), negative cumulative paths trigger a zero payoff. Gaussian additive shocks followed by that truncation change conditional means. The current referee review already records this distinction.

For a one-period illustration, \(G=1\) and \(e\sim N(0,4)\) give \(E[(G+e)_+]=1.395593>1\). This high-noise example proves the logical failure; it does not measure the error at the paper's calibration. A zero-noise test cannot detect it.

**Repair:** specify positive, conditionally mean-preserving shocks, or derive and implement the appropriate default/payoff normalization. Verify nonzero-noise conditional pricing as well as the timing identity. Also qualify the evidence table's last-cash-flow window claim: the empirical branch uses that endpoint, but the simulation branch of [the monthly constructor](</Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:184>) pads through December 2019.

### S6. [P2] The empirical narrative is stronger than the evidence

[Plan, opening empirical claim](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:77>) says the cited vintages deviate by several within-vintage standard errors each. Direct reproduction at beta 1.21 gives:

| Buyout vintage | Pooled residual | Within-vintage SE | Residual / SE |
|---|---:|---:|---:|
| 1996 | 0.1542 | 0.0865 | 1.78 |
| 2006 | -0.0247 | 0.0346 | -0.71 |
| 2009 | -0.0359 | 0.0252 | -1.42 |
| 2010 | -0.0200 | 0.0453 | -0.44 |

The sign runs are present, but the universal significance wording is false. There is also a positive 2005 vintage between the stated positive and negative runs. These within-vintage standard errors are descriptive and omit parameter-estimation uncertainty and cross-fund/common dependence.

**Repair:** describe persistent residual patterns as suggestive evidence. They do not identify the conditional covariance component \(c_v\) separately from conditional mean misspecification, selection, reporting or NAV effects. Likewise, a changed joint minimizer does not prove that alpha is unidentified; establish local rank/curvature and uncertainty before putting that conclusion in the abstract.

### S7. [P2] The proposed timing experiment does not ensure a wedge

[Plan, simulation design (c)](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:385>) assumes market-dependent exits violate (X) while preserving (P). That does not follow from the exit hazard alone. A payoff equal to realized benchmark growth at any benchmark-dependent exit has zero NPV pathwise and therefore still satisfies (X).

**Repair:** construct and check the joint payoff/timing model explicitly. A simple mean-wedge illustration is a bounded benchmark state \(S\in\{-1,1\}\), equally likely and independent of inception information, with discounted residual \(kS\), where \(0<k<1\) (payoff equals benchmark growth times \(1+kS\) after a unit call). Its inception mean is zero, while a state-dependent horizon weight has nonzero weighted mean when \(E[A(S)S]\ne0\). A timing-based version needs an equally explicit verification; market-timed exits alone are insufficient.

## Spec: coverage of the requested revision task

The independent spec review found the requested triage, revised assumptions, implementation agenda and minimal/effective comparison present. It found no material scope creep. Its findings remain separate from the standards review; overlap reflects two independent checks of the same requirement.

- **[P1] Sound mathematics: invalid learned-instrument guarantee.** The instruments paragraph promises orthogonality from leaving out one vintage despite persistent vintage dependence. See S1 for the counterexample and required conditions.
- **[P1] Sound mathematics: unsupported global exploding-alpha theorem.** T5 establishes a comparison with truth, not a global optimum. See S2.
- **[P1] Bug-free implementation: B1 omits the known default-truncation pricing defect.** The indexing repair and zero-noise test leave this stochastic failure in place. See S5.
- **[P2] Minimal revision: indispensable proof repairs are classified as optional.** [W5](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:326>) contains repaired asymptotics and is marked outside the minimal scope, while the minimal milestone is called a correct paper. Fixing the numerical sandwich in W1 does not repair the score-centering condition or the invalid common intermediate-Hessian argument. Move A8 and these proof repairs into the minimal scope, or remove the affected theorem and inference there.
- **[P2] Effective experiment: design (c) does not yet implement its claimed separation of (P) and (X).** Supply the DGP and its moment verification before using it to demonstrate the wedge. See S7.

The follow-up implementation check found another required repair:

- **[P1] B3 misses a separate gradient-vector bug.** [Current gradient code](</Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:507>) assigns `p.plus <- par[x] + delta` and likewise for minus. This replaces the whole vector with one scalar. The valuation function then drops all other parameter columns. Repair with `p.plus[x] <- par[x] + delta` and `p.minus[x] <- par[x] - delta`, after copying the full vector. Switching the loss to the unit squared error alone leaves this bug intact. With cash flows `(-1, 1.5)`, RF `(0, 0)`, market `(0, .1)`, alpha regressor `(1, 1)`, `(MKT=1, Alpha=.05)`, `max.month=1`, and zero penalty, the current absolute-error gradients are `(-.123967, -1.360544)`; correct full-vector perturbations give `(-.113422, -1.134216)`. No later definition overrides the function. B3's toy test must include multiple parameters, and preferably at least three to catch the separately identified mixed-Hessian bug.

## Additional reproducibility finding

**[P2] The Venture table uses an undocumented sample restriction.** The [appendix table](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan.tex:478>) does not reproduce from the cited script as written. The Buyout table, the product columns, the reported mean roots and joint grid minima reproduce. The Venture pooled and predetermined-origin squared-loss columns reproduce exactly only after imposing `minn=2`, which the script does not impose for those columns.

| Horizon | VC pooled LS, plan | Cited script | VC predetermined LS, plan | Cited script |
|---:|---:|---:|---:|---:|
| 0 | 1.93 | 1.94 | 1.97 | 1.98 |
| 5 | 2.20 | 2.11 | 2.16 | 2.13 |
| 10 | 1.78 | 1.74 | 1.80 | 1.78 |
| 15 | 1.32 | 1.23 | 1.34 | 1.28 |
| 20 | 1.24 | 1.17 | 1.23 | 0.90 |
| 30 | 0.90 | 0.83 | -0.48 | -0.50 |

**Repair:** use a stated common eligible-vintage sample for the direct comparison, or label the different sample rules explicitly, and generate the table from the saved output. The all-vintage predetermined VC result at 30 years is on the search-grid lower bound (-0.5); it is not an established interior optimum.

## Recommended adjustments to the roadmap

Keep W1–W7 and the gated status of W8, subject to the findings above. Make the minimum version include all essential proof repairs, multi-parameter derivative tests, exact pooling identities, and matched stochastic DGP checks. Treat origin choice, horizon endpoints and positivity/default conventions as a single valuation specification used by theory, R and the Python probes.

Write the valid theorem statements before promising a new estimator or committing to the target abstract. The common-shock decomposition and fixed-reference horizon moments can remain central. Learned instruments, a global exploding-alpha theorem, and a claim that common noise dominates the empirical data require additional work.

Choosing not to develop an optimal-H rule is a reasonable scope decision. The stated reason that a common population target leaves “nothing to tune” is too strong: weights can still change curvature and sampling variance. Similarly, a finite-sample horizon profile is a diagnostic, not by itself a test of (X); formal comparison requires joint uncertainty and attention to other sources of misspecification. Do not remove the multi-root discussion merely because a particular estimator extension is not adopted; it remains an identification caveat wherever mean moments are discussed.

## Verification and limits

- Ran `latexmk -pdf -interaction=nonstopmode -halt-on-error revision_plan.tex` and inspected rendered pages.
- Ran the complete existing `check_revision_probes.py` against the stored global files.
- Added [independent audit checks](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/check_revision_plan_audit.py>) and saved [their output](</Users/christausch/github/Fundwise_SDF/2025 sdf paper/jfinecon_referees/revision_plan_audit_checks.txt>). All deterministic assertions passed. These checks demonstrate counterexamples and reconcile the table sample; they do not validate the proposed estimators generally.
- Checked the relevant current R paths for loss, gradients, Hessians, Wald calculation, simulation timing/defaults and monthly endpoints. This is not a complete codebase audit or a regeneration of all manuscript tables.
- No raw-workbook rebuild, complete revised simulation campaign, estimator implementation, novelty determination or exhaustive literature review was performed. The global exploratory probes do not establish corrected North America estimates or empirical identification claims.

Standards axis: 7 findings, with P1 validity/target errors. Spec axis: 6 findings including the follow-up, with P1 mathematical and implementation omissions. The additional table-reproduction finding is P2; counts across axes overlap.
