# Primary-source literature checks for G1 R4 review

Checked 17 November 2026. This file contains literature verification and the reviewer’s inferences, not empirical reproductions. No licensed data were used or transmitted. Repository unchanged. Web PDF text was inspected at the sections/pages stated; screenshots of some remote PDFs were attempted but unavailable.

## 1. Delete-one jackknife, L1, and resampling

**Checked in source:** Bradley Efron (1979), *Bootstrap Methods: Another Look at the Jackknife*, Annals of Statistics 7(1), 1–26, especially Section 3 (pp. 5–7), Section 8 (p. 24). [Original article](https://blogs.helsinki.fi/bk-club/files/2012/05/Efron_Bootstrap_AS1979.pdf), [publisher record/DOI](https://doi.org/10.1214/aos/1176344552).

Efron explicitly demonstrates inconsistency of the standard jackknife variance estimator for the sample median, despite consistency of the ordinary bootstrap there under regular density assumptions. The discussion attributes failure to extrapolation from resamples too close to the empirical distribution and describes grouped deletion as a possible remedy. This verifies the substance of the brief’s recollection without relying on its exact Efron (1982) citation. I did not reread the 1982 book.

**Reviewer’s inference for this project:** This is a counterexample to a generic guarantee, not a theorem establishing inconsistency of this exact nonlinear horizon-weighted L1 estimator. The needed result would depend on the local representation of that estimator and vintage dependence. At an exact root, finite-difference second derivatives of the empirical absolute loss are not an estimator of LAD population curvature. No apparent “smooth” entry is certified merely because its SE is non-tiny.

A kink-located minimizer is itself routine for L1 and quantile estimators. The sample median is pinned by one observation. Therefore the economic/robustness usefulness of L1 cannot be rejected solely because the estimate is a vintage’s exact pricing root. T1 implies horizon invariance only when the same root remains selected; the theta- and horizon-dependent A weights can change selection or produce smooth minima.

**Checked in source:** Peter Bühlmann and Hans R. Künsch (1993 working paper; published 1995), *The Blockwise Bootstrap for General Parameters of a Stationary Time Series*, introduction (p. 2). [ETH author report](https://www.research-collection.ethz.ch/server/api/core/bitstreams/6ee4a421-2a0b-4d08-bcf8-d2ed766f6b27/content).

The introduction explains that ordinary independent-observation bootstrap resampling needs modification for time-series dependence and describes consecutive-block resampling. The paper also emphasizes that resampling objects must match the statistic’s structure. This supports treating dependence preservation as a requirement, without proving validity for the present estimator.

**Checked in primary publisher abstract:** Gregory, Lahiri, and Nordman (2018), *A Smooth Block Bootstrap for Quantile Regression with Time Series*, Annals of Statistics 46(3), 1138–1166. [Publisher issue/abstract](https://www.imstat.org/publications/aos/aos_46_3/aos_46_3.pdf), [DOI](https://doi.org/10.1214/17-AOS1580). The abstract describes smooth tapered block-bootstrap inference, improvements over unsmoothed versions, and corrections to earlier proofs. Full technical assumptions were not inspected; do not invoke it as a theorem covering this application.

**Reviewer’s inference:** CV fold dispersion is meaningful as a conditional measure of sensitivity to the training subsets and their composition, but is not a sampling SE: folds overlap, use different sample sizes, and may have structurally different vintage support. There is no automatic t-reference distribution for estimate/fold-SD. Show the fold estimates, ranges or SD, explicitly labeled stability, and omit significance language. A dependent bootstrap or subsampling method can be developed later, but V=25–27 leaves little room for block-size asymptotics. It needs Monte Carlo validation at these V and under common calendar shocks. More deletion alone does not automatically solve dependence or basin switches.

## 2. Kiefer–Vogelsang: source versions, values, and statistic

**Checked in source:** Kiefer and Vogelsang, January 2005 Cornell Working Paper 05-08, *A New Asymptotic Theory for Heteroskedasticity-Autocorrelation Robust Tests*, Section 3.3 (printed p. 11), Table I (printed p. 24; PDF page 26). [Cornell original](https://cae.economics.cornell.edu/05-08.pdf).

This version directly tabulates the one-sided 97.5% Bartlett t critical values **2.736 at b=.26** and **2.805 at b=.28**. Section 3.3 and the table note state **50,000 replications**, Brownian approximation with normalized partial sums of **1,000 independent standard normals**, and M=bT. They allow interpolation in b. Thus the brief’s reference numbers are verified for this January working-paper table. The 97.5% one-sided t quantile is the appropriate comparison to sqrt(the 95% quantile of t²) for a two-sided 5% test.

**Checked in source:** March 23, 2005 revised manuscript, Section 4.3 and Table I (printed p. 30). [Revised manuscript](https://users.ssc.wisc.edu/~bhansen/718/KieferVogelsang.pdf).

This Table I instead reports cubic response-function coefficients fitted to simulation quantiles, with the intercept constrained to the normal quantile. For Bartlett 97.5%: c(b)=1.9600+2.9694b+0.4160b²−0.5324b³. The paper states the same 50,000 and 1,000 simulation sizes. This is a fitted approximation, not the direct January grid. A comparison against it includes fit error. Numerically c(.26)=2.7508081376 and c(.28)=2.8123591552 (independently recomputed in Python).

**Source-version caveat:** The brief should name the January Cornell table explicitly rather than implying both references are the final published Table I. I did not obtain the journal typeset version; the March manuscript calls itself forthcoming. The numerical 2.736/2.805 comparison is still legitimate against a clearly identified primary working paper.

**Checked mathematical structure:** The revised paper’s Definition 1, Case (iii), gives the Bartlett Brownian-bridge denominator

Q_q(b)=(2/b) integral_0^1 Btilde_q(r) Btilde_q(r)' dr
        −(1/b) integral_0^(1−b) [Btilde_q(r+b) Btilde_q(r)' + Btilde_q(r) Btilde_q(r+b)'] dr.

Theorem 3 gives the scalar t limit B_1(1)/sqrt(Q_1(b)) and F limit B_q(1)'Q_q(b)^{-1}B_q(1)/q. Accordingly a Wald statistic W, rather than F, must use the latter numerator **without division by q**. This matches the code’s stated intent. The project’s partial-sum implementation is the discrete analogue; the parent is independently checking equality to direct Bartlett HAC.

**Reviewer’s arithmetic/inference:** If the project’s estimated critical-value Monte Carlo SD is 0.009 per 100,000 independent draws, and the 50,000-draw reference has the same quantile density, SD of the difference is sqrt(0.009²+(sqrt(2)×0.009)²)=0.0155884572681. A 0.028 difference is 1.7962 such SDs. This is plausible sampling variation, not positive evidence of a coding error, and not proof that all conventions agree. The relevant comparison must use the actual unrounded own value. Shared random seeds across methods alter covariance; simple combined-SE reasoning assumes independence.

For integer M and iid Gaussian y, the exact expectation of the discrete demeaned Bartlett denominator y_c'K y_c/T is 1−M/T+(M²−1)/(3T²). The Brownian expectation is 1−b+b²/3, also reported by KV (revised p. 10). This provides a cheap independent check of scale. It does not establish quantile accuracy. Check path resolution separately from replication count: more draws do not remove T=1000 discretization bias.

## 3. Fixed-b applicability and bandwidth rule

**Checked in source:** The revised KV manuscript is framed in general GMM/estimating equations, not solely the Gaussian location model. Assumptions 1–4 (pp. 5–6) require consistent estimation, a score FCLT with a constant long-run scale times Brownian motion, partial Jacobians converging uniformly to rG0, and a nonsingular limiting derivative/weight structure. Theorem 3 then derives pivotal fixed-b laws. Nonlinearity alone is not an objection.

**Reviewer’s inference for vintages:** The empirical vintage score sequence must justify these assumptions. Different maturity, vintage fund counts, deterministic horizon weights, early-vintage concentration, and shared calendar shocks can generate variance and derivative profiles different from the homogeneous Brownian bridge. A correct small-b sandwich formula and a perfect Gaussian critical-value table do not establish size. Simulate the complete estimator and reported test under the actual vintage/calendar structure, including constraints and optimizer selection. Also assess conditional-on-factor-path vs unconditional uncertainty: resampling vintages alone cannot create new common factor histories.

**Checked in source:** Lazarus, Lewis, Stock and Watson (2018), *HAR Inference: Recommendations for Practice*, JBES 36(4), 541–559, Eq. (3), p. 542 and footnote 3, [author’s journal PDF](https://stock.scholars.harvard.edu/sites/g/files/omnuum5911/files/stock/files/lazarus_lewis_stock_watson_har_inference_recommendations_for_practice_jbes_2018.pdf).

The rule is S=1.3 sqrt(T), rounding S upward for Newey–West and using fixed-b critical values. The constants trade size distortion against power loss. The paper restricts its proposed procedures to positive semidefinite long-run variance estimates with known fixed-b distributions and designs with persistence no greater than AR(1) rho=.7. Its reported Monte Carlo section uses T=200 with 50,000 repetitions. Thus ceil(1.3 sqrt(V)) is accurately attributed, but it is not a guarantee at 25–27 vintages or under strongly changing vintage precision.

**Recommendation inference:** Rebuilding the fixed-b table is sensible reproducibility work, but Step 4 size/coverage under the intended DGP is more scientifically consequential. A targeted high-precision rebuild at b=7/25 and 7/27, and the restriction dimensions actually reported, is a useful first step; no need to spend equal effort on irrelevant grid corners. Preserve seeds, simulation T, draws, statistic/kernel convention, table version/hash and quantile method. The table hash belongs in inference identity; ideally distinguish point-estimation identity from inference identity so pure calibration updates do not force reoptimization. Cached Wald p-values must be invalidated or recomputed from stored Wald statistics and b/q. 1e6 draws still gives only about 100 observations in a 1e−4 tail, so dense extreme-tail precision is not implied by precise 5% quantiles.

## 4. Economic restrictions, positivity, and local minima

**Checked in source:** Hansen and Jagannathan (1991), *Implications of Security Market Data for Models of Dynamic Economies*, JPE 99(2), 225–262, Restriction 2 (p. 229) and Section IV (pp. 240 ff.), [author-hosted original](https://larspeterhansen.org/wp-content/uploads/2016/11/ImplicationsofSecurityMarketData.pdf).

HJ’s restriction is positivity of the stochastic discount factor m. Section IV studies nonnegative/positive pricing kernels and corresponding volatility bounds. The argument does not endorse a universal lower floor such as .5 on a monthly benchmark gross return.

**Reviewer’s inference:** Here m_t=1/G_t under the reciprocal gross-return specification. Requiring G_t>=c>0 both guarantees positivity and caps the one-period discount multiplier at 1/c. c=.5 gives m_t<=2 and excludes a 98.3% one-month loss. That is substantive boundedness regularization, stronger than HJ positivity; the .017 candidate is already strictly positive, so positivity alone does not exclude it. A floor is economically interpretable and preferable to arbitrary coordinate boxes as the primary diagnostic, but its value needs prespecified sensitivity. It does not cap cumulative growth or rule out discount-away minima from very large positive gross returns.

Fix the public factor date domain once across full sample, folds, and jackknife refits. “Every sample month” must be explicit: only months actually entering the estimation through 2021, or an intentionally chosen external stress history extending to 2025. The latter is a different domain incorporating post-data-end information. Enforce identical domains in all refits; dropping a vintage must not silently relax the feasible set.

Boundary inference must be addressed if the new floor binds. An ordinary unconstrained Hessian sandwich/fixed-b t test is not automatically valid for constrained estimators. Report the active floor month and minimum gross return, as well as parameter bounds. A positive-by-construction alternative is valuable as model sensitivity but still can yield cumulative collapse/explosion; it is not an automatic identification repair.

A deterministic fixed-start local minimization algorithm **is** formally a well-defined estimator if all implementation choices are fixed. The objections to option (c) are dependence on the arbitrary basin and a different scientific target from the claimed global argmin, not a categorical impossibility of defining it. Likewise, a jackknife refit leaving the full-sample basin need not be “optimizer failure”: it can be evidence of actual sample sensitivity or switching among global minima. Log all basin switches and objective comparisons. Constraining and diagnosing are complementary.

The claim that 254 entries remain unchanged needs two conditions: the reported optima are true global optima and remain feasible on the chosen floor/date domain. A 30-start check is useful evidence, not a proof. The .69–.79 range in the brief describes three reported points, not necessarily every cell; scan all relevant models and refits. If those two conditions genuinely hold, minimizing on a subset containing a global optimum preserves its minimum value (though ties can alter the selected coefficient vector).

## Recommended wording to retain vs avoid

- Retain: “The current finite-difference L1 covariance estimator is invalid; we report L1 only as descriptive robustness until a suitable inference method is established.”
- Avoid: “All L1 estimators have no valid standard errors” or “a kink itself makes L1 invalid.”
- Retain: “The extreme basin is inadmissible under our stated bounded-gross-return domain.”
- Avoid: “HJ positivity proves that .5 is the correct floor” or “leaving the old basin is inherently an optimization failure.”
- Retain: “The .028 critical-value difference is compatible with the estimated Monte Carlo uncertainty, subject to matching simulation definitions.”
- Avoid: “Numerical agreement with KV establishes valid inference at V=25.”
