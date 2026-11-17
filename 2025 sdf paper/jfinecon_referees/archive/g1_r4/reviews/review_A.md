# Review of the G1 rerun and its three open decisions

17 November 2026. Reviewed checkout at af61927 on revision-2026, including the brief, the relevant state-of-evidence and execution-log entries, estimator/inference/cache/reporting code, archived diagnostics, and locally stored results and pooled cash flows.

**My recommendation:** adopt D1(a); adopt D2(a) only as an explicitly justified, sensitivity-tested regularization of the estimator, with boundary diagnostics and consistent refitting; adopt D3(a), with inference provenance, while giving the estimator-level size study higher scientific priority. Do not finalize the inference tables before resolving D2.

The most consequential new result is that a monthly gross-return floor of 0.5 would exclude **23 of the 254 otherwise confirmed reported Multi-table points**, all in VC. It cannot be presented as a repair confined to the six EG rows. The extreme EG solution also is not explained solely by December 2021 terminal NAV.

No repository files, configuration, paper, or stored results were changed. All experiments are under /private/tmp/g1_r4_review/. No licensed data were sent to a web service. The literature check used public papers only.

## F1. Current L1 standard errors

**Verdict: Verified by reproduction. Severity: critical for every current L1 inferential statement.** This verdict concerns the implemented covariance calculation, not the possibility of valid L1 inference in general.

### Evidence

I read the L1 loss in [estim_model_optimized.R](/Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:405), its inference contribution at line 556, and the finite-difference Hessian and sandwich in [inference.R](/Users/christausch/github/Fundwise_SDF/r_project/helper/inference.R:48). I independently recomputed the full-sample single-factor pricing errors, evaluated the criterion on both sides, and varied the differentiation step.

Command, from r_project/, with LANG=en_US.UTF-8 and SDF_MAX_WORKERS=4:

    Rscript /private/tmp/g1_r4_review/l1_checks.R

Results:

- The two L1 data CSVs contain 288 full-sample rows after excluding the H=1 alias. Exactly **250/288 market-coefficient fields** have SE<0.01 or |t|>50. The same criterion flags 268/288 “Coef” fields. Since the single-factor “Coef” field duplicates MKT, these are not independent counts: across 528 distinct model-parameter cells, 478 meet that numerical screen.
- In **40/48 single-factor cells**, the smallest absolute pricing error is below 1e-7. The largest of these smallest errors is **2.2508857916e-8**. The next-smallest absolute error ranges from **0.000587625491994 to 0.076283109463**.
- The loss is higher at both beta−0.01 and beta+0.01 in all 40 cases. The smallest increases are 8.45736251481322e-6 and 5.02586235354152e-6, respectively.
- The remaining eight single-factor SEs range from **0.505895878465115 to 11.5032237601079**.

The finite-difference dependence is decisive. For BO EW, H=0, beta=1.25225330121414:

| Differentiation step | Average numerical Hessian | Bartlett SE | Independent-units SE |
|---:|---:|---:|---:|
| 0.001 | 28.5611059052 | 0.00109480856606 | 0.00162132388996 |
| 0.0001 | 284.775855977 | 0.000109817121901 | 0.000162607865225 |
| 0.00001 | 2843.4854855 | 0.000011001186614 | 0.0000162852239377 |
| 0.000001 | 28086.7879544 | 0.00000111665110781 | 0.00000164871810302 |

Reducing the step tenfold makes the SE approximately ten times smaller. Changing the dependence estimator cannot repair this: both covariance estimators use the same defective curvature.

For the smooth VC FW H=60-month cell, the Hessian instead stays near 0.0110655 as the step changes from 0.001 to 0.0001, and its Bartlett SE stays near 11.50325. That numerical stability does **not** validate its LAD asymptotic interpretation.

The BO EW roots at H=0, 60, and 120 months are, respectively, **1.25225330121414, 1.252253323631, and 1.25225329410691**, with vintage 2010 closest to zero in all three. They agree to the optimizer’s numerical precision; they are not bit-identical.

### Why the diagnosis is right, and what it does not imply

Near a simple root, |e_v(theta)| behaves like |e'_v(theta_0)| |theta−theta_0|. Its centered second difference at the root grows as 2|e'_v|/delta. The inverse Hessian consequently shrinks with the arbitrary step. This is exactly the observed mechanism.

For regular LAD asymptotics, curvature belongs to the **population criterion**. Formally it includes both a smooth contribution involving sign(e_v) times the second derivative of e_v and a contribution from probability mass near zero involving the density at zero and products of first derivatives. Differentiating the observed absolute-loss sample twice, with a fixed numerical step, is not a justified estimator of that population curvature. Away from observed kinks it can simply omit the latter contribution.

Therefore **none of the current L1 sandwich SEs is established as usable**, including the eight smooth cells. Neither fixed-b critical values nor the independent-units sandwich rescues them.

T1 explains why a selected vintage root can recur across horizons. It does not imply that every L1 minimum must be a vintage root, or that the same root must remain optimal under different horizon weights.

## F2. Extreme lower MKT+EG minima

**Verdict: Verified by reproduction for the lower solutions; “254 are minima” and the economic diagnosis need qualification. Severity: high for the estimator’s definition and the affected estimates.**

### Original search and a broader independent search

I reran the archived 260-cell search. Resource auto-detection failed in this environment; a temporary copy changed only the worker count to 4. No objective, data, start, tolerance, or estimator setting was changed.

    Rscript /private/tmp/g1_r4_review/stall_check_4workers.R /private/tmp/g1_r4_review/stall.rds

All 260 rows matched their stored unit counts. Maximum relative objective-reproduction discrepancy: **4.33849892692e-15**.

Six rows have a lower solution by more than 1e-6 relative; these are three economically distinct cases, each repeated at H=1:

| Case, H=0 | Reported objective | Lower objective | Lower beta | Lower EG |
|---|---:|---:|---:|---:|
| BO FW | 0.0455868146129 | 0.0403469546722 | 3.29903865104 | 13.6494340758 |
| VC EW | 0.0530261771591 | 0.0458955950690 | 3.16644191850 | 13.6551535527 |
| VC FW | 0.0555857498172 | 0.0464461604770 | 3.19618673192 | 13.6486145346 |

Relative improvements are approximately 11.49%, 13.45%, and 16.44%. For the other 254 rows, the largest relative improvement is **2.0862156693e-7**.

The original random-start box is only a box for **initial values**. Nelder–Mead is unrestricted and demonstrably leaves it. Thus finding EG loadings near 13.65 from starts in [−4,4] is entirely consistent. Nevertheless the narrow initial box does not adequately probe every economically distinct basin.

I therefore ran a second search over all 260 cells, adding 16 starts per cell at 90% of the distance from (1,0) to the positivity boundary along evenly spaced directions. Alpha directions were scaled separately and respected its ±0.02 bounds. This used the actual estimation-period factor support; optimization remained the original unrestricted or alpha-bounded procedure.

    Rscript /private/tmp/g1_r4_review/wide_starts.R

It found **no additional improvement above 1e-6**. The largest improvement over the first search’s best value was **1.52582023145e-10**.

This is substantially stronger empirical evidence. It still is not a global-optimality proof. Appropriate wording is “no materially better solution found under the stated multistart checks,” rather than “all remaining points are global minima.”

### What the extreme solution exploits

I reconstructed each discounted cash-flow contribution and recovered the embedded terminal-NAV flows as the difference between the otherwise matched NAV100 and NAV0 pooled input files. This is preferable to treating all December 2021 cash flow as terminal NAV.

    Rscript /private/tmp/g1_r4_review/edge_checks.R
    Rscript /private/tmp/g1_r4_review/supplemental.R

The minimum monthly gross returns at the extreme points are:

- BO FW: **0.0173517186533**, January 2021.
- VC EW: **0.0169683471193**, January 2021.
- VC FW: **0.0174324843353**, January 2021.

A low gross benchmark return raises subsequent discounted values. The accompanying high growth in other months discounts much of the earlier/later stream down: over the actual estimation-period factor support, geometric annual growth factors are **2.42500736307, 2.55791509624, and 2.56904302482**. The corresponding reported points give **1.15927801699, 1.29334251701, and 1.34900269109**. These historical compounded benchmark paths are economically extreme.

The December 2021 hypothesis is too narrow:

| At the extreme point | BO FW | VC EW | VC FW |
|---|---:|---:|---:|
| RMS pricing error | 0.201734773361 | 0.238480507510 | 0.241341329288 |
| RMS discounted 2021 net cash flow | 0.0708002557317 | 0.0852786374976 | 0.0846227761827 |
| RMS all discounted terminal NAV | 0.0523588158453 | 0.0799174824948 | 0.0744297116833 |
| RMS discounted December 2021 net cash flow | 2.89332272918e-7 | 0.00306589904855 | 0.00262602487915 |

Terminal NAV matters, especially in younger vintages, but the year-end date itself is not the driver. Terminal NAV is embedded at funds’ own final reporting dates, not uniformly at the common data end.

For a transparent counterfactual, hold all other monthly gross returns at their extreme-point values and replace **only January 2021’s gross return with 1**. The sum of squared pricing errors changes:

- BO FW: **1.01742296958 → 1.13382161478**.
- VC EW: **1.53556971647 → 1.71640507034**.
- VC FW: **1.57263220501 → 1.74917928244**.

Holding those other months fixed, write each vintage error as a_v+c_v/G_Jan2021. The analytically optimal January return is **0.0176824915312, 0.0173693925131, and 0.0179094045126**, respectively: close to the actual extreme values. This directly supports a late-flow fitting mechanism; it is not just a visual inference from a low return.

But the extreme parameter vectors still fit better than the reported parameter vectors after removing **all** embedded terminal NAV from both:

| SSE after removing all terminal NAV, no reoptimization | Reported point | Extreme point |
|---|---:|---:|
| BO FW | 1.14552060237 | 1.09520120725 |
| VC EW | 2.41543962027 | 1.71690607312 |
| VC FW | 2.57232023969 | 1.73341089982 |

These are counterfactual evaluations at fixed parameter vectors, not newly estimated NAV0 minima. They rule out the simple explanation that the lower fit exists only because of terminal NAV.

**My reading:** the criterion rewards a combination of extreme cumulative growth and a strategically timed near-collapse in benchmark returns. The common factor path shapes the basin; actual cash-flow timing and sizes determine its benefit. Similar parameters across asset classes do not establish that cash flows are irrelevant.

## F3. Jackknife basin switch

**Verdict: Verified by reproduction. Severity: high for the reported robustness entry and the refitting protocol.**

I read [jackknife.se()](/Users/christausch/github/Fundwise_SDF/r_project/estim_model_optimized.R:580), reran the archived script, and then reran the calculation using the exact stored full-sample point instead of the rounded point hard-coded in that script.

    Rscript /private/tmp/g1_r4_review/jackknife_exact.R

Exact-start reproduction:

- Market SE: **2.633977635455374**, stored **2.6339776354554**; difference **−2.575717417130363e-14**.
- EG SE: **9.779749062129829**, stored **9.77974906213**; difference **−1.70530256582424e-13**.
- Deleting vintage 2003 gives **(2.86042348642455, 13.40247252022292)** under the pipeline’s default refit tolerances.
- Stored full-sample self-check: **0.0010853298818847**.

With tighter tolerances, deletion-2003 optimization from either the interior full-sample point or the extreme full-sample point reaches approximately **(2.860045, 13.402114)**, objective **0.0463040538163**. Thus this is not an obvious failed optimizer iterate.

The archived calculation excluding the one switching refit gives 0.1816 and 0.9439. Those are **not replacement jackknife SEs**; they illustrate the contribution of the switch.

The self-check asks whether refitting the **same full sample** moves the estimate. It says nothing decisive about a deletion sample’s competing basins. A small value cannot certify a jackknife.

“Contaminated” is reasonable relative to the intended interior-baseline interpretation. Under an actual unrestricted global-argmin estimator, a switch to a better basin could be genuine sample sensitivity. The deeper inconsistency is mixing a full-sample local solution with refits that sometimes select another basin.

## F4. Fixed-b numerical calibration

**Verdict: Verified by reproduction for the stored values and formula; the Monte Carlo explanation is compatible with the evidence, not a proof of exact published agreement. Severity: moderate for borderline markers; lower than the estimator and applicability issues.**

### Formula checks and exact stored numbers

I checked [fixedb.R](/Users/christausch/github/Fundwise_SDF/r_project/helper/fixedb.R:27), reconstructed the direct Bartlett quadratic form independently, and compared both scalar and bivariate versions.

    Rscript /private/tmp/g1_r4_review/fixedb_checks.R
    Rscript /private/tmp/g1_r4_review/fixedb_coupled.R

Stored two-sided 5% critical values:

- b=0.28: **2.776703970513**, versus the brief’s 2.777 and the reference’s 2.805.
- b=0.26: **2.729564199250**, versus 2.730 and 2.736.
- Actual VC b=7/27: **2.727215180881**. It must not be compared as if b were exactly 0.26.

For 30 Gaussian paths, T=1000 and M=280, the largest absolute difference between the partial-sum and direct-kernel denominators is **5.107025913276e-15**.

The independently computed exact denominator expectation is **0.746133**, agreeing with
1−M/T+(M²−1)/(3T²). The Brownian-limit expectation is **0.7461333333333**. This supports correct normalization.

For q=2, the covariance-matrix discrepancy is **5.551115123126e-16**; both representations yield Wald **9.58926614019**. The code simulates Wald W, not F=W/q, which is the appropriate scaling for its stored joint Wald statistic.

### Monte Carlo uncertainty and mesh error

Eight independent batches of 100,000 replications at b=0.28 gave critical values:

    2.796515327002  2.796782764027  2.774500466878  2.784687623374
    2.785970258108  2.776305492782  2.788012462717  2.762777794594

Their mean is **2.783194023685** and sample SD **0.01155804516783**. Eight batches estimate that SD imprecisely, but support the claimed order of Monte Carlo variation. This mean is a mean of quantile estimates, not a pooled 800,000-draw quantile.

Using the brief’s SD=0.009, the SD of the difference between an independent 100,000-draw quantile and 50,000-draw reference quantile is

    sqrt(0.009² + (sqrt(2) × 0.009)²) = 0.0155884572681.

The actual discrepancy **2.805−2.776703970513=0.028296029487** is approximately 1.82 combined SDs. That is compatible with Monte Carlo variation. It does not suggest a clear coding error.

I also coupled 200,000 paths on nested meshes: the T=1000 normals were pairwise aggregates of the T=2000 normals, preserving the Brownian path and numerator. The critical values are **2.788775057068** and **2.788963236004**, a difference of **0.0001881789359968**. This is much smaller than 0.028. It is evidence against mesh resolution explaining that gap in this experiment, not an exact bound on discretization error.

The claimed **16/864** near-threshold p-values is reproducible when the count comprises 288 p.MKT fields, 288 p.Coef fields, and 288 joint Wald fields: **10+2+4=16** lie in [0.045,0.055]. These are stored tests, not 864 distinct printed coefficient significance markers.

### Source-version qualification

**Checked in source:** the [January 2005 Cornell working paper](https://cae.economics.cornell.edu/05-08.pdf), Table I, printed p.24/PDF p.26, contains 2.736 at b=.26 and 2.805 at b=.28. Its simulation description states 50,000 replications and a 1,000-normal approximation.

The [March 2005 revised manuscript](https://kiefer.economics.cornell.edu/newasym2.pdf) instead reports cubic critical-value functions in Table I. Its Bartlett 97.5% polynomial gives **2.7508081376 at .26** and **2.8123591552 at .28**. That comparison also contains polynomial fit error. I did not verify the typeset journal’s Table I. The manuscript should identify which source/version supplies its numerical benchmark.

## D1. Recommended L1 decision

**Choose (a) now. Keep L1 point estimates as descriptive robustness, usually in an appendix; do not choose (d) merely because some estimates occur at kinks.**

Changes to the tables/paper:

1. Remove L1 SE_A, t, stars, inferential Wald p-values, and independent-units SEs wherever they might be interpreted as a valid alternative.
2. Remove the current L1 delete-one jackknife numbers from any table headed “inference robustness,” unless clearly presented as unvalidated deletion sensitivity rather than SEs.
3. Keep fold dispersion only as an explicitly descriptive stability measure. Show fold estimates/ranges and which vintage’s pricing root is selected. The existing signed stability ratio must carry no reference distribution or significance interpretation.
4. State that a valid L1 inferential procedure has not yet been established for this design.

A kink is routine for LAD: the sample median is observation-pinned and remains a legitimate estimator. T1 supplies a useful explanation of horizon plateaus, not a reason to discard the estimator.

**Checked in source:** [Efron (1979)](https://doi.org/10.1214/aos/1176344552), Section 3, demonstrates inconsistency of the usual delete-one jackknife for the sample median. This verifies the substance of the brief’s recollection; I did not reread the 1982 book. It is a counterexample to generic validity, not a theorem for this precise nonlinear, dependent-vintage estimator.

If L1 becomes scientifically central, (b) is the next development route I would investigate, with (c) as its theoretical benchmark. Dependence-preserving resampling is essential: iid resampling of vintage labels is not a remedy for the dependence motivating the HAC estimator. See the primary [Bühlmann–Künsch block-bootstrap report](https://www.research-collection.ethz.ch/server/api/core/bitstreams/6ee4a421-2a0b-4d08-bcf8-d2ed766f6b27/content), introduction. At 25–27 vintages, block/deletion-size sensitivity and Monte Carlo validation are necessary. No bootstrap label alone proves validity.

CV dispersion answers “how sensitive is this fit to these omitted blocks?” Folds overlap, omit different compositions, and have different effective information. It does not estimate a sampling SE by default, and division by sqrt(number of folds) would not fix it.

A further conceptual point: exact **mean** pricing of cash flows does not imply a zero **median/LAD score** at the same parameter. Before presenting L1 as an alternative consistent estimator of the same structural SDF parameter, establish its population target. Until then, “descriptive robustness” is the defensible interpretation.

## D2. Recommended admissibility and optimizer decision

**Choose (a) as the primary way to define a restricted estimator, with an economic rationale and broad sensitivity analysis. I do not recommend adopting 0.5 as an innocuous final cutoff on the evidence currently offered.**

### The 0.5 floor changes other results

I scanned the entire actual estimation-period factor support for all 260 reported Multi points.

    Rscript /private/tmp/g1_r4_review/supplemental.R

**23/260 reported points already violate 0.5**, all among the 254 points with no materially better solution found. All are VC; the binding historical month is October 1987. They comprise **15 MKT+Alpha, four MKT+ME, three MKT+EG, and one MKT+IA** points. None is an H=1 duplicate.

The smallest reported gross return is **0.394508405876584**. Every reported point is feasible at floors 0.05, 0.10, 0.25, 0.35, and 0.39; one fails 0.40; 23 fail 0.50. Feasibility does not establish that the constrained global solution remains the same.

Two direct constrained refits illustrate the effect:

| Otherwise confirmed point | Reported parameters | Parameters under G≥0.5 | Objective before → after |
|---|---|---|---|
| VC FW, MKT+ME, 5 years | (2.41179405232, 0.703211840272) | (2.03529186850, 0.458710391966) | 0.170494258432 → 0.173716788425 |
| VC FW, MKT+Alpha, 22.5 years | (2.53962394325, −0.011961581349) | (2.13279785339, −0.0115294897132) | 0.211604538774 → 0.217244615402 |

Both new points lie on the return-floor boundary. They are illustrative refits, not a completed constrained G1 rerun or a global-optimality certificate.

Conversely, I performed 13-start constrained searches for each of the three problematic EG H=0 cases. At 0.5, the best points return to the reported interior basin, up to small optimizer refinements:

- BO FW: (0.983211721642, 0.452896345378), Q=0.0455868144193.
- VC EW: (0.107714317025, 2.72807978515), Q=0.0530261767296.
- VC FW: (0.129643789787, 3.30510749351), Q=0.0555857494101.

But at floor **0.25**, competing boundary points still beat the reported objectives by **1.7850%, 1.8080%, and 5.4838%**, respectively. At 0.05 and 0.10 the gains are larger. Thus the restriction’s strength changes the selected estimator materially. A sensitivity table must show parameter paths, objective gaps, and active constraints; merely listing several floor values is insufficient.

### Economic interpretation and implementation requirements

For this reciprocal benchmark model, the one-period discount multiplier is m=1/G. A lower bound G≥c both guarantees positivity and places the upper bound m≤1/c. With c=.5, m≤2.

**Checked in source:** [Hansen–Jagannathan (1991)](https://larspeterhansen.org/wp-content/uploads/2016/11/ImplicationsofSecurityMarketData.pdf), Restriction 2, p.229 and Section IV, concerns positive/nonnegative pricing kernels. It does not prescribe a .5 monthly benchmark floor. The .017 solutions are already positive. Excluding them requires a **boundedness/economic admissibility restriction**, not positivity alone.

My preferred specification would make these choices explicit:

- An economically defended admissible set, fixed before interpreting significance. A return floor is more interpretable than separate arbitrary coordinate boxes, but cannot be justified merely by restoring preferred estimates.
- A fixed factor-date domain across the full sample, CV fits, and jackknife refits. Use estimation-period months for the baseline. A 1969–2025 stress history is a separate restriction using information beyond the cash-flow end; label it accordingly.
- A deterministic multistart rule, objective tolerance, and tie-breaking rule for **every** refit, including deletions and CV. A warm start can accelerate this but should not silently change which estimator the refit computes.
- Per-refit records of objective value, convergence, minimum G and its month, active constraints, and basin switches. Constraining and detecting switches are complementary.
- Suppression or replacement of ordinary interior SEs/tests whenever the new return constraint binds. Existing alpha-bound flags do not cover a binding return floor.
- Diagnostics of cumulative growth and reciprocal discount paths as well as min G. A lower floor alone does not rule out discount-away solutions produced by excessive growth in other months.

Option (b) can be useful as economically specified leverage/exposure regularization, rather than a substitute chosen solely for numerical convenience. Option (e) needs scaling and tuning rules: its lambda is not comparable across V or criterion normalizations without care. Option (d) is valuable as a model sensitivity check, but positive-by-construction does not itself rule out very small G or enormous cumulative growth.

I disagree with the statement that (c) is necessarily “not a well-defined estimator.” A completely specified deterministic fixed-start algorithm defines a statistic. The substantive objections are arbitrary basin dependence, a mismatch with a claimed global argmin, and inconsistent full-sample/refit selection. I would not choose it as the paper’s primary economic specification.

## D3. Recommended fixed-b precision/provenance decision

**Choose (a), but treat it as calibration and reproducibility work, not validation of inference.**

One million replications is reasonable, preferably concentrating first on b=7/25, b=7/27 and the restriction dimensions that the paper actually uses. That avoids interpolation noise at the principal empirical b values. Retain the broader grid for CV and other samples.

The larger scientific task is to establish size and coverage for the **complete estimator** at 25–27 vintages. Accurate Gaussian-reference quantiles cannot compensate for bias, heterogeneous vintage precision, weak information, active restrictions, or basin selection.

**Checked in source:** [Lazarus–Lewis–Stock–Watson (2018)](https://stock.scholars.harvard.edu/sites/g/files/omnuum5911/files/stock/files/lazarus_lewis_stock_watson_har_inference_recommendations_for_practice_jbes_2018.pdf), Eq.(3) and its rounding footnote, supports the chosen bandwidth rule and fixed-b calibration. [Kiefer–Vogelsang, March manuscript](https://kiefer.economics.cornell.edu/newasym2.pdf), Assumptions 2–3, requires the appropriate score functional CLT and partial Jacobians growing linearly in the sample fraction. Nonlinear estimating equations are within its scope; homogeneous partial-sum behavior in these vintages is the substantive question.

**Yes, the fixed-b table’s hash belongs in the inference identity.** The current result identity omits it, while cached joint p-values depend on it. A version bump is a safe short-term invalidation mechanism. The better design separates:

- Point-estimation identity: data, model, admissibility domain, objective, optimizer/selection rule.
- Covariance identity: point estimate, score/Hessian calculation, kernel/bandwidth, vintage geometry.
- Calibration identity: covariance/test statistic, restriction dimension, table content hash and simulation metadata.

A pure quantile update should recompute joint and marginal p-values from stored statistics without forcing optimization. Restart or invalidate the process-level fixed-b table cache as part of that workflow. Store seeds, RNG, T, replication count, quantile method, kernel convention and table source/version.

More draws do not automatically confer fine extreme-tail precision: at probability 1e-4, 100,000 draws provide about ten tail observations, and one million about 100. Report very small values as tail bounds or with appropriate precision.

## Additional points bearing on these decisions

### 1. The present N_eff is not parameter-specific

This matters because concentration is being used to qualify inference.

[inference.R](/Users/christausch/github/Fundwise_SDF/r_project/helper/inference.R:92) always uses the **first parameter’s raw score**, while the robustness tables often display the SE of EG, ME or Alpha beside that same N_eff.

I recomputed both raw-score and linearized-influence concentration:

| Cell | Stored N_eff (MKT score) | Second-parameter score concentration | Second-parameter influence concentration |
|---|---:|---:|---:|
| VC FW EG, H=0 | 7.57983927133 | 5.10582573580 | 5.87844769446 |
| VC EW ME, 20 years | 1.39577326122 | 2.00172686547 | 2.00277004844 |
| BO EW Alpha, 20 years | 4.99084555622 | 4.35871735562 | 6.27029269876 |

For a regular multivariate estimator, the coefficient-specific first-order contribution is proportional to the corresponding component of H^{-1}g_v, not generally to that parameter’s raw score. Its squared-share concentration is a more informative parameter-specific diagnostic. With dependence, even that is **not** degrees of freedom or an actual independent sample size. Label the existing quantity honestly or report coefficient-specific influence concentration.

### 2. Calendar gaps are another approximation, separate from Monte Carlo precision

The kernel uses calendar-year distance; the standard reference simulation uses equally spaced observations. BO lacks 1989; VC has 1982 followed by 1985. Replacing the sample span with V in b does not reproduce the entire observed kernel matrix.

As an intentionally narrow diagnostic, I simulated one million iid Gaussian **location-model** samples using either the contiguous or actual vintage-year Bartlett matrix, S=7, with the same random draws:

| Location-model design | Two-sided 5% critical value |
|---|---:|
| BO, 25 contiguous years | 2.81084468712 |
| BO, actual 25 vintage years | 2.77650180883 |
| VC, 27 contiguous years | 2.73946339050 |
| VC, actual 27 vintage years | 2.71226081218 |

The paired differences are roughly 0.0343 and 0.0272, comparable to the disputed calibration gap. This does **not** demonstrate that the SDF tests overreject; for BO the actual-calendar Gaussian critical value is fortuitously close to the stored reference. Nor is a location-model correction automatically appropriate for the nonlinear SDF estimator. It shows why exact vintage geometry belongs in the planned size study and why an ad hoc switch to b=S/calendar-span is not a complete solution.

### 3. Separate inference calibration from the criterion’s population target

The paper’s own scale/noise argument implies that the L2 population minimizer with finite pooling noise can differ from the true pricing parameter. If Step 4 records rejection of the true pricing parameter in that design, it combines criterion bias with inferential size distortion.

I recommend distinct checks:

1. A simple Gaussian/kernel benchmark, to isolate numerical calibration and vintage geometry.
2. Cash-flow designs where the criterion target is known or estimated to high precision, evaluating coverage of that target.
3. Economic-pricing designs evaluating bias, RMSE, coverage and rejection at the structural true parameter, explicitly identifying the contribution of target displacement.
4. The final constrained/multistart estimator, including boundary frequency, basin switching, failure frequency and the exact reported tests.

Use actual V, vintage gaps, fund-count/size profiles and shared calendar exposure. State whether inference is conditional on the realized factor path or also averages over new common factor histories. Resampling vintages alone cannot supply the latter.

### Suggested order of work

1. Remove unsupported L1 inference and specify the restricted estimator and its sampling target.
2. Document the admissible-set sensitivity, including the 23 already infeasible VC points, and apply one refitting/selection protocol everywhere.
3. Run the estimator-level size/coverage study with active boundaries and basin switches counted.
4. Increase fixed-b precision, record calibration hashes, and regenerate all dependent tables consistently.

The paper can become substantially stronger by making the criterion’s target, its admissible domain, and its uncertainty assessment explicit. Resolving these is more consequential than defending individual stars in the current tables.

## Reproduction record

All R commands ran from /Users/christausch/github/Fundwise_SDF/r_project/ with LANG=en_US.UTF-8 and SDF_MAX_WORKERS=4. Heavy jobs were run sequentially; the multistart jobs used at most four workers. R version 4.3.1; Rcpp 1.0.11; data.table 1.14.8.

Scripts, logs and derived diagnostic CSVs are in [the experiment directory](/private/tmp/g1_r4_review/). The main commands are listed in the corresponding sections. The additional concentration/calendar experiment is:

    Rscript /private/tmp/g1_r4_review/extra_inference.R

The original full G1 pipeline was not rerun. The two commits’ entire 60,736-row invariance claim was outside the requested recheck except insofar as the table/cell reproductions verify the rows used here. The NATRES extension was read in the archived script and plan; I did not independently rerun its 200-start search or its older CV fold. Those subsidiary NATRES numerical claims remain plausible but not reproduced in this review.

The separate [primary-source notes](/Users/christausch/.codex/visualizations/2026/11/17/01a1e97c-61a1-7f23-80a0-a1856988cc50/g1_r4_review/primary_sources.md) distinguish full-text checks from abstract-only checks. No decisive recommendation above relies on a citation recalled only from memory.

Final git status matched the initial status: only the user-provided, untracked review brief was listed.
