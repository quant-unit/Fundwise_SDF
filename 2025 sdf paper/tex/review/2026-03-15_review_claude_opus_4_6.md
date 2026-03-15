

# Referee Report

**Paper:** "Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows"
**Authors:** Christian Tausch and Alexander Bohnert
**Journal:** Journal of Financial and Quantitative Analysis

---

## Overall Assessment

This paper develops a semiparametric framework for estimating stochastic discount factors from pooled, non-traded cash flows, generalizing the inception-based approach of Driessen, Lin, and Phalippou (DLP, 2012) and clarifying its relation to Korteweg and Nagel (KN, 2016). The central innovation is "horizon averaging"—evaluating pricing errors at multiple compounding dates beyond fund inception—which the authors show can dramatically improve finite-sample performance even though it shifts the population pricing restriction. The paper applies the framework to Preqin buyout and venture capital data and concludes that only single-factor market models are stably identified in currently available panels.

The paper addresses a genuinely important problem and contains several valuable insights. The "Timing Risk Premium" concept is novel and the honest assessment of identification limits in PE data is a useful contribution to a literature that sometimes overstates the precision of its estimates. However, the paper suffers from excessive length, an asymptotic theory that is largely an application of known results to a specific setting, an empirical application whose main message is essentially negative, and several conceptual and expositional issues that limit its current suitability for JFQA. I recommend a **revise and resubmit** contingent on substantial revisions.

---

## Major Comments

### 1. Length and Focus

At 67 pages (including appendices), the paper is far too long for the novelty delivered. The methodology section (§3) could be compressed considerably. Much of §3.3–3.4 is standard extremum estimator theory (Newey & McFadden 1994) applied to a specific criterion, and the proofs in Appendix C, while correct, are largely mechanical verifications of textbook conditions. The paper would benefit from moving the asymptotic theory to an online appendix and focusing the main text on (i) the timing risk premium insight, (ii) the LMD formulation and its relation to DLP/KN, (iii) the simulation-based guidance on horizon choice, and (iv) the empirical application. I would target approximately 40–45 pages inclusive of tables and figures.

### 2. The Core Tradeoff Needs Sharper Characterization

The paper's central tension—that horizon averaging improves finite-sample performance but shifts the population target—is clearly stated but insufficiently resolved. Several issues:

- **The optimal horizon is entirely simulation-dependent.** The paper acknowledges this (p. 34: "determining the optimal horizon requires tailored simulation designs rather than a simple heuristic") but offers no practical guidance beyond "simulate your specific setting." For a methodology paper, this is unsatisfying. Can the authors provide even rough analytical guidance—e.g., conditions under which the MSE-minimizing horizon is bounded above or below, or a plug-in rule based on observable features of the cash-flow panel?

- **The MSE decomposition in Appendix A (Proposition 7) is exact but uninformative.** It separates four terms but says nothing about their relative magnitudes without further assumptions. The "heuristic link to second-order finite-sample bias" (p. 56) is left as a conjecture. Given that this is the paper's central methodological insight, a more rigorous treatment is needed—even if only for a simplified DGP.

- **The magnitude of the target shift is never quantified empirically.** The authors could estimate the TRP directly (or bound it) in their simulation DGP where the truth is known, which would make the bias-variance tradeoff concrete. How large is $\|\theta_H^\dagger - \theta_0\|$ at H = 15 years relative to the finite-sample bias reduction?

### 3. Relation to Korteweg and Nagel (2016)

The comparison to KN (§3.5.2) is somewhat superficial. Proposition 6 establishes the variance-decomposition distinction between loss-before-aggregation (LMD) and aggregation-before-loss (KN-style), but the economic implications are underexplored:

- When does the cross-sectional dispersion term matter quantitatively? Under what conditions would the two approaches yield materially different estimates?
- The claim that KN "is best viewed as a closely related aggregated-moment GMM analogue rather than as a special case" deserves more development. The reader would benefit from a simulation comparing the two estimators on the same DGP to understand when the distinction is economically meaningful.
- The paper notes that KN's empirical implementation identifies SDF parameters from benchmark moments alone, with PE pricing errors evaluated post-estimation. This is an important distinction that deserves more discussion—it means the two approaches answer somewhat different questions.

### 4. Asymptotic Framework: Standard Theory, Questionable Finite-Sample Relevance

The increasing-domain asymptotic framework (§3.3) is carefully stated but raises two concerns:

- **Practical relevance.** The authors themselves repeatedly note that asymptotic approximations are unreliable in their sample (e.g., asymptotic standard errors of 0.000 or >25 in Tables 8–9). If the asymptotics are so unreliable that the authors must rely primarily on cross-validation diagnostics, what is the value-added of the formal theory? The paper needs to either (a) demonstrate settings where the asymptotics work reasonably or (b) reframe the contribution as primarily about finite-sample methods, with asymptotics serving only as a benchmark.

- **Verification of primitive assumptions.** Assumptions 5 (uniform NED, stochastic equicontinuity) and 7 (local Hessian convergence) are stated as high-level primitives but never verified for the specific SDF specifications used in the empirical application. For a paper claiming to provide rigorous inference, some discussion of whether these conditions are plausible for the linear and exponential-affine SDFs applied to PE cash flows would be appropriate.

### 5. Empirical Application: Largely Negative Results

The empirical contribution is honest but essentially diagnostic rather than constructive. The main findings are:

- Single-factor MKT models are stable; multi-factor models are not.
- Alpha is weakly identified.
- Recent vintages with large NAVs inflate beta estimates.

While useful, these findings are not entirely surprising given the well-known data limitations in PE (sparse cross-sections, overlapping vintages, NAV contamination). The paper would be strengthened by:

- **A more substantive economic application.** For instance, using the stable single-factor estimates to compute risk-adjusted performance measures (PMEs, alphas) with proper uncertainty quantification, and comparing these to the existing literature.
- **Out-of-sample prediction.** The hv-block cross-validation is used as a stability diagnostic, but could it be used to formally compare models via out-of-sample pricing error criteria?
- **A broader data application.** The restriction to North American USD funds from Preqin limits generalizability. How do results change with other geographies or data providers?

### 6. Simulation Design and the DLP Comparison

The inability to replicate DLP's favorable simulation properties (Appendix B.1) is concerning and potentially important for the literature. However, the presentation is unsatisfying:

- The authors systematically vary features (error distribution, holding period, alpha bounds) but never identify the single dominant source of the discrepancy. Is it the error distribution? The holding period? The alpha bound? A clear decomposition would be valuable.
- Table 10 shows that even the closest approximation to DLP's design (S10) produces substantially more bias. The authors speculate about causes but do not resolve the discrepancy. This should either be resolved or the implications for the DLP literature should be stated more forcefully.

---

## Minor Comments

1. **Notation.** The paper introduces a large amount of notation (Tables 1–2), some of which is used only briefly. Streamlining would improve readability. For example, the deal-level notation ($\text{Inv}_{i,j}$, $\text{Div}_{i,j,k}$, $t^{\text{Inv}}_{i,j}$, $t^{\text{Div}}_{i,j,k}$, $K_{i,j}$) is introduced in Assumption 1 but serves only as a theoretical device; it could be condensed.

2. **Equation (25) vs. (26).** The simulation study finds that the exponential-affine and simple-linear SDFs perform similarly (Figure 9). This is an interesting negative result but the economic intuition is underdeveloped. Why doesn't the theoretical advantage of exponential-affine SDFs (guaranteed positivity of discount factors; Gourieroux and Monfort 2007) translate into finite-sample gains?

3. **Alpha bounds (Remark 6).** Bounding alpha at $\pm 2\%$/month is a consequential implementation choice. The sensitivity of results to this bound (partially explored in Appendix B.1) deserves more attention in the main text. A bound of $\pm 2\%$/month corresponds to $\pm 24\%$/year, which is very generous—yet even this does not prevent the "exploding alpha" problem in some specifications.

4. **NAV treatment (§4.3–4.4).** The systematic NAV discount analysis (Figure 5) is a valuable robustness exercise. However, the decision to treat the latest NAV as a final cash flow for non-liquidated funds is a strong assumption that deserves more discussion. The 50% NAV discount in Figure 4 is ad hoc—why 50%? Is there any empirical basis for this choice?

5. **Cross-validation implementation (Table 5).** The choice of h = 3 and v = 3 vintage years for the hv-block cross-validation is not formally justified. How sensitive are the results to these choices? The Racine (2000) framework provides conditions for consistency but not for optimal block size selection in finite samples.

6. **Figure quality.** Several figures (especially Figures 2–3, 12–15) contain a large number of overlapping lines that are difficult to distinguish. Consider using fewer lines with clearer differentiation, or using small-multiple panels.

7. **Table formatting.** Tables 8–9 contain many asymptotic standard errors that are clearly unreliable (0.000, 25.778, 40.420, 88.860). Rather than reporting these as if they were meaningful, consider flagging or suppressing them with a note explaining the numerical instability.

8. **p. 4:** "Horizon choice" is capitalized inconsistently throughout the paper (sometimes "Horizon," sometimes "horizon"). Choose one convention and apply it consistently.

9. **p. 11, Equation (12):** The role of $w_i$ is somewhat unclear. Is it a normalization factor or a weight? The text says "it is not an observation weight in the outer sample average," but its practical effect is similar. Clarify.

10. **p. 22, Proposition 6:** The variance decomposition is elementary and well-known. It would suffice to state it as a remark rather than a formal proposition.

11. **p. 26, Remark 7:** The definition of the Horizon parameter $H$ in years and its mapping to months via $12H$ is potentially confusing. State this mapping once and refer back.

12. **Sections 5.2–5.3:** The detailed horizon-by-horizon tables for two-factor models are extensive but repetitive. The main patterns (destabilization with alpha, sign reversals in q5 factors) could be conveyed more efficiently with summary statistics and fewer tables.

13. **Proofs (Appendix C).** These are correct but routine. The proof of Proposition 2 essentially says "apply Jenish and Prucha (2009, Theorem 6) to the transformed array," and the proofs of Propositions 3–4 follow Newey and McFadden (1994) step by step. Consider relegating these to an online appendix.

14. **Missing references.** The paper cites "Appendix ??" on pp. 11 and 35, suggesting unresolved cross-references. These should be fixed.

15. **Data description (§4.1).** The regional distribution of funds includes entries inconsistent with the North America filter (e.g., Asia, Europe, Japan for BO; China, Thailand, Australasia for VC). The paper mentions this but does not explain why these funds pass the filter. This should be cleaned up or explained more clearly.

16. **The q5 factors.** The choice to use q5 factors (Hou et al. 2015, 2021) rather than the more standard Fama-French five factors deserves justification. How different would results be with FF5?

17. **Remark 3 (p. 15):** The distinction between $\theta_H^\dagger$ and $\theta_0$ is important but could be made more intuitive with a simple example—e.g., what happens in a one-deal, one-period model where the TRP can be computed in closed form?

---

## Summary Recommendation

The paper addresses an important gap in the private equity asset pricing literature and contains several valuable contributions: the timing risk premium concept, the unified LMD framework, the systematic simulation study, and the honest assessment of identification limits. However, it needs substantial revision to reach the bar for JFQA publication:

1. **Shorten substantially** by compressing standard asymptotic material and relegating proofs to an online appendix.
2. **Provide sharper characterization** of the bias-variance tradeoff, ideally with analytical results for a simplified case.
3. **Strengthen the empirical application** with a more constructive economic contribution beyond stability diagnostics.
4. **Resolve or clearly state the implications** of the inability to replicate DLP's simulation results.
5. **Clean up the presentation** (unresolved references, inconsistent capitalization, overloaded figures, unreliable standard errors).

**Recommendation: Revise and Resubmit**

---



# Reframing the Paper as an Identification and Stability Contribution

## The Core Problem with the Current Framing

Your paper currently sits uncomfortably between three identities: a new estimator paper, an asymptotic theory paper, and an empirical application paper. The estimator generalizes DLP but the generalization is relatively straightforward. The asymptotic theory applies known results to your setting. And the empirical results are largely negative—multi-factor models fail, alpha is unidentified, and you conclude that current data support only a single market factor. Read as a factor-discovery paper, this feels like a null result. Read as an identification paper, it is your main finding.

The reframing I would suggest treats every element of the paper—the theory, the simulations, and the empirics—as evidence bearing on a single question: **What can semiparametric SDF estimation from pooled PE cash flows actually identify, and what can it not?**

## Concrete Suggestions

### 1. Restructure the Introduction Around an Identification Question

Your current introduction opens with the growth of PE as an asset class and proceeds to methodology. Instead, open with the identification problem directly. Something along these lines:

> A growing literature estimates stochastic discount factor parameters from pooled private-equity cash flows. These estimates are used to measure abnormal performance, assess systematic risk exposure, and evaluate factor structure. But how much of what these estimators report reflects genuine economic content, and how much reflects finite-sample instability, horizon choice, and NAV contamination? This paper provides a systematic answer.

This immediately signals that the contribution is diagnostic. The reader no longer expects you to deliver a new factor model or a definitive alpha estimate. Instead, they expect you to map the boundary of what is identifiable—and that is exactly what you do well.

### 2. Elevate the Timing Risk Premium from a Technical Byproduct to a Central Identification Result

Currently the TRP appears in Section 3.1 as a theoretical observation about discounting dates, and its implications trickle through the rest of the paper somewhat loosely. Reframe it as your first identification result:

**Result 1:** The choice of discounting date is not an innocuous implementation detail. It changes the population object being estimated. Any researcher using these methods must choose between (a) structural fidelity at the cost of severe finite-sample instability, and (b) finite-sample regularization at the cost of estimating a different population quantity.

This is a statement about identification, not about estimation technique. It tells the applied researcher that there is no free lunch in horizon choice—and that the DLP framework, by restricting to inception-only discounting, pays a high variance cost that is not visible in standard asymptotic standard errors.

You already have all the ingredients for this (Lemma 1, Proposition 1, the simulation evidence). What is missing is the explicit interpretive frame that says: this is a fundamental identification tradeoff, not a tuning-parameter choice.

### 3. Recast the Simulation Study as a Finite-Sample Identification Analysis

Your simulations currently read as a calibration exercise for implementation choices. Reframe them as a systematic investigation of identification strength. For each dimension you vary (portfolio formation, vintage span, number of funds, SDF specification, number of factors, horizon), state the finding as an identification result:

- **Portfolio formation** identifies market exposure more sharply than individual funds because it averages out the idiosyncratic noise that contaminates the nonlinear criterion. This is an identification result: the effective signal-to-noise ratio for the market factor is too low at the individual-fund level to support stable estimation.

- **Adding alpha destroys market-factor identification** even in correctly specified simulations with zero true alpha. This is not a statement about alpha being economically unimportant—it is a statement about the information content of 20 vintage-year portfolios being insufficient to separately identify two parameters in a nonlinear, multi-period pricing restriction.

- **Two-factor models are fragile** not because the second factor is economically irrelevant, but because the cross-section of vintage-year portfolios does not contain enough independent variation to separate correlated factor exposures. The sign reversals across horizons are a symptom of weak identification, not evidence against multi-factor pricing.

Each of these is currently in your paper, but buried in simulation commentary. Pull them out as explicit propositions about identification strength—informal propositions backed by simulation evidence rather than formal theorems, but stated with the same clarity.

### 4. Make the Empirical Section an Identification Map

Your empirical results currently proceed model by model (single factor, then alpha, then q5 factors). This reads as a search for the right specification. Instead, organize it as a systematic mapping of identification strength across the specification space:

**What is identified:**
- Broad market exposure in seasoned BO and VC portfolios. Stable across horizons, weighting schemes, and cross-validation folds. Market beta of roughly 1.0–1.3 for BO and 1.3–2.0 for VC depending on horizon.

**What is weakly identified:**
- The level of market beta (because it varies with horizon choice, and the horizon introduces a TRP wedge whose magnitude is not separately estimable).
- Size exposure for BO at short horizons (positive ME loading with moderate stability ratios, but decaying and reversing at longer horizons).

**What is not identified:**
- Alpha. The intercept and market loading are not separately identified in the available cross-section. Any alpha estimate from these data should be interpreted as a joint test of alpha and model specification, not as a clean measure of abnormal performance.
- Multi-factor structure. None of the q5 factors produces a stable loading across the horizon grid. This does not mean PE has no exposure to these factors—it means the data cannot tell us.

This organization turns your negative results into positive contributions. The reader walks away with a clear understanding of what current PE cash-flow panels can and cannot support.

### 5. Add an Explicit Identification Taxonomy

Consider adding a table or figure that summarizes your identification findings in a structured way. Something like:

| Object | Identified? | Evidence | Caveat |
|---|---|---|---|
| Market β (BO, seasoned) | Yes | Stable across H, EW/FW, CV folds | Level depends on H choice |
| Market β (VC, seasoned) | Yes | Higher than BO, stable | More sensitive to weighting |
| Market β (recent vintages) | Confounded | Inflated by NAV co-movement | Mechanical, not economic |
| Alpha (BO or VC) | No | Absorbs MKT variation, hits bounds | Weak separation from β |
| Second factor (any q5) | No | Sign reversals across H, unstable | Collinearity in small cross-section |
| Factor structure (>1 factor) | No | All two-factor models fragile | Information content too low |

This table would crystallize your paper's contribution in a way that is immediately useful to applied researchers.

### 6. Compress the Asymptotic Theory and Reposition It

The asymptotic framework is currently presented as if it were a primary contribution. In an identification paper, it serves a supporting role: it provides the benchmark against which finite-sample distortions are measured. Consider:

- Moving all proofs to an online appendix.
- Compressing Sections 3.3–3.4 into a single subsection that states the main results (consistency, asymptotic normality, SHAC inference) with references to the relevant theorems in Jenish and Prucha (2009, 2012), Newey and McFadden (1994), and Pötscher and Prucha (1997).
- Adding a short paragraph explaining why the asymptotic theory matters for interpretation even when it does not deliver reliable finite-sample inference: it defines the population objects, makes the target-shift distinction (θ₀ vs. θ†_H) precise, and motivates the SHAC covariance estimator.

### 7. Sharpen the Implications for the PE Performance Literature

Your conclusion currently states that your results "do not imply that private-equity alpha is zero, nor do they invalidate prior estimates of outperformance." This is appropriately cautious, but it could be sharpened into a more useful statement. Consider something like:

> Our findings imply that any semiparametric alpha estimate from pooled PE cash flows should be accompanied by (i) a specification of the discounting horizon and its implied target shift, (ii) evidence of stability across horizon grids and cross-validation folds, and (iii) an explicit assessment of whether the effective cross-section supports the dimensionality of the model. Studies that report alpha estimates without these diagnostics may be reporting identification artifacts rather than economic quantities.

This reframes your paper not as criticizing prior work, but as providing the diagnostic toolkit that should accompany any application of these methods.

### 8. Retitle the Paper

The current title, "Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows," signals a methodology paper. A title oriented toward identification would better set expectations. Possibilities:

- "What Can Pooled Cash Flows Identify? Horizon Choice, Finite-Sample Stability, and Factor Structure in Private Equity"
- "Identification Limits of Semiparametric SDF Estimation for Private Equity"
- "The Identification Frontier for Private-Equity Risk and Return"

Each of these signals that the paper's contribution is about what the data can tell us, not about a new estimation technique.

### 9. Rewrite the Abstract

Your abstract currently emphasizes the framework, the generalization of DLP, and the accommodation of dependence. In the reframed version, lead with the identification question and findings:

> We study what semiparametric stochastic discount factor estimation from pooled private-equity cash flows can identify in currently available data. We show that the choice of discounting horizon creates a fundamental tradeoff: inception-only pricing preserves the structural pricing restriction but produces highly unstable estimates, while horizon averaging stabilizes estimation but shifts the population target through a timing risk premium. In simulations calibrated to PE-like environments, single-factor market models are reliably identified, but alpha and multi-factor specifications are not—even when the true model is correctly specified. In Preqin buyout and venture capital data, market exposure is stably estimated for seasoned vintages, but abnormal performance and richer factor decompositions are weakly identified and highly sensitive to implementation choices. Our results imply that current PE cash-flow panels support inference about coarse market exposure but not about alpha or multi-factor structure under semiparametric methods.

This abstract tells the reader exactly what they will learn and why it matters—without requiring them to care about the LMD estimator per se.

---

## What This Reframing Preserves

Everything substantive in the current paper survives. The LMD framework remains as the estimation vehicle. The TRP remains as the theoretical insight. The simulations remain as the finite-sample evidence. The empirical application remains as the real-data assessment. What changes is the interpretive lens: instead of "here is a new estimator and what it finds," the paper says "here is what these methods can and cannot identify, and here is why."

This reframing also addresses the referee concern (from my earlier report) that the empirical results are "largely negative." In an identification paper, showing that something is not identified is a positive finding. The applied researcher who learns that alpha is not separately identified from market exposure in these data has learned something valuable—potentially more valuable than a specific alpha estimate that happens to be statistically significant in one sample.