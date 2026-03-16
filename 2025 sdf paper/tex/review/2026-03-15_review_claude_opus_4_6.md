

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

How to do this best in a minimal-invasive approach?

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

---

WHICH JOURNAL TO SUBMIT ????

This is an important question, and I want to give you an honest, strategically useful answer rather than a polite one.

Let me assess fit along several dimensions: scope, methodology balance, contribution type, competitive positioning, and realistic acceptance probability.

---

# How I classify the paper

Before recommending journals, it helps to be precise about what the paper is.

The paper sits at the intersection of three areas:

**Econometric methodology**
Semiparametric estimation, increasing-domain asymptotics, spatial dependence, nonlinear extremum estimators.

**Private-market asset pricing**
SDF estimation from non-traded cash flows, PE performance measurement, NAV contamination.

**Finite-sample identification and stability diagnostics**
Simulation-based evidence on bias-variance tradeoffs, cross-validation under dependence, model selection.

The paper is not primarily any one of these. It is a methods-for-applications paper with a strong simulation component and a diagnostic empirical exercise. That profile matters a lot for journal fit.

---

# Where the key reference papers were published

This anchors realistic expectations.

| Paper | Journal |
|-------|---------|
| Driessen, Lin, Phalippou (2012) | **JFQA** |
| Korteweg and Nagel (2016) | **JF** |
| Ang, Chen, Goetzmann, Phalippou (2018) | **JF** |
| Boyer, Nadauld, Vorkink, Weisbach (2023) | **JF** |
| Brown, Ghysels, Gredil (2023) | **RFS** |
| Franzoni, Nowak, Phalippou (2012) | **JF** |
| Harris, Jenkinson, Kaplan (2014) | **JF** |
| Cordell, Roberts, Schwert (2023) | **JF** |
| Jeffers, Lyu, Posenau (2024) | **JFE** |
| Li, Xu, Zhang (2016) | **JFQA** |
| Farnsworth, Ferson, Jackson, Todd (2002) | **JBus** |
| Gourieroux and Monfort (2007) | **JFnE** |
| Bertholon, Monfort, Pegoraro (2008) | **JFnE** |
| Kim and Sun (2011) | **J. Econometrics** |
| Conley (1999) | **J. Econometrics** |
| Jenish and Prucha (2009, 2012) | **J. Econometrics** |

Two observations:

First, the direct predecessor — Driessen et al. (2012) — was published in JFQA, not a top-3 journal. That is actually a useful signal: JFQA is a natural home for this type of contribution.

Second, the papers that reached JF or JFE on PE asset pricing typically either introduced a fundamentally new framework or delivered major new empirical findings. Your paper's contribution is more incremental and diagnostic, which is not a criticism but an honest assessment of where it sits.

---

# Journal-by-journal assessment

## Tier 1: Natural fits

### 1. Journal of Financial Econometrics (JFnE)

**This is probably the single best fit for the paper as currently written.**

Why:
- JFnE explicitly publishes "development of novel tools and techniques, and their applications to financial problems."
- The paper's mix of semiparametric estimation theory, spatial HAC inference, simulation evidence, and financial application is exactly what JFnE values.
- Reference papers like Gourieroux and Monfort (2007) and Bertholon et al. (2008) — which you cite — were published here.
- JFnE readers are comfortable with notation-heavy estimation papers.
- The asymptotic theory, which is a liability at more empirically oriented journals, is a strength here.
- The finite-sample simulation study is valued as a contribution in its own right at JFnE.
- JFnE does not require a big new empirical finding; methodological novelty with an application is enough.

Risks:
- JFnE will hold you to a high econometric standard on the asymptotics.
- A referee may ask whether the increasing-domain CLT adds genuinely new econometric content beyond applying Jenish and Prucha to a new setting.
- The PE application may be viewed as specialized unless you frame the methodology as applicable to other non-traded asset classes.

Probability of a favorable outcome: **Moderate to good**, perhaps 25-35% conditional on a clean submission. This is competitive but realistic for JFnE.

### 2. JFQA

**JFQA is a strong natural fit, especially given that Driessen et al. (2012) was published there.**

Why:
- Direct intellectual lineage: you generalize the foundational JFQA paper.
- JFQA publishes methodological contributions with applications in asset pricing.
- The paper's length and technical depth are consistent with JFQA norms.
- The identification/stability framing I recommended earlier would play well at JFQA.
- JFQA has published PE/alternative-asset papers before.

Risks:
- JFQA is now quite competitive, with acceptance rates around 5-7%.
- A referee may view the contribution as incremental over Driessen et al. (2012) plus Korteweg and Nagel (2016).
- The diagnostic/negative empirical findings may not feel like enough of a new result for JFQA unless reframed sharply.
- The asymptotic theory may be viewed as more elaborate than necessary for the empirical conclusions.

Probability: **Moderate**, perhaps 15-25%. Realistic but not easy. The paper needs to convince the editor that it is more than a technical refinement of Driessen et al.

### 3. Review of Finance (RF)

**RF is a strong option, particularly if you want a respected general-finance journal.**

Why:
- RF publishes solid methodological and empirical finance papers.
- PE asset pricing is well within scope.
- RF is more receptive than top-3 to papers that are methodologically interesting but not blockbuster in empirical findings.
- The identification/stability framing would resonate with RF's audience.

Risks:
- RF referees may be less technically sophisticated on the econometric side, so the asymptotic theory may need to be streamlined.
- RF wants clear economic implications, so the "what can we learn from PE cash flows" angle needs to be front and center.

Probability: **Moderate**, perhaps 20-30%. Good fit if the paper is well-positioned.

---

## Tier 2: Plausible with repositioning

### 4. Journal of Business and Economic Statistics (JBES)

**JBES could work if you reposition the paper toward the econometric contribution.**

Why:
- JBES publishes applied econometric methods papers with empirical illustrations.
- The finite-sample analysis, SHAC inference, cross-validation under dependence, and simulation design are valued at JBES.
- JBES readers appreciate careful Monte Carlo analysis and practical guidance for implementation.

Risks:
- The PE-specific framing may be too narrow for JBES unless you frame the methodology as applicable to any non-traded asset class.
- JBES will want the econometric contribution to be clearly separable from the application.
- You would likely need to deemphasize the PE-specific empirical results and emphasize the general estimation framework.

Probability: **15-20%** with repositioning. Less natural than JFnE or JFQA but plausible.

### 5. Journal of Empirical Finance

**JEmpFin could work if you shift emphasis toward the empirical findings.**

Why:
- JEmpFin is receptive to applied asset-pricing papers with methodological awareness.
- The NAV sensitivity analysis, identification diagnostics, and BO/VC comparison would be valued.
- JEmpFin readers care about practical guidance for empirical researchers.

Risks:
- The asymptotic theory is much more than JEmpFin typically publishes.
- You would need to move most of the theory to appendices and lead with empirics.
- JEmpFin may view the paper as too technical.

Probability: **20-25%** if streamlined toward empirics. Lower prestige than JFQA but higher acceptance probability.

### 6. Critical Finance Review (CFR)

**CFR could work if you frame the paper as challenging conventional PE alpha estimates.**

Why:
- CFR publishes methodologically rigorous papers that challenge conventional wisdom.
- The "identification limits of PE alpha" angle fits CFR's editorial identity.
- CFR is receptive to shorter, sharper papers with a clear point.

Risks:
- The paper is currently much too long for CFR.
- CFR would want you to lead with the challenge to existing PE alpha findings and treat the methodology as supporting evidence.
- You would need a very different presentation.

Probability: **10-15%** with major restructuring. High-risk, high-reward if the editor is interested.

---

## Tier 3: Possible but less natural

### 7. RFS (Review of Financial Studies)

Why it could work:
- RFS has published PE methodology papers, including Brown et al. (2023).
- If the paper's identification contribution is viewed as genuinely advancing how the field thinks about PE risk measurement, RFS could be interested.

Why it probably will not work in current form:
- RFS is extremely competitive.
- The empirical findings are diagnostic rather than revealing a new economic result.
- The paper would need to be significantly shorter and sharper.
- A desk reject is likely unless the editor sees a clear advance over Korteweg and Nagel (2016).

Probability: **5-10%**. High risk of desk reject.

### 8. Journal of Econometrics

Why it could work:
- The increasing-domain asymptotics, NED random-field arguments, and SHAC inference are within scope.
- Reference papers by Jenish and Prucha, Conley, and Kim and Sun were published here.

Why it probably will not work:
- JoE wants the econometric contribution to be the primary contribution.
- The PE application would need to be treated as an illustration, not the main event.
- The econometric novelty may be viewed as an application of existing tools to a new setting rather than a new econometric result.

Probability: **5-10%**. Only viable if you strip out most of the PE content and lead with the estimation theory.

### 9. Management Science

- Possible but unusual for this type of paper.
- MS publishes some PE papers but wants broader management/decision implications.
- The methodological depth would need to be reduced substantially.
- Probability: **5-10%**.

---

## Tier 4: Fallback options

### 10. Quantitative Finance
- Technically within scope.
- Less prestige but reasonable acceptance probability.
- Good fit for the technical content.

### 11. Journal of Financial Data Science
- You already published Tausch and Pietz (2024) here.
- The code-availability and reproducibility aspects would be valued.
- But this would undervalue the theoretical contribution.

### 12. Journal of Banking and Finance
- Broad scope, reasonable fit.
- Less prestige but higher acceptance probability.
- The paper would need to be shorter and more empirically focused.

---

# My honest ranking

If I had to advise you on where to submit, here is my ordering based on fit, realistic acceptance probability, and career value.

| Rank | Journal | Fit | Acceptance probability | Notes |
|------|---------|-----|----------------------|-------|
| 1 | **JFnE** | Excellent | 25-35% | Best methodological fit |
| 2 | **JFQA** | Very good | 15-25% | Natural home given Driessen et al. |
| 3 | **Review of Finance** | Good | 20-30% | Good fallback, respected |
| 4 | **JBES** | Good with repositioning | 15-20% | Needs econometric emphasis |
| 5 | **J. Empirical Finance** | Good with repositioning | 20-25% | Needs empirical emphasis |
| 6 | **CFR** | Niche fit | 10-15% | Needs major restructuring |
| 7 | **RFS** | Stretch | 5-10% | High desk-reject risk |

---

# Strategic recommendation

My recommended submission strategy would be:

## Option A: Aim high, prepare for revision
Submit to **JFQA** first, with the identification/stability reframing.
- If R&R: great.
- If reject with constructive feedback: revise and send to JFnE or RF.

## Option B: Maximize expected value
Submit to **JFnE** directly.
- Best fit for the paper as written.
- Highest probability of a constructive review process.
- JFnE is well-respected in the financial econometrics community.
- Your paper's technical depth is a strength there, not a liability.

## Option C: Dual-track strategy
If the paper can be split:
- A shorter, sharper identification/stability paper for **JFQA** or **RF**, focused on the horizon tradeoff and empirical identification limits.
- A companion technical paper on the increasing-domain asymptotics and SHAC inference for **JFnE** or **JBES**.

This is aggressive and only worthwhile if the two papers can stand alone.

---

# What I would actually do

If this were my paper, I would:

1. **Reframe around identification/stability** as discussed in the previous conversation.
2. **Tighten to around 40-45 pages** including appendices.
3. **Submit to JFQA first.**
4. **If rejected, revise for JFnE.**

The reasoning: JFQA is the natural home given Driessen et al. (2012), and the identification/stability framing is strong enough to be competitive. But JFnE is the safety option where the paper's technical profile is a genuine advantage rather than a potential liability.

---

# One more thing

Whichever journal you target, the single most important revision for journal fit is the one we already discussed: **reframe the paper as an identification paper, not a factor-discovery paper.** That reframing helps at every journal on this list.

At JFQA, it avoids the "what is the true PE factor structure" trap.
At JFnE, it positions the horizon tradeoff as a methodological insight.
At RF, it frames diagnostic empirical findings as a contribution.
At JBES, it emphasizes the econometric content of identification analysis.

The substance of the paper barely needs to change. What needs to change is how the paper presents itself.