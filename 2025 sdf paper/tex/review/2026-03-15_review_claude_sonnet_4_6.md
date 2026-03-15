# Referee Report
## Journal of Financial and Quantitative Analysis
**Manuscript:** "Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows"
**Recommendation:** **Major Revision**

---

## I. Summary

The paper develops a generalized Least-Mean-Distance (LMD) framework for estimating stochastic discount factors from pooled non-traded cash flows. The central conceptual contribution is the formalization of a **Timing Risk Premium (TRP)**: when pricing errors are averaged over discounting dates beyond fund inception, the estimator targets a horizon-specific pseudo-true parameter rather than the structural SDF parameter. The authors show this creates a finite-sample bias-variance tradeoff — averaging over longer horizons can dramatically improve estimator stability, at the cost of shifting the population target. An increasing-domain asymptotic framework accommodating spatial dependence across vintage years is developed, and the framework is applied to Preqin buyout and venture-capital data.

The paper addresses a genuine and important methodological gap in the private equity asset-pricing literature. The identification of the TRP as a structural feature of cash-flow-based SDF estimation is original and clearly articulated. That said, I have several major concerns that I believe must be addressed before publication.

---

## II. Major Comments

### 1. The Horizon Selection Problem Is Characterized but Not Solved

The paper's central finding — that the MSE-minimizing horizon is interior, DGP-dependent, and approximately 12.5–15 years for PE-like environments — is important. However, the prescription offered to applied researchers ("tailored simulation designs") is operationally unsatisfying. The key questions left unanswered are:

- **How should a researcher select H in practice when the DGP is unknown?** The paper provides a grid of estimates across horizons, but without a principled selection criterion, the researcher faces a model-selection problem that could itself introduce data-snooping biases.
- The paper mentions hv-block cross-validation as a stability diagnostic but does not formally propose it as a horizon-selection criterion. If CV fold dispersion (SD$_{CV}$) is already being computed, why not formalize a CV-based horizon selector analogous to standard bandwidth selection in nonparametric estimation?

I would strongly encourage the authors to either (a) propose and evaluate a formal data-driven horizon selector, or (b) be more explicit about the inferential implications of reporting results over a grid of horizons without selection.

---

### 2. Conceptual Tension Between Theory and Practice

There is a fundamental tension that the paper acknowledges but does not fully resolve. The theoretical consistency result (Corollary 1) applies only at H=0 — the inception-only NPV criterion. Yet the simulation study strongly recommends H>0 in practice, since inception-only estimates are severely biased in finite samples. This means:

- Practitioners following the paper's empirical recommendations are estimating a **pseudo-true parameter θ†_H**, not the structural θ₀.
- The structural interpretation of the estimated alpha and beta coefficients becomes ambiguous.

The paper partially addresses this in Appendix A, but the treatment is asymptotic and does not yield operational guidance. I suggest the authors provide a more direct discussion — perhaps a dedicated subsection in Section 3 — on the **economic interpretation of θ†_H** relative to θ₀, including conditions under which the target shift is likely to be small enough that structural inference remains approximately valid.

---

### 3. The Relationship to Korteweg and Nagel (2016) Deserves More Formal Treatment

Section 3.5.2 argues that KN (2016) is "a closely related time-series GMM analogue rather than a literal special case." This characterization, while probably correct, is insufficiently precise. The argument that the two criteria differ due to "order of aggregation" (loss applied before vs. after aggregation) deserves formal development. Specifically:

- Under what conditions do the two estimators agree asymptotically? Do they share the same pseudo-true parameter?
- The claim that "different aggregation and weighting choices imply different sample objectives and estimators" should be formalized with at least one concrete illustrative result.

Given that KN (2016) is one of the two foundational papers the authors claim to generalize or contextualize, this relationship deserves treatment comparable in rigor to Proposition 5 (the Driessen et al. connection).

---

### 4. Characterization of the TRP Magnitude Is Incomplete

Lemma 2 and the surrounding analysis establish that the averaged population moment at θ₀ is shifted by post-investment covariance terms, and Remark 10 provides a Cauchy-Schwarz bound. However:

- The sign of the TRP is not characterized. In which direction does horizon averaging bias the estimated beta? The simulation evidence consistently shows downward attenuation of β̂_MKT with increasing H, but the theoretical mechanism producing this direction is not made explicit.
- There is no analytic approximation of the magnitude of the target shift as a function of observable DGP features (e.g., fund lifetime, idiosyncratic volatility, factor loadings). Such an approximation — even a rough one — would considerably strengthen the paper's practical guidance.

---

### 5. Treatment of Non-Liquidated Funds and NAV Endogeneity

The paper treats the latest reported NAV as the final cash flow for non-liquidated funds and then conducts a sensitivity analysis under various NAV haircuts. This is a reasonable pragmatic approach but raises a concern that deserves more explicit treatment:

- Reported NAVs are known to be managed and smooth relative to true valuations (Brown et al., 2019; Jenkinson et al., 2020). Using them as final cash flows without adjustment **mechanically reduces estimated volatility** of the cash-flow stream, which could bias not only the beta estimates but also the standard errors.
- The NAV sensitivity analysis (Section 4.4) is informative but treats all NAV discounts as uniform across funds. A more realistic approach might discount NAVs by fund age, vintage year, or investment style. The paper should at least discuss why uniform discounts are a reasonable approximation.

---

### 6. Simulation Calibration and External Validity

The baseline simulation uses 20 vintage years and 20 funds per vintage (400 total fund-years). This is sparse even relative to the authors' own empirical sample. Several concerns:

- The comparison to Driessen et al. (2012) in Appendix B is useful but reveals a **large unexplained discrepancy**: even under the most favorable scenario (S10), the authors cannot reproduce the Driessen et al. benchmark (β̂ = 0.93 vs. 0.98, alpha bias 0.27% vs. 0.05%). The authors attribute this to differences in the error distribution and deal timing, but the underlying source of the discrepancy is never definitively identified. This should be resolved or the implications for the paper's simulation conclusions should be discussed more carefully.
- The simulation uses i.i.d. cross-sectional error terms, abstracting from the very spatial dependence that motivates the asymptotic framework. A simulation that introduces vintage-year dependence would both stress-test the SHAC estimator and more closely reflect the data environment the paper describes.

---

## III. Minor Comments

**3.1 Notation and Presentation**
The paper is notation-heavy. Tables 1 and 2 help, but the simultaneous use of T, T_i, T_{i,H}, T^max_i, and τ in overlapping roles creates avoidable confusion. I recommend consolidating the notation and providing a concrete worked example (beyond Figure 1) illustrating the full estimation pipeline for a small hypothetical dataset.

**3.2 Alpha Bounds**
The alpha bound of ±2%/month is tight relative to the DGP range. Table 10 shows that 32% of draws hit the bound under Scenario S6, and the bound qualitatively affects simulation conclusions. The choice of ±2% should be better justified — either by reference to the empirical distribution of plausible PE alphas or by showing robustness to the bound across a range of values.

**3.3 hv-Block Parameters**
The hv-block implementation uses h=3 and v=3 vintage years, motivated informally by the typical duration of cash-flow overlap between adjacent vintages. A more rigorous justification — or at least a sensitivity check to h∈{2,4} — would strengthen confidence in the cross-validation diagnostics.

**3.4 Asymptotic Standard Errors and Their Reliability**
Tables 6–9 show that asymptotic standard errors are highly erratic, sometimes rounding to zero and sometimes exceeding 25, within the same table. The paper correctly attributes this to the limited cross-section and notes that CV diagnostics are more informative. However, if the asymptotic standard errors are this unreliable in the empirical application, one wonders about their reliability in the asymptotic theory sections. The authors should be clearer about the effective sample sizes at which the asymptotic approximation becomes adequate — perhaps via a simulation specifically designed to address this question.

**3.5 Bandwidth Choice for SHAC**
The choice of bandwidth b_V = 12 years for the Bartlett kernel (Section 4.6) is stated without derivation. Given that the total vintage span in the seasoned sample is roughly 25 years (1985–2010), a bandwidth of 12 years covers nearly half the data. The paper should report sensitivity of standard errors to bandwidth choices of, say, 6 and 18 years, and discuss how the theoretical requirement b_V/V → 0 maps to the finite-sample choice.

**3.6 Geographic Inconsistencies in Preqin**
Section 4.1 notes that some funds pass the North America filter despite geographic labels inconsistent with that region (Asia: 16, Canada: 1, Europe: 6, etc. in the BO sample). The authors dismiss this by noting all retained funds are USD-denominated. This is not fully convincing — a USD-denominated fund managed in Asia or Europe may have different risk exposures than a North American fund. At minimum, the authors should report robustness results excluding these geographically ambiguous observations.

**3.7 Discussion of the Literature on PE Outperformance**
The paper's empirical conclusion — that alpha is not reliably estimable with current data — has direct implications for the large literature on PE outperformance. The discussion in Section 5.4 and the conclusion touches on this, but a more explicit engagement with what the paper implies for prior estimates of PE alpha (e.g., Kaplan and Schoar 2005; Harris et al. 2014; Korteweg 2019) would strengthen the paper's contribution to the broader literature.

---

## IV. Summary Assessment

| Dimension | Assessment |
|---|---|
| Originality | Strong — TRP concept and LMD formulation are genuine contributions |
| Theoretical rigor | Strong — proofs are careful and the asymptotic framework is appropriate |
| Empirical execution | Adequate — honest about limitations but some robustness checks are missing |
| Practical guidance | Weak — the horizon selection problem is not resolved |
| Exposition | Adequate — heavy notation reduces accessibility |

The paper makes meaningful theoretical and empirical contributions to an active literature, and I encourage the authors to revise and resubmit. The core conceptual insight — that horizon averaging trades structural validity for finite-sample stability — is original and important. Addressing the concerns above, particularly regarding horizon selection, the KN (2016) relationship, and the TRP magnitude, would substantially strengthen the paper's impact.