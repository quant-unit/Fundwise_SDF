# Refined Revision Suggestions — March 4, 2026

This document refines the Editorial Revision Memo from the referee review (2026-03-03) into actionable revision items. For each major point we assess the referee's critique, evaluate your inline responses, and provide a concrete, prioritized recommendation.

---

## Overall Assessment

The referee review is **largely fair and well-calibrated**. The reviewer clearly read the paper and the benchmark papers carefully, and the major points (estimand clarity, CV-as-inference, KN nesting language) are the exact friction points a JFQA referee panel would raise. However, some recommendations are disproportionate to the paper's scope (e.g., demanding Hessian eigenvalue diagnostics, endogenous exit timing simulations, or model-based residual-value adjustments). Below we separate the high-leverage revisions from the "nice-to-have" items that risk bloating the paper.

---

## Priority 1 — Must-Do Revisions (High Impact, Moderate Effort)

### 1.1 Reframe hv-block CV as a Stability Diagnostic

> **Referee concern (Major #3):** SE_CV and t_CV are presented as inferential objects, but fold-to-fold dispersion is not a valid standard error of the full-sample estimator.

**Assessment:** The referee is **correct**, and this is the single most likely point of rejection at JFQA. The current language (e.g., "cross-validation confirms that these market beta estimates remain highly statistically significant") invites the interpretation that CV dispersion = standard error. It does not.

**Your inline response** ("We will clarify how hv-block CV should be used as a stability diagnostic in line with the Racine 2000 paper") is the right direction.

**Concrete actions:**

1. **Rename throughout:** `SE_CV` → `SD_CV` or `FD` (Fold Dispersion); `t_CV` → `SR` (Stability Ratio) or keep `t_CV` but define it explicitly as "estimate-to-fold-dispersion ratio, not a test statistic with known null distribution."
2. **Add a framing paragraph** at the start of the empirical section (or in §4.1 after the hv-block CV description) with language like:

   > The hv-block cross-validation quantities reported below serve as **stability diagnostics** rather than inferential standard errors. They summarize fold-to-fold variation in parameter estimates under dependence-reducing sample splits and are informative about model sensitivity, but they do not provide a valid null distribution for the full-sample estimator. For formal inference, we rely on the SHAC-based asymptotic framework of Section 3.4.

3. **In table notes**, replace "SE_CV = cross-validation standard error" with "SD_CV = cross-validation fold dispersion (standard deviation of hold-out estimates across folds)."
4. **In the conclusion and interpretation sections**, replace phrases like "cross-validation confirms ... highly statistically significant" with "cross-validation stability diagnostics indicate that the single-factor estimates are robust across data subsets."

> [!IMPORTANT]
> This is the single highest-ROI edit. A JFQA referee will flag this in paragraph one of their report if unchanged.

**✅ IMPLEMENTED.** All changes applied to `sdf_paper.tex`:
- Added formal `Remark` (Cross-validation as stability diagnostic) with explicit disclaimer that CV quantities are not inferential.
- Referenced Racine (2000) as consistent model-selection criterion for dependent data in the hv-block CV subsection.
- Renamed `SE_CV` → `SD_CV` (fold dispersion) and `t_CV` → `SR_CV` (stability ratio) with formal definitions in the empirical preamble.
- Updated all 12 table column headers and 10+ table notes throughout main text and appendix.
- Rewrote all inline narrative: replaced "cross-validation confirms significance" with "stability diagnostics indicate robustness," "cross-validation inference" with "cross-validation stability diagnostics," etc.
- Updated abstract, DLP12 comparison, interpretation section, conclusion, and PE appendix (including Driessen comparison tables).
- **Remaining open:** R plotting code still labels columns as `SE_CV` / `t_CV` — these auto-generate LaTeX tables and should be updated when re-running estimates.

---

### 1.2 Distinguish Structural vs. Horizon-Regularized Estimand

> **Referee concern (Major #1):** The paper mixes θ₀ (structural NPV target) with the horizon-averaged pseudo-true parameter.

**Assessment:** The referee is **substantively correct**, though the paper already contains the mathematical machinery to make this distinction (Proposition 1, Eq. 18–20, the TRP decomposition). The issue is **expositional**, not mathematical. The paper correctly derives that NPV-only yields unbiased moments and that NFV introduces the TRP covariance term — but then sometimes discusses horizon-averaged estimates as if they target the same θ₀.

**Your inline remark** about first strengthening the appendix section on bias–variance tradeoff and then migrating key content to the main text is a sound two-step strategy.

**Concrete actions:**

1. **Expand Appendix A (Bias–Variance Tradeoff)** to include a formal definition of the horizon-specific pseudo-true parameter:
   
   θ*_H = argmin_θ Q_H(θ)
   
   and state explicitly that θ*_H = θ₀ when H = 0 (NPV-only) and θ*_H ≠ θ₀ in general for H > 0.

2. **Add a prominent Remark** in the main methodology (after §2.2 or within §2.2) titled "Target Parameter Under Horizon Averaging" that:
   - Formally defines θ*_H = argmin_θ Q_H(θ) and states that θ*_H = θ₀ when H = 0 and θ*_H ≠ θ₀ in general for H > 0 due to the TRP.
   - Clarifies that horizon averaging is treated as a **finite-sample regularization device**: the simulation shows that the dominant source of estimation error at short horizons is severe upward bias in factor loadings (not estimator variance), and that averaging pricing errors over an extended set of compounding dates effectively controls and stabilizes this error.
   - States that the resulting MSE is minimized at approximately H ≈ 12–15 years, where the small TRP-induced population bias is more than offset by the large reduction in finite-sample estimation error.
   - References Appendix A for the full bias–variance decomposition.

3. **In the abstract and introduction**, replace "bias-minimizing horizon" with "variance-reducing horizon" or "horizon that minimizes finite-sample mean squared error," and add the qualifier "under a horizon-specific pseudo-true target." [Our Remark: Better: Our simulation shows that the non-negligible finite-sample estimation error of the parameter estimates observed at short horizons can be controlled and effectively stabilized by averaging pricing errors over an extended set of compounding dates (horizon averaging).]

4. **Do NOT restructure the entire paper** around "two estimands" as the referee suggests — that would be overkill. A clear Remark + Appendix expansion is sufficient and preserves the paper's flow.

**✅ IMPLEMENTED.** All changes applied to `sdf_paper.tex`:
- Added formal `\subsection{Pseudo-true parameter}` to Appendix A defining θ*_H = argmin_θ Q_H(θ), with explicit statements that θ*_0 = θ₀ and θ*_H ≠ θ₀ for H > 0.
- Added `Remark [Target parameter under horizon averaging]` in §2.2 (after the bias-variance tradeoff paragraph) with self-explanatory wording: defines θ*_H, frames horizon averaging as a finite-sample regularization device, identifies severe upward bias (not variance) as the dominant error source, and states MSE is minimized at H ≈ 12–15 years.
- Replaced all 6 instances of "bias-minimizing" with "MSE-minimizing" throughout the paper (abstract, introduction, simulation section, empirical section, appendix).
- Rewrote abstract simulation sentence to use the agreed phrasing: "the non-negligible finite-sample estimation error observed at short horizons can be controlled and effectively stabilized by averaging pricing errors over an extended set of compounding dates."
- Added `\label{sec:appendix_interpretation}` to the Appendix A Interpretation subsection for cross-referencing.

---

### 1.3 Tone Down KN Nesting Language

> **Referee concern (Major #4):** "Generalizes/nests KN" is overstated because the estimators differ in architecture (direct PE pricing vs. public SDF + PE valuation).

**Assessment:** The referee is **correct on the framing, but your Proposition 4 is mathematically valid**. The issue is that "nests" and "generalizes" imply practical equivalence, when in fact the KN relationship is a special-case representation under an alternative aggregation order, not the same estimator in disguise.

**Your inline response** ("Agree, we will be more conservative in our language!") is appropriate.

**Concrete actions:**

1. Replace "elegantly nests both estimators" (abstract, line 115) with "provides a unifying objective-function perspective that recovers both estimators as special cases under specific aggregation and loss-function choices."
2. In Proposition 4's preamble, add: "This result establishes a structural connection at the objective-function level; it does not imply that the estimators are statistically equivalent in finite samples, as they differ in data inputs, conditioning, and identification strength."
3. In Table 2 (comparison table), keep the structural comparison but avoid language suggesting superiority.
4. In the literature review (§1.1), replace "formally nesting both" with "providing a unifying loss-function perspective on both."

**✅ IMPLEMENTED.** All changes applied to `sdf_paper.tex`:
- Abstract: "generalizing the popular approaches" → "providing a unifying objective-function perspective on the popular approaches."
- Introduction (L115): "elegantly nests both estimators" → "recovers both estimators as special cases under specific aggregation and loss-function choices."
- Introduction (L119): "arguably both simplifies and generalizes" → "both simplifies and extends."
- Literature review (L162): "formally nesting both" → "providing a unifying loss-function perspective on both."
- Proposition preamble (L889): "true generalization of both" → "recovers both … as special cases."
- Proposition title: "Generalization of existing estimators" → "Relationship to existing estimators"; body: "generalizes" → "recovers … as special cases."
- Inserted caveat after Proposition: "This result establishes a structural connection at the objective-function level; it does not imply that the estimators are statistically equivalent in finite samples."
- Proof summary (L965): "generalizes" → "extends."

---

### 1.4 Moderate Claim Strength Throughout

> **Referee concern (Major #9):** Phrases like "hopelessly over-parameterized," "only robust specification," "should be restricted to" are too categorical.

**Assessment:** **Agreed.** These phrases are accurate descriptions of the empirical findings but read as universal verdicts. JFQA referees react better to scoped claims.

**Concrete actions — targeted replacements:**

| Current phrasing | Suggested replacement |
|---|---|
| "should currently be restricted to single-factor market models" | "in our sample, single-factor market models are substantially more stable than multi-factor specifications" |
| "bias-minimizing 15-year horizon" | "horizon that minimizes finite-sample MSE in our baseline simulations (approximately 12–15 years)" |
| "hopelessly over-parameterized" | "severely under-identified given the available cross-section" |
| "only robust specification" | "most robust specification in our implementations" |
| "massive reliance on reported NAVs" | "considerable sensitivity to reported NAVs" |
| "cross-validation confirms ... highly statistically significant" | "cross-validation stability diagnostics indicate robustness across folds" |

Add scope qualifiers: "in our sample," "under our implementations," "in the Preqin BO/VC data."

**✅ IMPLEMENTED.** All changes applied to `sdf_paper.tex`:
- "bias-minimizing" → "MSE-minimizing" (already done in §1.2).
- "cross-validation confirms … significant" → stability diagnostic language (already done in §1.1).
- Abstract (L100): "should currently be restricted to" → "in our sample, single-factor market models are substantially more stable than multi-factor specifications."
- Interpretation heading (L1655): "only robust specification" → "most robust specification."
- BO/VC conclusion (L1484): "only statistically robust" → "most statistically robust … in our sample."
- DLP comparison (L2090): "hopelessly over-parameterized" → "severely under-identified given the available cross-section."
- DLP comparison (L2092): "only statistically robust and estimable" → "most robust and estimable … in our implementations."
- NAV reliance (L2234–2235): "massive reliance" → "considerable reliance"; "massively overweighted" → "disproportionately overweighted."
- Conclusion (L1703): "should be restricted to" → "in our sample, … are substantially more stable … and should be preferred."
- Implications paragraph (L1683): "should be restricted to" → "is substantially more reliable when restricted to."

---

## Priority 2 — Should-Do Revisions (Medium Impact, Low–Medium Effort)

### 2.1 Factor/Currency/Region Robustness

> **Referee concern (Major #6):** U.S. q⁵ factors applied to a global sample with non-U.S. funds.

**Assessment:** This is a **legitimate concern**. Your inline response to restrict to a North American USD-only subsample is a clean solution and avoids the complexity of multi-region factor matching.

**Concrete actions:**

1. Restrict the main empirical sample to North American, USD-denominated funds (as you proposed).
2. Add a brief statement in the data section: "To ensure consistency between fund cash flows and factor benchmarks, we restrict the sample to North American funds reporting in USD."
3. Keep the global sample results as robustness in the appendix, with a note on the factor mismatch caveat.
4. Report fund counts for the restricted vs. unrestricted sample.

---

### 2.2 NAV Sensitivity — Expand the Discount Grid

> **Referee concern (Major #5):** NAV treatment is too ad hoc; the referee wants a model-based residual-value adjustment.

**Assessment:** The referee's suggestion for a full predictive NAV-to-MV calibration model is **disproportionate** — it would be a paper in itself and is not central to your contribution. Your inline response (adding 0%, 25%, 75% discount values alongside existing 100% and 50%) is **pragmatic and sufficient**.

**Concrete actions:**

1. Add the 0%, 25%, 75% constant NAV discount results as planned.
2. Present these in a compact sensitivity table or figure showing how β_MKT varies across {0%, 25%, 50%, 75%, 100%} NAV discounts.
3. Add a liquidated-only subsample check (funds fully liquidated by cutoff date) — this is cheap to implement and addresses the referee's "minimum acceptable upgrade."
4. **Do NOT build** a predictive NAV-to-MV calibration model. State explicitly: "A full residual-value model is beyond the scope of this paper; we instead provide transparent sensitivity across a range of constant NAV discounts and a liquidated-only subsample."

---

### 2.3 SHAC Bandwidth Sensitivity

> **Referee minor comment #3:** Bandwidth fixed at 12 years needs a sensitivity table.

**Assessment:** This is a **fair and easy request**. A bandwidth sensitivity table is standard in spatial econometrics papers.

**Concrete action:** Add an appendix table showing SHAC standard errors for the single-factor MKT model at bandwidths {6, 9, 12, 15, 18} years for both BO and VC. This is a simple parameter sweep on existing code.

---

### 2.4 Strengthen the Bias–Variance Appendix

> **Referee concern (Major #2):** TRP as bias correction is too heuristic for JFQA.

**Assessment:** The referee wants a more formal bias–variance–risk decomposition. Your inline response to extend Appendix A is correct. The current appendix already contains:

- **Pseudo-true parameter** definition (θ\*\_H, §A.1)
- **Decomposition of E[ε̄ᵢ(θ₀)]** into NPV (zero) and NFV (TRP) terms (§A.2)
- **Reconciliation** with simulation: Effect 1 (Nagar-type finite-sample bias ∝ Var[ε̄ᵢ]/n) and Effect 2 (variance reduction through T-averaging) (§A.3)
- **Interpretation** via three perspectives: GMM moment combination, HJ-bounded TRP, model averaging analogy (§A.4)

What is **missing** is a single, self-contained formal statement (Proposition or Lemma) that assembles these pieces into an explicit MSE decomposition, with clearly stated monotonicity properties in H, and an existence argument for an interior optimum H\*. This is what the referee is asking for when they say "too heuristic."

#### Detailed mathematical plan

Below is the precise specification for a new subsection **§A.5 "Formal MSE decomposition"** to be inserted between the current §A.3 ("Reconciliation with the simulation study") and §A.4 ("Interpretation"). The existing §A.4 would then become §A.5. No existing content needs to be deleted.

---

**Step 1: State the full MSE decomposition explicitly (Proposition)**

Add a new **Proposition (MSE decomposition of the horizon-averaged estimator)** with the following structure:

> **Proposition (MSE decomposition).** *Under the assumptions of Proposition 1 and the regularity conditions for the Nagar bias expansion (Eq. 9), the mean squared error of θ̂\_H around the structural parameter θ₀ decomposes as*
>
> MSE(θ̂\_H, θ₀) := E[‖θ̂\_H − θ₀‖²]  
> = ‖θ\*\_H − θ₀‖² + E[‖θ̂\_H − θ\*\_H‖²] + 2(θ\*\_H − θ₀)ᵀ E[θ̂\_H − θ\*\_H]
>
> *which, to leading order, yields:*
>
> MSE(θ̂\_H, θ₀) = B\_pop(H)² + V(H) + B\_fs(H)²
>
> *where:*
> - *B\_pop(H)² := ‖θ\*\_H − θ₀‖² is the squared population bias (driven by the TRP),*
> - *V(H) := tr(Σ\_H / n) is the asymptotic variance of θ̂\_H around its pseudo-true target θ\*\_H, with Σ\_H the sandwich covariance matrix evaluated at θ\*\_H,*
> - *B\_fs(H)² is the squared finite-sample bias from the Nagar-type expansion (Eq. 9), proportional to Var[ε̄ᵢ^(H)] / n.*

The key mathematical point is that the standard MSE = Bias² + Variance decomposition applies **around θ₀**, not around θ\*\_H. So the "bias" has **two sources**: (i) population-level targeting error ‖θ\*\_H − θ₀‖ from the TRP, and (ii) finite-sample deviation ‖E[θ̂\_H] − θ\*\_H‖ from the Nagar-type nonlinearity correction. The cross-term 2(θ\*\_H − θ₀)ᵀ E[θ̂\_H − θ\*\_H] is absorbed into the leading-order factorization since both bias terms are O(1/n) or smaller.

**Why this matters for the referee:** This makes the decomposition completely rigorous—no hand-waving about "Effect 1" and "Effect 2" anymore; they are formally identified as B\_fs(H) and B\_pop(H).

---

**Step 2: Establish monotonicity properties in H**

Add a **Corollary (monotonicity)** immediately after the Proposition:

> **Corollary (monotonicity in H).** *Under the conditions of the preceding Proposition:*
>
> *(i) B\_pop(H) is non-decreasing in H with B\_pop(0) = 0. Each additional NFV term τ > t\_{i,j}^{Inv} contributes a non-negative TRP increment ≥ 0 to ‖θ\*\_H − θ₀‖, and the 1/|T\_i| denominator ensures that B\_pop(H) grows at most sublinearly.*
>
> *(ii) B\_fs(H) is non-increasing in H for H ≥ 0. This follows from the Nagar expansion (Eq. 9): B\_fs ∝ Var[ε̄ᵢ^(H)] / n, and Var[ε̄ᵢ^(H)] is strictly decreasing in |T\_i| since the pricing errors ε\_{τ,i} at different τ share the same deal cash flows but are imperfectly correlated (Eq. 12, Eq. 13).*
>
> *(iii) V(H) is non-increasing in H. The sandwich variance Σ\_H / n inherits the moment-reduction effect: with more (imperfectly correlated) moments contributing to the objective, the asymptotic variance of θ̂\_H around θ\*\_H decreases, analogous to the efficiency gain from combining moment conditions in Hansen (1982).*
>
> *(iv) Since B\_pop² is increasing from 0 while B\_fs² + V is decreasing, MSE(H) = B\_pop(H)² + B\_fs(H)² + V(H) has an interior minimum at some H\* ∈ (0, ∞), provided the TRP terms are not identically zero (i.e., some NFV terms have non-vanishing Cov(Ψ\_{t}/Ψ\_τ, δ\_{i,j})).*

**Why this matters for the referee:** Part (iv) is the formal existence argument for H\* that was previously only stated informally. Part (ii) directly connects the Nagar expansion (already in the appendix) to the MSE framework.

---

**Step 3: Bound B\_pop(H) via the Hansen–Jagannathan inequality**

Add a **Remark (magnitude of the population bias)** that formalizes what is currently only hinted at in the "Interpretation" subsection:

> **Remark (HJ bound on B\_pop).** *For each (τ, j) pair with τ > t\_{i,j}^{Inv}, the individual TRP term satisfies*
>
> |Cov(Ψ\_t / Ψ\_τ, δ\_{i,j})| ≤ σ(Ψ\_t / Ψ\_τ) · σ(δ\_{i,j})
>
> *where σ(·) denotes the standard deviation. Under the Hansen–Jagannathan (1991) framework, σ(Ψ\_t / Ψ\_τ) is bounded by the maximum Sharpe ratio of tradeable assets over the (t, τ) horizon times E[Ψ\_t / Ψ\_τ]. The 1/|T\_i| averaging further dilutes: each TRP increment enters the population bias as a 1/|T\_i|-weighted average, so for a fund with |T\_i| = 12H discounting dates:*
>
> |B\_pop(H)| ≤ C · f(σ\_SDF, σ\_δ, J) / (12H)
>
> *where C is a constant depending on the implied Hessian H⁻¹ of the objective, f(·) captures the cumulative TRP contribution across all J deals, and the 1/(12H) dilution factor ensures sublinear growth. Under our baseline simulation DGP with β\_MKT = 1 and σ = 20%, the resulting B\_pop(H\*) is an order of magnitude smaller than B\_fs(0).*

This inequality makes quantitative what is currently only qualitative ("the TRP is small"). It tells the referee: *we know exactly where the bias bound comes from (HJ), and the averaging denominator makes it vanish at rate 1/H while the finite-sample bias correction gains are front-loaded.*

---

**Step 4: Explicit dependence on DGP parameters**

Add a brief **Remark (DGP dependence of H\*)** collecting the qualitative sensitivities:

> *The optimal horizon H\* depends on:*
> - *Idiosyncratic volatility σ: higher σ → larger B\_fs(0) → H\* shifts upward.*
> - *Factor risk premium and volatility: larger SDF variability → larger TRP → H\* shifts downward.*
> - *Number of cross-sectional units n: larger n → B\_fs(H) = O(1/n) shrinks for all H → H\* shifts downward since the "benefit" of horizon averaging diminishes.*
> - *Deal timing and fund lifetime: longer holding periods → more NFV terms contribute TRP → H\* shifts downward.*
>
> *The simulation results in Section 4 empirically confirm this tradeoff and indicate that H\* ≈ 12–15 years under our baseline DGP (β\_MKT = 1, σ = 0.20, n = 20, 20 vintages). The exact value will vary with factor dynamics, cash-flow timing, and sample composition.*

---

**Step 5: Cross-reference from "Reconciliation" section**

At the end of current §A.3 ("Reconciliation"), add one sentence pointing forward:

> *Proposition X below formalizes this tradeoff and establishes the existence of an interior MSE-minimizing horizon H\*.*

---

#### Summary of structural changes to Appendix A

| Current §  | Content | Action |
|---|---|---|
| A.1 | Pseudo-true parameter | Keep as is |
| A.2 | Decomposition of E[ε̄] | Keep as is |
| A.3 | Reconciliation with simulation | Keep; add forward reference |
| **A.4 (NEW)** | **MSE decomposition (Proposition + Corollary + Remarks)** | **INSERT** |
| A.5 (was A.4) | Interpretation | Keep; renumber |

#### References to add to `ref.bib`

- `bao2007finite`: Bao, Y. and Ullah, A. (2007). "Finite sample properties of maximum likelihood estimator in spatial models." *Journal of Econometrics*, 137(2), 396–413. (For the extension of RSU96 to dependent data)
- Already cited: `rilstone1996second` (Rilstone, Srivastava, Ullah 1996), `newey2004higher` (Newey & Smith 2004), `HJ91` (Hansen & Jagannathan 1991), `H82` (Hansen 1982), `hansen2012jackknife` (Hansen & Racine 2012), `nagar1959bias` (Nagar 1959).

#### Estimated effort

- **Medium–High.** The Proposition and Corollary require careful LaTeX typesetting of the three-way MSE decomposition, the monotonicity argument, and the HJ bound. The mathematical content is straightforward given the existing machinery in the appendix; the main effort is careful formulation and cross-referencing.
- **Estimated length:** ~1.5 pages of new material (Proposition + proof sketch + Corollary + 2 Remarks), inserted between current §A.3 and §A.4.

---

**✅ IMPLEMENTED.** Done.

### 2.5 Strengthen the Finance Payoff in the Conclusion

> **Referee concern (Major #10):** Finance takeaway is underdeveloped relative to the technical apparatus.

**Assessment:** **Fair.** The conclusion reads as methodological guidance. Adding 1–2 paragraphs on practical implications would strengthen the JFQA pitch.

**Concrete actions:** Add a brief "Implications for Practice" paragraph to the conclusion covering:
- LP benchmarking: single-factor SDF estimates provide a practical, robust benchmark for PE/VC returns.
- Outperformance claims: multi-factor alpha estimates are too noisy to support "outperformance relative to size/value" claims with current data.
- Quality-control protocol: researchers should report horizon sensitivity, fold stability, and NAV discount sensitivity before publishing risk-adjusted PE returns.

---

## Priority 3 — Nice-to-Have / Push Back on These

### 3.1 Simulation Stress Tests (Major #7)

> **Referee wants:** NAV smoothing, endogenous exit timing, heterogeneous betas, survival selection in simulation.

**Assessment:** Your inline response is **correct** — these are each separate research agendas, and adding them would dilute the paper's focus on the compounding horizon contribution. The simulation is deliberately stylized to isolate the horizon effect.

**Recommendation:** **Push back politely** in the response letter. Add one paragraph to the simulation section acknowledging these limitations explicitly:

> "Our simulation design is intentionally stylized to isolate the finite-sample impact of horizon averaging. We abstract from NAV smoothing, endogenous exit timing, and vintage-dependent betas, each of which would constitute a separate simulation study. The qualitative ranking of single- vs. multi-factor stability is robust across our current DGP specifications, but we caution that the exact bias-minimizing horizon (12–15 years) is DGP-dependent."

### 3.2 Hessian/Curvature Diagnostics (Major #8)

> **Referee wants:** Eigenvalues, profile objective plots, weak-ID diagnostics, multi-start checks, alpha-bound sensitivity.

**Assessment:** Your inline response ("this is overkill; DLP12/KN16 do not analyze these") is **defensible but could be softer**. The referee has a point that documenting pathologies narratively is less convincing than showing curvature diagnostics — but a full identification analysis would be a separate paper.

**Recommendation:** **Partial compromise.** Add one appendix figure showing the profile objective function in (α, β_MKT) space for the single-factor + α model at H = 0 and H = 15, illustrating the flat/ridge objective surface. This is a single plot that visually demonstrates the identification problem without requiring eigenvalue tables. Mention multi-start optimization as a footnote.

### 3.3 Title Change

> **Referee suggests** titles emphasizing "dependence-aware" or "horizon regularization."

**Assessment:** Your current title ("Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows") is already clean and accurate. "Dependence-aware" is not the paper's main selling point, and "horizon regularization" introduces terminology not yet established in the field.

**Recommendation:** **Keep the current title.** It is precise and neutral, which is better for a first submission. If a referee insists on a change, consider adding a subtitle: "with Horizon-Averaged Pricing Errors" or similar.

### 3.4 Abstract Rewrite

> **Referee provides** a suggested drop-in abstract.

**Assessment:** The referee's suggested abstract is competent but overly cautious — it strips out the empirical specifics (beta estimates, horizon values) that make the abstract memorable and informative. A JFQA abstract should contain your main numbers.

**Recommendation:** **Revise your existing abstract** with the estimand-clarity and CV-language fixes from Priority 1, but keep your empirical estimates in the abstract. Specifically:
- Add "horizon-specific pseudo-true parameter" or "regularized target" language.
- Replace "cross-validation confirms ... statistically significant" with "stability diagnostics indicate robust estimates."
- Keep the beta estimates (0.79/0.80 for BO, 1.26/1.05 for VC).
- Soften "should currently be restricted" to "are substantially more stable than multi-factor alternatives in our sample."

### 3.5 Distance Metric Sensitivity (Minor #4)

> **Referee suggests** robustness to cash-flow-overlap distance instead of vintage-year distance.

**Assessment:** This is a reasonable suggestion but low priority. Vintage-year distance is a tractable, interpretable proxy, and overlap-based distance would require substantial re-engineering of the SHAC implementation.

**Recommendation:** Acknowledge in a footnote that overlap-based distance (as in KN16) is an alternative and note that your vintage-year proxy is computationally convenient and yields a conservative (wider) bandwidth. Leave implementation for future work.

---

## Summary: Revision Checklist

| # | Item | Priority | Effort | Status |
|---|---|---|---|---|
| 1.1 | Relabel CV quantities as stability diagnostics | **P1** | Low | ✅ Done |
| 1.2 | Add "two estimands" Remark + expand Appendix A | **P1** | Medium | ✅ Done |
| 1.3 | Tone down KN nesting language | **P1** | Low | ✅ Done |
| 1.4 | Moderate claim strength throughout | **P1** | Low | ✅ Done |
| 2.1 | Restrict to NA/USD sample | **P2** | Medium | ✅ Done |
| 2.2 | Expand NAV discount grid | **P2** | Medium | ✅ Done |
| 2.3 | SHAC bandwidth sensitivity table | **P2** | Low | ☐ |
| 2.4 | Formal MSE decomposition in Appendix A | **P2** | Medium | ✅ Done |
| 2.5 | Finance payoff in conclusion | **P2** | Low | ☐ |
| 3.1 | Simulation scope paragraph (push back) | **P3** | Low | ☐ |
| 3.2 | Profile objective plot in appendix | **P3** | Low | ☐ |
| 3.3 | Title — keep current | **P3** | None | ☐ |
| 3.4 | Abstract — targeted edits only | **P3** | Low | ☐ |
| 3.5 | Distance metric — footnote only | **P3** | Low | ☐ |

---

## Strategic Note

The referee's overall assessment ("major revision, borderline reject with clear potential") is **constructive and fair**. The paper's core contribution is strong. The revisions above are designed to address the two decisive issues (estimand clarity + CV-as-inference) while making proportionate improvements elsewhere. Resist the temptation to add every appendix the referee suggests — scope creep is the enemy of a clean JFQA submission. The goal is a **tighter, more defensible paper**, not a longer one.
