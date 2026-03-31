# Methodological Review: Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows

## I. Formal/Technical Errors and Typos

### 1. Region filter inconsistency in the data section (L939, L966–968, L976)

Resolved!

---

### 2. Notation: $CF_t$ vs. $CF_{i,t}$ in Proposition 4 proof (L736)

Resolved!

---

### 3. "0.1 yr" Horizon label in Tables (L1407, L1439, etc.)

Resolved!

---

### 4. Panel B header says "$\beta_{\text{Alpha}}$" (L1520, L1572)

Resolved!

---

### 5. Capitalization inconsistency: "Horizon"

Resolved!

---

### 6. Mixed labeling: "FW" vs "value-weighted"
The paper introduces "FW" as "fund-size-weighted" (L942) but some captions and text refer to "value weighting" and the tables column header says "Value-Weighted (FW)". The acronym derivation should be made explicit once: *FW = fund-size-weighted (value-weighted)*.

---

### 7. Equation 11 time index
In Eq. 11 (L497), the pricing error sums from $t=0$ to $T$. But for funds with inception $t_i^0 > 0$, this includes periods before inception where $CF_{i,t}=0$ by Assumption 1 (part 4). While mathematically harmless (the terms vanish), it could be cleaner to sum from $t_i^0$ to $T_i^{\max}$, consistent with the fund-specific notation used elsewhere.

---

## II. Internal Inconsistencies

### 8. Abstract/conclusion claim vs. empirical evidence on "seasoned portfolios"
The abstract (L59) claims: *"the framework yields more stable evidence on broad market exposure **in seasoned portfolios**."*
However, the paper never estimates on *unseasoned* funds and compares stability — it only reports results for max vintage 2010 (seasoned). The vintage cutoff sensitivity (Section 4.3) shows *point estimates shift* when unseasoned vintages are included, but this is not a stability-vs-instability comparison for the full specification grid. The claim is slightly stronger than the evidence supports.

---

### 9. Simulation DGP uses $\alpha=0$ but the two-factor section estimates $\alpha$
The simulation calibrates $\alpha=0$ (L1115: *"we just use the MKT factor with $\beta_{\text{MKT}}=1$"*). In the two-factor $\alpha+$MKT simulations (L1182), $\alpha=0$ is the true value. Hence the simulations test whether the estimator can *correctly recover zero alpha*, but they never test whether a **non-zero** alpha can be recovered. An adversarial referee could argue that the paper concludes $\alpha$ is weakly identified partly because the DGP never forces the estimator to find one.

> [!NOTE]
> Consider adding a supplementary simulation with a true $\alpha \neq 0$ (e.g., $\alpha=0.002$/month) to demonstrate that even when alpha is present, the estimator remains unstable/weakly identified. This would strengthen the identification narrative.

---

### 10. Cross-validation fold structure and the 2010 cutoff

Resloved! 

---

### 11. Claim about DLP12 asymptotic device
Line 789 states: *"Their Theorem 1 fixes the number of fund-of-funds portfolios and lets the number of underlying projects or funds within each portfolio diverge."*

Resolved!

---

## III. Conceptual/Methodological Weaknesses

### 12. No formal characterization of $\theta_H^\dagger$ or the wedge magnitude
The paper's central contribution — the Compounding-Horizon Wedge — is stated qualitatively but never bounded. The reader knows that $\theta_H^\dagger \neq \theta_0$ for $H>0$ but has no sense of *how far* they can diverge. Even a simple special case (e.g., one-deal fund with one factor, Gaussian SDF) computing $\theta_H^\dagger - \theta_0$ in closed form would make the wedge tangible.

> [!IMPORTANT]
> This is the most likely point a top referee will press. The paper demonstrates the *existence* of the wedge but not its *magnitude*, sign, or dependence on model parameters. Adding even an approximate analytical expression would significantly strengthen the contribution.

---

### 13. The empirical analysis has no formal identification test
The paper argues extensively that $\alpha$ and two-factor models are weakly identified, but never applies a formal weak-identification test (e.g., Stock-Yogo-type statistic, rank test of the Hessian, or moment condition redundancy test). The evidence is entirely visual (Horizon plots) and descriptive (SR$_{CV}$ ratios). Consider at minimum reporting the condition number of $\hat{\mathcal{H}}_H$ or the smallest eigenvalue of the Hessian across specifications to give a numerical measure of near-singularity.

---

### 14. Quadratic loss only — no robustness to loss specification
Equation 14 defines a general loss $L(\cdot)$, but every simulation and empirical result uses $L(x) = x^2$. The paper does not discuss robustness to alternative loss functions (e.g., absolute deviation, Huber loss). Given that PE cash flows are heavy-tailed and often have extreme outliers, the quadratic loss may give excessive weight to outlier pricing errors. Acknowledge this limitation or test one alternative.

---

### 15. Treatment of NAV as "final cash flow" conflates two distinct issues

Resolved!

---

### 16. Vintage-year portfolio formation introduces a selection/weighting choice not fully explored
The paper emphasizes VYP formation for variance reduction but does not explore sensitivity to the portfolio construction rule within VYPs. For instance, equal-weighting vs. value-weighting *within* a vintage-year portfolio (before the outer EW/FW weighting across VYPs) is not discussed. This is a second-order weighting choice that can matter when fund sizes are heterogeneous within a vintage.

---

## IV. Exposition and Presentation Issues

### 17. Introduction length
The introduction is ~140 lines (L68–144), which is very long for a finance paper. The three-paragraph structure of "our results suggest..." (L122–133) could be shortened since the same messages are repeated in the interpretation section. Consider trimming by ~30%.

---

### 18. Literature review placement
Section 2 (literature review) appears before the methodology, which is standard in some traditions but unusual for top finance journals where the contribution is expected to be clear before the literature. Consider whether integrating the literature comparisons into the methodology discussion (as is partly done in Section 3.4) and shortening Section 2 would improve flow.

---

### 19. Simulation section organization
The simulation study (Section 4.4) mixes five distinct experiments in a single subsection using `\paragraph{}` separators. Each experiment has its own figure but shares a common preamble. Converting these into numbered subsections (4.4.1–4.4.5) would improve navigability and make cross-referencing easier.

---

### 20. Missing discussion: why $q^5$ factors and not Fama-French?

Resolved.

---

## V. Missing Elements

### 21. No goodness-of-fit measure
The paper reports parameter estimates but never reports any measure of model fit (e.g., average pricing error, R² analogue, cross-validated MSE). Including the value of the minimized objective function $Q_n(\hat{\theta})$ across models and horizons would give readers a sense of how well each specification fits.

---

### 22. No bootstrap comparison
The paper adopts SHAC inference over bootstrap (DLP12's approach). But if the asymptotic SEs are "highly erratic" as acknowledged, a natural question is whether bootstrap SEs are more reliable. Even if computationally expensive, reporting a single bootstrap comparison would address the obvious referee question.

---

### 23. Online Appendix not reviewed
The paper references Online Appendix sections (e.g., `\ref{sec:proofs_asymptotics_H}`, `\ref{sec:appendix_bias_variance}`, `\ref{sec:dlp_simulation}`) which are in `sdf_paper_online_appendicies.tex`. Key regularity conditions and proofs are deferred there. Ensure the regularity conditions are stated precisely and verifiably.

---

## Summary Assessment

| Category | Count | Severity |
|---|---|---|
| Typos/notation errors | 4 | Low–Medium |
| Internal inconsistencies | 4 | Medium |
| Methodological weaknesses | 5 | Medium–High |
| Exposition issues | 4 | Low |
| Missing elements | 3 | Medium |

The paper's core theoretical contribution (the Compounding-Horizon Wedge) and empirical message (weak identification of alpha and multi-factor models) are clear and well-argued. The most significant vulnerability is the absence of any analytical or numerical characterization of the wedge magnitude (point 12), which a referee will likely request. The data inconsistency in the region filter (point 1) is a factual issue that should be addressed before submission.
