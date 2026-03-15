## Referee Report: "Semiparametric SDF Estimators for Pooled, Non-Traded Cash Flows"

### 1. Summary of the Paper

The manuscript develops a generalized semiparametric framework for estimating Stochastic Discount Factors (SDF) using pooled, non-traded cash flows, specifically tailored for private equity (PE). The authors introduce a Least-Mean-Distance (LMD) estimator that endogenizes the choice of the discounting date. The paper effectively bridges the inception-only NPV approach of Driessen et al. (2012) and the time-series approach of Korteweg and Nagel (2016).

A central theoretical contribution is demonstrating that evaluating cash flows beyond fund inception introduces a "Timing Risk Premium" (TRP) which shifts the population target parameter. Through simulation and an empirical application using Preqin data (Buyout and Venture Capital), the authors show that extending the compounding horizon creates a bias-variance tradeoff: it stabilizes finite-sample bias significantly, but shifts the target estimand. Empirically, they find that single-factor market models are robust, while multi-factor models (including alpha terms or $q^5$ factors) severely destabilize and collapse.

### 2. Overall Assessment

This is a strong, rigorously executed paper that tackles a pervasive econometric issue in the private markets literature: the fragility of semiparametric SDF estimates in finite, overlapping samples. The formalization of the TRP and the exact MSE decomposition (Proposition 6)  are excellent theoretical additions to the literature. Furthermore, the application of vintage-year increasing-domain asymptotics  is highly appropriate for the data generation process of PE funds. The "negative" finding—that current cross-sectional PE data cannot reliably support multi-factor identification or unconstrained alpha estimation —is a candid and valuable contribution that JFQA readers would appreciate.

However, to cross the threshold for publication, the manuscript needs to provide more constructive guidance on *how* applied researchers should navigate the tradeoffs it identifies.

### 3. Major Comments

* 
**Navigating the Bias-Variance Tradeoff in Practice:** You demonstrate that setting the horizon $H > 0$ stabilizes the estimator but targets a pseudo-true parameter $\theta_{H}^{\dagger}$ rather than the structural parameter $\theta_{0}$. You note that the MSE-minimizing horizon in your simulations is around 12.5 to 15 years.


* *Critique:* For an empiricist bringing a new dataset to this framework, how should they endogenously select the optimal $H$? You mention in the conclusion that determining this optimal configuration is an objective for future research. However, for this paper to be highly cited as a methodological benchmark, providing a data-driven heuristic (e.g., based on cross-validation or fund duration) for selecting $H$ would significantly elevate the paper's impact.

* 
**The Alpha Problem:** The empirical finding that introducing an alpha intercept causes the MKT loading to collapse is striking (e.g., MKT drops to 0.14 for BO under FW).

* 
*Critique:* Given that evaluating abnormal performance (alpha) is the primary goal of much of the PE literature , the conclusion that single-factor models without alpha are the only robust specification  leaves a gap. Can you suggest or test regularization techniques (like Ridge penalties on the alpha parameter) to prevent this collapse, rather than just bounding alpha between -2% and 2%?

* 
**Truncation of Recent Vintages:** You correctly identify that unseasoned funds with large reported NAVs mechanically inflate market beta. Your solution is to truncate the sample to vintages 2010 and earlier.


* *Critique:* While methodologically sound, this means your framework cannot evaluate the last 15 years of PE activity. You briefly mention not imposing a model-based residual-value adjustment. Acknowledging this limitation more prominently in the introduction, or briefly discussing how recent NAV-smoothing models could be integrated into the $Q_{n,H}(\theta)$ objective, would strengthen the paper.


* **Replication of Driessen et al. (2012):** In Appendix B, you show that even after aligning your simulation design with Driessen et al. (using shifted lognormal errors, etc.), you still observe substantial beta attenuation (0.93 vs. 0.98) and alpha bias.


* 
*Critique:* Claiming that a seminal paper's small-sample properties cannot be reproduced  is a strong statement. Ensure you are absolutely confident that no minor parameterization differences remain, and consider softening the language to focus on the *sensitivity* of their model rather than an outright failure to replicate.


### 4. Minor Comments

* 
**Notation Clarity:** The use of $Q_{n,H}(\theta):=-\frac{1}{n}\sum_{i=1}^{n}L(\overline{\epsilon}_{i}^{(H)}(\theta))$  is clean and standard. However, ensure that when you transition to empirical results (Tables 6-9), it is explicitly clear to the reader that the asymptotic standard errors $SE_{A}$ relate to the pseudo-true parameter $\theta_{H}^{\dagger}$, not $\theta_{0}$.


* 
**Cross-Validation Metric:** Using hv-block cross-validation to account for spatial dependence  is a highly practical addition. Emphasizing this as a core robustness tool for PE empiricists in the abstract could attract more readers.
