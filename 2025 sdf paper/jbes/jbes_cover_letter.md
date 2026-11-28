# Cover Letter for Submission to the *Journal of Business & Economic Statistics*

**Date:** November 28, 2026  
**To:** Joint Editors (Professors Yingying Fan, Michal Kolesár, and Dacheng Xiu), *Journal of Business & Economic Statistics*  
**From:** Christian Tausch (AssetMetrix GmbH, Munich; corresponding author: `christian.tausch@quant-unit.com`) and Alexander Bohnert (Munich University of Applied Sciences)  
**Manuscript Title:** *Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference*

---

Dear Editors,

We are pleased to submit our manuscript, *"Cash-Flow Pricing of Non-Traded Assets: Finite-Pool Identification, Horizon Selection, and Inference,"* for consideration for publication in the *Journal of Business & Economic Statistics*. The manuscript is original, has not been published elsewhere, and is not under consideration by any other journal.

### Scope and Econometric Contributions

Across the non-traded economy—from private companies (themselves pools of untraded corporate investment projects) and closed-end private capital funds (\$13+ trillion across buyout, venture capital, private real estate, private credit, infrastructure, natural resources, and impact) to illiquid structured-credit residuals (such as CLO equity tranches), specialized claim pools, and cohort-grouped direct real assets, whole loans, or human-capital contracts—periodic market prices are unobserved or noisy appraisals, and irregular cash flows are the primary record. Whether units are natively pooled or grouped into origination cohorts to dampen noise, estimating nonlinear discounting parameters directly from cash flows creates a finite-pool M-estimation problem fundamentally distinct from linear traded-return factor regressions. Our paper makes four methodological contributions tailored to JBES's econometrics and financial econometrics readership:

1. **Finite-Pool Identification Failure (Theorem 1, T1–T6):** We prove that minimizing squared pooled net present values (including the widely used inception criterion of Driessen, Lin, and Phalippou 2012) fails consistency at finite cohort sizes whenever parameter-dependent payoff variance shifts the population target even at inception ($H = 0$), and establish an exact factorization showing how post-inception horizon compounding (Tausch and Pietz 2024, eq. 4) creates an opposing Scale-Jacobian pull.
2. **Cross-Unit Pair-Product Estimator $U_{v,H}$ (Theorem 1, T7–T8):** Connecting to leave-one-out quadratic-form debiasing, U-statistic M-estimation, and conditional-moment methods (Angrist, Imbens, and Krueger 1999; Domínguez and Lobato 2004; Honoré and Powell 1994; Kline, Saggio, and Sølvsten 2020; Mikusheva and Sun 2022; Muandet, Jitkrittum, and Kübler 2020), we introduce the within-cohort cross-unit pair-product criterion $U_{v,H}$ ($q_v = 1$), a tuning-free $O(n_v)$ upgrade eliminating the idiosyncratic diagonal contribution at $H = 0$ with as few as two units per cohort when strict factor-path exogeneity, sign-aligned conditional means, and absent shared cohort shocks ($c_v = 0$) hold.
3. **Cross-Cohort Rank-1 Impossibility and Duration-Matched Horizon Selection (Proposition 2):** Under shared vintage or calendar shocks (which persist even as $n_v \to \infty$), we prove that extending *uninstrumented* leave-out pairing across distant cohorts fails structurally via rank-1 collapse and finite-$V$ indefiniteness, and show how forward-horizon compounding becomes a constructive lever: matching compounding sensitivity to payoff duration yields an exact duration-matching score root $H_0^\ast$ that cancels first-order shared-shock bias and restores coercivity, with global bounds, marginal-law AR(1) invariance, calendar-overlap conditions, and finite-$V$ MSE expansions.
4. **Primitive 1D Increasing-Domain and Identification-Robust Inference (Proposition 3, Corollary 1, and Monte Carlo Study):** Modeling overlapping vintages directly on the one-dimensional calendar-time lattice $y_v \in \mathbb{Z}$, we establish mixing and uniform moment inheritance from primitive bounded-lifetime DGP conditions (Assumptions A1–A5, Lemma D.1), establish valid small-sample (fixed-$b$) Bartlett HAC and identification-robust Anderson–Rubin (AR) asymptotics, and validate finite-sample behavior across 72 calibrated Monte Carlo designs plus an increasing-$V$ ablation study through $V = 1{,}600$ that resolves compounding-kurtosis studentization dynamics. An empirical application to North American Buyout and Venture Capital cash flows (383 BO and 432 VC funds) illustrates all theoretical mechanisms cleanly.

Given the paper's focus on nonlinear financial econometrics, asset-pricing identification, horizon selection, and identification-robust inference, we respectfully believe it aligns particularly well with the financial econometrics editorial track (e.g., Professors Dacheng Xiu or Andrew Patton), though we of course leave editorial assignment entirely to your discretion.

### Formal ASA / JBES Data-Exemption Request (Proprietary Preqin Data)

In accordance with the ASA/JBES Data and Code Availability Policy, we respectfully request an **exemption from publicly redistributing the raw fund-level private-equity cash-flow files**, which are proprietary commercial data licensed from Preqin (April 14, 2022 export; North American USD Buyout and Venture Capital funds, vintages through 2010) and are strictly restricted against redistribution under Preqin's licensing agreement.

To ensure full transparency and end-to-end reproducibility despite this commercial data restriction, our accompanying **anonymized Replication Code Archive** (`jbes_replication.zip`) provides:
- All **public U.S. $q^5$ factor data** (`global-q.org`) and **aggregate-only** vintage fund-count tables;
- Complete R scripts reproducing **100% of our 72-design Monte Carlo study**, Proposition 2 analytical lookup tables, and the $V = 25\text{--}1{,}600$ studentization ablation study (together with pre-computed public simulation summary CSVs);
- A standalone **synthetic Preqin-calibrated verification script** (`run_synthetic_verification.R`, runtime ~15 seconds) that generates calibrated Buyout and Venture Capital cash-flow panels and executes our entire single-factor, joint $(\alpha,\beta)$, calendar-vs.-index fixed-$b$ Bartlett inference, $N_{\mathrm{eff}}$, and Anderson–Rubin confidence-set pipeline end-to-end without requiring a Preqin license; and
- Exact Preqin data-ingestion and cleaning scripts (`prep_preqin.R`, `run_empirical_study.R`, and exhibit scripts) enabling any licensed Preqin user to replicate every empirical table and figure directly.

### Submitted Files

1. `jbes_main_blind.pdf`: Anonymized double-spaced main manuscript (42 pages inclusive of 3 figures, 5 tables, and references; 140-word abstract).
2. `jbes_supplement_blind.pdf`: Anonymized standalone Online Supplementary Material (Appendices A–F, 68 pages; all mathematical proofs, primitive verification lemmas, and extended Monte Carlo and empirical exhibits).
3. `jbes_titlepage.pdf`: Separate non-blind title page with author affiliations, acknowledgments, funding, and conflict-of-interest disclosure.
4. `jbes_replication.zip`: Anonymized, self-contained R replication archive with detailed `README.md` and synthetic end-to-end verification script.

Thank you very much for considering our manuscript for publication in the *Journal of Business & Economic Statistics*.

Sincerely,

**Christian Tausch** (Corresponding Author)  
AssetMetrix GmbH, Theresienhöhe 13, D-80339 Munich, Germany  
`christian.tausch@quant-unit.com` | +49 89 5432880333

**Alexander Bohnert**  
Munich University of Applied Sciences, Am Stadtpark 20, D-81243 Munich, Germany  
`alexander.bohnert@hm.edu`
