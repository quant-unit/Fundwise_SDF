# JBES Methodological Referee Revision Plan & Checklist

*Created:* 2026-11-27  
*Source report:* `2025 sdf paper/jbes/jbes_methodological_referee_report.pdf` (8 pp., dated 27 November 2026)  
*Baseline commit:* `bd6be4b` (branch `revision-2026`)  
*Companion log:* `2025 sdf paper/jbes/jbes_revision_plan.md` (Sections 0–8)

---

## 0. Executive Summary & Scope Boundaries

### Verdict on the Referee Report
- **Almost every mathematical critique, counterexample, and literature citation in the 8-page report has been independently verified and holds up.**
- **Zero re-estimation is required.** We audited `r_project/empirical/prep_preqin.R` (lines 294–305), `r_project/helper/pool_vintage.R` (lines 37–66), and `r_project/helper/objective.R` (lines 137–146): the R code does **not** double-weight fund size in FW estimation. The double-sizing error flagged on p. 2 is purely a notational contradiction in Online Supplement Appendix A Table A.1.
- **No new estimator or new Monte Carlo suite is required.** Comparators raised on pp. 6 & 8 (instrumented linear moments $\mathbb{E}[z_v PV_v(\theta_0)] = 0$, Domínguez & Lobato 2004, Muandet, Jitkrittum & Kübler 2020) are addressed theoretically in prose under matched assumptions.
- **Hard page budget:** Keep `jbes_main_blind.pdf` strictly **$\le 40$ double-spaced pages** (currently 40 pp.) and abstract **$\le 200$ words** (currently 198 words). Trimming the Brownian volatility special case, held-out $H=0$ squared-loss CV mention, and repetitive prose frees space for the new 2-fund lead-in, estimator-comparison table, and matched-vintage VC highlight.

### Master Progress Tracker

| Priority | Workstream | Primary Targets | Status |
| :---: | :--- | :--- | :---: |
| **P0** | Must-fix math, theorems, proofs & formal contradictions (6 items) | `sdf_paper.tex`, `sdf_paper_online_appendices.tex`, `mc_paper.R`, `step5_diagnostics.R` | `[x]` |
| **P1** | Attribution, units/scaling, inference qualifiers & empirical framing (7 items) | `sdf_paper.tex`, `sdf_paper_online_appendices.tex`, `sdf_paper_refs.bib`, `main_exhibits_u.R`, `step5_diagnostics.R` | `[x]` |
| **P2** | Structural hierarchy, 2-fund lead-in, summary table & abstract polish (4 items) | `sdf_paper.tex`, `sdf_paper_online_appendices.tex` | `[x]` |
| **P3** | Table-note regeneration, LaTeX compile, page-count, blind-review & archive sync (5 items) | `jbes_main_blind.pdf`, `jbes_supplement_blind.pdf`, `jbes_cover_letter.tex`, `jbes_revision_plan.md`, `build_jbes_replication.sh` | `[x]` |

---

## 1. Priority P0 — Must-Fix Mathematical & Formal Errors

### P0.1 — Delete False Proposition 2(iii) Efficiency Claim & Repair Eq. (17) MSE Expansion
*Referee report p. 1, p. 5 (two-fund counterexample).*

**Verified error:**
- Differentiating $U_{v,H} = A_{v,H}^2 U_{v,0}$ gives $s_{v,H} = A_{v,H}^2 s_{v,0} + 2A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$. Appendix C.11(iii) dropped the second term, which is random for finite $n_v$ and can reduce variance (in the two-fund, one-period example, $V_0 = \sigma^2/2 + \sigma^4 > V_{\text{later}} = \sigma^2/2 + \sigma^4/4$). Homoscedasticity alone also does not factor $\mathbb{E}[A_{v,H}^2 h_{v,0}]$ into $\mathbb{E}[A_{v,H}^2]\mathbb{E}[h_{v,0}]$.
- Eq. (17) optimizes the *asymptotic proxy* $B_0(H)^2 + V^{-1}\mathcal{V}(H)$. Actual finite-$V$ MSE contains $2B_0(H)b_1(H)/V$, whose derivative at $H_0^\ast$ adds $2B_0'(H_0^\ast)b_1(H_0^\ast)/V$ of the same $O(V^{-1})$ order and can change the sign (explaining why $H_{\mathrm{RMSE}}^\ast > H_{\mathrm{zb}}$ in Table AS6 Panel A Regimes 3a–3b).

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§3.1, Prop. 2(iii), lines 399–406): Delete the claim that $H_{\mathrm{MSE}}^\ast(V) = 0$ for every $V$ when $\sigma_{\mathrm{eff},v}^2 \equiv 0$ and $q_v = 1$. State explicitly that under $c_v = 0$, $H$ trades compounding-scale variance against the $2A(\nabla A)U_0$ score term (with $H=0$ minimizing RMSE in our stochastic-return core simulations).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Table 1, lines 435–438): Remove "$H_{\mathrm{MSE}}^\ast = 0$" and "$H=0$ minimizes variance under $q_v=1$" from Regime 1 cells; replace with accurate wording (zero population bias for all $H \ge 0$; $H=0$ minimizes RMSE in the core simulations).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Prop. 2(iii), eq. (17), lines 400–406): Label eq. (17) explicitly as the minimizer of the asymptotic proxy $B_0(H)^2 + V^{-1}\mathcal{V}(H)$, and add the exact finite-$V$ second-order bias term: numerator $\mathcal{V}'(H_0^\ast) + 2B_0'(H_0^\ast)b_1(H_0^\ast)$, noting that the sign can flip and that discrete grid optima require discrete comparison.
- [x] `2025 sdf paper/tex/sdf_paper_online_appendices.tex` (§C.11(iii), lines 455–476): Delete the false Cauchy–Schwarz proof completely; write out the exact score identity $s_{v,H} = A_{v,H}^2 s_{v,0} + 2A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$; derive both the asymptotic-proxy expansion and the $2B_0' b_1$ adjustment for actual MSE.
- [x] `r_project/simulation/mc_paper.R` (line 803, generating `table_AS6_optimal_H_regimes.tex`): Remove "confirming Proposition~\ref{prop:optimal_H}(iii)" from the Table AS6 note; state neutrally that Regime 1 $U_{v,H}$ achieves its empirical grid minimum at $H=0$.

---

### P0.2 — Restrict and Harden Proposition 2(i) (Duration Matching, Caps, Score Root & Brownian Model)
*Referee report p. 1, p. 4.*

**Verified gaps:**
- $D^\top \Sigma D > 0$ does not guarantee $T_{\mathrm{eff}} > 0$ (counterexample: dates $(1,2)$, innovation $(Z, -\frac{3}{4}Z)$ gives $D^\top \Sigma D = \mathrm{Var}(Z)/16 > 0$ and $T_{\mathrm{eff}} = -2$).
- For pair products, an individual cross-fund matrix $\operatorname{Cov}(e_{i,s}, e_{j,t}\mid \mathcal{G}_v)$ need not be symmetric in $(s,t)$; symmetry holds for the **pair-weighted average over ordered pairs** $i \neq j$.
- Under eq. (7), windows are capped at $T_i^{\max}$ (possibly vintage-specific), so continuous $2T_{\mathrm{eff}}$ can be unreachable, and integer months can bracket the root without hitting it.
- Differentiating $H_0^\ast$ w.r.t. $r$ holds $T_{\mathrm{eff}}$ fixed although its weights depend on $r$.
- $\nabla_\theta Q_H(\theta_0) = 0$ establishes a **score root (stationary point)**, not $\theta_{H_0^\ast}^\dagger = \theta_0$, without a global minimization / shape argument.
- The Brownian volatility proof switches from additive $G_u = 1 + r_{f,u} + \alpha + \beta^\top F_u$ to multiplicative $R_t(\alpha) = \exp(\alpha t + X_t)$ and factors out a fixed-duration noise weight; the deterministic $\beta$ proof requires constant $\bar{F} \neq 0$.

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§3.1, lines 361–392) & `sdf_paper_online_appendices.tex` (§C.11(i), lines 316–398):
  - Define $\Sigma_v^{(s,t)}$ for $U_{v,H}$ explicitly as the pair-weighted average $\frac{\sum_{i \neq j} \omega_{v,i}\omega_{v,j}\operatorname{Cov}(e_{i,s}, e_{j,t}\mid\mathcal{G}_v)}{1 - \sum_k \omega_{v,k}^2}$ (manifestly symmetric in $(s,t)$).
  - Add an explicit sufficient condition for $T_{\mathrm{eff}} > 0$ (e.g., nonnegative discounted covariance row sums $\sum_s D_{v,s}\Sigma_v^{(s,t)} \ge 0$ with strict positivity on a future date $t > t_v^0$, or state $T_{\mathrm{eff}} > 0$ as an explicit condition; note payment-timing shifts can produce negative off-diagonal covariances).
  - Separate **(a) General aggregate score identity** eq. (14) (retaining vintage-specific caps $\mathcal{T}_{v,H}$ and requiring an attainable-range condition for a root to exist) from **(b) Constant-return, uncapped continuous-time corollary** eqs. (15)–(16) (requiring constant $G = e^r$, fixed $T_{\mathrm{eff}} > 0$, and constant nonzero excess return $\bar{F} \neq 0$ for $\beta$).
  - State $\partial H_0^\ast / \partial r < 0$ explicitly **at fixed $T_{\mathrm{eff}}$** (and optionally note in Appendix C that when $\Sigma^{(s,t)} \ge 0$, $\partial T_{\mathrm{eff}}/\partial r = -\frac{1}{2}\mathrm{Var}_w(s+t) \le 0$ reinforces the decrease).
  - Call $H_0^\ast$ a **zero-score horizon** ($\nabla_\theta Q_{H_0^\ast}(\theta_0) = 0$) and state explicitly that $\theta_{H_0^\ast}^\dagger = \theta_0$ requires $\theta_0$ to be the global minimizer (e.g., unique stationary point / unimodality on $\Theta$), not merely A1 + A5.
  - Move the Brownian volatility-shrinkage result out of the main Proposition 2(i) statement into a clearly labeled **stylized multiplicative-shift, fixed-duration special case** (Remark/Corollary in Appendix C, brief 1-sentence pointer in main text), noting that empirical $\beta$ uses the actual additive sensitivities $\sum_u F_u/G_u$ in eq. (14).

---

### P0.3 — Fix Proposition 2(ii) (Calendar-Overlap Counterexample, Level vs. Slope & Sign Overclaims)
*Referee report p. 6.*

**Verified gaps:**
- Unconditional calendar overlap can *shorten* duration (counterexample: deals on $(0,4]$ and $(5,6]$ give static $T_{\mathrm{eff}} = 5.0$ vs. overlap-weighted $4.4$).
- Positive semidefiniteness gives $c_v c_w \ge c_{v,w}^2 \implies c_v > 0$, a statement about **levels**, not $\nabla c_v(\theta_0) \neq 0$ or its sign.
- Universal "upward inception tilt" / "downward Scale-Jacobian pull" language holds only under named models (scalar $\alpha$, or scalar $\beta$ with positive excess returns and $T_{\mathrm{eff}} > 0$), not universally or in coupled vector models.
- AR(1) target invariance requires holding the *entire within-vintage marginal law* fixed (not just scalar $c_v(\theta_0)$).
- Positive $c_{v,w}$ need not raise $\Lambda_H^{\mathrm{cross}}$ when $\mu_v'(\theta_0)\mu_w'(\theta_0) < 0$ (possible for $\beta$); the $O(1/V)$ bias has two signed terms; and a zero average score at $H_0^\ast$ does not force every vintage's quadratic score $s_{v,H}^{\mathrm{quad}}$ to vanish pointwise.

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Abstract lines 75–76; Intro lines 104, 106; §2.2; Prop. 2(ii) lines 393–398; Conclusion line 675): Replace universal "upward"/"downward" claims with exact derivative statements ($\nabla_\theta \sigma_{\mathrm{eff},v}^2(\theta_0)$ and $\nabla_\theta \log A_{v,H}(\theta_0)$), specifying the upward/downward direction explicitly under named models (scalar $\alpha$, positive-excess-return single-factor $\beta$ with $T_{\mathrm{eff}} > 0$) and in the empirical application.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Prop. 2(ii)(a)–(c), lines 394–398) & `sdf_paper_online_appendices.tex` (§C.11(ii), lines 400–454):
  - (a) State that PSD implies $c_v(\theta_0) > 0$, and that a nonzero inception score additionally requires $\mathbb{E}[q_v \nabla_\theta c_v(\theta_0)] \neq 0$.
  - (b) State AR(1) invariance under fixed **within-vintage marginal law**; replace the flawed unconditional Chebyshev claim for calendar shocks with the exact sufficient condition: **positive weighted covariance between exit date $b_2$ and expected active-life overlap** (noting this holds at fixed entry date and lengthens $T_{\mathrm{eff}}$ by ~2 yr in our Preqin deal calibration).
  - (c) Qualify $\Lambda_H^{\mathrm{cross}}$ sign by the sign of $\mathbb{E}[q_v q_w A_v^2 A_w^2 \mu_v'(\theta_0)\mu_w'(\theta_0) c_{v,w}(\theta_0)]$; note that $c_{v,w}$ enters both signed terms ($C_{hs}$ and $\Lambda_H$) of the $O(1/V)$ expansion; and restrict pointwise vanishing of $s_{v,H}^{\mathrm{quad}}$ at $H_0^\ast$ to homogeneous constant-return / common-duration settings.

---

### P0.4 — Fix Lemma C.1(a)–(b), Assumption A2 Singleton Contradiction & Merton (1980) Error
*Referee report p. 5, p. 7, p. 8.*

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper_online_appendices.tex` (Lemma C.1, lines 251–261):
  - Part (a): Add the first-cash-flow-origin edge case where only one fund has a nonzero cash flow at $t_v^0$: then every distinct-fund pair product vanishes at $t_v^0$, so $U_{v,0}(\alpha,\beta) \to 0$ as $\alpha \to +\infty$ while $Q_0^{\mathrm{L2}}$ retains a strictly positive limiting floor.
  - Part (b): Replace "If $H > 0$, so that $\mathcal{T}_{v,H}$ contains a date $\tau^\ast$ after $t_v^{\mathrm{first}}$" with the exact requirement that $\mathcal{T}_{v,H}$ contains a valuation date $\tau^\ast$ strictly after $t_v^{\mathrm{first}}$ (for L2) and sufficiently late relative to contributing funds' first outlays — specifically $2\tau^\ast > t_{v,i}^{\mathrm{first}} + t_{v,j}^{\mathrm{first}}$ for at least one same-sign pair $(i,j)$ (for $U_{v,H}$). Add one sentence noting that tail divergence as $\alpha \to +\infty$ does not by itself force an interior minimizer inside a finite box $\alpha \in [-0.02, 0.02]$.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Assumption A2, line 485) & `sdf_paper_online_appendices.tex` (A2, line 581): Restrict $\inf_{v:\,n_v\ge 2}(1 - \sum_{i\in\mathcal{I}_v}\omega_{v,i}^2) \ge \underline{d} > 0$ explicitly to pair-eligible ($n_v \ge 2$) vintages for $U_{v,H}$, noting that L2 and the linear AR moment test have no pair denominator and therefore accommodate singleton vintages (all 27 VC vintages).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§5.3, line 657) & `r_project/results/step5_diagnostics.R` (lines 267–268, generating `table_alpha_beta_identification.tex`): Remove the false claim that Merton (1980) imposes an "irreducible" $2.54\%$/yr ($0.21\%$/mo) lower bound on $\hat\alpha$ precision (which contradicts Table 3 Panel A's $\mathrm{SD}(\hat\alpha) = 0.066\%$/mo). Either drop Merton (1980) or reframe it explicitly as a heuristic for *unconditional* mean-return estimation, not a precision bound conditional on realized factor paths.

---

### P0.5 — Tighten Section 2.1–2.2 Conditioning, Measurability, T7 Identification & Kernel Meaning
*Referee report p. 3.*

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.1, Assumption 1 / Lemma 1 / Prop. 1, lines 161–220) & `sdf_paper_online_appendices.tex` (§C.1, lines 98–113): State explicitly that random deal arrival dates $t_{i,j}^{\mathrm{Inv}}$, outlays $\mathrm{Inv}_{i,j}$, and deal participation/count are $\mathcal{F}_{t_{i,j}^{\mathrm{Inv}}}$-measurable (adapted stopping/participation structure with integrability) so term-by-term iterated expectations hold. State explicitly that the sample-end window cap $T^{\max}$ is administrative and predetermined (not outcome-dependent).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2, Conditions (P) and (X), lines 260–265):
  - For Condition (X) in the multiplicative deal example, condition on deal entry date $a$ and outlay $I_a$ as well ($\mathbb{E}[N \mid \mathcal{G}_v, a, I_a] = 1$, or state weighted orthogonality $\mathbb{E}[D_{t_v^0,a} I_a (N-1)\mid \mathcal{G}_v] = 0$ directly), since $\mathbb{E}[N\mid \mathcal{G}_v] = 1$ alone does not zero $\mathbb{E}[D_{t_v^0,a} I_a (N-1)\mid \mathcal{G}_v]$ if $(a, I_a)$ correlates with $N$.
  - For Condition (P) transferred to the empirical pool, state the pricing restriction explicitly for the **selected and weighted empirical population** (including per-commitment normalization, complete-history reporting eligibility, and FW fund sizes), noting non-informative selection/weighting or conditioning on the selected design.
  - For persistent cohort shocks, distinguish unconditional $\mathbb{E}[PV_v(\theta_0)] = 0$ from conditional $\mathbb{E}[PV_v(\theta_0)\mid \mathcal{F}_{t_v^0}] = 0$ and specify what belongs to the investor's inception information $\mathcal{F}_{t_v^0}$ (noting prior cohorts' cash flows are not yet realized at $t_v^0$).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Abstract line 77; Intro line 110; §2.2 motivating block line 256; Theorem 1(T7) lines 319–323) & `sdf_paper_online_appendices.tex` (§C.8, lines 235–241):
  - Clarify that dropping the diagonal proves an expectation identity (not automatically global minimization under heterogeneous $\mu_{v,i}(\theta)$, where opposing slopes can make $\theta_0$ a local maximum).
  - State the generalized sufficient condition $\mu_{v,i}(\theta) = a_{v,i}(\theta)\mu_v(\theta)$ with $a_{v,i} > 0$, write out the local Hessian $2\mathbb{E}[q_v A_{v,H}(\theta_0)^2 J_v J_v^\top]$ ($J_v = \nabla_\theta \mu_v(\theta_0)$) requiring cross-vintage spanning of $J_v$, and explain the economic source of identifying variation (vintages exposed to distinct realized factor paths separate $\alpha$ from $\beta$, illustrated via a crisp two-state one-period note).
  - Note clearly that an unbiased population objective (zero score at $\theta_0$) can still have a biased finite-$V$ minimizer.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.1, line 182): Narrow the incomplete-markets sentence: $\Psi_t(\theta)$ is a private-claim pricing kernel *only on claims for which the pricing restriction holds* (incompleteness alone does not establish those conditions); refer to $\beta$ as the **loading of the fitted cash-flow benchmark** (noting that equating it to the beta of a traded total-return series requires an additional return-generating model).

---

## 2. Priority P1 — Attribution, Units/Scaling, Inference & Empirical Framing

### P1.1 — Fix Appendix A Table A.1 Units, DLP Normalization & Pair-Weighting Economics
*Referee report p. 2.*

**Verified against R code:** `prep_preqin.R` scales FW fund cash flows by $s_{v,i}$ (USD mn), and `pool_vintage.R` / `objective.R` divide by $\sum_k s_{v,k}$ and $(\sum_k s_{v,k})^2 - \sum_k s_{v,k}^2$ — so the mathematical inputs to both L2 and $U_{v,H}$ are **per-commitment** cash flows weighted by $\omega_{v,i} = s_{v,i}/S_v$ (no double sizing in code).

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper_online_appendices.tex` (Appendix A Table A.1, lines 40–45): Fix the unit contradiction! Define $CF_{v,i,t}$ and $PV_{v,i}(\theta)$ uniformly as **per unit of commitment** in the mathematical expressions (both EW and FW), and state explicitly how dollar fund cash flows $C_{v,i,t} = s_{v,i} CF_{v,i,t}$ map to per-commitment vintage-pool NPV $PV_v^{\mathrm{FW}} = S_v^{-1}\sum_i C_{v,i} = \sum_i \omega_{v,i} PV_{v,i}$ ($S_v = \sum_i s_{v,i}$) and to the vintage criteria.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2 line 238, §3.4 line 535) & `r_project/results/main_exhibits_u.R` (line 100): Qualify the DLP (2012) nesting statement precisely — DLP's benchmark (p. 533) sums dollar cash flows within each vintage FoF, so squaring DLP's dollar NPV scales each vintage term by $S_v^2$ relative to our per-commitment FW L2; our per-commitment normalization recovers DLP's inverse-FoF-size weighting (DLP Table 7, p. 533) or equal vintage commitments $S_v$.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2, after eq. (11), around line 280): Add a concise remark on the economic meaning of FW pair weighting in $U_{v,H}$: with $n_v = 2$, size weights cancel completely ($U_{v,H} = \bar\epsilon_{1,H}\bar\epsilon_{2,H}$); for $n_v > 2$, fund $i$'s marginal pair weight is proportional to $\omega_{v,i}(1-\omega_{v,i})$ rather than $\omega_{v,i}$, so under heterogeneous exposures "fund-size weighted" differs between $U_{v,H}$ and pooled L2.

---

### P1.2 — Cite Tausch & Pietz (2024) Eq. (4), Credit DLP (2012), & Narrow Prop. 2(iv) Comparator
*Referee report p. 2, p. 6, p. 8.*

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2, right where eq. (8) defines horizon averaging, lines 233–238, and §1 Related Literature): Cite `\citet{TP24}` (strictly third-person to preserve double-blind review: e.g., *"Prior empirical work by \citet{TP24} (eq.~4) averages pricing errors over compounding dates before squaring to mitigate small-sample bias and variance; Theorem~\ref{thm:finite_pool} and Proposition~\ref{prop:optimal_H} establish what that horizon-averaged criterion targets and when it removes or merely shifts bias."*). Double-check author first name in `sdf_paper_refs.bib` line 161 (`Pietz, Marcus` — verified correct).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§1 Intro / §2.2): Explicitly credit \citet[p.~517 fn.~8; Table~1 and p.~522]{DLP12} in the main text for flagging the large-$\alpha$ attenuation pathology and upward-$\alpha$/downward-$\beta$ bias in high-volatility joint designs, clearly delineating our contribution (finite-pool population target, pair-product debiasing, and failure/horizon analysis under shared shocks).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (Prop. 2(iv) lines 408–410, Table 1 caption line 452, §1 Related Literature lines 128–131) & `sdf_paper_online_appendices.tex` (§C.11(iv), lines 478–506):
  - Narrow the rank-1 collapse and impossibility claim strictly to the **uninstrumented** cross-vintage pair-product construction $\widetilde{U}_{H,m}^{\mathrm{cross}}$.
  - Note two formal precision details: (i) an AR(1) covariance remains nonzero at every finite lag $m$ (its normalized contribution vanishes asymptotically under summability); (ii) a scalar restriction with rank-1 Hessian destroys *regular local identification* for $p \ge 2$ (though an isolated singular root can theoretically exist).
  - Add a matched-assumptions comparator discussion of **instrumented linear moments** $\mathbb{E}[z_v PV_v(\theta_0)] = 0$ (inception-measurable $z_v$ under (P); public-path functions $z_v = h(\mathcal{F}^{\mathrm{pub}})$ under (X) that retain validity with mean-zero shared shocks and can supply full-rank Jacobians, alongside their practical trade-offs: ad hoc instrument choice, multiple/empty roots in nonlinear discounting, and weak relevance at $V \approx 25$) and **conditional-moment / kernel U-statistic methods**.
  - Update Table 1 caption (line 452) so it no longer says AR is the "only" valid $H=0$ check under Regimes 2–3 (say *"uninstrumented $H=0$ check under Condition~(P) alone"*).
- [x] `2025 sdf paper/tex/sdf_paper_refs.bib`: Add verified BibTeX entries for:
  - Domínguez, M. A. & Lobato, I. N. (2004), "Consistent Estimation of Models Defined by Conditional Moment Restrictions," *Econometrica*, 72(5), 1601–1615.
  - Muandet, K., Jitkrittum, W., & Kübler, J. M. (2020), "Kernel Conditional Moment Test via Maximum Moment Restriction," *Proceedings of the 36th Conference on Uncertainty in Artificial Intelligence (UAI)*, PMLR 124, 41–50.

---

### P1.3 — Tighten Inference Scope, Bandwidth Qualifier & Horizon-Invariance Diagnostic Framing
*Referee report p. 7.*

**Checklist:**
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§3.3 lines 502–530, §5.2 line 639) & `sdf_paper_online_appendices.tex` (§D.1–D.4, Remark D.1 lines 691–696, Remark D.2 lines 763–768):
  - **Per-vintage centering under T8:** Explicitly note in §3.3, Remark D.1, and §5.2 that per-vintage centering $\mathbb{E}[s_{v,H}(\theta_H^\dagger)] = 0$ can fail under shared shocks (T8) just as under heterogeneous L2 (average score vanishes at $H_0^\ast$ while individual vintage expected scores remain nonzero). Qualify/label any $p$-values attached to the $H=12.5$~yr VC estimates accordingly (in addition to persistent-shock over-rejection).
  - **Fixed-$b$ finite-sample approximation:** State explicitly right after eq. (20) (line 511) that $S = \lceil 1.3\sqrt{V}\rceil$ implies $S/V \to 0$ asymptotically, so evaluating Corollary 1's fixed-$b$ critical values at the realized ratio $b = S/V$ is a finite-sample approximation.
  - **Trending $n_v$ wording:** Change "which our pools violate" (line 515) to "which trending pool sizes threaten" (matching the more accurate phrasing in Appendix D line 606).
  - **Gaussian-factor MC domain:** In A1 ([sdf_paper_online_appendices.tex](file:///Users/christausch/github/Fundwise_SDF/2025%20sdf%20paper/tex/sdf_paper_online_appendices.tex#L579)), replace "so the theory covers them only conditionally on that domain" with a clean statement that Gaussian-factor MC designs operate outside A1's bounded-support domain (conditioning on a realized positivity floor does not prove an unconditional increasing-domain limit).
  - **Horizon-invariance test:** Present $t_{\mathrm{HI}}(H_1)$ explicitly as a **specification diagnostic** testing the necessary implication $\beta_0^\dagger = \beta_{H_1}^\dagger$ of T7 (noting common bias at both horizons is invisible, consistency requires separated targets $\beta_0^\dagger \neq \beta_{H_1}^\dagger$ and a shrinking contrast SE, and power is only 33.2% even at VC-scale shocks).

---

### P1.4 — Reframe Empirical Interpretation (VC Sensitivity Pair, Matched Gaps, AR vs. Wald, $(\alpha,\beta)$, Fix-$A$, CV)
*Referee report p. 2, p. 7, p. 8.*

**Checklist:**
- [x] **Mark all single-factor results "$\alpha$ fixed at zero":** Update Abstract, Intro, §5.2, Conclusion, and table captions/notes in `r_project/results/main_exhibits_u.R` (Table E1), `r_project/results/horizon_invariance_tables.R` (Table D.1), and `r_project/results/table_published_estimates.R` (Table F.12).
- [x] **Drop VC "band" / bracketing language:** Across Abstract (line 79), Intro (line 124), §5.2 (line 639), and Conclusion (line 678), present VC inception estimates ($1.493$/$1.742$) and calibrated $H=12.5$~yr estimates ($1.382$/$1.583$) side by side as **sensitivity estimates under two identifying regimes** (one static shared-shock calibration), stating explicitly that heterogeneity or path-exogeneity failure could move both estimates to the same side of truth — neither brackets the structural $\beta$.
- [x] **Highlight matched 23-vintage VC L2–$U$ contrasts & qualify DLP $2.73$ comparison:** In §5.2 and §5.3 (lines 637, 666), place the matched 23-vintage VC L2–$U$ contrasts ($+0.435$ EW, $+0.551$ FW under identical $\alpha=0$, data, and weighting) next to the actual $U$ estimates to carry the diagonal-debiasing evidence, and qualify the DLP $2.73$ reference explicitly as coming from a joint $(\alpha, \beta)$ + jointly-estimated final-MV specification on TVE 1980–1993 data (so the difference vs. $1.49/1.74$ cannot be attributed to diagonal correction).
- [x] **Clarify 27 vs. 23 VC vintages, AR sets vs. Wald tests, and $N_{\mathrm{eff}}$:**
  - State clearly in §5.2 and Table 4 notes that all-fund VC AR sets and full-sample L2 use 27 vintages while $U_{v,H}$ uses the 23 vintages with $n_v \ge 2$.
  - State plainly that under the weak inception-mean restriction (Condition P alone), Buyout AR sets ($[0.949, 3.994]$ EW, $[0.805, 3.994]$ FW) rule out low loadings while extending to the search ceiling, whereas VC AR sets ($[-0.244, 3.994]$ EW, $[-0.549, 3.994]$ FW) include both $0$ and $1$ and give little direction — explicitly distinguishing this from the stronger-model (T7) VC Wald rejections of $\beta=0$.
  - Clarify that $N_{\mathrm{eff}} = 4.5\text{--}8.9$ measures **score concentration**, not independent degrees of freedom, and keep "large point estimate", "non-rejection of $\beta=1$", and "precise estimate" distinct.
- [x] **Tighten $(\alpha,\beta)$ framing ("ridge / poor identification"):** In Abstract (line 79), Intro (line 126), §5.3 (lines 648–665), and Conclusion (line 679), describe joint $(\alpha,\beta)$ and second-factor profiles as showing a **ridge, collinearity, and sensitivity to bounds and cohorts** (rather than asserting formal weak-identification asymptotic sequences from sampling correlation alone); keep formal "weak-identification-robust" strictly for AR moment inversion.
- [x] **Qualify Fix-$A$ counterfactuals & held-out $H=0$ squared-loss CV:**
  - `2025 sdf paper/tex/sdf_paper.tex` (§5.2, line 633) & `r_project/results/step5_diagnostics.R` (line 190, `table_fixA_jacobian.tex`): State that the exact T5 decomposition applies to the score at a common $\theta$, whereas differences between separately reoptimized estimators ($\hat\beta_H - \hat\beta_H^{\mathrm{Fix}A}$) are reference-dependent counterfactuals (evaluated at $\bar\theta = \hat\beta_{U,0}$).
  - `2025 sdf paper/tex/sdf_paper.tex` (§5.1, line 616) & `sdf_paper_online_appendices.tex` (§F.1, line 857): Either remove the held-out $H=0$ squared-loss CV clause from main §5.1 (saving lines!) and qualify it in Appendix F as a purely descriptive predictive-fit check for L2 inception loss (explicitly noting that L2 inception loss has a shifted population target $\theta_0^\dagger \neq \theta_0$).

---

## 3. Priority P2 — Structural Hierarchy & Exposition Polish (JBES Positioning)

### P2.1 — Lead with the 2-Fund, 1-Vintage Choice & Add Estimator-Comparison Table
*Referee report p. 8 (Priority 1).*

- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2, lines 240–257): Upgrade the motivating one-deal algebra to a **one-vintage, two-fund contrast** comparing the squared average NPV $(\frac{1}{2}(PV_1 + PV_2))^2$ against the pair product $PV_1 PV_2$, showing transparently what each targets as vintages accumulate at fixed $n=2$ (connecting directly to $R^\dagger = R_0(1+\sigma^2/n)$ vs. exact cancellation).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2 or §3.2): Add a compact summary table comparing **L2 (squared pooled NPV)**, **$U_{v,H}$ (cross-fund pair product)**, and **Inception Mean / AR** across four columns: *Moment/criterion*, *Conditioning information*, *Separation/identification restriction*, and *Population target* (reconciling cleanly with Table 1 so total float count stays within JBES budget).

### P2.2 — Usable Hierarchy of Results & Inputs Determining $H$
*Referee report p. 8 (Priority 2).*

- [x] `2025 sdf paper/tex/sdf_paper.tex` (§2.2 & §3.1): Emphasize T1 and T7–T8 as the central algebra/identification results; present eq. (14) as the general horizon score and eqs. (15)–(16) as a restricted constant-return corollary; state explicitly which unobserved inputs determine $H_0^\ast$ (latent covariance $\Sigma_v$, deal timing, return process) — noting a researcher observing only fund totals needs a model or estimates for those inputs, and that an estimated $\hat H$ introduces pre-estimation uncertainty if plugged into inference.

### P2.3 — Empirical 3-Column Readout, Simulation Question-Headings & Direct Novelty Positioning
*Referee report p. 8 (Priorities 3 & 5).*

- [x] `2025 sdf paper/tex/sdf_paper.tex` (§5.2): Present single-factor empirical findings cleanly in the sequence **inception estimate ($\alpha=0$) / calibrated $H=12.5$~yr sensitivity / inception AR set**, with matched 23-vintage VC L2–$U$ gaps prominently beside the $U$ estimates.
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§4.1–4.4): Give simulation subsection titles crisp economic/statistical questions (e.g., *Does $U_{v,0}$ remove finite-pool bias without inflating size?*, *When does moderate $H>0$ rescue joint $(\alpha,\beta)$ under shared shocks?*, *Why does finite-$V$ studentization converge slowly?*).
- [x] `2025 sdf paper/tex/sdf_paper.tex` (§1 Related Literature, lines 128–131): Position the statistical novelty directly around the three distinctive features: repeated fund observations under a **common factor path**, **bounded pool sizes**, and **parameter-dependent horizon revaluation**, contrasting cleanly with Tausch & Pietz (2024) eq. (4) and instrumented conditional-moment U-statistics (Domínguez & Lobato 2004; Muandet, Jitkrittum & Kübler 2020).

### P2.4 — Rewrite Abstract Conclusion ($\le 200$ Words)
*Referee report p. 8 (Suggested abstract conclusion).*

- [x] `2025 sdf paper/tex/sdf_paper.tex` (Abstract, lines 72–80): Tighten the abstract so it (a) avoids universal "upward/downward" claims and clarifies T7's separation/common-mean condition, (b) removes "only weakly identified" in favor of "imprecisely identified", and (c) adopts a crisp conclusion aligned with p. 8:
  > *"With alpha fixed at zero, inception estimates of market exposure are 1.15--1.20 for Buyout and 1.49--1.74 for Venture Capital; at a 12.5-year horizon motivated by a shared-shock calibration, Venture Capital estimates are 1.38--1.58. Inception mean-moment tests rule out low Buyout exposures but include both zero and unity for Venture Capital, while abnormal performance and second-factor loadings remain poorly identified."*
  Verify total word count stays strictly $\le 200$ words.

---

## 4. Priority P3 — Regeneration, Compilation, QA & Sync

- [x] **P3.1 — Regenerate affected table `.tex` files (notes/captions only, zero estimation reruns):**
  - Run `Rscript results/main_exhibits_u.R`, `Rscript results/step5_diagnostics.R` (or targeted note update if `step5_diagnostics.R` runs bootstrap — check runtime first so we don't accidentally trigger slow recomputation!), `Rscript results/horizon_invariance_tables.R`, `Rscript results/table_published_estimates.R`, and `Rscript simulation/mc_paper.R` (note-only update on Table AS6). Confirm all numerical entries in CSVs/tables are bit-identical before and after.
- [x] **P3.2 — Recompile both PDFs cleanly:**
  - Build `jbes_main_blind.pdf` and `jbes_supplement_blind.pdf` via full `pdflatex` + `bibtex` cycle.
  - Verify **0 undefined references (`??`)**, **0 BibTeX errors**, **0 broken cross-document `xr-hyper` refs**, and **main blind PDF strictly $\le 40$ double-spaced pages**.
- [x] **P3.3 — Double-blind anonymization check:**
  - Run anonymization grep across `.tex` and `.pdf` files, verifying zero identifying leaks while explicitly allowing the standard third-person bibliographic citation to `\citet{TP24}` (Tausch & Pietz 2024).
- [x] **P3.4 — Sync cover letter, revision log, and replication archive:**
  - Update `2025 sdf paper/jbes/jbes_cover_letter.tex` & `.md` (final page count, abstract word count, refined methodological summary).
  - Log all commits and decisions in `2025 sdf paper/jbes/jbes_revision_plan.md` Section 8 (and link to this file).
  - Run `2025 sdf paper/jbes/build_jbes_replication.sh` so updated R table-note scripts propagate into `jbes_replication.zip` and pass its leak scan.
- [x] **P3.5 — Final visual & numerical spot-check:**
  - Visually inspect modified main-text pages, new tables, and revised Proposition 2 / Appendix C–D displays in the compiled PDFs before committing.

---

## 5. Round-2 Focused Revision Pass (2026-11-27, after `a4d444e`)

Addressing the **Second-Round Referee Report** (`2025 sdf paper/jbes/archive/jbes_second_round_referee_report.pdf`, verdict downgraded to *"A focused revision before submission"* with all major Round-1 objections signed off; zero re-estimation required; full follow/pushback rationale documented in `2025 sdf paper/jbes/archive/jbes_referee_response_decisions.md`):

- [x] **R2.1 — Population targets in Table 1, Intro, Conclusion, Lemma C.2, Prop. 2(iv), and Appendix B:**
  - **L2 row & Lemma C.2:** Replaced "whenever $n_v^{-1}\sigma_{\nu,v}^2+c_v$ varies with $\theta$" in Table 1 with the exact sufficient displacement condition for interior $\theta_0$ ($\nabla_\theta\mathbb{E}[A_{v,H}^2\sigma_{\mathrm{eff},v}^2]|_{\theta_0}\neq 0$), revised T3 accordingly, changed Intro L104 from "requires $n_v\to\infty$" to diffuse within-vintage infill ($\sum_i \omega_{v,i}^2\to 0$) as a *sufficient* route absent shared shocks, and updated Lemma C.2 (and Table B.1) to require Herfindahl weight diffusion $h_v = \sum_i w_{v,i}^2 \to 0$ (noting $n_v\to\infty$ alone can leave $h_v \ge 1/4$ under FW, and $O(1/n_v)$ requires a rate condition $h_v = O(1/n_v)$).
  - **Pair products $U_{v,H}$ row & Appendix B:** Replaced "when $c_v > 0$" in Table 1 with nonzero vintage-average shared-noise score $\nabla_\theta\mathbb{E}[q_v A_{v,H}^2 c_v]|_{\theta_0} \neq 0$, carried T7's global separation condition ($\mu_v(\theta)\neq 0$ on positive weight for $\theta\neq\theta_0$) into Table 1, Intro contribution 2, and Conclusion, and updated Supplement Appendix B (`sdf_paper_appendix_sdf.tex`) to cite Theorem 1(T7) directly.
  - **Inception mean / AR row & Prop. 2(iv):** Stated the population target in Table 1 as the set of roots of $\mathbb{E}[PV_v(\theta)]=0$ (one or several for scalar $\beta$) and eq. (19) as providing nominal pointwise coverage of $\theta_0$; softened Prop. 2(iv) (main & supplement & Table C.1) to "can have multiple roots", noting that $A_{v,H}(\beta)$ at $H>0$ can shift or eliminate the one-period growth-optimal root $\beta^\ast$.
- [x] **R2.2 — Two-fund worked example (§2.2) & T3 proof (§C.4):**
  - Named both return conventions in §2.2: per-period $R^\dagger = R_0(1+\sigma_N^2/n)^{1/T}$ via $z=(R_0/R)^T$, versus cumulative holding-period $R^\dagger = R_0(1+\sigma_N^2/n)$ (or $T=1$).
  - Conditioned the worked example on the full public path $\mathcal{F}^{\mathrm{pub}}$ (or $\mathcal{F}_{\max(T,K-1)}^{\mathrm{pub}}$) so $A_H(\theta)$ exits the expectation even when $K=2T+1 > T$.
  - Updated Supplement §C.4 proof of T3 to cite Equation (9) with the $n=2$ coefficient $\sigma_N^2$ (and $(2\sigma_N^2/n)$ for general $n$) rather than the old single-deal $2\sigma_N^2$.
- [x] **R2.3 — Horizon $H=12.5$~yr sources, Figure 1B axis, Tables C.2–C.3, and Prop. 2(ii) polish:**
  - Relabeled Figure 1B vertical axis (`r_project/simulation/optimal_H_theory/figure_optimal_H_theory.tex` L126) from *"Zero-bias horizon"* to *"Zero-score horizon $H_0^\ast$ (yr)"*, recompiled `figure_optimal_H_theory.pdf`, and updated Supplement Tables C.2–C.3 captions (and Appendix D L827) to *"continuous-time zero-score horizon"*.
  - Distinguished the two sources behind $H=12.5$~yr in §3.1 and §5.2: the stylized static-shock duration calculation ($12.65$~yr) serves as motivation, while the VC shared-shock Monte Carlo design ($-0.008$ bias for additive-model $U_{v,12.5}$) supplies calibration evidence. Left the approved Abstract wording unchanged.
  - Tightened Prop. 2(ii) (main & Appendix C.11(ii)): restricted raw overlap length $|[a_1,b_1]\cap[a_2,b_2]|$ to independent calendar innovations $e_u$ and equal outlays (with the double autocovariance sum $\Gamma_e$ under serial correlation), clarified common-entry monotonicity gives nonnegative covariance (strict positivity requires exit-date variation altering expected overlap), and qualified §C.11(ii)(c) so positive linear-score covariance contributes positively to the *linear component* of $\Lambda_H^{\mathrm{cross}}$ (signing total $\Lambda_H$ requires sign restrictions on quadratic and mixed terms).
- [x] **R2.4 — Prop. 2(iii) local vs. global optimum, missing $q_v$ in eq. (C.3), and normalized criterion expectation:**
  - Stated eq. (17) (main & Appendix C.11(iii)) as the **local** MSE minimizer near a designated simple zero $H_0^\ast$ (global if $H_0^\ast$ is the unique globally separated zero on a compact horizon domain; otherwise competing roots are ranked by $\mathcal{V}$) and carried the uniform differentiability requirement of the $O(V^{-1})$ expansion directly into Prop. 2(iii) for actual finite-$V$ MSE.
  - Restored the missing $q_v$ factor in $2q_v A_{v,H}(\nabla_\theta A_{v,H})U_{v,0}$ in both Prop. 2(iii) and eq. (C.3) (noting the two-fund counterexample uses baseline $q_v = 1$).
  - Replaced the unconditional pull-out of $(\sum_v q_v)^{-1}$ in Appendix C.11(ii)(a) with $\mathbb{E}[\widehat{U}_H(\theta)\mid\mathcal{W}]$ and unconditional $\sum_v \mathbb{E}[\tilde{q}_v A_{v,H}^2(\mu_v^2 + c_v)]$ (noting deterministic $q_v \equiv 1$ exits directly).
- [x] **R2.5 — Matched 23-vintage VC L2–$U_{v,0}$ attribution qualification (partially pushed back; see `archive/jbes_referee_response_decisions.md` §2.1):**
  - Added an exact remark right after eq. (11) noting that moving from $L_v$ to $U_{v,H}$ drops own-fund diagonal squares and rescales the surviving off-diagonal products by $(1-h_v)^{-1}$ ($2$ at $n_v=2$ vs. $10/9$ at $n_v=10$ under EW, analogous size-concentration scaling under FW)—preserving identical unit weight on the common vintage signal $\mu_v^2 + c_v$ relative to L2 while scaling realized cross-fund noise and within-vintage heterogeneity in small/concentrated cohorts. Phrased $+0.435$/$+0.551$ in §5.2 and §5.3 as the *"matched 23-vintage difference between pooled L2 and the normalized pair estimator $U_{v,0}$."*
- [x] **R2.6 — Remaining theorem-alignment qualifications:**
  - **Conditional means:** Defined $\mu_v(\theta) := \mathbb{E}[PV_v(\theta)\mid\mathcal{G}_v]$ as the conditional mean of the weighted pool in T4 and Appendix C.5 (making T4 exact without equal fund means away from $\theta_0$); phrased $\mathbb{E}[N\mid\mathcal{G}_v,a,I_a]=1$ as *sufficient* for Condition~(X); made explicit in Condition~(P) that eligibility based on future reporting imposes an additional pricing restriction on the selected cohort.
  - **Coercivity & Remark C.1(iii):** Qualified §5.3 and Conclusion coercivity via Lemma C.1(b), reported factually that interior VC $(\alpha,\beta)$ fits occur at $H = 5, 15$, and $30$~yr (tail growth alone cannot force an interior minimum in $[-0.02,0.02]$), updated Fig. C.1 and Fig. F.2 (`step5_diagnostics.R`) notes, and clarified in Remark C.1(iii) and §3.2 that eq. (18)'s fixed reference-parameter instrument $z_{v,H} = A_{v,H}(1)/A_{v,0}(1)$ does not grow with candidate $\alpha$ and still yields a vanishing unstudentized moment under a predetermined zero-CF origin.
  - **Instruments & DLP FW scope:** Added finite-range ($m$-dependent) / explicit orthogonality condition for distant-vintage instruments after Lemma C.3 (noting AR(1) retains correlation across every fixed gap); limited the exact DLP (2012) inverse-FoF-size dollar-portfolio equivalence strictly to FW in §2.2, §3.4, Table A.1, and `main_exhibits_u.R` notes (noting EW averages per-dollar fund NPVs equally and defines a different portfolio when sizes differ).
- [x] **R2.7 — Rebuild, QA & archive refresh:**
  - Rebuilt `figure_optimal_H_theory.pdf`, re-ran `main_exhibits_u.R` (zero diff in `exhibits_data.csv`), recompiled `jbes_main_blind.pdf` (**40 pages**, 3 figures, 5 tables, 0 warnings, 0 >10pt overfull hboxes, 184-word abstract, clean blind scan) and `jbes_supplement_blind.pdf` (**67 pages**, 0 >10pt overfull hboxes, clean blind scan), and refreshed `jbes_replication.zip`.

---

## 6. Round-3 Pre-Submission Polish Pass (Completed, 2026-11-28, after `971680e`)

Addressing the **Third-Round Referee Report** (`2025 sdf paper/jbes/archive/jbes_third_round_referee_report.pdf`, 3 pp., GPT 6.2 Astra; closes the main objection to $U_{v,H}$, retains the VC application as written, and signs off on keeping T1–T8 unified) plus the author's **fixed-$b$ positioning decision** (zero re-estimation or R code changes required; full independent algebraic verification and decision rationale archived in `2025 sdf paper/jbes/archive/jbes_referee_response_decisions.md` §5):

- [x] **R3.1 — Fixed-$b$ positioning in Abstract, Keywords, and Cover Letter (Author TODO):**
  - Kept fixed-$b$ Bartlett inference intact throughout §3.3, §4, §5, and Online Supplement Appendix D (essential finite-$V$ craft for $V = 23\text{--}27$ overlapping vintages).
  - `2025 sdf paper/tex/sdf_paper.tex` (Keywords, L82): Dropped `Fixed-$b$ inference` and replaced with `Horizon selection` (matching the title and Proposition 2 while avoiding HAR-asymptotics referee misrouting).
  - `2025 sdf paper/tex/sdf_paper.tex` (Abstract, L77): Replaced *"and develop one-dimensional fixed-$b$ and identification-robust inference"* with *"and establish valid small-sample (fixed-$b$) and identification-robust inference."* Verified total abstract word count is **186 words** ($\le 200$).
  - `2025 sdf paper/jbes/jbes_cover_letter.tex` & `jbes_cover_letter.md`: Aligned inference phrasing, editorial-track sentence (focusing on nonlinear financial econometrics, asset-pricing identification, horizon selection, and identification-robust inference), date (November 28, 2026), and final page/word counts.

- [x] **R3.2 — Table 1 L2 variance conditioning & Abstract/Conclusion wording (Astra R3, p. 2 — verified exact):**
  - `2025 sdf paper/tex/sdf_paper.tex` (Table 1, L331): Split the L2 population-target condition cleanly by information set and vintage averaging: at $H=0$ under (P), vintage-avg. $\nabla_\theta\mathbb{E}[\operatorname{Var}(PV_v\mid\mathcal{F}_{t_v^0})]\big|_{\theta_0} \neq 0$ (T3); at any $H \ge 0$ under (X) and T4–T6, vintage-avg. $\nabla_\theta\mathbb{E}[A_{v,H}^2\sigma_{\mathrm{eff},v}^2]\big|_{\theta_0} \neq 0$. Updated T3 (L290) and Supplement §C.4 (L162) to distinguish $\mathcal{F}_{t_v^0}$- vs. $\mathcal{G}_v$-conditioning explicitly; kept $U_{v,H}$ and AR rows unchanged.
  - `2025 sdf paper/tex/sdf_paper.tex` (Abstract L75, Conclusion L657) & `jbes_cover_letter.tex/.md`: Changed *"whose derivative tilts/shifts the population target"* to *"whose derivative can shift the population target"*.

- [x] **R3.3 — Proposition 2(ii)(b) nonnegative calendar autocovariance condition (Astra R3, p. 2 — verified exact):**
  - `2025 sdf paper/tex/sdf_paper.tex` (Prop. 2(ii)(b), L408): Inserted `"with nonnegative calendar autocovariances"` into the common-entry parenthesis (matching Supplement §C.11(ii)(b) L433). Retained reported calibration shift unchanged.

- [x] **R3.4 — Finite-vintage $O(1/V)$ bias centering, competing-root ties, and two bookkeeping fixes (Astra R3, p. 3 — verified exact):**
  - **Bias centering & root ties:** Added `bao2007timeseries` (Bao & Ullah 2007, *JoE* 140(2), 650–669) to `sdf_paper_refs.bib`. Centered §C.11(ii)(c)'s stochastic $O(1/V)$ expansion on the finite-$V$ expected-criterion minimizer $\theta_{V,H}^\dagger$ plus $(\theta_{V,H}^\dagger - \theta_H^\dagger)$, stated the explicit matching rate condition $V^{-1}\sum_v \mathbb{E}[s_{v,H}(\theta_H^\dagger)] = \gamma_1(H)/V + o(V^{-1})$ (adding $-\mathcal{H}_H^{-1}\gamma_1(H)$ to $b_1(H)$, vanishing under $o(V^{-1})$) plus higher-order expansion regularity across Prop. 2(ii)(c), Prop. 2(iii), and §C.11(ii)(c)–(iii), and noted competing-root ties in $\mathcal{V}$ remain unresolved at $O(V^{-1})$.
  - **Bookkeeping 1 — Normalization consistency:** Harmonized Main L371/L463 and Supplement §C.11(ii)(a) (L412–425), using $\mu_q := \lim_{V\to\infty} V^{-1}\sum_v\mathbb{E}[q_v] > 0$ (avoiding notation collision with Assumption A2's upper bound $\bar{q}$) so $\widehat{U}_H \xrightarrow{p} Q_H/\mu_q$ with the same minimizer, aligned $q_v$ in L463, and moved the eligibility indicator $E_v = \mathbf{1}_{\{n_v \ge 2\}}$ and denominator inside the unconditional expectation in §C.11(ii)(a) (noting $\tilde{q}_v$ exits only when both $q_v$ and eligibility are deterministic/conditioned on).
  - **Bookkeeping 2 — Fund-size weighting & "can fail":** Replaced `"$O(1/n_v)$ idiosyncratic noise"` / `"$\sigma_{\nu,v}^2/n_v$"` with `"idiosyncratic diagonal contribution"` ($h_v\sigma_{\nu,v}^2$, equal to $\sigma_{\nu,v}^2/n_v$ under EW) at Main L109, L319, L407 and Supplement Table C.1; changed `"fails"` to `"can fail"` at Main L484 and softened Remark D.1 L745 (`"are generally nonzero and need only average to zero"`).

- [x] **R3.5 — Compile, QA, Cover Letter Sync & Archive Update:**
  - Recompiled `jbes_main_blind.pdf` (**40 pages**, 32 lines on p. 40, 3 figures, 5 tables, **186-word abstract**, **0 LaTeX/BibTeX warnings**, **0 overfull hboxes >10pt**, **0 `??`**, clean blind scan) and `jbes_supplement_blind.pdf` (**68 pages**, **0 overfull hboxes >10pt**, **0 `??`**, clean blind scan; Bao & Ullah 2007 verified in bibliography).
  - Recompiled `jbes_cover_letter.pdf` (**2 pages**) with synced counts (40 pp. main, 186-word abstract, 68 pp. supplement).
  - Archived `jbes_third_round_referee_report.pdf` into `2025 sdf paper/jbes/archive/`, updated `jbes_referee_response_decisions.md` and `jbes_revision_plan.md` §8.12. (Zero R scripts or CSVs modified; `jbes_replication.zip` unaffected.)
